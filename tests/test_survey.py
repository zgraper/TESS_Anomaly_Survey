import csv
import json
from pathlib import Path
import sqlite3
import numpy as np
import pytest
from tess_anomaly_survey import core
from tess_anomaly_survey.cli import main


def make_file(root,name,tic,depth=0):
    path=root/'processed'/f'{name}.npz'
    path.parent.mkdir(parents=True,exist_ok=True)
    n=1100
    t=1000+np.arange(n,dtype=float)/720
    flux=np.random.default_rng(int(tic)).normal(0,.0005,n)
    if depth:
        flux[625:680]-=depth
    good=np.ones(n,bool)
    np.savez_compressed(path,time_btjd=t,pdc_flux_rel=flux,
                        pdc_flux_err_rel=np.full(n,.0005),good_pdc=good,
                        quality=np.zeros(n,dtype=np.int32))
    return f'processed/{name}.npz'


def write_catalog(path,records):
    with path.open('w',encoding='utf-8',newline='') as h:
        w=csv.DictWriter(h,fieldnames=['tic_id','sector','processed_file','status','split'])
        w.writeheader();w.writerows(records)


def test_frozen_sources():
    core.frozen_check()


def test_quantile_strict():
    assert core.empirical_quantile([1,2,3,4,5],.8)==4.0
    assert core.empirical_quantile([1,2,3,4,5],.99)==5.0


def test_gap_preservation(tmp_path):
    p=tmp_path/'gap.npz'
    t=1000+np.arange(1300)/720
    t[700:]+=1
    np.savez_compressed(p,time_btjd=t,pdc_flux_rel=np.zeros(1300),
                        good_pdc=np.ones(1300,bool))
    windows=list(core.iter_windows(p,stride=256))
    assert windows and all((w['time'][-1]-w['time'][0])<1 for w in windows)
    assert not any(w['start']<700 and w['start']+512>700 for w in windows)


def test_full_pipeline(tmp_path):
    root=tmp_path/'data';root.mkdir()
    cal_records=[]; eval_records=[]
    for i in range(7):
        ident=str(112200+i)
        file=make_file(root,'cal_'+str(i),ident)
        cal_records.append({'tic_id':ident,'sector':'1','processed_file':file,'status':'ok','split':'train'})
    for i in range(5):
        ident=str(113300+i)
        file=make_file(root,'scan_'+str(i),ident,depth=.003)
        eval_records.append({'tic_id':ident,'sector':'11','processed_file':file,'status':'ok','split':'test'})
    calib=tmp_path/'cal.csv'; survey=tmp_path/'survey.csv'
    write_catalog(calib,cal_records);write_catalog(survey,eval_records)
    artifact=tmp_path/'thresholds.json'
    main(['calibrate','--root',str(root),'--catalog',str(calib),
          '--target-alert-fraction','.1','--grouping','noise_tertiles',
          '--min-stars','1','--min-windows','1','--allow-small-calibration',
          '--output',str(artifact)])
    cal=json.loads(artifact.read_text())
    assert cal['n_stars']==7 and cal['method']=='fixed_box_2h'
    assert all(g in cal['thresholds'] for g in ('0','1','2'))
    validate=tmp_path/'validation.json'
    main(['validate','--root',str(root),'--catalog',str(survey),
          '--calibration',str(artifact),'--bootstrap','20','--output',str(validate)])
    assert json.loads(validate.read_text())['distinct_stars']==5
    db=tmp_path/'anomalies.sqlite'
    command=['scan','--root',str(root),'--catalog',str(survey),
             '--calibration',str(artifact),'--db',str(db)]
    main(command+['--limit-files','2'])
    main(command)
    main(command)
    cx=sqlite3.connect(db)
    assert cx.execute('SELECT COUNT(*) FROM processed_files').fetchone()[0]==5
    assert cx.execute('SELECT COUNT(*) FROM windows').fetchone()[0]>5
    assert cx.execute('SELECT COUNT(*) FROM candidates').fetchone()[0]>=1
    cx.close()
    out=tmp_path/'candidates.csv'
    main(['export','--db',str(db),'--output',str(out),'--limit','2'])
    assert out.read_text().startswith('id,run_id,tic')
    cid=int(list(csv.DictReader(out.open()))[0]['id'])
    main(['review','--db',str(db),'--candidate-id',str(cid),'--status','interesting','--note','Test'])
    assert sqlite3.connect(db).execute('SELECT review_status FROM candidates WHERE id=?',(cid,)).fetchone()[0]=='unreviewed'
    from tess_anomaly_survey.web import SurveyDatabaseAdapter
    overlay = SurveyDatabaseAdapter(db).review_db_path
    assert sqlite3.connect(overlay).execute('SELECT status FROM candidate_review_overlay WHERE candidate_id=?',(cid,)).fetchone()[0]=='interesting'
    with pytest.raises(SystemExit,match='data leak'):
        main(['validate','--root',str(root),'--catalog',str(calib),
              '--calibration',str(artifact),'--output',str(tmp_path/'bad.json')])


def test_calibration_requires_enough_stars(tmp_path):
    root=tmp_path/'data';root.mkdir()
    file=make_file(root,'f','7711')
    catalog=tmp_path/'c.csv'
    write_catalog(catalog,[{'tic_id':'7711','sector':'1','processed_file':file,'status':'ok','split':'train'}])
    with pytest.raises(SystemExit,match='Calibration insufficient'):
        main(['calibrate','--root',str(root),'--catalog',str(catalog),
              '--output',str(tmp_path/'cal.json')])


def test_split_catalog_is_star_disjoint(tmp_path):
    root=tmp_path/'data';root.mkdir()
    rows=[]
    for i in range(12):
        tid=str(44000+i)
        f=make_file(root,'file'+str(i),tid)
        rows.append({'tic_id':tid,'sector':'1','processed_file':f,'status':'ok','split':'train'})
    cat=tmp_path/'all.csv';write_catalog(cat,rows)
    out=tmp_path/'splits'
    main(['split-catalog','--root',str(root),'--catalog',str(cat),
          '--output-dir',str(out),'--calibration-fraction','.4',
          '--validation-fraction','.25'])
    sets=[]
    for group in ['calibration','validation','survey']:
        sets.append({r['tic_id'] for r in csv.DictReader((out/(group+'_catalog.csv')).open())})
    assert sum(map(len,sets))==12
    assert not (sets[0]&sets[1] or sets[1]&sets[2] or sets[0]&sets[2])
    with pytest.raises(SystemExit,match='refusing overwrite'):
        main(['split-catalog','--root',str(root),'--catalog',str(cat),'--output-dir',str(out)])
