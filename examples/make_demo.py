"""Generate fabricated demo NPZs and two disjoint catalogs; NOT real astronomy data."""
from pathlib import Path
import csv
import numpy as np

root=Path(__file__).resolve().parent/'demo_data'
root.mkdir(parents=True,exist_ok=True)

def write(name, records):
    with (root/name).open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=['tic_id','sector','processed_file','status','split'])
        w.writeheader();w.writerows(records)

rng=np.random.default_rng(42)
cal=[];survey=[]
for i in range(40):
    tic=950000+i
    t=1600+np.arange(2000)/720
    flux=rng.normal(0,.0006,len(t))
    if i>=20:
        flux[900:950]-=.008 if i%3==0 else 0
    path=root/f'lc_{tic}.npz'
    np.savez_compressed(path,time_btjd=t,pdc_flux_rel=flux,
                        pdc_flux_err_rel=np.full(len(t),.0006),
                        good_pdc=np.ones(len(t),bool),quality=np.zeros(len(t),int))
    item={'tic_id':tic,'sector':1 if i<20 else 11,
          'processed_file':path.name,'status':'ok','split':'train' if i<20 else 'test'}
    (cal if i<20 else survey).append(item)
write('calibration.csv',cal)
write('survey.csv',survey)
print('Created TEST-ONLY demo data in',root)
print('Calibrate with --allow-small-calibration; never report as scientific evidence.')
