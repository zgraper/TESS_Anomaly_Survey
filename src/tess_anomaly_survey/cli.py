"""Command-line interface for the exploratory TESS survey and SQLite candidate library."""
from __future__ import annotations
import argparse
from collections import defaultdict
import csv
import hashlib
import json
import math
import random
from pathlib import Path
import sys
import numpy as np
import experiment_b0 as b0
import experiment_b1 as b1
from . import __version__
from . import core, db, streaming


def options(parser,settings=True):
    parser.add_argument('--root',required=True,help='Root containing processed NPZs')
    parser.add_argument('--catalog',required=True,help='CSV with tic_id,sector,processed_file')
    if settings:
        parser.add_argument('--length',type=int,default=512)
        parser.add_argument('--stride',type=int,default=256)
        parser.add_argument('--min-good',type=float,default=.97)
        parser.add_argument('--max-windows-per-star',type=int,default=0,
                            help='0=all valid overlapping windows; use same value for calibration and scan')


def get_settings(args):
    values={'length':args.length,'stride':args.stride,'min_good':args.min_good,
            'max_windows':args.max_windows_per_star}
    if args.length!=512:
        raise ValueError('Frozen B1 transit-template library is validated at 512 cadences only')
    if args.stride<1 or args.max_windows_per_star<0 or not(.9<=args.min_good<=1.0):
        raise ValueError('Invalid scanning settings')
    return values


def require_calibration(cal_path,args):
    cal=core.read_json(cal_path)
    if cal.get('format_version')!='SURVEY_CALIBRATION_V1':
        raise RuntimeError('Not a survey calibration (B1.2 5% thresholds are not valid survey 1% thresholds)')
    if cal['b0_sha256']!=core.ORIGINAL_B0_SHA256 or cal['b1_sha256']!=core.ORIGINAL_B1_SHA256:
        raise RuntimeError('Detector source hash mismatch')
    for k,v in cal['settings'].items():
        # Scan settings come from calibration artifact to enforce identical window definitions.
        if k not in {'length','stride','min_good','max_windows'}:
            raise ValueError('Unexpected scanner settings')
    return cal




def do_split_catalog(args):
    if not (0 < args.calibration_fraction < 1 and 0 < args.validation_fraction < 1 and
            args.calibration_fraction + args.validation_fraction < 1):
        raise ValueError('Fractions must be positive, and calibration + validation < 1')
    groups=defaultdict(list)
    for row in core.eligible_catalog_rows(args.catalog,args.root):
        if row.get('status','ok')!='ok':
            continue
        groups[b0.star_key(row)].append(row)
    if len(groups)<3:
        raise RuntimeError('Need at least 3 distinct TICs for independent splits')
    keys=sorted(groups)
    random.Random(args.seed).shuffle(keys)
    a=round(len(keys)*args.calibration_fraction)
    b=round(len(keys)*args.validation_fraction)
    if min(a,b,len(keys)-a-b)<1:
        raise RuntimeError('Insufficient distinct stars for requested split proportions')
    subsets={'calibration':keys[:a],'validation':keys[a:a+b],'survey':keys[a+b:]}
    dest=Path(args.output_dir)
    dest.mkdir(parents=True,exist_ok=True)
    filenames={k:dest/(k+'_catalog.csv') for k in subsets}
    manifest=dest/'catalog_split_manifest.json'
    if manifest.exists() or any(p.exists() for p in filenames.values()):
        raise RuntimeError('Output split catalog or manifest exists; refusing overwrite')
    summary={}
    for group,ids in subsets.items():
        rows=[r for key in ids for r in groups[key]]
        fields=list(rows[0])
        with filenames[group].open('w',newline='',encoding='utf-8') as h:
            writer=csv.DictWriter(h,fieldnames=fields)
            writer.writeheader();writer.writerows(rows)
        summary[group]={'unique_tics':len(ids),'observations':len(rows),
                        'sha256':core.sha256_file(filenames[group])}
    core.write_json(manifest,{'created_utc':core.now(),'source_sha256':core.sha256_file(args.catalog),
       'seed':args.seed,'star_disjoint':True,'subsets':summary,
       'warning':'Exploratory random split. Stars previously used in historical development or final tests are NOT magically new independent confirmations.'})
    print('Created star-disjoint exploratory catalogs:',json.dumps(summary,indent=2))

def do_calibrate(args):
    settings=get_settings(args)
    if args.method=='planet_template':
        core.template_cache(args.method,1/720)  # Fail before committing partial results.
    # Aim for at least ~100 samples in the empirical upper tail.
    args.min_windows=max(args.min_windows, int(math.ceil(100/args.target_alert_fraction))) if args.target_alert_fraction>0 else args.min_windows
    rows=list(core.eligible_catalog_rows(args.catalog,args.root))
    if not rows:
        raise RuntimeError('No processed source files found')
    bootstrap = getattr(args, 'bootstrap', 200)
    seed = getattr(args, 'seed', 2026)
    cal = core.build_calibration(rows, args.root, args.method, settings,
                                args.target_alert_fraction, args.grouping,
                                args.min_stars, args.min_windows,
                                args.allow_small_calibration,
                                bootstrap=bootstrap, seed=seed)
    cal['calibration_catalog_sha256'] = core.sha256_file(args.catalog)
    cal['calibration_root'] = str(Path(args.root).resolve())
    core.write_json(args.output, cal)
    print(f"Calibrated {cal['n_windows']} windows on {cal['n_stars']} unique stars")
    print('Thresholds:', cal['thresholds'])
    if cal.get('uncertainty') and 'pooled' in cal['uncertainty']:
        u = cal['uncertainty']['pooled']
        print(f"Bootstrap 95% CI: [{u['ci_95'][0]}, {u['ci_95'][1]}] (std err: {u['std_err']}, relative: {u['relative_stderr']*100:.2f}%)")
    print('Saved:',args.output)
    if cal['small_calibration_test_only']:
        print('WARNING: TEST-ONLY calibration; not evidence of a controlled alert rate',file=sys.stderr)


def do_validate(args):
    cal=require_calibration(args.calibration,args)
    if cal['method']=='planet_template': core.template_cache(cal['method'],1/720)
    forbidden=set(cal['stars'])
    per_star=defaultdict(lambda:[0,0])
    skipped=0
    cache={}
    for row in core.eligible_catalog_rows(args.catalog,args.root):
        tic=b0.star_key(row)
        if tic in forbidden:
            raise RuntimeError(f'Validation data leak: {tic} used during calibration')
        try:
            for w in core.iter_windows(Path(args.root)/row['processed_file'],**cal['settings']):
                _,noise,score=core.score_window(w,cal['method'],cache)
                threshold,_=core.threshold_for(cal,noise)
                per_star[tic][0]+=1
                per_star[tic][1]+=int(score>threshold)
        except Exception as e:
            skipped+=1
    n=sum(x[0] for x in per_star.values()); k=sum(x[1] for x in per_star.values())
    if n<1:
        raise RuntimeError('No validation windows')
    rng=np.random.default_rng(args.seed)
    matrix=np.asarray(list(per_star.values()),dtype=float)
    boot=[]
    for _ in range(args.bootstrap):
        batch=matrix[rng.integers(len(matrix),size=len(matrix))].sum(axis=0)
        boot.append(batch[1]/batch[0] if batch[0] else 0.0)
    result={'format_version':'SURVEY_VALIDATION_V1','calibration_sha256':core.sha256_file(args.calibration),
            'validation_catalog_sha256':core.sha256_file(args.catalog),
            'target_alert_fraction':cal['target_alert_fraction'],'observed_flagged_fraction':k/n,
            'flagged_windows':k,'total_windows':n,'distinct_stars':len(per_star),
            'star_bootstrap_95_interval':list(map(float,np.quantile(boot,[.025,.975]))),
            'skipped_files':skipped,'created_utc':core.now(),
            'warning':'Window-level alert fraction on PRESUMED backgrounds; multiple windows per TIC correlated. Independent validation does not prove astrophysical false-discovery rates.'}
    core.write_json(args.output,result)
    print(f"Star-disjoint observed alert-window fraction: {k}/{n} = {100*k/n:.3f}%")
    print('Star-bootstrap interval:',result['star_bootstrap_95_interval'])
    print('Saved:',args.output)


def do_scan(args):
    cal=require_calibration(args.calibration,args)
    if cal['method']=='planet_template': core.template_cache(cal['method'],1/720)
    calhash=core.sha256_file(args.calibration)
    cathash=core.sha256_file(args.catalog)
    run_id=hashlib.sha256((calhash+'|'+cathash+'|'+str(Path(args.root).resolve())+'|'+
             str(args.candidates_per_window)).encode()).hexdigest()[:16]
    con=db.connect(args.db)
    db.begin_run(con,run_id,calhash,cal,dict(cal['settings'],
                 candidates_per_window=args.candidates_per_window),args.catalog,core.now())
    cache={}; bad=0;completed=0
    allowed=set(cal['stars'])
    for row in core.eligible_catalog_rows(args.catalog,args.root):
        source=row['processed_file']; tic=b0.star_key(row)
        if tic in allowed:
            raise RuntimeError(f'Survey data leak: {tic} belongs to calibration cohort')
        if db.is_done(con,run_id,source):
            if args.retry_errors:
                with con:
                    con.execute('DELETE FROM processed_files WHERE run_id=? AND source_file=? AND status=?', (run_id,source,'error'))
                if db.is_done(con,run_id,source): continue
            else:
                continue
        if args.limit_files and completed>=args.limit_files:
            break
        count=0; flagged=0
        rowdata={'tic':tic,'sector':row['sector']}
        try:
            with con:
                for w in core.iter_windows(Path(args.root)/source,**cal['settings']):
                    raw,noise,top=core.score_window(w,cal['method'],cache)
                    cut,g=core.threshold_for(cal,noise)
                    alert=int(top>cut)
                    con.execute('''INSERT INTO windows(run_id,source_file,tic,start_index,
                       start_btjd,noise,max_score,cutoff,noise_bin,flagged)
                       VALUES (?,?,?,?,?,?,?,?,?,?)''',(
                       run_id,source,tic,int(w['start']),w['window_start_btjd'],
                       noise,top,cut,g,alert))
                    count+=1;flagged+=alert
                    if alert:
                        candidates=b1.deduplicate(raw,w['cadence_days'],cut)
                        for candidate in candidates[:args.candidates_per_window]:
                            db.upsert_candidate(con,run_id,rowdata,source,candidate,w['cadence_days'])
                con.execute('INSERT INTO processed_files VALUES (?,?,?,?,?,?,?,?)',(
                    run_id,source,tic,int(row['sector']),count,flagged,'ok',None))
            completed+=1
        except Exception as e:
            with con:
                # Do not silently commit partial observations from a failed star.
                con.execute('INSERT INTO processed_files VALUES (?,?,?,?,?,?,?,?)',(
                    run_id,source,tic,int(row['sector']),0,0,'error',str(e)[:900]))
            bad+=1;completed+=1
        if completed%args.log_every==0:
            print(f'New files processed: {completed}; new errors: {bad}; run {run_id}',flush=True)
    print('Resume-safe run ID:',run_id)
    for stat in db.statistics(con):
        if stat['run_id']==run_id:
            print(json.dumps(stat,indent=2))
    con.close()


def do_stats(args):
    cx=db.connect(args.db)
    for s in db.statistics(cx):
        print(json.dumps(s,indent=2))
    cx.close()


def do_export(args):
    from . import ranking
    cx = db.connect(args.db)
    group_by = getattr(args, 'group_by', 'none')
    if getattr(args, 'group_by_star', False):
        group_by = 'star'
    elif getattr(args, 'group_by_system', False):
        group_by = 'system'

    ranking_method = getattr(args, 'ranking', 'raw')

    if group_by == 'system':
        systems = ranking.export_periodic_systems(
            cx, run_id=args.run_id, method=ranking_method, limit=args.limit
        )
        cx.close()
        cols = ['system_id', 'run_id', 'tic', 'is_periodic', 'classification', 'period_days',
                'period_uncertainty_days', 'epoch_btjd', 'epoch_uncertainty_days',
                'alternate_period_days', 'n_events', 'period_snr', 'odd_even_sigma',
                'note', 'created_utc']
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        with Path(args.output).open('w', newline='', encoding='utf-8') as fh:
            writer = csv.writer(fh)
            writer.writerow(cols)
            for s in systems:
                writer.writerow([s.get(c) for c in cols])
        print(f'Exported {len(systems)} candidate systems to {args.output} (grouped by periodic system)')
        return

    cond = 'WHERE run_id=?' if args.run_id else ''
    vals = [args.run_id] if args.run_id else []
    if args.status:
        cond += (' AND ' if cond else ' WHERE ') + 'review_status=?'
        vals.append(args.status)
    cur = cx.execute('''SELECT id,run_id,tic,sector,source_file,center_btjd,score,width_hr,
                     depth_ppm,n_support_windows,review_status,note
                     FROM candidates ''' + cond + ' ORDER BY score DESC', vals)
    raw_cands = []
    for r in cur:
        raw_cands.append({
            'id': r[0], 'run_id': r[1], 'tic': r[2], 'sector': r[3],
            'source_file': r[4], 'center_btjd': r[5], 'score': r[6],
            'width_hr': r[7], 'depth_ppm': r[8], 'n_support_windows': r[9],
            'review_status': r[10], 'note': r[11]
        })
    cx.close()

    if group_by == 'star':
        star_recs = ranking.group_candidates_by_star(
            raw_cands, method=ranking_method, limit=args.limit
        )
        cols = ['tic', 'sector', 'candidate_count', 'max_score', 'mean_score',
                'max_normalized_score', 'first_center_btjd', 'last_center_btjd',
                'median_depth_ppm', 'median_width_hr', 'top_candidate_id', 'review_status']
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        with Path(args.output).open('w', newline='', encoding='utf-8') as fh:
            writer = csv.writer(fh)
            writer.writerow(cols)
            for st in star_recs:
                writer.writerow([st.get(c) for c in cols])
        print(f'Exported {len(star_recs)} unique stars to {args.output} (grouped by star, ranked by {ranking_method})')
        return

    # Default: candidate-level export with preserved raw scores and optional normalized ranking
    max_per_star = getattr(args, 'max_per_star', 2)
    ranked_cands = ranking.compute_candidate_ranking(
        raw_cands, method=ranking_method, max_per_star=max_per_star, limit=args.limit
    )
    cols = ['id', 'run_id', 'tic', 'sector', 'source_file', 'center_btjd', 'score',
            'normalized_score', 'width_hr', 'depth_ppm', 'n_support_windows',
            'review_status', 'note']
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with Path(args.output).open('w', newline='', encoding='utf-8') as fh:
        writer = csv.writer(fh)
        writer.writerow(cols)
        for c in ranked_cands:
            writer.writerow([c.get(col) for col in cols])
    diversity_msg = f' (max {max_per_star} per TIC star)' if max_per_star > 0 else ' (unrestricted)'
    print(f'Exported {len(ranked_cands)} ranked candidates to {args.output}{diversity_msg} (ranking={ranking_method}); this is a review cap, not a calibrated alert rate')


def do_review(args):
    from .web import SurveyDatabaseAdapter
    SurveyDatabaseAdapter(args.db).record_event_review(args.candidate_id, args.status, args.note)
    print(f'Candidate {args.candidate_id} review updated to {args.status}')


def do_plot(args):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    cx=db.connect(args.db)
    row=cx.execute('SELECT source_file,center_btjd,tic,score FROM candidates WHERE id=?',(
        args.candidate_id,)).fetchone()
    if not row: raise ValueError('Unknown candidate ID')
    source,center,tic,score=row
    t,f,good,_=core.series(Path(args.root)/source)
    delta=abs(t-center)*24
    idx=(delta<=args.hours)&good
    if not np.any(idx):raise RuntimeError('No valid data around candidate center')
    fig,ax=plt.subplots(figsize=(9,4))
    ax.scatter(t[idx],f[idx]*1e6,s=7,alpha=.65)
    ax.axvline(center,color='tab:orange',linestyle='--',label='Candidate center')
    ax.set(xlabel='TESS BTJD',ylabel='Relative flux (ppm)',
           title=f'{tic} | candidate {args.candidate_id} | score {score:.2f}')
    ax.legend(loc='best');fig.tight_layout()
    Path(args.output).parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(args.output,dpi=160);plt.close(fig)
    print('Plot saved:',args.output)


def do_analyze_periods(args):
    from . import periodicity
    if args.min_events < 3:
        raise ValueError("Periodicity analysis requires at least 3 candidate events (min_events >= 3)")
    cx = db.connect(args.db)
    cond = "WHERE run_id=?" if args.run_id else ""
    vals = [args.run_id] if args.run_id else []

    cur = cx.execute(f'''SELECT id, run_id, tic, center_btjd, score, width_hr, depth_ppm
        FROM candidates {cond} ORDER BY run_id, tic, center_btjd''', vals)

    stars = defaultdict(list)
    for r in cur:
        stars[(r[1], r[2])].append({
            'id': r[0], 'run_id': r[1], 'tic': r[2], 'center_btjd': r[3],
            'score': r[4], 'width_hr': r[5], 'depth_ppm': r[6]
        })

    n_analyzed = 0
    n_periodic = 0
    for (run_id, tic), events in stars.items():
        if len(events) < args.min_events:
            continue
        n_analyzed += 1
        cand_ids = [e['id'] for e in events]
        result = periodicity.find_periodicity(
            events, min_events=args.min_events,
            max_period_days=args.max_period, min_period_days=args.min_period
        )
        sys_id = db.save_candidate_system(cx, run_id, tic, result, cand_ids)
        if result.is_periodic:
            n_periodic += 1
            p_err = f" +/- {result.period_uncertainty_days:.4f}" if result.period_uncertainty_days is not None else ""
            t0_err = f" +/- {result.epoch_uncertainty_days:.4f}" if result.epoch_uncertainty_days is not None else ""
            alt_str = f", P_alt={result.alternate_period_days:.4f} d" if result.alternate_period_days is not None else ""
            tag = result.classification.upper() if getattr(result, 'classification', None) else "PERIODIC"
            print(f"TIC {tic} (run {run_id}): [{tag}] P={result.period_days:.4f}{p_err} d, T0={result.epoch_btjd:.4f}{t0_err}, N={result.n_events}, SNR={result.period_snr:.1f}{alt_str}{' ' + result.note if result.note else ''}")

    cx.close()
    print(f"Periodicity analysis complete: {n_analyzed} stars analyzed (>= {args.min_events} events), {n_periodic} periodic candidate systems detected.")


def do_review_ui(args):
    from . import web
    web.serve_review_ui(
        args.db,
        root_path=args.root,
        host=args.host,
        port=args.port,
        cache_dir=getattr(args, 'cache_dir', None),
        review_db=getattr(args, 'review_db', None),
        read_only=getattr(args, 'read_only', False),
        experiments_dir=getattr(args, 'experiments_dir', 'experiments')
    )


def do_benchmark(args):
    from . import benchmark
    depths = [float(x.strip()) for x in args.depths.split(',')] if args.depths else None
    durations = [float(x.strip()) for x in args.durations.split(',')] if args.durations else None
    shapes = [x.strip() for x in args.shapes.split(',')] if args.shapes else None

    print(f"Starting injection-recovery benchmark on {args.catalog}...")
    manifest = benchmark.run_benchmark(
        catalog_path=args.catalog,
        root=args.root,
        calibration_path=args.calibration,
        depths_ppm=depths,
        durations_hr=durations,
        shapes=shapes,
        n_injections_per_cell=args.n_injections,
        seed=args.seed,
        output_path=args.output
    )
    print(f"Benchmark complete: {manifest['total_evaluations']} injections evaluated across {len(manifest['grid_results'])} parameter cells.")
    print(f"Results saved to {args.output} and {Path(args.output).with_suffix('.csv')}")


def do_stream(args):
    import sqlite3
    db_path = args.db

    # Status check operation
    if getattr(args, 'status', False):
        if not Path(db_path).is_file():
            print(f"Streaming results database not found: {db_path}")
            return
        con = sqlite3.connect(db_path)
        obs_total = con.execute("SELECT count(*) FROM observations").fetchone()[0]
        obs_comp = con.execute("SELECT count(*) FROM observations WHERE status='completed'").fetchone()[0]
        obs_zero = con.execute("SELECT count(*) FROM observations WHERE status='no_anomalies'").fetchone()[0]
        obs_fail = con.execute("SELECT count(*) FROM observations WHERE status='failed'").fetchone()[0]
        cands_total = con.execute("SELECT count(*) FROM candidates").fetchone()[0]
        cands_ver = con.execute("SELECT count(*) FROM candidates WHERE scientific_classification='verified_localized_dimming'").fetchone()[0]
        
        bench = con.execute("SELECT total_wall_sec, avg_obs_per_hour, avg_download_mb_per_sec, peak_ram_mb, avg_cpu_pct FROM performance_benchmarks ORDER BY created_utc DESC LIMIT 1").fetchone()
        
        print("="*60)
        print(f"PHASE 7 STREAMING SURVEY STATUS: {db_path}")
        print("="*60)
        print(f"Total Observations Recorded : {obs_total}")
        print(f" - Completed with Anomalies : {obs_comp}")
        print(f" - Completed (No Anomalies) : {obs_zero}")
        print(f" - Failed Observations      : {obs_fail}")
        print(f"Total Candidate Events      : {cands_total}")
        print(f" - Verified Localized Dips  : {cands_ver}")
        if bench:
            print(f"Latest Benchmark Telemetry  :")
            print(f" - Total Wall Time          : {bench[0]:.1f} s")
            print(f" - Processing Throughput    : {bench[1]:.1f} obs/hr")
            print(f" - Network Throughput       : {bench[2]:.2f} MB/s")
            print(f" - Peak RAM Usage           : {bench[3]:.1f} MB")
            print(f" - Avg CPU Utilization      : {bench[4]:.1f}%")
        print("="*60)
        con.close()
        return

    # Export operation
    if getattr(args, 'export', None):
        if not Path(db_path).is_file():
            raise FileNotFoundError(f"Database not found: {db_path}")
        con = sqlite3.connect(db_path)
        import pandas as pd
        df = pd.read_sql_query("SELECT * FROM candidates", con)
        df.to_csv(args.export, index=False)
        print(f"Exported {len(df)} candidate records to {args.export}")
        con.close()
        return

    # Targets gathering
    targets = []
    if getattr(args, 'manifest', None):
        with open(args.manifest, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                targets.append(row)
    elif getattr(args, 'tic', None):
        targets.append({'tic': args.tic, 'sector': getattr(args, 'sector', None)})
    elif getattr(args, 'sector', None):
        raise NotImplementedError("For sector scans, provide a target manifest CSV.")
    else:
        raise ValueError("Must specify --manifest, --tic, --status, or --export")

    if getattr(args, 'dry_run', False):
        dev_tics = streaming.load_development_universe() if not getattr(args, 'allow_any_target', False) else set()
        print(f"[DRY-RUN] Loaded {len(targets)} targets from input.")
        if dev_tics:
            in_dev = [t for t in targets if str(t.get('tic_id') or t.get('tic')).replace('TIC:', '').strip() in dev_tics]
            print(f"[DRY-RUN] Targets within 8,984 development universe: {len(in_dev)} / {len(targets)}")
            if len(in_dev) < len(targets):
                print(f"[DRY-RUN] WARNING: {len(targets) - len(in_dev)} targets outside development universe!")
        return

    run_id = args.run_id or f"stream_{int(time.time())}"
    metrics = streaming.stream_pipeline(
        targets=targets,
        db_path=db_path,
        run_id=run_id,
        temp_dir=args.temp_dir,
        retain_fits=args.retain_fits,
        resume=not args.force,
        check_dev_universe=not getattr(args, 'allow_any_target', False),
        rate_limit_delay_sec=args.rate_limit_delay,
        max_observations=args.limit or 0,
        max_workers=getattr(args, 'workers', 1),
        log_every=args.log_every
    )
    print("\n" + "="*60)
    print(f"STREAMING SURVEY RUN COMPLETE: {run_id}")
    print("="*60)
    print(f"Processed Observations : {metrics['n_processed']}")
    print(f"Successful             : {metrics['n_successful']}")
    print(f"Failed                 : {metrics['n_failed']}")
    print(f"Skipped (Resumed)      : {metrics['n_skipped']}")
    print(f"Wall Clock Time        : {metrics['total_wall_sec']:.1f} s")
    print(f"Processing Throughput  : {metrics['avg_obs_per_hour']:.1f} obs/hr")
    print(f"Download Throughput    : {metrics['avg_download_mb_sec']:.2f} MB/s")
    print(f"Peak RAM Usage         : {metrics['peak_ram_mb']:.1f} MB")
    print(f"Average CPU Load       : {metrics['avg_cpu_pct']:.1f}%")
    print("="*60)


def build_parser():
    p=argparse.ArgumentParser(description='TESS anomaly survey v'+__version__)
    sub=p.add_subparsers(dest='command',required=True)
    sp=sub.add_parser('split-catalog',help='Reproducible star-disjoint exploratory catalog split')
    options(sp,settings=False)
    sp.add_argument('--output-dir',required=True)
    sp.add_argument('--calibration-fraction',type=float,default=.45)
    sp.add_argument('--validation-fraction',type=float,default=.20)
    sp.add_argument('--seed',type=int,default=2026)
    sp.set_defaults(action=do_split_catalog)
    ca=sub.add_parser('calibrate',help='Fit empirical alert cutoffs using unmodified calibration stars')
    options(ca)
    ca.add_argument('--method',choices=core.METHODS,default='fixed_box_2h')
    ca.add_argument('--target-alert-fraction',type=float,default=.01)
    ca.add_argument('--grouping',choices=['pooled','noise_tertiles'],default='pooled')
    ca.add_argument('--min-stars',type=int,default=1000)
    ca.add_argument('--min-windows',type=int,default=10000)
    ca.add_argument('--allow-small-calibration',action='store_true',help='For software smoke tests only')
    ca.add_argument('--bootstrap',type=int,default=200,help='Star-level bootstrap iterations for threshold uncertainty (0 to disable)')
    ca.add_argument('--seed',type=int,default=2026,help='Random seed for bootstrap resampling')
    ca.add_argument('--output',required=True)
    ca.set_defaults(action=do_calibrate)
    va=sub.add_parser('validate',help='Test window alert incidence on independent TICs')
    options(va,settings=False)
    va.add_argument('--calibration',required=True)
    va.add_argument('--bootstrap',type=int,default=300)
    va.add_argument('--seed',type=int,default=2026)
    va.add_argument('--output',required=True)
    va.set_defaults(action=do_validate)
    sc=sub.add_parser('scan',help='Scan stars and build resumable SQLite anomaly library')
    options(sc,settings=False)
    sc.add_argument('--calibration',required=True)
    sc.add_argument('--db',required=True)
    sc.add_argument('--limit-files',type=int,default=0)
    sc.add_argument('--retry-errors',action='store_true',help='Retry previously failed files without affecting completed files')
    sc.add_argument('--candidates-per-window',type=int,default=3)
    sc.add_argument('--log-every',type=int,default=50)
    sc.set_defaults(action=do_scan)
    st=sub.add_parser('stats',help='Summarize window and candidate counts')
    st.add_argument('--db',required=True);st.set_defaults(action=do_stats)
    ex=sub.add_parser('export',help='Export ranked candidates from library')
    ex.add_argument('--db',required=True);ex.add_argument('--output',required=True)
    ex.add_argument('--limit',type=int,default=100)
    ex.add_argument('--run-id')
    ex.add_argument('--status',choices=['unreviewed','interesting','known','artifact','dismissed'])
    ex.add_argument('--max-per-star',type=int,default=2,
                    help='Maximum candidates exported per TIC star (default 2; set 0 for unrestricted)')
    ex.add_argument('--ranking',choices=['raw','normalized'],default='raw',
                    help='Candidate ranking: raw matched-filter score (default) or star-normalized score')
    ex.add_argument('--group-by',choices=['none','star','system'],default='none',
                    help='Grouping for export: none (candidate level), star (star level), or system (periodic systems)')
    ex.add_argument('--group-by-star',action='store_true',
                    help='Convenience flag to group export by star')
    ex.add_argument('--group-by-system',action='store_true',
                    help='Convenience flag to group export by periodic system')
    ex.set_defaults(action=do_export)
    rv=sub.add_parser('review',help='Annotate a candidate after visual inspection')
    rv.add_argument('--db',required=True);rv.add_argument('--candidate-id',type=int,required=True)
    rv.add_argument('--status',required=True,choices=['unreviewed','interesting','known','artifact','dismissed'])
    rv.add_argument('--note',default='');rv.set_defaults(action=do_review)
    pl=sub.add_parser('plot',help='Create candidate inspection figure')
    pl.add_argument('--db',required=True);pl.add_argument('--root',required=True)
    pl.add_argument('--candidate-id',type=int,required=True)
    pl.add_argument('--hours',type=float,default=8)
    pl.add_argument('--output',required=True);pl.set_defaults(action=do_plot)
    ap=sub.add_parser('analyze-periods',help='Detect periodic recurring candidate events per TIC star')
    ap.add_argument('--db',required=True,help='Path to SQLite anomaly library')
    ap.add_argument('--run-id',help='Optional run ID filter')
    ap.add_argument('--min-events',type=int,default=3,help='Minimum events per star (strictly >= 3)')
    ap.add_argument('--max-period',type=float,default=20.0,help='Maximum search period in days (default: 20.0)')
    ap.add_argument('--min-period',type=float,default=0.2,help='Minimum search period in days (default: 0.2)')
    ap.set_defaults(action=do_analyze_periods)
    ui=sub.add_parser('review-ui',help='Launch local web browser review interface')
    ui.add_argument('--db',default=None,help='Path to SQLite anomaly library (optional, defaults to active experiment)')
    ui.add_argument('--root',help='Root directory containing processed TESS NPZs for light curve display')
    ui.add_argument('--host',default='127.0.0.1',help='Host binding address (default: 127.0.0.1)')
    ui.add_argument('--port',type=int,default=8000,help='Port to serve (default: 8000)')
    ui.add_argument('--cache-dir',help='Local cache directory for on-demand downloaded FITS products')
    ui.add_argument('--review-db',help='Path to external overlay SQLite database for non-destructive review storage')
    ui.add_argument('--read-only',action='store_true',help='Open survey database in strictly read-only mode')
    ui.add_argument('--experiments-dir',default='experiments',help='Directory for experiment metadata, configs, and catalogs (default: experiments)')
    ui.set_defaults(action=do_review_ui)
    bm=sub.add_parser('benchmark',help='Synthetic signal injection-recovery benchmark')
    bm.add_argument('--catalog',required=True,help='Path to validation/benchmark catalog CSV')
    bm.add_argument('--root',required=True,help='Root directory for TESS processed files')
    bm.add_argument('--calibration',required=True,help='Path to calibration JSON file')
    bm.add_argument('--depths',default='200,500,1000,2000,5000,10000',help='Comma-separated depths in ppm')
    bm.add_argument('--durations',default='0.5,1.0,2.0,4.0,8.0',help='Comma-separated durations in hours')
    bm.add_argument('--shapes',default='box,asymmetric',help='Comma-separated signal shapes (box, asymmetric)')
    bm.add_argument('--n-injections',type=int,default=20,help='Injections per parameter grid cell (default: 20)')
    bm.add_argument('--seed',type=int,default=2026,help='Random seed for reproducibility')
    bm.add_argument('--output',required=True,help='Path to output results JSON file')
    bm.set_defaults(action=do_benchmark)

    sm = sub.add_parser('stream', help='Streaming ingest and real-time processing of TESS light curves from MAST')
    sm.add_argument('--db', default='data/streaming_pilot.sqlite', help='Path to dedicated Phase 7 SQLite database')
    sm.add_argument('--manifest', help='Path to target manifest CSV with tic_id and optional sector')
    sm.add_argument('--tic', help='Single target TIC identifier')
    sm.add_argument('--sector', type=int, help='Optional observing sector filter')
    sm.add_argument('--run-id', help='Unique survey run identifier')
    sm.add_argument('--temp-dir', help='Custom temporary cache directory for light curve files')
    sm.add_argument('--retain-fits', action='store_true', help='Retain downloaded FITS files on disk (default: delete after processing)')
    sm.add_argument('--force', action='store_true', help='Force re-processing of already completed targets')
    sm.add_argument('--limit', type=int, default=0, help='Maximum number of observations to process (0 for unlimited)')
    sm.add_argument('--workers', type=int, default=2, help='Number of parallel archive acquisition and processing workers (default: 2)')
    sm.add_argument('--rate-limit-delay', type=float, default=0.05, help='Delay in seconds between remote downloads')
    sm.add_argument('--log-every', type=int, default=10, help='Progress reporting frequency')
    sm.add_argument('--allow-any-target', action='store_true', help='Bypass 8,984 development universe guard (CAUTION: only for pre-registered non-Phase 6B targets)')
    sm.add_argument('--dry-run', action='store_true', help='Validate targets and check development universe without downloading')
    sm.add_argument('--status', action='store_true', help='Display survey progress and benchmark telemetry for the specified database')
    sm.add_argument('--export', help='Export all candidate records from database to specified CSV file')
    sm.set_defaults(action=do_stream)

    return p


def main(argv=None):
    args=build_parser().parse_args(argv)
    core.frozen_check()
    if getattr(args,'bootstrap',1)<1 or getattr(args,'candidates_per_window',1)<1 or getattr(args,'log_every',1)<1:
        raise SystemExit('Counts must be positive')
    try:
        args.action(args)
    except (RuntimeError,ValueError,FileNotFoundError) as e:
        raise SystemExit('TESS SURVEY ERROR: '+str(e)) from e

if __name__=='__main__': main()
