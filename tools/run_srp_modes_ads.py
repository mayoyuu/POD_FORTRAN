#!/usr/bin/env python3
"""Audit, execute and atomically publish a separate C2 modes ADS run."""
import argparse
import csv
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT=Path(__file__).resolve().parents[1]


class RunError(RuntimeError):
    pass


def require(ok,message):
    if not ok: raise RunError(message)


def load_tool(name):
    spec=importlib.util.spec_from_file_location(name,ROOT/'tools'/f'{name}.py')
    module=importlib.util.module_from_spec(spec)
    sys.modules[name]=module
    spec.loader.exec_module(module)
    return module


AUDITOR=load_tool('audit_boxwing_mode_table')
RECONSTRUCTOR=load_tool('reconstruct_srp_ads_snapshot')


def file_sha256(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream,'sha256').hexdigest()


def read_config(path):
    values={}
    for line in Path(path).read_text(encoding='utf-8').splitlines():
        line=line.split('#',1)[0].strip()
        if not line: continue
        require('=' in line,'invalid config record')
        key,value=map(str.strip,line.split('=',1))
        require(key not in values,f'duplicate config key {key}')
        values[key]=value
    return values


def config_path(value):
    path=Path(value)
    return (path if path.is_absolute() else ROOT/path).resolve()


def preflight(opm,config):
    opm=Path(opm).resolve(); config=Path(config).resolve()
    cfg=read_config(config)
    require(cfg.get('srp_model')=='time_varying_modes','C2 config must select time_varying_modes')
    require(cfg.get('use_srp','').lower()=='true','C2 requires active SRP')
    require(cfg.get('srp_mode_table_file') and cfg.get('srp_mode_metadata_file'),'missing mode paths')
    csv_path=config_path(cfg['srp_mode_table_file']); meta_path=config_path(cfg['srp_mode_metadata_file'])
    try:
        audit=AUDITOR.audit_table(csv_path,meta_path)
    except (AUDITOR.AuditError,OSError) as exc:
        raise RunError(f'mode preflight: {exc}') from exc
    require(audit.get('audit_complete') is True and audit.get('gate_passed') is True,'mode data are not finalized')
    require(file_sha256(opm)==audit['opm_sha256'],'runtime OPM hash does not match mode identity')
    for key in ('srp_cannonball_cr','srp_cannonball_area_m2','srp_mass_kg','srp_pressure_1au_n_m2'):
        require(key in cfg,f'missing common physical parameter {key}')
        actual=float(cfg[key])
        if key=='srp_pressure_1au_n_m2' and actual==-1: actual=1367/299792458
        expected=audit['physical_parameters'][key]
        require(math.isfinite(actual) and math.isclose(actual,expected,rel_tol=1e-12),f'common physical parameter mismatch: {key}')
    source=read_config(audit['config_file'])
    keys=('LU','TU','MU','VU','AccU','use_earth_nspheric','use_moon_nspheric','use_third_body',
          'use_srp','use_drag','use_relativity','earth_degree','moon_degree','active_planets','reference_frame','central_body')
    for key in keys:
        require(key in cfg and key in source,f'non-SRP identity missing {key}')
        require(cfg[key].replace(',',' ').split()==source[key].replace(',',' ').split(),f'non-SRP dynamics mismatch: {key}')
    return audit,csv_path,meta_path


def unique_object(pairs):
    value={}
    for key,item in pairs:
        require(key not in value,f'duplicate JSON key {key}')
        value[key]=item
    return value


def load_json(path):
    def invalid(value): raise RunError(f'non-finite JSON constant {value}')
    return json.loads(Path(path).read_text(),object_pairs_hook=unique_object,parse_constant=invalid)


def rows(path):
    with Path(path).open(newline='',encoding='utf-8') as stream:
        reader=csv.DictReader(stream)
        require(reader.fieldnames and len(reader.fieldnames)==len(set(reader.fieldnames)),f'invalid CSV header: {path}')
        result=list(reader)
    for row in result:
        require(None not in row and all(value is not None and value.strip() for value in row.values()),f'invalid CSV fields: {path}')
        require(all(math.isfinite(float(value)) for value in row.values()),f'nonfinite CSV values: {path}')
    return result


def close(a,b,label,absolute=1e-13):
    require(math.isclose(float(a),float(b),rel_tol=1e-12,abs_tol=absolute),f'{label} mismatch: {a} versus {b}')


def bind_run_settings(report,expected):
    for key,value in expected.items():
        if key in ('rel_tol','abs_tol'):
            require(math.isclose(float(report.get(key,float('nan'))),value,rel_tol=1e-12,abs_tol=0.0),
                    f'run setting {key} mismatch')
        else:
            close(report.get(key,float('nan')),value,f'run setting {key}')


def audit_outputs(staging,report,csv_hash,meta_hash):
    staging=Path(staging).resolve()
    require(report.get('schema_version')=='srp_modes_ads_v1' and report.get('formal_result') is False,'invalid raw report contract')
    require(report.get('n_variables')==8 and report.get('parameterization')=='SRP_MODES','report is not 8D modes')
    require(report.get('variable_roles')==['x','y','z','vx','vy','vz','s','u'],'invalid variable roles')
    require(report.get('mode_csv_sha256')==csv_hash and report.get('mode_metadata_sha256')==meta_hash,'raw mode hashes mismatch')
    require(report.get('da_handles_final')==report.get('da_handles_baseline')==0,'DA handle lifecycle failed')
    require(report.get('validation_failed') is False,'ADS validation failed')
    duration=float(report['duration_hours']); cadence=float(report['save_hours'])
    require(math.isfinite(duration) and math.isfinite(cadence) and duration>0 and cadence>0,'invalid report time limits')
    intervals=round(duration/cadence)
    close(intervals*cadence,duration,'report cadence')
    history=rows(staging/'c2_history.csv')
    require(len(history)==intervals+1,'history checkpoint count mismatch')
    for step,row in enumerate(history):
        close(row['hour'],step*cadence,'history time')
        close(row['epoch_et'],float(report['epoch0_tdb_s'])+step*cadence*3600,'absolute TDB time',1e-6)
        for i in range(1,9): require(float(row[f'split_dim{i}'])>=0,'missing eight-direction split counts')
    initial=rows(staging/'c2_initial_points.csv')
    require(len(initial)==report.get('shape_count')==512,'invalid study point count')
    points={}
    validation_ids=set()
    for i,row in enumerate(initial,1):
        require(float(row['point_id'])==i,'initial point IDs are not unique/sequential')
        q=[float(row[f'u{k}']) for k in range(1,9)]
        require(all(abs(x)<=1+1e-13 for x in q),'initial point outside 8D box')
        close(row['s'],0.5*(1+q[6]),'physical s'); close(row['u'],0.5*(1+q[7]),'physical u')
        require(float(row['validation']) in (0,1),'invalid validation marker')
        points[i]=row
        if float(row['validation'])==1: validation_ids.add(i)
    require(len(validation_ids)==report.get('validation_count')==337,'invalid 8D validation grid')
    shape=rows(staging/'c2_shape.csv'); probes=rows(staging/'c2_probes.csv')
    require(len(shape)==512*len(history) and len(probes)==337*len(history),'point/checkpoint record counts mismatch')
    for data,ids,label in ((shape,set(points),'shape'),(probes,validation_ids,'probe')):
        seen=set()
        for row in data:
            hour=float(row['hour']); point_id=int(row['point_id']); step=round(hour/cadence)
            require(0<=step<=intervals and math.isclose(hour,step*cadence,abs_tol=1e-12),f'{label} off-grid time')
            require(point_id in ids and (step,point_id) not in seen and float(row['point_id'])==point_id,f'{label} invalid point identity')
            seen.add((step,point_id))
        require(len(seen)==len(history)*len(ids),f'{label} records are incomplete')
    per_hour={step:[0.0,0.0] for step in range(len(history))}
    for row in probes:
        point_id=int(row['point_id']); step=round(float(row['hour'])/cadence)
        close(row['s'],points[point_id]['s'],'probe s'); close(row['u'],points[point_id]['u'],'probe u')
        errors=[math.hypot(*(float(row[f'ads_{axis}'])-float(row[f'real_{axis}']) for axis in axes))
                for axes in (('x','y','z'),('vx','vy','vz'))]
        close(row['ads_numerical_error_pos_km'],errors[0],'probe position error')
        close(row['ads_numerical_error_vel_kms'],errors[1],'probe velocity error',1e-15)
        per_hour[step]=[max(a,b) for a,b in zip(per_hour[step],errors)]
    max_errors=[max(values[i] for values in per_hour.values()) for i in (0,1)]
    for step,row in enumerate(history):
        close(row['ads_numerical_error_pos_km'],per_hour[step][0],'history probe position')
        close(row['ads_numerical_error_vel_kms'],per_hour[step][1],'history probe velocity',1e-15)
    close(report['ads_numerical_error']['max_position_km'],max_errors[0],'report probe position')
    close(report['ads_numerical_error']['max_velocity_kms'],max_errors[1],'report probe velocity',1e-15)
    require(max_errors[0]<=report['position_tolerance_km'] and max_errors[1]<=report['velocity_tolerance_kms'],'ADS numerical gate failed')
    zero=[max(float(row[key]) for row in history) for key in ('zero_mode_pos_km','zero_mode_vel_kms')]
    close(report['zero_mode_error']['max_position_km'],zero[0],'zero position')
    close(report['zero_mode_error']['max_velocity_kms'],zero[1],'zero velocity',1e-15)
    require(zero[0]<=1e-6 and zero[1]<=1e-9,'zero-mode control gate failed')
    events=rows(staging/'c2_split_events.csv')
    require(len(events)==report['split_event_count']==sum(sum(float(row[f'split_dim{i}']) for i in range(1,9)) for row in history),'split event totals mismatch')
    for event_id,row in enumerate(events,1):
        require(float(row['event_id'])==event_id and 1<=float(row['direction'])<=8,'invalid 8D event direction/id')
        for i in range(1,9):
            require(float(row[f'parent_width{i}'])>0 and float(row[f'direction_score{i}'])>=0,'invalid event patch/score')
    depth=sum(float(row['depth_limited']) for row in history)
    require(report['depth_limited_patches']==depth and report['precision_unconfirmed']==(depth>0),'depth-limited precision flag mismatch')
    prefixes=report.get('snapshot_prefixes',[])
    require(prefixes and len(prefixes)==len(set(prefixes)),'missing or duplicate snapshots')
    snapshot_audits=[]; snapshot_times=set()
    for prefix in prefixes:
        path=Path(prefix)
        if not path.is_absolute(): path=staging/path
        require(path.resolve().parent==staging,'snapshot path escapes staging')
        snapshot=RECONSTRUCTOR.load_snapshot(path)
        require(snapshot.schema_version==2 and snapshot.n_variables==8,'snapshot is not v2 8D')
        require(snapshot.mode_hashes=={'mode_csv_sha256':csv_hash,'mode_metadata_sha256':meta_hash},'snapshot mode identity mismatch')
        step=round(snapshot.time_seconds/(cadence*3600))
        require(0<=step<=intervals,'snapshot time outside run')
        close(snapshot.time_seconds,step*cadence*3600,'snapshot checkpoint')
        nominal=[float(history[step][key]) for key in ('nominal_x_km','nominal_y_km','nominal_z_km','nominal_vx_kms','nominal_vy_kms','nominal_vz_kms')]
        _,state=snapshot.evaluate_global((0,)*6+(-1,-1))
        close(math.dist(state[:3],nominal[:3]),history[step]['zero_mode_pos_km'],'snapshot zero position',1e-10)
        close(math.dist(state[3:],nominal[3:]),history[step]['zero_mode_vel_kms'],'snapshot zero velocity',1e-13)
        reconstructed=RECONSTRUCTOR.reconstruct_snapshot(path,points_per_patch=512,verify=True)
        snapshot_audits.append(reconstructed); snapshot_times.add(step)
    require({0,intervals}<=snapshot_times,'initial/final snapshots missing')
    return snapshot_audits


def environment():
    env=os.environ.copy()
    libraries=[ROOT/'lib',Path.home()/'local_dace/lib',Path.home()/'local_openblas/lib']
    for key in ('LD_LIBRARY_PATH','LIBRARY_PATH'):
        env[key]=':'.join(map(str,libraries))+':'+env.get(key,'')
    return env


def source_hashes():
    suffixes={'.f','.for','.f90','.f95','.f03','.f08','.fpp','.c','.cpp','.h','.hpp','.inc'}
    files=[path for path in (ROOT/'src').rglob('*') if path.is_file() and path.suffix.lower() in suffixes]
    files += [ROOT/'app/run_srp_modes_ads.f90',ROOT/'fpm.toml']
    files += [ROOT/'tools'/name for name in ('run_srp_modes_ads.py','audit_boxwing_mode_table.py','reconstruct_srp_ads_snapshot.py',
        'generate_srp_modes_probability_samples.py','pack_srp_modes_probability_hdf5.py','requirements-srp-modes-probability.txt')]
    return {str(path.relative_to(ROOT)):file_sha256(path) for path in sorted(files)}


def build_executable():
    result=subprocess.run(['fpm','run','--target','run_srp_modes_ads','--','--help'],cwd=ROOT,env=environment(),capture_output=True,text=True)
    require(result.returncode==0 and '8D' in result.stdout and '--hours' in result.stdout,'C2 target was not built: '+result.stdout[-3000:]+result.stderr[-3000:])
    paths=list((ROOT/'build').glob('gfortran_*/app/run_srp_modes_ads'))
    require(len(paths)==1,'expected one bound C2 executable; pass --executable for another compiler/profile')
    return paths[0].resolve()


def relative_paths(value,staging):
    if isinstance(value,dict): return {k:relative_paths(v,staging) for k,v in value.items()}
    if isinstance(value,list): return [relative_paths(v,staging) for v in value]
    if isinstance(value,str) and value.startswith(str(staging)+os.sep):
        return Path(value).relative_to(staging).as_posix()
    return value


def probability_disk_preflight(parent,hours,cadence,n,v,da_order,max_depth):
    while not parent.exists(): parent=parent.parent
    checkpoints=round(hours/cadence)+1
    fixed={0.0,hours}|{h for h in (0,24,72,168,240,360) if h<=hours and math.isclose(h/cadence,round(h/cadence),abs_tol=1e-10)}
    # At most two extra frames for each of the three first-event classes.
    keyframes=min(checkpoints,len(fixed)+6)
    point_records=checkpoints*v+keyframes*n
    csv_bytes=point_records*500
    hdf5_bytes=point_records*100+n*200
    # Conservative all-patch full polynomial upper bound, not only saved points.
    snapshot_bytes=checkpoints*(2**max_depth)*6*math.comb(da_order+8,8)*160
    required=csv_bytes+hdf5_bytes+snapshot_bytes+n*800+2_000_000_000
    free=shutil.disk_usage(parent).free
    require(free>=required,f'insufficient probability output disk space: need {required}, available {free}')
    return {'estimated_csv_bytes':csv_bytes,'estimated_hdf5_bytes':hdf5_bytes,
            'conservative_snapshot_bytes':snapshot_bytes,'required_bytes_with_buffer':required,'available_bytes':free}


def run_probability_native(args,staging,log):
    with subprocess.Popen(args,cwd=ROOT,env=environment(),stdout=log,stderr=subprocess.STDOUT) as process:
        try:
            while True:
                require(shutil.disk_usage(staging).free>=2_000_000_000,'probability disk space buffer exhausted')
                try: return process.wait(timeout=1)
                except subprocess.TimeoutExpired: pass
        except BaseException:
            process.terminate()
            try: process.wait(timeout=10)
            except subprocess.TimeoutExpired: process.kill(); process.wait()
            raise


def run_modes(opm,config,output_dir,hours=360.0,save_hours=1.0,da_order=4,max_depth=8,
              pos_tol_km=0.1,vel_tol_kms=1e-6,executable=None,rel_tol=1e-12,abs_tol=1e-12,
              probability_samples=None,probability_manifest=None,probability_count=50000,
              probability_visualization_count=2000,probability_validation_count=64,probability_keep_csv=False):
    output=Path(output_dir).absolute()
    require(not output.exists(),f'output directory already exists: {output}')
    require(all(math.isfinite(x) and x>0 for x in (hours,save_hours,pos_tol_km,vel_tol_kms,rel_tol,abs_tol))
            and da_order>=1 and max_depth>=0,'invalid run limits')
    require(bool(probability_samples)==bool(probability_manifest),'probability samples and manifest must be supplied together')
    probability_enabled=probability_samples is not None
    probability_tool=None; disk_estimate=None
    if probability_enabled:
        require(type(probability_count) is int and probability_count>=2 and
                type(probability_visualization_count) is int and probability_visualization_count>0 and
                type(probability_validation_count) is int and probability_validation_count>0,'invalid probability counts')
        try:
            probability_tool=load_tool('pack_srp_modes_probability_hdf5')
            population,probability_identity=probability_tool.SAMPLER.load_sample_bundle(
                probability_samples,probability_manifest,probability_count)
        except (ImportError,OSError,ValueError) as exc:
            raise RunError(f'probability preflight (use isolated NumPy/h5py environment): {exc}') from exc
        disk_estimate=probability_disk_preflight(output.parent,hours,save_hours,probability_count,
            min(probability_visualization_count,probability_count),da_order,max_depth)
    opm=Path(opm).resolve(); config=Path(config).resolve()
    audit,csv_path,meta_path=preflight(opm,config)
    inputs={'opm':opm,'config':config,'mode_csv':csv_path,'mode_metadata':meta_path}
    if probability_enabled:
        inputs.update(probability_samples=Path(probability_samples).resolve(),probability_manifest=Path(probability_manifest).resolve())
    hashes={key:file_sha256(path) for key,path in inputs.items()}
    sources=source_hashes()
    exe=Path(executable).resolve() if executable else build_executable()
    require(source_hashes()==sources,'sources changed during build')
    exe_hash=file_sha256(exe)
    revision=subprocess.run(['git','rev-parse','HEAD'],cwd=ROOT,check=True,text=True,capture_output=True).stdout.strip()
    dirty=bool(subprocess.run(['git','status','--porcelain'],cwd=ROOT,check=True,text=True,capture_output=True).stdout.strip())
    output.parent.mkdir(parents=True,exist_ok=True)
    lock=output.parent/f'.{output.name}.publish.lock'
    try:
        lock_fd=os.open(lock,os.O_CREAT|os.O_EXCL|os.O_WRONLY,0o600)
    except FileExistsError as exc:
        raise RunError(f'publication already locked: {lock}') from exc
    staging=None
    try:
        require(not output.exists(),f'output directory already exists: {output}')
        staging=Path(tempfile.mkdtemp(prefix=f'.{output.name}.staging-',dir=output.parent)).resolve()
        bound=staging/'runner'; shutil.copy2(exe,bound)
        require(file_sha256(bound)==exe_hash,'executable changed during freezing')
        args=[str(bound),'-opm',str(opm),'-cfg',str(config),'-o',str(staging/'c2'),
              '--hours',str(hours),'--save-hours',str(save_hours),'--da-order',str(da_order),
              '--max-depth',str(max_depth),'--pos-tol-km',str(pos_tol_km),'--vel-tol-kms',str(vel_tol_kms),
              '--rel-tol',str(rel_tol),'--abs-tol',str(abs_tol),
              '--mode-meta-sha256',hashes['mode_metadata'],'--source-revision',revision]
        if probability_enabled:
            sample_copy=staging/'probability_samples.csv'; manifest_copy=staging/'probability_samples.json'
            shutil.copy2(inputs['probability_samples'],sample_copy); shutil.copy2(inputs['probability_manifest'],manifest_copy)
            require(file_sha256(sample_copy)==hashes['probability_samples'] and
                    file_sha256(manifest_copy)==hashes['probability_manifest'],'probability input changed while freezing')
            args += ['--probability-samples',str(sample_copy),'--probability-count',str(probability_count),
                     '--probability-visualization-count',str(probability_visualization_count),
                     '--probability-validation-count',str(probability_validation_count)]
        with (staging/'run.log').open('w') as log:
            if probability_enabled: returncode=run_probability_native(args,staging,log)
            else: returncode=subprocess.run(args,cwd=ROOT,env=environment(),stdout=log,stderr=subprocess.STDOUT).returncode
        require(returncode==0,f'C2 executable failed ({returncode})')
        report=load_json(staging/'c2_report.json')
        expected={'duration_hours':hours,'save_hours':save_hours,'da_order':da_order,'max_depth':max_depth,
                  'position_tolerance_km':pos_tol_km,'velocity_tolerance_kms':vel_tol_kms,
                  'rel_tol':rel_tol,'abs_tol':abs_tol}
        bind_run_settings(report,expected)
        report['snapshot_audits']=audit_outputs(staging,report,hashes['mode_csv'],hashes['mode_metadata'])
        if probability_enabled:
            prob_audit=probability_tool.audit_probability_run(staging,report,population,probability_validation_count,
                RECONSTRUCTOR,rows(staging/'c2_history.csv'),rows)
            require(report['probability_outputs']['visualization_sample_count']==min(probability_visualization_count,probability_count),
                    'native probability visualization setting mismatch')
            require(prob_audit['numerical_validation_passed'],'probability numerical/precision gate failed')
            opm_data=load_json(opm)
            covariance_scale=[math.sqrt(opm_data[key]) for key in
                ('CX_X','CY_Y','CZ_Z','CX_DOT_X_DOT','CY_DOT_Y_DOT','CZ_DOT_Z_DOT')]
            cfg=read_config(config)
            probability_outputs=probability_tool.pack_probability(staging/'c2',sample_copy,manifest_copy,metadata={
                'frame':cfg['reference_frame'],'center':cfg['central_body'],'time_scale':'TDB_seconds',
                'covariance_scale':covariance_scale,'input_sha256':hashes,'source_sha256':sources,
                'executable_sha256':exe_hash,'epoch0_tdb_s':report['epoch0_tdb_s']})
            probability_outputs.update(prob_audit,disk_preflight=disk_estimate)
            report['probability_outputs']=probability_outputs
            require(file_sha256(sample_copy)==hashes['probability_samples'] and
                    file_sha256(manifest_copy)==hashes['probability_manifest'],'frozen probability identity changed')
        require({key:file_sha256(path) for key,path in inputs.items()}==hashes,'inputs changed during run')
        require(source_hashes()==sources and file_sha256(exe)==exe_hash and file_sha256(bound)==exe_hash,'source/executable identity changed during run')
        if probability_enabled and not probability_keep_csv:
            # Only owned staging artifacts, already audited byte-for-value against recoverable HDF5.
            for name in ('c2_probability_visualization.csv','c2_probability_keyframes.csv'):
                (staging/name).unlink()
            report['probability_outputs']['staging_points_csv_removed_after_verified_hdf5']=True
        report.update(formal_result=True,source_revision=revision)
        report['provenance']={'inputs':{key:{'path':str(path),'sha256':hashes[key]} for key,path in inputs.items()},
            'executable':{'path':'runner','source_path':str(exe),'sha256':exe_hash},
            'source_files_sha256':sources,'working_tree_dirty':dirty,'mode_audit_gate_passed':audit['gate_passed']}
        report['mode_csv_file']=str(csv_path); report['mode_metadata_file']=str(meta_path)
        report=relative_paths(report,staging)
        for path in staging.glob('*_reconstruction.json'):
            path.write_text(json.dumps(relative_paths(load_json(path),staging),indent=2,allow_nan=False)+'\n')
        report['artifacts_sha256']={path.name:file_sha256(path) for path in sorted(staging.iterdir()) if path.is_file() and path.name!='c2_report.json'}
        (staging/'c2_report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
        require(not output.exists(),f'output directory already exists: {output}')
        staging.rename(output)
        return report
    except (OSError,ValueError,KeyError,subprocess.SubprocessError,RECONSTRUCTOR.SnapshotError,RunError) as exc:
        raise RunError(f'{exc}; diagnostics preserved at {staging}') from exc
    finally:
        os.close(lock_fd); lock.unlink()


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--opm',required=True,type=Path); parser.add_argument('--config',required=True,type=Path)
    parser.add_argument('--output-dir',required=True,type=Path); parser.add_argument('--executable',type=Path)
    parser.add_argument('--hours',type=float,default=360); parser.add_argument('--save-hours',type=float,default=1)
    parser.add_argument('--da-order',type=int,default=4); parser.add_argument('--max-depth',type=int,default=8)
    parser.add_argument('--pos-tol-km',type=float,default=0.1); parser.add_argument('--vel-tol-kms',type=float,default=1e-6)
    parser.add_argument('--rel-tol',type=float,default=1e-12); parser.add_argument('--abs-tol',type=float,default=1e-12)
    parser.add_argument('--probability-samples',type=Path); parser.add_argument('--probability-manifest',type=Path)
    parser.add_argument('--probability-count',type=int,default=50000)
    parser.add_argument('--probability-visualization-count',type=int,default=2000)
    parser.add_argument('--probability-validation-count',type=int,default=64)
    parser.add_argument('--probability-keep-csv',action='store_true')
    args=vars(parser.parse_args(argv))
    try:
        report=run_modes(**args)
    except (RunError,OSError,ValueError) as exc:
        parser.exit(2,f'C2 audit/run failed: {exc}\n')
    print(json.dumps({'output_dir':str(args['output_dir']),'ads_numerical_error':report['ads_numerical_error'],
                      'zero_mode_error':report['zero_mode_error'],'precision_unconfirmed':report['precision_unconfirmed']},indent=2))
    return 0


if __name__=='__main__': raise SystemExit(main())
