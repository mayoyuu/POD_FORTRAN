#!/usr/bin/env python3
"""Publish independent nominal Real SRP endpoint truth, not a control-domain rerun."""
import argparse
import csv
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('srp_truth_common',ROOT/'tools/run_srp_modes_ads.py')
COMMON=importlib.util.module_from_spec(spec); spec.loader.exec_module(COMMON)
RunError=COMMON.RunError
require=COMMON.require
STATE_KEYS=('x_km','y_km','z_km','vx_kms','vy_kms','vz_kms')
SCENARIOS={'Deta','proxy_earth','proxy_moon','boxwing_earth','boxwing_moon'}


def preflight(opm,config):
    audit,csv_path,meta_path=COMMON.preflight(opm,config)
    cfg=COMMON.read_config(config)
    for key,expected in audit['physical_parameters'].items():
        require(key in cfg,f'truth physical parameter missing: {key}')
        if isinstance(expected,str):
            require(cfg[key]==expected,f'truth physical parameter mismatch: {key}')
        else:
            actual=[float(x) for x in cfg[key].replace(',',' ').split()]
            if key=='srp_pressure_1au_n_m2' and actual==[-1.0]: actual=[1367/299792458]
            values=expected if isinstance(expected,list) else [expected]
            require(len(actual)==len(values) and all(math.isfinite(a) and math.isclose(a,b,rel_tol=1e-12,abs_tol=1e-15)
                    for a,b in zip(actual,values)),f'truth physical parameter mismatch: {key}')
    return audit,csv_path,meta_path


def source_hashes():
    result=COMMON.source_hashes()
    for name in ('app/run_srp_modes_truth.f90','tools/run_srp_modes_truth.py'):
        result[name]=COMMON.file_sha256(ROOT/name)
    return result


def validate_truth_output(folder,report,hours,save_hours,epoch0,initial):
    require(report.get('schema_version')=='srp_modes_truth_v1' and
            report.get('state_source')=='independent_real_nominal','invalid Real truth schema')
    require(report.get('da_handles_baseline')==report.get('da_handles_final')==0,'truth DA handle leak')
    fallback=report.get('geometry_fallback_count')
    require(type(fallback) is int and fallback>=0,'invalid geometry fallback diagnostic')
    COMMON.close(report.get('epoch0_tdb_s',math.nan),epoch0,'truth epoch',1e-6)
    with (Path(folder)/'truth_states.csv').open(newline='') as stream:
        reader=csv.DictReader(stream)
        require(reader.fieldnames==['scenario','hour','epoch_et',*STATE_KEYS],'invalid truth CSV schema')
        data=list(reader)
    intervals=round(hours/save_hours)
    require(len(data)==report.get('row_count')==len(SCENARIOS)*(intervals+1),'incomplete truth records')
    seen=set()
    for row in data:
        require(None not in row and all(v is not None for v in row.values()),'invalid truth CSV fields')
        values=[float(row[k]) for k in ('hour','epoch_et',*STATE_KEYS)]
        require(all(math.isfinite(x) for x in values),'nonfinite truth state')
        hour,epoch,*state=values
        step=round(hour/save_hours); key=(row['scenario'],step)
        require(row['scenario'] in SCENARIOS and 0<=step<=intervals and key not in seen,'invalid truth identity')
        require(math.isclose(hour,step*save_hours,abs_tol=1e-12),'off-grid truth hour')
        require(math.isclose(epoch,epoch0+hour*3600,rel_tol=0,abs_tol=1e-6),'truth epoch mismatch')
        if step==0:
            require(all(math.isclose(a,b,rel_tol=0,abs_tol=1e-9) for a,b in zip(state,initial)),'truth initial state mismatch')
        seen.add(key)
    return data


def build_executable():
    result=subprocess.run(['fpm','run','--target','run_srp_modes_truth','--','--help'],cwd=ROOT,
        env=COMMON.environment(),text=True,capture_output=True)
    require(result.returncode==0 and 'Real-only' in result.stdout,'Real truth target unavailable: '+result.stdout[-3000:]+result.stderr[-1000:])
    binaries=list((ROOT/'build').glob('gfortran_*/app/run_srp_modes_truth'))
    require(len(binaries)==1,'expected one Real truth executable')
    return binaries[0].resolve()


def run_truth(opm,config,output_dir,hours=360.0,save_hours=1.0,rel_tol=1e-13,abs_tol=1e-13,executable=None):
    output=Path(output_dir).absolute()
    require(not output.exists(),f'output directory exists: {output}')
    require(all(math.isfinite(v) and v>0 for v in (hours,save_hours,rel_tol,abs_tol)) and hours<=360 and
            hours/save_hours<=10000 and math.isclose(round(hours/save_hours)*save_hours,hours,abs_tol=1e-10),'invalid truth limits')
    opm=Path(opm).resolve(); config=Path(config).resolve()
    audit,mode_csv,metadata=preflight(opm,config)
    inputs={'opm':opm,'config':config,'mode_csv':mode_csv,'mode_metadata':metadata}
    identities={k:COMMON.file_sha256(p) for k,p in inputs.items()}
    initial_opm=COMMON.load_json(opm)
    initial=[float(initial_opm[k]) for k in ('X','Y','Z','X_DOT','Y_DOT','Z_DOT')]
    sources=source_hashes()
    exe=Path(executable).resolve() if executable else build_executable()
    require(source_hashes()==sources,'truth sources changed during build')
    exe_hash=COMMON.file_sha256(exe)
    output.parent.mkdir(parents=True,exist_ok=True)
    lock=output.parent/f'.{output.name}.publish.lock'
    try: lock_fd=os.open(lock,os.O_CREAT|os.O_EXCL|os.O_WRONLY,0o600)
    except FileExistsError as exc: raise RunError(f'truth publication already locked: {lock}') from exc
    staging=None
    try:
        require(not output.exists(),f'output directory exists: {output}')
        staging=Path(tempfile.mkdtemp(prefix=f'.{output.name}.staging-',dir=output.parent)).resolve()
        frozen=staging/'runner'; shutil.copy2(exe,frozen)
        require(COMMON.file_sha256(frozen)==exe_hash,'truth executable changed during freezing')
        args=[str(frozen),'-opm',str(opm),'-cfg',str(config),'-o',str(staging/'truth'),
              '--hours',str(hours),'--save-hours',str(save_hours),'--rel-tol',str(rel_tol),'--abs-tol',str(abs_tol)]
        with (staging/'run.log').open('w') as log:
            result=subprocess.run(args,cwd=ROOT,env=COMMON.environment(),stdout=log,stderr=subprocess.STDOUT)
        require(result.returncode==0,f'truth executable failed ({result.returncode})')
        report=COMMON.load_json(staging/'truth_report.json')
        require(report.get('formal_result') is False,'invalid raw truth report')
        COMMON.bind_run_settings(report,{'duration_hours':hours,'save_hours':save_hours,'rel_tol':rel_tol,'abs_tol':abs_tol})
        validate_truth_output(staging,report,hours,save_hours,audit['epoch0_tdb_s'],initial)
        require({k:COMMON.file_sha256(p) for k,p in inputs.items()}==identities,'truth inputs changed during run')
        require(source_hashes()==sources and COMMON.file_sha256(exe)==exe_hash and
                COMMON.file_sha256(frozen)==exe_hash,'truth source/executable identity changed during run')
        report.update(formal_result=True,scenarios=sorted(SCENARIOS),
            frame='J2000',center='EARTH',timescale='TDB',illumination='full_sun',
            scenario_parameters={'Deta':{'eta':0.3},'proxy_earth':{'s':1,'u':0.5},'proxy_moon':{'s':0,'u':1},
                                 'boxwing_earth':{'attitude':'earth'},'boxwing_moon':{'attitude':'moon'}})
        report['provenance']={'inputs':{k:{'path':str(p),'sha256':identities[k]} for k,p in inputs.items()},
            'executable':{'path':'runner','sha256':exe_hash},'source_files_sha256':sources,
            'physical_parameters':audit['physical_parameters']}
        report['artifacts_sha256']={p.name:COMMON.file_sha256(p) for p in sorted(staging.iterdir())
                                    if p.is_file() and p.name!='truth_report.json'}
        (staging/'truth_report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
        require(not output.exists(),f'output directory exists: {output}')
        staging.rename(output)
        return report
    except (OSError,ValueError,KeyError,subprocess.SubprocessError,RunError) as exc:
        raise RunError(f'{exc}; diagnostics preserved at {staging}') from exc
    finally:
        os.close(lock_fd); lock.unlink()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for flag in ('opm','config','output-dir'): parser.add_argument('--'+flag,type=Path,required=True)
    parser.add_argument('--hours',type=float,default=360); parser.add_argument('--save-hours',type=float,default=1)
    parser.add_argument('--rel-tol',type=float,default=1e-13); parser.add_argument('--abs-tol',type=float,default=1e-13)
    parser.add_argument('--executable',type=Path)
    try: report=run_truth(**vars(parser.parse_args()))
    except (RunError,OSError,ValueError) as exc: parser.exit(2,str(exc)+'\n')
    print(json.dumps({'formal_result':report['formal_result'],'hours':report['duration_hours'],'row_count':report['row_count']}))


if __name__=='__main__': main()
