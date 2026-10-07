#!/usr/bin/env python3
"""Run and audit the 15-orbit paired SRP ADS Monte Carlo campaign."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
import time
import uuid
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
import h5py
import numpy as np

from pack_srp_mc_points_hdf5 import pack_points
from summarize_srp_mc_pair import summarize_files,write_flat_csv

ORBITS=(
    'DRO-1','DRO-2','DRO-3',
    'L1Halo-1','L1Halo-2','L1Halo-3',
    'L1PL-1','L1PL-2','L1PL-3',
    'L2Halo-1','L2Halo-2','L2Halo-3',
    'L2PL-1','L2PL-2','L2PL-3',
)
OUTPUT_SCHEMA='srp_mc_campaign_v1'
COVARIANCE_KEYS=(
    ('CX_X',0,0),('CY_X',1,0),('CY_Y',1,1),
    ('CZ_X',2,0),('CZ_Y',2,1),('CZ_Z',2,2),
    ('CX_DOT_X',3,0),('CX_DOT_Y',3,1),('CX_DOT_Z',3,2),
    ('CY_DOT_X',4,0),('CY_DOT_Y',4,1),('CY_DOT_Z',4,2),
    ('CZ_DOT_X',5,0),('CZ_DOT_Y',5,1),('CZ_DOT_Z',5,2),
    ('CX_DOT_X_DOT',3,3),('CY_DOT_X_DOT',4,3),('CY_DOT_Y_DOT',4,4),
    ('CZ_DOT_X_DOT',5,3),('CZ_DOT_Y_DOT',5,4),('CZ_DOT_Z_DOT',5,5),
)


class CampaignError(RuntimeError):
    pass


def file_sha256(path: Path) -> str:
    digest=hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda:handle.read(1024*1024),b''):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path,data: dict) -> None:
    path.parent.mkdir(parents=True,exist_ok=True)
    fd,name=tempfile.mkstemp(prefix=path.name+'.',suffix='.tmp',dir=path.parent)
    try:
        with os.fdopen(fd,'w',encoding='utf-8') as handle:
            json.dump(data,handle,indent=2,ensure_ascii=False,allow_nan=False)
            handle.write('\n')
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(name,path)
    except BaseException:
        try: os.unlink(name)
        except FileNotFoundError: pass
        raise


def atomic_csv(path: Path,rows: list[dict],fields: list[str]) -> None:
    path.parent.mkdir(parents=True,exist_ok=True)
    fd,name=tempfile.mkstemp(prefix=path.name+'.',suffix='.tmp',dir=path.parent)
    try:
        with os.fdopen(fd,'w',newline='',encoding='utf-8') as handle:
            writer=csv.DictWriter(handle,fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(name,path)
    except BaseException:
        try: os.unlink(name)
        except FileNotFoundError: pass
        raise


def discover_opms(opm_root: str|Path) -> dict[str,Path]:
    root=Path(opm_root)
    result={}
    problems=[]
    for orbit in ORBITS:
        directory=root/orbit
        candidates=[directory/f'{orbit}_init.opm',
                    directory/f'{orbit}_init.opm.json']
        present=[path for path in candidates if path.is_file()]
        if len(present)!=1:
            problems.append(f'{orbit}: expected one direct init OPM, found {len(present)}')
        else:
            result[orbit]=present[0]
    if problems:
        raise CampaignError('; '.join(problems))
    return result


def read_covariance(path: Path) -> np.ndarray:
    try:
        data=json.loads(path.read_text(encoding='utf-8'))
    except (OSError,json.JSONDecodeError) as exc:
        raise CampaignError(f'invalid JSON OPM {path}: {exc}') from exc
    covariance=np.zeros((6,6),dtype=np.float64)
    for key,row,column in COVARIANCE_KEYS:
        try:
            value=float(data[key])
        except (KeyError,TypeError,ValueError) as exc:
            raise CampaignError(f'{path}: invalid covariance element {key}') from exc
        if not math.isfinite(value):
            raise CampaignError(f'{path}: nonfinite covariance element {key}')
        covariance[row,column]=value
        covariance[column,row]=value
    try:
        np.linalg.cholesky(covariance)
    except np.linalg.LinAlgError as exc:
        raise CampaignError(f'{path}: covariance is not positive definite') from exc
    return covariance


def audit_opm(path: Path,expected_pos_km: float=10.0,
              expected_vel_mps: float=0.03,
              relative_tolerance: float=1e-9) -> dict:
    covariance=read_covariance(path)
    pos=float(np.sqrt(np.trace(covariance[:3,:3])))
    vel=float(1000.0*np.sqrt(np.trace(covariance[3:,3:])))
    accepted=(math.isclose(pos,expected_pos_km,rel_tol=relative_tolerance,abs_tol=1e-12)
              and math.isclose(vel,expected_vel_mps,rel_tol=relative_tolerance,
                               abs_tol=1e-12))
    return {'opm':str(path),'opm_sha256':file_sha256(path),
            'position_rss_km':pos,'velocity_rss_mps':vel,
            'positive_definite':True,
            'status':'eligible' if accepted else 'initial_error_mismatch',
            'reason':'' if accepted else
                f'expected RSS {expected_pos_km:g} km and {expected_vel_mps:g} m/s'}


def publish_directory(staging: str|Path,final: str|Path,validator) -> None:
    staging=Path(staging)
    final=Path(final)
    valid,reasons=validator(staging)
    if not valid:
        raise CampaignError('staged output validation failed: '+';'.join(reasons))
    final.parent.mkdir(parents=True,exist_ok=True)
    backup=final.with_name(final.name+'.__old__'+uuid.uuid4().hex)
    if final.exists():
        final.rename(backup)
    try:
        staging.rename(final)
    except BaseException:
        if backup.exists() and not final.exists():
            backup.rename(final)
        raise
    if backup.exists():
        shutil.rmtree(backup)


def find_executable(project_root: Path) -> Path:
    candidates=[path for path in (project_root/'build').glob(
        '*/app/run_srp_ads_mc_pair') if path.is_file() and os.access(path,os.X_OK)]
    if not candidates:
        raise CampaignError('cannot locate build/*/app/run_srp_ads_mc_pair')
    return max(candidates,key=lambda path:path.stat().st_mtime)


def stage_executable(executable: Path,experiment_root: Path) -> Path:
    digest=file_sha256(executable)
    runtime=experiment_root/'.runtime'
    runtime.mkdir(parents=True,exist_ok=True)
    target=runtime/f'run_srp_ads_mc_pair_{digest[:16]}'
    if not target.exists() or file_sha256(target)!=digest:
        temporary=target.with_suffix('.tmp')
        shutil.copy2(executable,temporary)
        temporary.chmod(temporary.stat().st_mode|0o111)
        temporary.replace(target)
    return target


def inject_provenance(report_path: Path,orbit: str,opm: Path,config: Path,
                      samples: Path,executable: Path,args) -> dict:
    report=json.loads(report_path.read_text(encoding='utf-8'))
    report.update({
        'orbit_name':orbit,'source_opm':str(opm.resolve()),
        'opm_sha256':file_sha256(opm),'source_config':str(config.resolve()),
        'config_sha256':file_sha256(config),'sample_file':str(samples.resolve()),
        'sample_sha256':file_sha256(samples),'sample_seed':20261005,
        'executable_sha256':file_sha256(executable),
        'campaign_schema_version':OUTPUT_SCHEMA,
        'expected_hour_count':int(round(args.days*24.0/args.save_hours))+1,
    })
    atomic_json(report_path,report)
    return report


def _csv_rows(path: Path) -> list[dict[str,str]]:
    with path.open(newline='',encoding='utf-8') as handle:
        return list(csv.DictReader(handle))


def _all_numeric_finite(rows: list[dict[str,str]],text_fields=()) -> bool:
    for row in rows:
        for key,value in row.items():
            if key in text_fields:
                continue
            try:
                if not math.isfinite(float(value)):
                    return False
            except (TypeError,ValueError):
                return False
    return True


def organize_outputs(stage: Path,prefix: Path,orbit: str,args) -> None:
    visualization=Path(str(prefix)+'_mc_visualization_points.csv')
    keyframes=Path(str(prefix)+'_mc_keyframe_points.csv')
    report=Path(str(prefix)+'_mc_report.json')
    pack_points(visualization,stage/f'{orbit}_visualization_points.h5',
                min(args.visualization_count,args.mc_count),'visualization',report)
    pack_points(keyframes,stage/f'{orbit}_keyframe_points.h5',
                args.mc_count,'keyframe',report)
    expected_hours=int(round(args.days*24.0/args.save_hours))+1
    summary=summarize_files(Path(str(prefix)+'_mc_hourly_stats.csv'),
                            Path(str(prefix)+'_mc_paired_stats.csv'),report,
                            expected_hours=expected_hours,save_hours=args.save_hours)
    atomic_json(stage/f'{orbit}_summary.json',summary)
    write_flat_csv(summary,stage/f'{orbit}_summary.csv')

    mappings={
        Path(str(prefix)+'_mc_hourly_stats.csv'):stage/f'{orbit}_hourly_mc_stats.csv',
        Path(str(prefix)+'_mc_paired_stats.csv'):stage/f'{orbit}_paired_effect_stats.csv',
        report:stage/f'{orbit}_events.json',
        Path(str(prefix)+'_treatment_split_events.csv'):
            stage/f'{orbit}_treatment_split_events.csv',
        Path(str(prefix)+'_control_split_events.csv'):
            stage/f'{orbit}_control_split_events.csv',
        Path(str(prefix)+'_run.log'):stage/f'{orbit}_run.log',
    }
    for source,target in mappings.items():
        if not source.is_file():
            raise CampaignError(f'missing staged artifact: {source.name}')
        source.replace(target)
    visualization.unlink()
    keyframes.unlink()
    snapshot_dir=stage/'patch_snapshots'
    snapshot_dir.mkdir()
    for source in sorted(stage.glob(f'{orbit}_*_patch_s*')):
        source.replace(snapshot_dir/source.name)
    if not any(snapshot_dir.iterdir()):
        raise CampaignError('no patch snapshots were produced')


def _split_criterion_matches(row: dict[str,str],criterion: object) -> bool:
    if not isinstance(criterion,dict):
        return False
    scalar_fields={
        'component':'component','direction':'direction',
        'parent_depth':'parent_depth','bfs_iteration':'bfs_iteration',
    }
    try:
        for report_key,csv_key in scalar_fields.items():
            if int(criterion[report_key])!=int(row[csv_key]):
                return False
        arrays=(
            ('truncation_estimates','error',6),
            ('component_tolerances','tolerance',6),
            ('positive_excess','excess',6),
            ('direction_scores','direction_score',7),
            ('parent_unit_center','parent_center',7),
            ('parent_unit_width','parent_width',7),
        )
        for report_key,csv_prefix,count in arrays:
            values=criterion[report_key]
            if not isinstance(values,list) or len(values)!=count:
                return False
            for index,value in enumerate(values,1):
                if not math.isclose(float(value),float(row[f'{csv_prefix}{index}']),
                                    rel_tol=1e-12,abs_tol=1e-14):
                    return False
    except (KeyError,TypeError,ValueError,OverflowError):
        return False
    return True


def audit_split_group(group: str,group_report: object,event_rows: list[dict[str,str]],
                      hourly_rows: list[dict[str,str]]) -> list[str]:
    reasons=[]
    if group not in {'treatment','control'} or not isinstance(group_report,dict):
        return [f'{group}_event_report_missing']
    try:
        group_hourly=[row for row in hourly_rows if row['group']==group]
        finest=max(int(row['sample_count']) for row in group_hourly)
        finest_rows=[row for row in group_hourly
                     if int(row['sample_count'])==finest]
        expected_events=sum(sum(int(row[f'split_dim{i}']) for i in range(1,8))
                            for row in finest_rows)
    except (KeyError,ValueError,TypeError):
        return [f'{group}_split_hourly_parse']
    if expected_events!=len(event_rows):
        reasons.append(f'{group}_split_event_count')
    try:
        if [int(row['event_id']) for row in event_rows]!=list(
                range(1,len(event_rows)+1)):
            reasons.append(f'{group}_split_event_ids')
        if not _all_numeric_finite(event_rows):
            reasons.append(f'{group}_split_event_nonfinite')
    except (KeyError,ValueError,TypeError):
        reasons.append(f'{group}_split_event_parse')

    first=group_report.get('first_split')
    if not isinstance(first,dict) or bool(first.get('observed'))!=bool(event_rows):
        reasons.append(f'{group}_first_split_presence')
    elif event_rows:
        row=event_rows[0]
        try:
            matches=(math.isclose(float(first['interval_start_hour']),
                                  float(row['interval_start_hour']),
                                  rel_tol=0.0,abs_tol=1e-12) and
                     math.isclose(float(first['hour']),float(row['interval_end_hour']),
                                  rel_tol=0.0,abs_tol=1e-12) and
                     _split_criterion_matches(row,first.get('full_criterion')))
        except (KeyError,TypeError,ValueError,OverflowError):
            matches=False
        if not matches:
            reasons.append(f'{group}_first_split_mismatch')

    eta_rows=[]
    try:
        eta_rows=[row for row in event_rows if int(row['direction'])==7]
    except (KeyError,ValueError,TypeError):
        pass
    if group=='control':
        if eta_rows:
            reasons.append('control_eta_split')
    else:
        first_eta=group_report.get('first_eta_split')
        if not isinstance(first_eta,dict) or bool(first_eta.get('observed'))!=bool(eta_rows):
            reasons.append('treatment_first_eta_presence')
        elif eta_rows:
            row=eta_rows[0]
            try:
                matches=(math.isclose(float(first_eta['interval_start_hour']),
                                      float(row['interval_start_hour']),
                                      rel_tol=0.0,abs_tol=1e-12) and
                         math.isclose(float(first_eta['hour']),
                                      float(row['interval_end_hour']),
                                      rel_tol=0.0,abs_tol=1e-12) and
                         _split_criterion_matches(row,first_eta.get('full_criterion')))
            except (KeyError,TypeError,ValueError,OverflowError):
                matches=False
            if not matches:
                reasons.append('treatment_first_eta_mismatch')
    return reasons


def audit_result(directory: str|Path,orbit: str,args) -> tuple[bool,list[str]]:
    directory=Path(directory)
    reasons=[]
    hourly=[]
    paired=[]
    report=None
    expected_hours=int(round(args.days*24.0/args.save_hours))+1
    required={
        'hourly':directory/f'{orbit}_hourly_mc_stats.csv',
        'paired':directory/f'{orbit}_paired_effect_stats.csv',
        'events':directory/f'{orbit}_events.json',
        'summary':directory/f'{orbit}_summary.json',
        'visual':directory/f'{orbit}_visualization_points.h5',
        'keyframe':directory/f'{orbit}_keyframe_points.h5',
        'treatment_events':directory/f'{orbit}_treatment_split_events.csv',
        'control_events':directory/f'{orbit}_control_split_events.csv',
        'log':directory/f'{orbit}_run.log',
    }
    for key,path in required.items():
        if not path.is_file():
            reasons.append(f'missing_{key}')
    if reasons:
        return False,reasons
    try:
        hourly=_csv_rows(required['hourly'])
        paired=_csv_rows(required['paired'])
        if len(hourly)!=expected_hours*2*3:
            reasons.append('hourly_row_count')
        if len(paired)!=expected_hours*3:
            reasons.append('paired_row_count')
        if {float(row['hour']) for row in hourly}!=set(
                float(i)*args.save_hours for i in range(expected_hours)):
            reasons.append('hourly_time_grid')
        if {float(row['hour']) for row in paired}!=set(
                float(i)*args.save_hours for i in range(expected_hours)):
            reasons.append('paired_time_grid')
        if {row['group'] for row in hourly}!={'treatment','control'}:
            reasons.append('hourly_groups')
        if {int(row['sample_count']) for row in hourly}!={10000,20000,50000}:
            reasons.append('hourly_prefixes')
        if not _all_numeric_finite(hourly,('group',)):
            reasons.append('hourly_nonfinite')
        if not _all_numeric_finite(paired):
            reasons.append('paired_nonfinite')
    except (OSError,ValueError,KeyError):
        reasons.append('csv_parse')

    try:
        report=json.loads(required['events'].read_text(encoding='utf-8'))
        if report.get('output_schema_version')!='srp_ads_mc_pair_v1':
            reasons.append('report_schema')
        if report.get('pairing_complete') is not True:
            reasons.append('pairing_incomplete')
        if int(report.get('sample_count',-1))!=args.mc_count:
            reasons.append('report_sample_count')
        if report.get('orbit_name')!=orbit:
            reasons.append('report_orbit')
        for key in ('opm_sha256','config_sha256','sample_sha256','executable_sha256'):
            value=report.get(key)
            if not isinstance(value,str) or len(value)!=64 or any(
                    character not in '0123456789abcdef' for character in value.lower()):
                reasons.append(f'report_{key}')
        treatment_rows=_csv_rows(required['treatment_events'])
        control_rows=_csv_rows(required['control_events'])
        reasons.extend(audit_split_group('treatment',
            report.get('treatment_events'),treatment_rows,hourly))
        reasons.extend(audit_split_group('control',
            report.get('control_events'),control_rows,hourly))
    except (OSError,json.JSONDecodeError,TypeError,ValueError,KeyError):
        reasons.append('report_parse')

    try:
        summary=json.loads(required['summary'].read_text(encoding='utf-8'))
        if summary.get('data_complete') is not True:
            reasons.append('summary_incomplete')
        if 'converged' not in summary.get('sampling_convergence',{}):
            reasons.append('convergence_missing')
    except (OSError,json.JSONDecodeError):
        reasons.append('summary_parse')

    try:
        with h5py.File(required['visual'],'r') as h5:
            if h5.attrs.get('output_schema_version')!='srp_mc_points_hdf5_v1':
                reasons.append('visual_schema')
            if h5['treatment_error'].shape!=(expected_hours,
                    min(args.visualization_count,args.mc_count),6):
                reasons.append('visual_shape')
        with h5py.File(required['keyframe'],'r') as h5:
            if h5['treatment_error'].shape[1:]!=(args.mc_count,6):
                reasons.append('keyframe_shape')
            duration_hours=args.days*24.0
            fixed_hours={0.0,duration_hours}
            fixed_hours.update(hour for hour in (24.0,72.0,168.0,240.0,360.0)
                               if hour<=duration_hours)
            if h5['treatment_error'].shape[0]<len(fixed_hours):
                reasons.append('keyframe_fixed_hours')
    except (OSError,KeyError):
        reasons.append('hdf5_parse')

    if not (directory/'patch_snapshots').is_dir() or not any(
            (directory/'patch_snapshots').iterdir()):
        reasons.append('patch_snapshots')
    return not reasons,reasons


def write_file_manifest(directory: Path,orbit: str) -> None:
    files=[]
    for path in sorted(directory.rglob('*')):
        if path.is_file() and path.name!=f'{orbit}_file_manifest.json':
            files.append({'relative_path':str(path.relative_to(directory)),
                          'size_bytes':path.stat().st_size,
                          'sha256':file_sha256(path)})
    atomic_json(directory/f'{orbit}_file_manifest.json',{
        'output_schema_version':'srp_mc_orbit_manifest_v1',
        'orbit_name':orbit,'files':files,
    })


def run_one(orbit: str,opm: Path,executable: Path,project_root: Path,
            experiment_root: Path,samples: Path,config: Path,args,
            environment: dict[str,str]) -> dict:
    results_root=experiment_root/'results'
    final=results_root/orbit
    if final.is_dir() and not args.force:
        valid,reasons=audit_result(final,orbit,args)
        if valid:
            return {'orbit':orbit,'status':'skipped_complete','elapsed_seconds':0.0,
                    'reason':'','result_directory':str(final)}
        if not args.refresh_incompatible:
            return {'orbit':orbit,'status':'incompatible_existing',
                    'elapsed_seconds':0.0,'reason':';'.join(reasons),
                    'result_directory':str(final)}

    staging_root=results_root/'.staging'
    staging_root.mkdir(parents=True,exist_ok=True)
    stage=staging_root/f'{orbit}.{uuid.uuid4().hex}'
    stage.mkdir()
    prefix=stage/orbit
    log=Path(str(prefix)+'_run.log')
    command=[
        str(executable),'-opm',str(opm.resolve()),'-o',str(prefix.resolve()),
        '--mc-samples',str(samples.resolve()),'--mc-count',str(args.mc_count),
        '--visualization-count',str(args.visualization_count),
        '--eta-sigma',str(args.eta_sigma),'--days',str(args.days),
        '--save-hours',str(args.save_hours),'--da-order',str(args.da_order),
        '--max-depth',str(args.max_depth),'--pos-tol-km',str(args.pos_tol_km),
        '--vel-tol-kms',str(args.vel_tol_kms),'-cfg',str(config.resolve()),
    ]
    started=time.monotonic()
    with log.open('w',encoding='utf-8') as handle:
        handle.write('COMMAND: '+' '.join(command)+'\n')
        handle.flush()
        completed=subprocess.run(command,cwd=project_root,env=environment,
                                 stdout=handle,stderr=subprocess.STDOUT,check=False)
    elapsed=round(time.monotonic()-started,3)
    try:
        report_path=Path(str(prefix)+'_mc_report.json')
        if completed.returncode!=0 or not report_path.is_file():
            raise CampaignError(f'propagator failed with return code {completed.returncode}')
        inject_provenance(report_path,orbit,opm,config,samples,executable,args)
        organize_outputs(stage,prefix,orbit,args)
        valid,reasons=audit_result(stage,orbit,args)
        if not valid:
            raise CampaignError(';'.join(reasons))
        write_file_manifest(stage,orbit)
        publish_directory(stage,final,lambda path:audit_result(path,orbit,args))
        return {'orbit':orbit,'status':'succeeded','elapsed_seconds':elapsed,
                'reason':'','result_directory':str(final)}
    except Exception as exc:
        with log.open('a',encoding='utf-8') as handle:
            handle.write(f'CAMPAIGN ERROR: {exc}\n')
        return {'orbit':orbit,'status':'failed','elapsed_seconds':elapsed,
                'reason':str(exc),'result_directory':str(stage)}


def write_population_outputs(experiment_root: Path,args) -> dict:
    completed=[]
    for orbit in ORBITS:
        directory=experiment_root/'results'/orbit
        valid,_=audit_result(directory,orbit,args) if directory.is_dir() else (False,[])
        if valid:
            completed.append(orbit)
    summary_rows=[]
    metric_names=(
        'delta_position_rms_km','delta_position_p95_km','delta_position_p99_km',
        'delta_position_p999_km','delta_position_max_km',
        'delta_velocity_rms_kms','delta_velocity_p95_kms','delta_velocity_p99_kms',
        'delta_velocity_p999_kms','delta_velocity_max_kms',
        'delta_position_max_to_treatment_axis1',
        'delta_position_max_to_control_axis1',
    )
    for orbit in completed:
        directory=experiment_root/'results'/orbit
        data=json.loads((directory/f'{orbit}_summary.json').read_text())
        row={'orbit':orbit,'data_complete':data['data_complete'],
             'sampling_converged':data['sampling_convergence']['converged'],
             'first_validation_failure_hour':data['first_validation_failure_hour'],
             'last_continuously_valid_hour':data['last_continuously_valid_hour']}
        for metric in metric_names:
            if metric not in data['metrics']:
                continue
            for window in ('valid','all'):
                values=data['metrics'][metric][window]
                for key in ('endpoint_value','peak_value','peak_hour',
                            'time_mean','trapezoid_integral'):
                    row[f'{metric}.{window}.{key}']=values.get(key)
        summary_rows.append(row)
    population=experiment_root/'population'
    population.mkdir(parents=True,exist_ok=True)
    if summary_rows:
        fields=list(summary_rows[0])
        for row in summary_rows[1:]:
            for key in row:
                if key not in fields:
                    fields.append(key)
        atomic_csv(population/'orbit_summary.csv',summary_rows,fields)

    hourly_data={metric:{} for metric in metric_names}
    for orbit in completed:
        path=experiment_root/'results'/orbit/f'{orbit}_paired_effect_stats.csv'
        for row in _csv_rows(path):
            if int(row['sample_count'])!=50000:
                continue
            hour=float(row['hour'])
            for metric in metric_names:
                if metric in row:
                    hourly_data[metric].setdefault(hour,[]).append(float(row[metric]))
    quantile_rows=[]
    for metric,hours in hourly_data.items():
        for hour,values in sorted(hours.items()):
            array=np.asarray(values,dtype=np.float64)
            quantile_rows.append({
                'hour':hour,'metric':metric,'orbit_count':len(values),
                'minimum':float(np.min(array)),
                'q25':float(np.quantile(array,0.25)),
                'median':float(np.quantile(array,0.5)),
                'q75':float(np.quantile(array,0.75)),
                'maximum':float(np.max(array)),
            })
    if quantile_rows:
        atomic_csv(population/'population_hourly_quantiles.csv',quantile_rows,
                   list(quantile_rows[0]))
    result={'output_schema_version':'srp_mc_population_v1',
            'completed_orbits':completed,'completed_count':len(completed),
            'expected_count':len(ORBITS),
            'complete':len(completed)==len(ORBITS)}
    atomic_json(population/'population_summary.json',result)
    return result


def make_parser() -> argparse.ArgumentParser:
    project_root=Path(__file__).resolve().parents[1]
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-root',type=Path,default=project_root)
    parser.add_argument('--opm-root',type=Path,default=project_root/'OPM')
    parser.add_argument('--experiment-root',type=Path,
                        default=project_root/'SRP'/'20261005SRP')
    parser.add_argument('--samples',type=Path)
    parser.add_argument('--config',type=Path,default=project_root/'config'/'config.txt')
    parser.add_argument('--executable',type=Path)
    parser.add_argument('--jobs',type=int,default=2)
    parser.add_argument('--mc-count',type=int,default=50000)
    parser.add_argument('--visualization-count',type=int,default=2048)
    parser.add_argument('--eta-sigma',type=float,default=0.3)
    parser.add_argument('--days',type=float,default=15.0)
    parser.add_argument('--save-hours',type=float,default=1.0)
    parser.add_argument('--da-order',type=int,default=4)
    parser.add_argument('--max-depth',type=int,default=8)
    parser.add_argument('--pos-tol-km',type=float,default=0.1)
    parser.add_argument('--vel-tol-kms',type=float,default=1e-6)
    parser.add_argument('--expected-position-rss-km',type=float,default=10.0)
    parser.add_argument('--expected-velocity-rss-mps',type=float,default=0.03)
    parser.add_argument('--error-relative-tolerance',type=float,default=1e-9)
    parser.add_argument('--no-build',action='store_true')
    parser.add_argument('--dry-run',action='store_true')
    parser.add_argument('--audit-only',action='store_true')
    parser.add_argument('--refresh-incompatible',action='store_true')
    parser.add_argument('--force',action='store_true')
    return parser


def main(argv=None) -> int:
    args=make_parser().parse_args(argv)
    args.project_root=args.project_root.resolve()
    args.opm_root=args.opm_root.resolve()
    args.experiment_root=args.experiment_root.resolve()
    args.config=args.config.resolve()
    if args.jobs<1 or args.mc_count<2 or args.visualization_count<1:
        raise SystemExit('jobs/counts must be positive and MC count at least two')
    if args.days<=0.0 or args.save_hours<=0.0 or args.da_order<1 or args.max_depth<0:
        raise SystemExit('invalid propagation duration, cadence, order, or depth')
    if not math.isclose((args.days*24.0)/args.save_hours,
                        round((args.days*24.0)/args.save_hours),abs_tol=1e-12):
        raise SystemExit('duration must contain an integral number of save intervals')
    if not args.config.is_file():
        raise SystemExit(f'configuration file does not exist: {args.config}')
    try:
        opms=discover_opms(args.opm_root)
    except CampaignError as exc:
        print(f'ERROR: {exc}',file=sys.stderr)
        return 2

    for name in ('config','samples','results','population','logs','audit','figures'):
        (args.experiment_root/name).mkdir(parents=True,exist_ok=True)
    audits=[]
    for orbit,opm in opms.items():
        try:
            row=audit_opm(opm,args.expected_position_rss_km,
                          args.expected_velocity_rss_mps,
                          args.error_relative_tolerance)
        except CampaignError as exc:
            row={'opm':str(opm),'opm_sha256':'','position_rss_km':'',
                 'velocity_rss_mps':'','positive_definite':False,
                 'status':'invalid_opm','reason':str(exc)}
        row['orbit']=orbit
        audits.append(row)
    atomic_csv(args.experiment_root/'audit'/'input_audit.csv',audits,
               ['orbit','opm','opm_sha256','position_rss_km','velocity_rss_mps',
                'positive_definite','status','reason'])
    rejected=[row for row in audits if row['status']!='eligible']
    if rejected:
        print(f'ERROR: {len(rejected)} whitelist OPM inputs failed preflight',
              file=sys.stderr)
        return 2

    config_data={
        'output_schema_version':OUTPUT_SCHEMA,'orbits':list(ORBITS),
        'mc_count':args.mc_count,'visualization_count':args.visualization_count,
        'sample_seed':20261005,'eta_sigma':args.eta_sigma,'days':args.days,
        'save_hours':args.save_hours,'da_order':args.da_order,
        'max_split_depth':args.max_depth,'position_tolerance_km':args.pos_tol_km,
        'velocity_tolerance_kms':args.vel_tol_kms,
        'expected_position_rss_km':args.expected_position_rss_km,
        'expected_velocity_rss_mps':args.expected_velocity_rss_mps,
        'config_file':str(args.config),'config_sha256':file_sha256(args.config),
    }
    atomic_json(args.experiment_root/'config'/'experiment_config.json',config_data)
    atomic_csv(args.experiment_root/'config'/'orbit_manifest.csv',audits,
               ['orbit','opm','opm_sha256','position_rss_km','velocity_rss_mps',
                'positive_definite','status','reason'])
    if args.dry_run:
        print(json.dumps({'status':'dry_run','eligible_orbits':len(audits)}))
        return 0

    canonical_samples=args.samples is None
    if canonical_samples:
        samples=args.experiment_root/'samples'/f'mc_standard_normal_{args.mc_count}.csv'
    else:
        samples=args.samples.resolve()
    sample_manifest=args.experiment_root/'samples'/'sample_manifest.json'
    sample_npz=args.experiment_root/'samples'/f'mc_standard_normal_{args.mc_count}.npz'
    sample_gzip=args.experiment_root/'samples'/f'mc_standard_normal_{args.mc_count}.csv.gz'
    bundle_missing=(not samples.is_file() or
                    (canonical_samples and
                     (not sample_npz.is_file() or not sample_gzip.is_file())))
    if bundle_missing:
        command=[sys.executable,str(args.project_root/'tools'/'generate_srp_mc_samples.py'),
                 '--csv',str(samples),'--manifest',str(sample_manifest),
                 '--count',str(args.mc_count),'--seed','20261005',
                 '--eta-sigma',str(args.eta_sigma)]
        if canonical_samples:
            command.extend(['--npz',str(sample_npz),'--csv-gz',str(sample_gzip)])
        generated=subprocess.run(command,cwd=args.project_root,check=False)
        if generated.returncode!=0:
            return 2
    if not sample_manifest.is_file():
        atomic_json(sample_manifest,{
            'format_version':'external-srp-mc-samples-v1',
            'sample_count':args.mc_count,'seed':20261005,
            'csv_file':str(samples),'csv_sha256':file_sha256(samples),
        })

    if args.audit_only:
        rows=[]
        for orbit in ORBITS:
            directory=args.experiment_root/'results'/orbit
            valid,reasons=audit_result(directory,orbit,args) if directory.is_dir() else (
                False,['missing_result_directory'])
            rows.append({'orbit':orbit,'complete':valid,'reasons':';'.join(reasons)})
        atomic_csv(args.experiment_root/'audit'/'output_completeness.csv',rows,
                   ['orbit','complete','reasons'])
        population=write_population_outputs(args.experiment_root,args)
        return 0 if population['complete'] else 1

    if args.executable:
        executable=args.executable.resolve()
    else:
        if not args.no_build:
            build_log=args.experiment_root/'logs'/'build.log'
            with build_log.open('w',encoding='utf-8') as handle:
                built=subprocess.run(['fpm','build'],cwd=args.project_root,
                                     stdout=handle,stderr=subprocess.STDOUT,check=False)
            if built.returncode!=0:
                print(f'ERROR: build failed; see {build_log}',file=sys.stderr)
                return 2
        try:
            executable=find_executable(args.project_root)
        except CampaignError as exc:
            print(f'ERROR: {exc}',file=sys.stderr)
            return 2
    executable=stage_executable(executable,args.experiment_root)
    environment=os.environ.copy()
    environment['OMP_NUM_THREADS']='1'
    environment['OPENBLAS_NUM_THREADS']='1'
    started=time.monotonic()
    results=[]
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures={
            pool.submit(run_one,orbit,opms[orbit],executable,args.project_root,
                        args.experiment_root,samples,args.config,args,environment):orbit
            for orbit in ORBITS
        }
        for index,future in enumerate(as_completed(futures),1):
            orbit=futures[future]
            try:
                result=future.result()
            except Exception as exc:
                result={'orbit':orbit,'status':'worker_exception',
                        'elapsed_seconds':0.0,'reason':str(exc),
                        'result_directory':''}
            results.append(result)
            print(f'[{index}/{len(ORBITS)}] {orbit}: {result["status"]}',flush=True)
    results.sort(key=lambda row:ORBITS.index(row['orbit']))
    atomic_csv(args.experiment_root/'audit'/'campaign_results.csv',results,
               ['orbit','status','elapsed_seconds','reason','result_directory'])
    population=write_population_outputs(args.experiment_root,args)
    failed=[row for row in results if row['status'] not in
            ('succeeded','skipped_complete')]
    batch={
        'output_schema_version':OUTPUT_SCHEMA,
        'status':'complete' if not failed and population['complete']
                 else 'completed_with_failures',
        'orbits_requested':len(ORBITS),'failed_count':len(failed),
        'completed_count':population['completed_count'],
        'elapsed_seconds':round(time.monotonic()-started,3),
        'executable':str(executable),'executable_sha256':file_sha256(executable),
        'sample_file':str(samples),'sample_sha256':file_sha256(samples),
        'campaign_results_csv':str(
            args.experiment_root/'audit'/'campaign_results.csv'),
    }
    atomic_json(args.experiment_root/'population'/'batch_summary.json',batch)
    return 1 if failed or not population['complete'] else 0


if __name__=='__main__':
    raise SystemExit(main())
