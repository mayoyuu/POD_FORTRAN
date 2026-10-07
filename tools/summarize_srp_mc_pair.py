#!/usr/bin/env python3
"""Validate and summarize one paired SRP ADS Monte Carlo history."""
from __future__ import annotations
import argparse
import csv
import json
import math
import os
import tempfile
from pathlib import Path
from typing import Iterable

PREFIX_COUNTS=(10000,20000,50000)
GROUPS=('treatment','control')
MEAN_POS=('mean_dx_km','mean_dy_km','mean_dz_km')
MEAN_VEL=('mean_dvx_kms','mean_dvy_kms','mean_dvz_kms')
SKEW=('skew_dx','skew_dy','skew_dz','skew_dvx','skew_dvy','skew_dvz')
KURT=('excess_kurt_dx','excess_kurt_dy','excess_kurt_dz',
      'excess_kurt_dvx','excess_kurt_dvy','excess_kurt_dvz')
PAIRED_METRICS=(
    'center_position_difference_km','center_velocity_difference_kms',
    'delta_position_rms_km','delta_velocity_rms_kms',
    'delta_position_p95_km','delta_velocity_p95_kms',
    'delta_position_p99_km','delta_velocity_p99_kms',
    'delta_position_p999_km','delta_velocity_p999_kms',
    'delta_position_max_km','delta_velocity_max_kms',
    'position_max_increment_km','velocity_max_increment_kms',
    'axis1_position_increment_km','axis1_velocity_increment_kms',
    'position_max_amplification','velocity_max_amplification',
    'axis1_position_amplification','axis1_velocity_amplification',
    'delta_position_max_to_treatment_axis1',
    'delta_position_max_to_control_axis1',
    'delta_velocity_max_to_treatment_axis1',
    'delta_velocity_max_to_control_axis1',
)
GROUP_METRICS=(
    'position_axis1_rms_km','velocity_axis1_rms_kms',
    'position_norm_rms_km','velocity_norm_rms_kms',
    'position_p95_km','velocity_p95_kms','position_p99_km','velocity_p99_kms',
    'position_p999_km','velocity_p999_kms','position_max_km','velocity_max_kms',
    'position_max_to_axis1','velocity_max_to_axis1','patch_count',
    'depth_limited_patches','validation_position_error_km',
    'validation_velocity_error_kms',
)


def _read_csv(path: Path) -> list[dict[str,str]]:
    with path.open(newline='') as handle:
        rows=list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f'empty CSV: {path}')
    return rows


def _number(row: dict[str,str],key: str) -> float:
    try:
        return float(row[key])
    except (KeyError,TypeError,ValueError) as exc:
        raise ValueError(f'invalid or missing numeric field {key!r}') from exc


def _ratio_difference(a: float,b: float) -> float:
    scale=abs(b)
    if scale <= 100.0*math.ulp(1.0):
        return 0.0 if abs(a-b)<=100.0*math.ulp(1.0) else math.inf
    return abs(a-b)/scale


def _trapz(hours: list[float],values: list[float]) -> float:
    return sum((values[i]+values[i+1])*0.5*(hours[i+1]-hours[i])
               for i in range(len(hours)-1))


def scalar_summary(hours: Iterable[float],values: Iterable[float],
                   structurally_valid: bool=True) -> dict:
    hours=list(hours)
    values=list(values)
    valid=(structurally_valid and bool(hours) and len(hours)==len(values)
           and all(math.isfinite(v) for v in hours+values)
           and all(hours[i+1]>hours[i] for i in range(len(hours)-1)))
    if not valid:
        return {'valid':False,'endpoint_hour':None,'endpoint_value':None,
                'peak_value':None,'peak_hour':None,'time_mean':None,
                'trapezoid_integral':None,'minimum':None,'minimum_hour':None,
                'maximum':None,'maximum_hour':None,'maximum_absolute':None,
                'maximum_absolute_hour':None}
    peak_index=max(range(len(values)),key=lambda i:(values[i],-i))
    min_index=min(range(len(values)),key=lambda i:(values[i],i))
    max_abs_index=max(range(len(values)),key=lambda i:(abs(values[i]),-i))
    integral=_trapz(hours,values)
    duration=hours[-1]-hours[0]
    mean=integral/duration if duration>0.0 else values[-1]
    return {
        'valid':True,'start_hour':hours[0],'endpoint_hour':hours[-1],
        'duration_hours':duration,'endpoint_value':values[-1],
        'peak_value':values[peak_index],'peak_hour':hours[peak_index],
        'time_mean':mean,'trapezoid_integral':integral,
        'minimum':values[min_index],'minimum_hour':hours[min_index],
        'maximum':values[peak_index],'maximum_hour':hours[peak_index],
        'maximum_absolute':values[max_abs_index],
        'maximum_absolute_hour':hours[max_abs_index],
    }


def _first_issue(rows: list[dict[str,str]],expected: list[float],
                 key_fields: tuple[str,...]) -> tuple[float|None,dict]:
    seen={}
    nonfinite=[]
    duplicate=[]
    for row in rows:
        try:
            hour=_number(row,'hour')
        except ValueError:
            return 0.0,{}
        key=tuple(row.get(name,'') for name in key_fields)+(hour,)
        if key in seen:
            duplicate.append(hour)
        seen[key]=row
        for name,value in row.items():
            if name in ('group',) or value is None or value=='':
                continue
            try:
                numeric=float(value)
            except ValueError:
                return hour,seen
            if not math.isfinite(numeric):
                nonfinite.append(hour)
    present={float(row['hour']) for row in rows if math.isfinite(float(row['hour']))}
    missing=[hour for hour in expected if hour not in present]
    candidates=missing+nonfinite+duplicate
    return (min(candidates) if candidates else None),seen


def _convergence(group_index: dict,valid_hours: list[float]) -> dict:
    thresholds={
        'mean_over_axis':0.005,'longest_axis_relative':0.01,
        'max_abs_skew_difference':0.02,'max_abs_kurtosis_difference':0.05,
        'p99_relative':0.02,
    }
    worst={name:0.0 for name in thresholds}
    missing=False
    for hour in valid_hours:
        for group in GROUPS:
            r20=group_index.get((group,'20000',hour))
            r50=group_index.get((group,'50000',hour))
            if r20 is None or r50 is None:
                missing=True
                continue
            for components,axis in ((MEAN_POS,'position_axis1_rms_km'),
                                    (MEAN_VEL,'velocity_axis1_rms_kms')):
                difference=math.sqrt(sum((_number(r20,c)-_number(r50,c))**2
                                         for c in components))
                denom=abs(_number(r50,axis))
                value=difference/denom if denom>0.0 else (0.0 if difference==0.0 else math.inf)
                worst['mean_over_axis']=max(worst['mean_over_axis'],value)
            for axis in ('position_axis1_rms_km','velocity_axis1_rms_kms'):
                worst['longest_axis_relative']=max(
                    worst['longest_axis_relative'],
                    _ratio_difference(_number(r20,axis),_number(r50,axis)))
            skew20=max(abs(_number(r20,c)) for c in SKEW)
            skew50=max(abs(_number(r50,c)) for c in SKEW)
            kurt20=max(abs(_number(r20,c)) for c in KURT)
            kurt50=max(abs(_number(r50,c)) for c in KURT)
            worst['max_abs_skew_difference']=max(
                worst['max_abs_skew_difference'],abs(skew20-skew50))
            worst['max_abs_kurtosis_difference']=max(
                worst['max_abs_kurtosis_difference'],abs(kurt20-kurt50))
            for metric in ('position_p99_km','velocity_p99_kms'):
                worst['p99_relative']=max(
                    worst['p99_relative'],
                    _ratio_difference(_number(r20,metric),_number(r50,metric)))
    checks={name:(math.isfinite(worst[name]) and worst[name]<=limit)
            for name,limit in thresholds.items()}
    return {'converged':not missing and bool(valid_hours) and all(checks.values()),
            'missing_required_rows':missing,'thresholds':thresholds,
            'worst_observed':worst,'checks':checks}


def summarize_files(hourly_path: str|Path,paired_path: str|Path,
                    report_path: str|Path,expected_hours: int=361,
                    save_hours: float|None=None) -> dict:
    hourly_path=Path(hourly_path)
    paired_path=Path(paired_path)
    report_path=Path(report_path)
    hourly_rows=_read_csv(hourly_path)
    paired_rows=_read_csv(paired_path)
    with report_path.open() as handle:
        report=json.load(handle)
    if save_hours is None:
        save_hours=float(report.get('save_hours',1.0))
    if not math.isfinite(save_hours) or save_hours<=0.0:
        raise DataError('save_hours must be finite and positive')
    expected=[float(i)*save_hours for i in range(expected_hours)]
    issue_group,group_index=_first_issue(
        hourly_rows,expected,('group','sample_count'))
    issue_pair,pair_index=_first_issue(
        paired_rows,expected,('sample_count',))
    required_group={(group,str(count),hour) for group in GROUPS
                    for count in PREFIX_COUNTS for hour in expected}
    required_pair={(str(count),hour) for count in PREFIX_COUNTS for hour in expected}
    # _first_issue uses raw string key fields.
    missing_group=required_group-set(group_index)
    missing_pair=required_pair-set(pair_index)
    issue_candidates=[x for x in (issue_group,issue_pair) if x is not None]
    if missing_group:
        issue_candidates.extend(key[-1] for key in missing_group)
    if missing_pair:
        issue_candidates.extend(key[-1] for key in missing_pair)
    first_issue=min(issue_candidates) if issue_candidates else None
    complete=(first_issue is None and report.get('output_schema_version')==
              'srp_ads_mc_pair_v1' and report.get('pairing_complete') is True)

    failure_hours=[]
    for hour in expected:
        for group in GROUPS:
            row=group_index.get((group,'50000',hour))
            if row is not None and _number(row,'validation_failed')!=0.0:
                failure_hours.append(hour)
    first_failure=min(failure_hours) if failure_hours else None
    valid_hours=[hour for hour in expected
                 if first_failure is None or hour<first_failure]
    if first_issue is not None:
        valid_hours=[hour for hour in valid_hours if hour<first_issue]

    metrics={}
    paired50={hour:pair_index.get(('50000',hour)) for hour in expected}
    for metric in PAIRED_METRICS:
        if not any(row is not None and metric in row for row in paired50.values()):
            continue
        all_values=[(_number(paired50[hour],metric)
                     if paired50[hour] is not None and metric in paired50[hour]
                     else math.nan) for hour in expected]
        valid_values=[(_number(paired50[hour],metric)
                       if paired50[hour] is not None and metric in paired50[hour]
                       else math.nan) for hour in valid_hours]
        metrics[metric]={
            'all':scalar_summary(expected,all_values,complete),
            'valid':scalar_summary(valid_hours,valid_values,
                                   complete or (first_issue is None)),
        }
    for group in GROUPS:
        rows50={hour:group_index.get((group,'50000',hour)) for hour in expected}
        for metric in GROUP_METRICS:
            if not any(row is not None and metric in row for row in rows50.values()):
                continue
            all_values=[(_number(rows50[hour],metric)
                         if rows50[hour] is not None and metric in rows50[hour]
                         else math.nan) for hour in expected]
            valid_values=[(_number(rows50[hour],metric)
                           if rows50[hour] is not None and metric in rows50[hour]
                           else math.nan) for hour in valid_hours]
            metrics[f'{group}.{metric}']={
                'all':scalar_summary(expected,all_values,complete),
                'valid':scalar_summary(valid_hours,valid_values,
                                       complete or (first_issue is None)),
            }
    return {
        'output_schema_version':'srp_ads_mc_pair_summary_v1',
        'source_report_schema':report.get('output_schema_version'),
        'data_complete':complete,'first_data_issue_hour':first_issue,
        'expected_hour_count':expected_hours,'save_hours':save_hours,
        'first_validation_failure_hour':first_failure,
        'last_continuously_valid_hour':valid_hours[-1] if valid_hours else None,
        'sampling_convergence':_convergence(group_index,valid_hours),
        'metrics':metrics,
    }


def _atomic_text(path: Path,text: str) -> None:
    path.parent.mkdir(parents=True,exist_ok=True)
    fd,temp=tempfile.mkstemp(prefix=path.name+'.',suffix='.tmp',dir=path.parent)
    try:
        with os.fdopen(fd,'w',newline='') as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temp,path)
    except BaseException:
        try: os.unlink(temp)
        except FileNotFoundError: pass
        raise


def write_flat_csv(result: dict,path: Path) -> None:
    row={'data_complete':result['data_complete'],
         'first_data_issue_hour':result['first_data_issue_hour'],
         'first_validation_failure_hour':result['first_validation_failure_hour'],
         'last_continuously_valid_hour':result['last_continuously_valid_hour'],
         'sampling_converged':result['sampling_convergence']['converged']}
    for metric,windows in result['metrics'].items():
        for window,values in windows.items():
            for key in ('endpoint_value','peak_value','peak_hour',
                        'time_mean','trapezoid_integral'):
                row[f'{metric}.{window}.{key}']=values.get(key)
    output=[]
    from io import StringIO
    buffer=StringIO()
    writer=csv.DictWriter(buffer,fieldnames=list(row))
    writer.writeheader()
    writer.writerow(row)
    _atomic_text(path,buffer.getvalue())


def main() -> int:
    parser=argparse.ArgumentParser()
    parser.add_argument('--prefix',required=True,type=Path)
    parser.add_argument('--output-json',required=True,type=Path)
    parser.add_argument('--output-csv',type=Path)
    parser.add_argument('--expected-hours',type=int,default=361)
    parser.add_argument('--save-hours',type=float)
    args=parser.parse_args()
    prefix=str(args.prefix)
    result=summarize_files(prefix+'_mc_hourly_stats.csv',
                           prefix+'_mc_paired_stats.csv',
                           prefix+'_mc_report.json',args.expected_hours,args.save_hours)
    _atomic_text(args.output_json,json.dumps(result,indent=2,allow_nan=False)+'\n')
    if args.output_csv:
        write_flat_csv(result,args.output_csv)
    print(json.dumps({'data_complete':result['data_complete'],
                      'sampling_converged':result['sampling_convergence']['converged']}))
    return 0 if result['data_complete'] else 2


if __name__=='__main__':
    raise SystemExit(main())
