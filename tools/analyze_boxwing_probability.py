#!/usr/bin/env python3
"""Portable analysis of the original fifteen saved SRP ADS probability cases."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import numpy as np

from srp_probability_analysis_data import (COVARIANCE_UNITS, ERROR_FIELDS, LIMITS, UNITS,
                                           discover_cases, load_case, load_cloud)


def write_csv(path, rows):
    with path.open('w',newline='',encoding='utf-8-sig') as stream:
        writer = csv.DictWriter(stream,fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def save_figure(fig,path):
    fig.tight_layout(); fig.savefig(path,dpi=140); plt.close(fig)


def plot_case(name, case, output):
    hours = case['hours']; mask = case['recorded_prefix']; selected = hours[mask]
    label = f'{name}: recorded numerical-check prefix ends at {case["last_passing_hour"]} h'
    fig,axes = plt.subplots(2,1,figsize=(11,7),sharex=True)
    for axis,kind,unit in zip(axes,('position','velocity'),('km','km/s')):
        for source in ('original_probe','probability_probe','zero_mode'):
            key = source+'_'+kind; field = ERROR_FIELDS[key]
            values = [r[field] for r in case['diagnostics']]
            line, = axis.plot(hours,np.maximum(values,1e-18),label=source)
            axis.axhline(LIMITS[key],color=line.get_color(),linestyle='--',alpha=.75,label=source+' tolerance')
        if case['first_invalid_hour'] is not None: axis.axvline(case['first_invalid_hour'],color='black',linestyle=':')
        axis.set(yscale='log',ylabel=f'ADS numerical error ({unit})'); axis.legend(fontsize=8,ncol=2); axis.grid(alpha=.2)
    axes[0].set_title(name+': same-surrogate ADS vs Real numerical checks; all requested hours')
    axes[-1].set_xlabel('Elapsed time (h)'); save_figure(fig,output/'numerical_errors.png')
    fig,axes = plt.subplots(2,2,figsize=(12,7),sharex=True)
    for block,(kind,unit) in enumerate((('position','km'),('velocity','km/s'))):
        for i in range(3):
            axes[block,0].plot(selected,case['mean_error'][mask,block*3+i],label=('x','y','z')[i])
            axes[block,1].plot(selected,case['std'][mask,block*3+i],label=('x','y','z')[i])
        axes[block,0].set_ylabel(f'{kind} mean offset ({unit})'); axes[block,1].set_ylabel(f'{kind} std ({unit})')
    for axis in axes.flat: axis.legend(); axis.grid(alpha=.2); axis.set_xlabel('Elapsed time (h)')
    fig.suptitle(label+'\nOffsets from same-run cannonball Real nominal; 50,000 ADS evaluations')
    save_figure(fig,output/'mean_std.png')
    fig,axes = plt.subplots(2,1,figsize=(10,7),sharex=True)
    for axis,kind,unit in zip(axes,('position','velocity'),('km','km/s')):
        for q in ('p50','p95','p99','p999','max'):
            axis.plot(selected,[r[f'{kind}_{q}'] for r,m in zip(case['rows'],mask) if m],label=q)
        axis.set_ylabel(f'{kind} offset norm ({unit})'); axis.legend(ncol=5); axis.grid(alpha=.2)
    axes[0].set_title(label+'\nEmpirical ADS-map norm quantiles; p999 is not a rare-event guarantee')
    axes[-1].set_xlabel('Elapsed time (h)'); save_figure(fig,output/'norm_quantiles.png')
    fig,axes = plt.subplots(2,3,figsize=(13,8))
    if selected.size:
        hour = float(selected[-1]); ids,values = load_cloud(case['folder'],hour)
        index = int(np.flatnonzero(mask)[-1])
        for block,unit in enumerate(('km','km/s')):
            for axis,(a,b) in zip(axes[block],((0,1),(0,2),(1,2))):
                a += 3*block; b += 3*block
                axis.scatter(values[:,a],values[:,b],s=2,alpha=.25,label='fixed 2,000-point subset')
                cov = case['covariance'][index][np.ix_([a,b],[a,b])]
                eigen,vectors = np.linalg.eigh(cov)
                angle = np.degrees(np.arctan2(vectors[1,1],vectors[0,1]))
                center = case['mean_error'][index,[a,b]]
                for sigma in (1,3):
                    axis.add_patch(Ellipse(center,2*sigma*np.sqrt(max(0,eigen[1])),2*sigma*np.sqrt(max(0,eigen[0])),
                                           angle=angle,fill=False,color='red',linestyle='-' if sigma == 1 else '--',label=f'{sigma}-sigma covariance'))
                axis.set_xlabel(f'{("dx","dy","dz","dvx","dvy","dvz")[a]} ({unit})')
                axis.set_ylabel(f'{("dx","dy","dz","dvx","dvy","dvz")[b]} ({unit})'); axis.grid(alpha=.2)
        axes[0,0].legend(fontsize=6)
    else:
        for axis in axes.flat: axis.text(.5,.5,'No recorded passing prefix',ha='center',transform=axis.transAxes)
    fig.suptitle(label+'\nOffsets from same-run cannonball Real nominal; 50,000-sample covariance contours\nContours are not Gaussian probability guarantees')
    save_figure(fig,output/'last_prefix_cloud.png')


def analyze(root,output,max_hour):
    cases = discover_cases(root); summaries = []; provenance = []
    output.mkdir(parents=True,exist_ok=True)
    for name,folder in cases.items():
        try:
            case = load_case(folder,max_hour)
            target = output/name; target.mkdir(exist_ok=True)
            write_csv(target/'hourly_diagnostics.csv',case['diagnostics'])
            np.savez_compressed(target/'moments.npz',hours=case['hours'],mean_state=case['mean_state'],
                                mean_error=case['mean_error'],covariance=case['covariance'],std=case['std'],
                                recorded_prefix=case['recorded_prefix'],component_order=np.array(['x','y','z','vx','vy','vz']),
                                state_units=UNITS,covariance_units=COVARIANCE_UNITS,
                                reference='same_run_cannonball_real_nominal',statistics_sample_count=50000)
            summary = dict(orbit=name,source_folder=str(folder.resolve()),observed_saved_last_hour=case['observed_saved_last_hour'],
                           requested_window_endpoint=max_hour,recorded_passing_prefix_endpoint=case['last_passing_hour'],
                           first_invalid_hour=case['first_invalid_hour'],first_failure_causes=';'.join(case['first_failure_causes']),
                           final_report_state=case['final_report_state'],formal_result=case['formal_result'],warnings=';'.join(case['warnings']))
            write_csv(target/'summary.csv',[summary]); summaries.append(summary)
            plot_case(name,case,target)
            hashes = {}
            for filename in ('c2_probability_hourly_stats.csv','c2_history.csv','c2_report.json'):
                path = folder/filename
                if path.exists(): hashes[filename] = hashlib.sha256(path.read_bytes()).hexdigest()
            provenance.append(dict(orbit=name,source_folder=str(folder.resolve()),source_revision=case['source_revision'],
                                   input_sha256=hashes,cloud_source='c2_probability_points.h5' if (folder/'c2_probability_points.h5').exists()
                                   else 'c2_probability_visualization.csv'))
            print(f'{name}: passing prefix {case["last_passing_hour"]} h; {case["final_report_state"]}',flush=True)
        except (ValueError,OSError,KeyError) as error:
            raise ValueError(f'{name} ({folder}): {error}') from error
    write_csv(output/'campaign_summary.csv',summaries)
    manifest = dict(schema='srp_probability_analysis_v1',root=str(root.resolve()),max_hour=max_hour,
                    state_units=UNITS,covariance_units=COVARIANCE_UNITS,numerical_tolerances=LIMITS,
                    reference='same_run_cannonball_real_nominal',statistics_sample_count=50000,visualization_sample_count=2000,
                    interpretation='Recorded numerical-check prefixes, not new formal acceptance; ADS vs Real of the same surrogate, not live Box-Wing truth.',
                    numpy_version=np.__version__,matplotlib_version=matplotlib.__version__,cases=provenance,
                    copy_verification='See package copy manifest; this analysis does not recertify the copy.')
    (output/'provenance.json').write_text(json.dumps(manifest,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parent)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--max-hour',type=float,default=300)
    args = parser.parse_args(argv)
    try: analyze(args.root,args.output or args.root/'analysis',args.max_hour)
    except (ValueError,OSError,KeyError) as error:
        print(f'Analysis failed: {error}',file=sys.stderr); return 1
    return 0


if __name__ == '__main__': raise SystemExit(main())
