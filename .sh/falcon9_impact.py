#!/usr/bin/env python3
"""Four Falcon9 OPM MC populations at one UTC epoch, using run_srp_unified."""
from __future__ import annotations
import argparse
import concurrent.futures
import copy
import csv
from datetime import datetime, timezone
from pathlib import Path
import sys
import shutil
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
import run_srp_unified as runner

IMPACT_UTC = '2026-08-05T06:35:42.450Z'
SOURCE = 'https://www.projectpluto.com/pluto/mpecs/25010d_new.htm'
ORBITS = ('FAC_2604_1', 'FAC_2605_1', 'FAC_2606_1', 'FAC_2607_1')
CASES = ('C0', 'C1', 'T1')


def utc(text):
    value = datetime.fromisoformat(text.replace('Z', '+00:00'))
    if value.tzinfo is None or value.utcoffset().total_seconds() != 0:
        raise ValueError('epoch must explicitly use UTC (Z or +00:00)')
    return value


def export_snapshot(directory, end_utc, expected):
    """Export the full last saved cloud; these are states, not surface intersections."""
    import h5py
    with h5py.File(directory/'results.h5') as f:
        hour = float(f['hour'][-1])
        end_et = float(f.attrs['end_et'])
    count = 0
    with (directory/'points.csv').open(newline='') as source, (directory/'impact_epoch_points.csv').open('w', newline='') as dest:
        reader = csv.DictReader(source)
        writer = csv.DictWriter(dest, ['epoch_utc', 'epoch_et', *reader.fieldnames])
        writer.writeheader()
        for row in reader:
            if abs(float(row['hour'])-hour) <= 1e-10:
                writer.writerow(dict(epoch_utc=end_utc, epoch_et=format(end_et, '.17g'), **row))
                count += 1
    if count != expected:
        raise ValueError(f'final snapshot has {count} points; expected {expected}')


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--count', type=int, default=50000)
    p.add_argument('--jobs', type=int, default=8)
    p.add_argument('--seed', type=int, default=20261005)
    p.add_argument('--cases', nargs='+', choices=CASES, default=list(CASES))
    p.add_argument('--final-only', action='store_true', help='store only terminal MC distributions')
    p.add_argument('--impact-utc', default=IMPACT_UTC,
                   help='common reference epoch, not a per-particle impact event')
    p.add_argument('--output', type=Path)
    p.add_argument('--executable', type=Path)
    p.add_argument('--opm-root', type=Path, default=ROOT/'OPM/Falcon9')
    p.add_argument('--kernel-dir', type=Path,
                   default=Path(runner.os.environ.get('POD_KERNEL_DIR', ROOT/'kernels')))
    p.add_argument('--prepare-only', action='store_true')
    a = p.parse_args(argv)
    if a.count < 2 or a.jobs < 1:
        raise ValueError('count>=2 and jobs>=1 are required')
    if len(set(a.cases)) != len(a.cases): raise ValueError('duplicate cases')
    end = utc(a.impact_utc)
    a.impact_utc = end.isoformat(timespec='microseconds').replace('+00:00', 'Z')
    a.output = (a.output or a.opm_root/('mc_impact_'+datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S_%f'))).resolve()
    if a.output.exists():
        raise FileExistsError('output directory already exists; use a new path')
    a.opm_root = a.opm_root.resolve()
    a.kernel_dir = a.kernel_dir.resolve()
    if not a.prepare_only and not a.kernel_dir.is_dir():
        raise ValueError('missing --kernel-dir / POD_KERNEL_DIR')
    config = ROOT/'config/config_falcon9_cannonball.txt'
    # Six-sigma truncation effectively retains the OPM Gaussian covariance.
    # eta uses an independent uniform [-1,1] draw.
    args = runner.parser().parse_args(['--methods', 'MC', '--cases', *a.cases, '--count', str(a.count),
        '--seed', str(a.seed), '--eta-prior', 'uniform', '--eta-sigma', str(1/6),
        '--domain-sigma', '6', '--state-scale', '1', '--saved-points', '0',
        '--save-hours', '1000000000', '--end-utc', a.impact_utc, '--config', str(config),
        '--kernel-dir', str(a.kernel_dir)])
    args.executable = a.executable
    args.final_only = a.final_only
    inputs = []
    for name in ORBITS:
        opm = a.opm_root/(name+'.opm')
        mean, cov, metadata = runner.read_opm(opm)
        if metadata.get('REF_FRAME') not in ('GCRF', 'J2000', 'EME2000'):
            raise ValueError(f'{name}: expected geocentric inertial OPM')
        duration = (end-utc(metadata['EPOCH'])).total_seconds()/3600
        if duration <= 0:
            raise ValueError(f'{name}: impact epoch is not after OPM epoch')
        inputs.append((name, opm, metadata, duration))
    a.output.mkdir(parents=True)
    print('Output: '+str(a.output), flush=True)
    design, physical, ids = runner.load_samples(None, a.count, a.seed, 6., 1/6, 'uniform')
    np.savez_compressed(a.output/'input_design.npz', design=design, sample_id=ids)
    manifest = dict(schema='falcon9_srp_impact_epoch_v1', impact_utc=a.impact_utc,
        impact_epoch_kind='published prediction; common snapshot epoch, not collision detection',
        source=SOURCE, effective_cr_area_mass_m2_kg=.0079025,
        scale_prior='T1: uniform [0,2]; C0: exactly zero; C1: exactly one',
        state_prior='top-level OPM covariance, Gaussian truncated at six whitened sigmas',
        count=a.count, seed=a.seed, cases=a.cases, final_only=a.final_only,
        driver_pid=runner.os.getpid(), worker_limit=a.jobs,
        started_utc=datetime.now(timezone.utc).isoformat(),
        command=[sys.executable, str(Path(__file__).resolve()), *(argv if argv is not None else sys.argv[1:])],
        jobs=[], orbits={})
    jobs = []
    for name, opm, metadata, duration in inputs:
        job_args = copy.copy(args)
        job_args.hours = duration
        manifest['orbits'][name] = dict(epoch_utc=metadata['EPOCH'], utc_duration_hours=duration,
                                       opm_sha256=runner.sha(opm))
        print(f'{name}: {metadata["EPOCH"]} -> {a.impact_utc}; UTC duration {duration:.12f} h', flush=True)
        for case_name in a.cases:
            base = a.output/name/case_name
            pop, cfg, points, provenance = runner.prepare_case(job_args, name, case_name, base, design, physical, ids, opm)
            job = (name, case_name, 'MC', 'MC', 4, False, pop, cfg, points, provenance, opm, base)
            jobs.append((job, job_args))
    if a.prepare_only:
        manifest['jobs'] = [dict(orbit=j[0], case=j[1], method='MC', status='prepared', saved_points=0) for j, _ in jobs]
        runner.write_json(a.output/'campaign.json', manifest)
        print(f'Prepared {len(jobs)} MC jobs; no propagation.', flush=True)
        return 0
    manifest['planned_jobs'] = len(jobs)
    runner.write_json(a.output/'campaign.json', manifest)
    env = runner.environment()
    exe, build_seconds = runner.build_executable(args, env)
    # Pin this campaign's binary so a later local rebuild cannot change queued jobs.
    runtime=a.output/'.runtime'
    runtime.mkdir()
    pinned=runtime/'run_srp_unified'
    shutil.copy2(exe,pinned)
    exe=pinned
    manifest.update(executable_path=str(exe), executable_sha256=runner.sha(exe), build_seconds=build_seconds)
    runner.write_json(a.output/'campaign.json', manifest)
    with concurrent.futures.ThreadPoolExecutor(max_workers=a.jobs) as pool:
        pending = {pool.submit(runner.run_one, job, jargs, exe, env): (job, jargs) for job, jargs in jobs}
        for future in concurrent.futures.as_completed(pending):
            job, jargs = pending[future]
            try:
                report = future.result()
                if report['status'] == 'completed':
                    export_snapshot(job[-1]/'MC', a.impact_utc, a.count)
                    report['final_snapshot'] = str(job[-1]/'MC/impact_epoch_points.csv')
                    runner.write_json(job[-1]/'MC/report.json', report)
            except Exception as exc:
                report = dict(orbit=job[0], case=job[1], method='MC', status='driver_failed', reason=str(exc))
            manifest['jobs'].append(report)
            runner.write_json(a.output/'campaign.json', manifest)
    fields = ('orbit', 'case', 'method', 'status', 'saved_sample_count', 'final_snapshot', 'reason')
    with (a.output/'campaign_summary.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fields, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(manifest['jobs'])
    failed = sum(j['status'] != 'completed' for j in manifest['jobs'])
    print(f'{len(jobs)-failed}/{len(jobs)} jobs completed. Output: {a.output}', flush=True)
    return int(failed > 0)


if __name__ == '__main__':
    try:
        sys.exit(main())
    except Exception as exc:
        print('ERROR: '+str(exc), file=sys.stderr)
        sys.exit(2)
