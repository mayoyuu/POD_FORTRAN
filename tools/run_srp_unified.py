#!/usr/bin/env python3
"""One portable campaign interface; physics remains in the shared Fortran driver."""
from __future__ import annotations
import argparse
import concurrent.futures
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shlex
import subprocess
import sys
import time
import numpy as np
from srp_unified_math import shared_design, joint_sigma_points, henze_zirkler

ROOT = Path(__file__).resolve().parents[1]
CASES = {
    'C0': dict(model='off', attitude='sun', error='none'),
    'C1': dict(model='cannonball', attitude='sun', error='none'),
    'C2': dict(model='box_wing', attitude='moon', error='none'),
    'C3': dict(model='box_wing', attitude='earth', error='none'),
    'T1': dict(model='cannonball', attitude='sun', error='eta'),
    'T2': dict(model='time_varying_modes', attitude='earth', error='attitude'),
}
STATE_KEYS = ('X', 'Y', 'Z', 'X_DOT', 'Y_DOT', 'Z_DOT')
COV_KEYS = (
    ('CX_X',), ('CY_X', 'CY_Y'), ('CZ_X', 'CZ_Y', 'CZ_Z'),
    ('CX_DOT_X', 'CX_DOT_Y', 'CX_DOT_Z', 'CX_DOT_X_DOT'),
    ('CY_DOT_X', 'CY_DOT_Y', 'CY_DOT_Z', 'CY_DOT_X_DOT', 'CY_DOT_Y_DOT'),
    ('CZ_DOT_X', 'CZ_DOT_Y', 'CZ_DOT_Z', 'CZ_DOT_X_DOT', 'CZ_DOT_Y_DOT', 'CZ_DOT_Z_DOT'))
POINT_HEADER = 'sample_id,x,y,z,vx,vy,vz,eta,s,u,wm,wc'


def resolve_case(name):
    if name not in CASES: raise ValueError(f'unknown case: {name}')
    return CASES[name].copy()


def parse_method(name, order):
    value = name.upper()
    if value in ('ADS', 'UT', 'MC', 'DAMC'): return value, order, False
    if value == 'DAMC-AUTO': return 'DAMC', order, True
    if value.startswith('DAMC-') and value[5:].isdigit() and int(value[5:]) > 0:
        return 'DAMC', int(value[5:]), False
    raise ValueError(f'unknown method: {name}')


def read_opm(path, scale=1.):
    data = json.loads(Path(path).read_text(encoding='utf-8-sig'))
    mean = np.array([float(data[k]) for k in STATE_KEYS])
    cov = np.zeros((6, 6))
    for i, row in enumerate(COV_KEYS):
        for j, key in enumerate(row): cov[i, j] = cov[j, i] = float(data[key])*scale*scale
    if not np.all(np.isfinite(mean)) or not np.all(np.isfinite(cov)) or scale <= 0:
        raise ValueError('OPM state/covariance/scale is invalid')
    np.linalg.cholesky(cov)
    return mean, cov, data


def population_for_case(design, mean, cov, eta_sigma, case):
    p = np.zeros((len(design), 9))
    p[:, :6] = mean+design[:, :6]@np.linalg.cholesky(cov).T
    if case['error'] == 'eta': p[:, 6] = eta_sigma*design[:, 6]
    if case['error'] == 'attitude': p[:, 7:9] = design[:, 7:9]
    return p


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''): h.update(block)
    return h.hexdigest()


def array_sha(a):
    return hashlib.sha256(np.asarray(a, dtype='<f8').tobytes()).hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False)+'\n', encoding='utf-8')


def write_points(path, pop, ids, wm, wc):
    data = np.column_stack((ids, pop, wm, wc))
    np.savetxt(path, data, delimiter=',', header=POINT_HEADER, comments='', fmt=['%d']+['%.17e']*11)


def load_samples(path, n, seed, sigma_limit, eta_sigma=.3, eta_prior='normal'):
    if eta_prior == 'uniform' and sigma_limit*eta_sigma > 1+1e-12:
        raise ValueError('uniform eta prior support exceeds the physical SRP bound [-1,1]')
    if path is None:
        return shared_design(n, seed, sigma_limit, eta_prior), None, np.arange(1, n+1)
    with Path(path).open(encoding='utf-8-sig') as f:
        first = f.readline()
    rows = np.genfromtxt(path, delimiter=',', names=True, dtype=float,
                        skip_header=int(first.startswith('schema=')), encoding='utf-8-sig')
    rows = np.atleast_1d(rows)
    if len(rows) != n: raise ValueError(f'sample CSV has {len(rows)} rows; --count is {n}')
    fields = set(rows.dtype.names)
    ids = rows['sample_id'] if 'sample_id' in fields else np.arange(1, n+1)
    if not np.all(np.isfinite(ids)) or not np.all(ids == ids.astype(np.int64)) or len(np.unique(ids)) != n:
        raise ValueError('sample ids must be unique finite integers')
    design = shared_design(n, seed, sigma_limit, eta_prior)
    physical = None
    if all(f'z{i}' in fields for i in range(1, 7)):
        for i in range(6): design[:, i] = rows[f'z{i+1}']
        if 'z7' in fields: design[:, 6] = rows['z7']
        for i, key in enumerate(('s', 'u'), 7):
            if key in fields: design[:, i] = rows[key]
    elif set(('x', 'y', 'z', 'vx', 'vy', 'vz')).issubset(fields):
        physical = np.column_stack([rows[k] for k in ('x', 'y', 'z', 'vx', 'vy', 'vz')])
        for i, key in enumerate(('s', 'u'), 7):
            if key in fields: design[:, i] = rows[key]
        if 'eta' in fields: design[:, 6] = rows['eta']/eta_sigma
    else: raise ValueError('samples require z1..z6 or physical x,y,z,vx,vy,vz columns')
    if not np.all(np.isfinite(design)) or (physical is not None and not np.all(np.isfinite(physical))):
        raise ValueError('sample coordinates must be finite')
    if np.any(np.abs(design[:, :6]) > sigma_limit) and physical is None:
        raise ValueError('standardized samples exceed --domain-sigma')
    if np.any((design[:, 7:] < 0) | (design[:, 7:] > 1)):
        raise ValueError('mode samples s,u must lie in [0,1]')
    return design, physical, ids.astype(np.int64)


def configured_text(path, changes):
    values = {}
    for line in Path(path).read_text(encoding='utf-8-sig').splitlines():
        text = line.split('#', 1)[0].strip()
        if '=' in text:
            key, value = text.split('=', 1)
            values[key.strip()] = value.strip()
    values.update({k: str(v) for k, v in changes.items()})
    return '\n'.join(f'{k} = {v}' for k, v in values.items())+'\n'


def prepare_case(args, orbit, name, directory, design, physical, ids, opm):
    case = resolve_case(name)
    for arg, key in ((args.srp_model, 'model'), (args.attitude, 'attitude'), (args.srp_error, 'error')):
        if arg is not None: case[key] = arg
    if case['model'] == 'off' and case['error'] != 'none': raise ValueError('SRP off requires error=none')
    if (case['model'] == 'time_varying_modes') != (case['error'] == 'attitude'):
        raise ValueError('time_varying_modes and attitude error must be selected together')
    mean, cov, metadata = read_opm(opm, args.state_scale)
    pop = population_for_case(design, mean, cov, args.eta_sigma, case)
    if physical is not None:
        pop[:, :6] = physical
    if not np.all(np.isfinite(pop)): raise ValueError('sample population is not finite')
    q = np.linalg.solve(np.linalg.cholesky(cov), (pop[:, :6]-mean).T).T
    if np.any(np.abs(q) > args.domain_sigma+1e-9): raise ValueError('physical input exceeds OPM uncertainty domain')
    if case['error'] == 'eta' and (np.min(pop[:, 6]) < -1 or np.max(np.abs(pop[:, 6])) > args.domain_sigma*args.eta_sigma+1e-12):
        raise ValueError('eta samples exceed physical or ADS bounds')
    cfg = args.config or ROOT/'config/config_box_wing_earth.txt'
    changes = dict(use_srp='false' if case['model'] == 'off' else 'true',
                   srp_model='cannonball' if case['model'] == 'off' else case['model'],
                   srp_attitude_mode=case['attitude'], srp_scale_da_span='0.0',
                   srp_attitude_bias_span_arcsec='0.0 0.0 0.0', srp_array_angle_span_deg='0.0',
                   max_threads='1', parallel_computation='false', verbose_output='false')
    provenance = dict(case=case, orbit=orbit, epoch=metadata['EPOCH'], opm_sha256=sha(opm),
                      input_state_sha256=array_sha(pop[:, :6]), joint_population_sha256=array_sha(pop),
                      population_count=len(pop), eta_prior=args.eta_prior,
                      eta_support=args.domain_sigma*args.eta_sigma)
    if case['error'] == 'attitude':
        mode_dir = args.mode_root/orbit
        table, meta = mode_dir/'boxwing_modes.csv', mode_dir/'boxwing_modes_meta.json'
        m = json.loads(meta.read_text(encoding='utf-8-sig'))
        if m.get('csv_sha256') != sha(table) or m.get('opm_sha256') != sha(opm):
            raise ValueError('mode table hash or OPM identity mismatch')
        if m.get('gate_passed') is not True or m.get('audit_complete') is not True:
            raise ValueError('mode table has not passed its original audit')
        if args.hours > m['duration_hours']: raise ValueError('mode table does not cover requested duration')
        changes.update(srp_mode_table_file=table.resolve(), srp_mode_metadata_file=meta.resolve())
        provenance.update(mode_table_sha256=sha(table), mode_metadata_sha256=sha(meta),
                          physical_model='existing time-varying surrogate; not live Box-Wing truth')
    directory.mkdir(parents=True, exist_ok=False)
    config_path = directory/'resolved_config.txt'
    config_path.write_text(configured_text(cfg, changes), encoding='utf-8')
    full_points = directory/'initial_points.csv'
    write_points(full_points, pop, ids, np.full(len(pop), 1/len(pop)), np.full(len(pop), 1/(len(pop)-1)))
    provenance['resolved_config_sha256'] = sha(config_path)
    write_json(directory/'population.json', provenance)
    return pop, config_path, full_points, provenance


def ut_input(pop, case, args, path, mean, cov):
    active = list(range(6))+([6] if case['error'] == 'eta' else [7, 8] if case['error'] == 'attitude' else [])
    points, wm, wc = joint_sigma_points(pop[:, active], args.ut_alpha, args.ut_beta, args.ut_kappa)
    physical = np.zeros((len(points), 9))
    physical[:, active] = points
    q = np.linalg.solve(np.linalg.cholesky(cov), (physical[:, :6]-mean).T).T
    if np.max(np.abs(q)) > args.domain_sigma+1e-9:
        raise ValueError('UT sigma points exceed state domain: reduce --ut-alpha')
    if case['error'] == 'eta' and np.min(physical[:, 6]) < -1:
        raise ValueError('UT eta below -1: reduce --ut-alpha')
    if case['error'] == 'attitude' and np.any((physical[:, 7:] < 0) | (physical[:, 7:] > 1)):
        raise ValueError('UT mode sigma points exceed [0,1]: reduce --ut-alpha')
    write_points(path, physical, np.arange(1, len(points)+1), wm, wc)
    return len(points)


def build_executable(args, env):
    if args.executable:
        path = args.executable.resolve(strict=True)
        if not os.access(path, os.X_OK): raise ValueError('executable is not executable')
        return path, 0.
    begin = time.perf_counter()
    build_dir = ROOT/'build'
    build_dir.mkdir(exist_ok=True)
    log_path = build_dir/'unified-build.log'
    with log_path.open('w') as log:
        result = subprocess.run(['fpm', 'run', '--target', 'run_srp_unified', '--profile', 'release',
                                 '--flag', '-ffree-line-length-none', '--', '--help'],
                                cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError('Fortran build failed: '+log_path.read_text()[-4000:])
    binaries = list((ROOT/'build').glob('gfortran_*/app/run_srp_unified'))
    if not binaries: raise RuntimeError('fpm did not produce run_srp_unified')
    return max(binaries, key=lambda p: p.stat().st_mtime).resolve(), time.perf_counter()-begin


def environment():
    env = os.environ.copy()
    dirs = [ROOT/'lib']
    for key in ('DACE_ROOT', 'OPENBLAS_ROOT'):
        if env.get(key): dirs.append(Path(env[key])/'lib')
    if env.get('SPICE_LIB_DIR'): dirs.append(Path(env['SPICE_LIB_DIR']))
    for key in ('LIBRARY_PATH', 'LD_LIBRARY_PATH'):
        env[key] = ':'.join(map(str, dirs))+':'+env.get(key, '')
    return env


def validate_and_diagnose(directory, method, args, expected_points):
    begin = time.perf_counter()
    moments = np.atleast_1d(np.genfromtxt(directory/'hourly_moments.csv', delimiter=',', names=True))
    points = np.atleast_1d(np.genfromtxt(directory/'points.csv', delimiter=',', names=True))
    grid = np.arange(math.floor(args.hours/args.save_hours+1e-10)+1)*args.save_hours
    if len(grid) == 0 or not np.isclose(grid[-1], args.hours): grid = np.append(grid, args.hours)
    if len(moments) != len(grid) or not np.allclose(moments['hour'], grid, rtol=0, atol=1e-8):
        raise ValueError('hourly output time grid is incomplete or duplicated')
    for table in (moments, points):
        if not all(np.all(np.isfinite(table[k])) for k in table.dtype.names):
            raise ValueError('output contains nonfinite numeric data')
    cov_cols = [f'cov{i}{j}' for i in range(1, 7) for j in range(1, 7)]
    covariance = np.column_stack([moments[k] for k in cov_cols]).reshape(-1, 6, 6)
    means = np.column_stack([moments[f'mean{i}'] for i in range(1, 7)])
    if not np.allclose(covariance, covariance.transpose(0, 2, 1), rtol=1e-10, atol=1e-15):
        raise ValueError('output covariance is not symmetric')
    hz, cloud, reference, ids0 = [], [], [], None
    hz_seconds = 0.
    for h in grid:
        p = points[np.isclose(points['hour'], h, rtol=0, atol=1e-8)]
        if len(p) != expected_points or len(np.unique(p['sample_id'])) != expected_points:
            raise ValueError(f'hour {h}: expected {expected_points} unique saved points, got {len(p)}')
        if ids0 is None: ids0 = p['sample_id']
        elif not np.array_equal(ids0, p['sample_id']): raise ValueError('saved sample ids changed over time')
        x = np.column_stack([p[k] for k in ('x', 'y', 'z', 'vx', 'vy', 'vz')])
        cloud.append(x)
        m = moments[np.isclose(moments['hour'], h, rtol=0, atol=1e-8)][0]
        ref = np.array([m[f'reference{i}'] for i in range(1, 7)])
        reference.append(ref)
        mark = time.perf_counter()
        if method == 'UT':
            diag = dict(status='not_applicable_weighted_sigma_points', n=len(p), dimension=6,
                        alpha=args.hz_alpha, statistic=None, p_value=None, reject_normal=None, rank=None)
        elif len(p) <= 7:
            diag = dict(status='insufficient_samples', n=len(p), dimension=6,
                        alpha=args.hz_alpha, statistic=None, p_value=None, reject_normal=None, rank=None)
        else: diag = henze_zirkler(x, args.hz_alpha)
        hz_seconds += time.perf_counter()-mark
        hz.append(dict(hour=float(h), **diag))
    load_check_seconds = time.perf_counter()-begin-hz_seconds
    io_start = time.perf_counter()
    with (directory/'hz.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, list(hz[0])); writer.writeheader(); writer.writerows(hz)
    import h5py
    with h5py.File(directory/'results.h5', 'x') as f:
        f.attrs.update(schema='srp_unified_v1', method=method, state_units='km,km/s',
                       covariance_units='rr:km^2;rv:km^2/s;vv:km^2/s^2',
                       reference='same initial OPM cannonball Real nominal',
                       hz_sample_count=expected_points, hz_is_full_population=False)
        f.create_dataset('hour', data=grid)
        f.create_dataset('mean', data=means)
        f.create_dataset('covariance', data=covariance)
        f.create_dataset('reference', data=np.array(reference))
        f.create_dataset('sample_id', data=ids0.astype(np.int64))
        f.create_dataset('samples', data=np.array(cloud), compression='gzip', shuffle=True)
        for key in ('accuracy_valid', 'precision_confirmed', 'order', 'patches', 'depth_limited'):
            f.create_dataset(key, data=moments[key])
        if method == 'UT':
            f.create_dataset('mean_weights', data=p['wm'])
            f.create_dataset('covariance_weights', data=p['wc'])
    post_io_seconds = time.perf_counter()-io_start
    valid = (moments['accuracy_valid'] == 1) & (moments['precision_confirmed'] == 1)
    valid = np.logical_and.accumulate(valid)
    timing = next(csv.DictReader((directory/'timing.csv').open()))
    timing = {k: float(v) for k, v in timing.items()}
    timing.update(hz_seconds=hz_seconds, postprocess_seconds=load_check_seconds,
                  postprocess_io_seconds=post_io_seconds,
                  compute_with_diagnostics_excluding_io_seconds=timing['compute_excluding_io_seconds']+hz_seconds)
    return dict(hour_count=len(grid), saved_sample_count=expected_points,
                last_continuously_valid_hour=float(grid[valid][-1]) if np.any(valid) else None,
                all_hours_numerically_valid=bool(np.all(valid)), timing=timing)


def run_one(job, args, exe, env):
    orbit, case_name, label, method, order, adaptive, pop, config_path, full_points, provenance, opm, base = job
    directory = base/label
    directory.mkdir(exist_ok=False)
    point_file = full_points
    expected = min(2000, args.count)
    if method == 'UT':
        point_file = directory/'initial_sigma_points.csv'
        mean, cov, _ = read_opm(opm, args.state_scale)
        expected = ut_input(pop, provenance['case'], args, point_file, mean, cov)
    cmd = [str(exe), '--method', method, '--opm', str(opm), '--config', str(config_path),
           '--points', str(point_file), '--outdir', str(directory), '--error', provenance['case']['error']]
    for key, value in (('hours', args.hours), ('save-hours', args.save_hours), ('eta-sigma', args.eta_sigma),
                       ('state-scale', args.state_scale), ('domain-sigma', args.domain_sigma), ('order', order),
                       ('max-order', max(order, args.max_order)), ('adaptive', int(adaptive)), ('probes', args.probes),
                       ('pos-tol-km', args.pos_tol_km), ('vel-tol-kms', args.vel_tol_kms),
                       ('kernel-dir', str(args.kernel_dir.resolve())+'/')):
        cmd += ['--'+key, str(value)]
    (directory/'command.sh').write_text('#!/usr/bin/env bash\nset -euo pipefail\n'+shlex.join(cmd)+'\n', encoding='utf-8')
    started = time.perf_counter()
    result = {**provenance, 'orbit': orbit, 'case': case_name,
              'case_parameters': provenance['case'], 'method': label}
    with (directory/'run.log').open('w') as log:
        proc = subprocess.run(cmd, cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
    if proc.returncode:
        result.update(status='process_failed', returncode=proc.returncode, elapsed_seconds=time.perf_counter()-started)
    else:
        try:
            result.update(validate_and_diagnose(directory, method, args, expected))
            result.update(status='completed' if result['all_hours_numerically_valid'] else 'numerically_invalid',
                          elapsed_seconds=time.perf_counter()-started)
        except Exception as exc:
            result.update(status='output_validation_failed', reason=str(exc), elapsed_seconds=time.perf_counter()-started)
    write_json(directory/'report.json', result)
    print(f'{orbit} {case_name} {label}: {result["status"]}; {result["elapsed_seconds"]:.2f}s', flush=True)
    return result


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--methods', '--method', nargs='+', default=['UT', 'ADS'])
    p.add_argument('--cases', nargs='+', choices=list(CASES), default=list(CASES))
    p.add_argument('--orbits', nargs='+', default=['L1Halo-1'])
    p.add_argument('--opm', type=Path)
    p.add_argument('--opm-root', type=Path, default=ROOT/'OPM')
    p.add_argument('--mode-root', type=Path, default=ROOT/'modes')
    p.add_argument('--config', type=Path)
    p.add_argument('--samples', type=Path)
    p.add_argument('--output', type=Path)
    p.add_argument('--build-only', action='store_true', help='build once, print executable path, do not propagate')
    p.add_argument('--kernel-dir', type=Path, default=Path(os.environ.get('POD_KERNEL_DIR', ROOT/'kernels')))
    p.add_argument('--executable', type=Path)
    p.add_argument('--srp-model', choices=('off', 'cannonball', 'box_wing', 'time_varying_modes'))
    p.add_argument('--attitude', choices=('sun', 'moon', 'earth'))
    p.add_argument('--srp-error', choices=('none', 'eta', 'attitude'))
    p.add_argument('--eta-prior', choices=('normal', 'uniform'), default='normal',
                   help="T1 eta prior: truncated normal (legacy) or uniform on the ADS domain")
    for key, default in (('hours', 24.), ('save-hours', 1.), ('eta-sigma', .3), ('state-scale', 1.),
                         ('domain-sigma', 3.), ('pos-tol-km', .1), ('vel-tol-kms', 1e-6),
                         ('ut-alpha', .5), ('ut-beta', 2.), ('ut-kappa', 0.), ('hz-alpha', .05)):
        p.add_argument('--'+key, type=float, default=default)
    for key, default in (('count', 50000), ('seed', 20261005), ('jobs', 1), ('probes', 64), ('order', 4), ('max-order', 6)):
        p.add_argument('--'+key, type=int, default=default)
    p.add_argument('--prepare-only', action='store_true', help='validate input and save resolved files without propagation')
    return p


def main(argv=None):
    args = parser().parse_args(argv)
    if args.build_only:
        exe, _ = build_executable(args, environment())
        print(exe, flush=True)
        return 0
    if args.output is None: raise ValueError('--output is required unless --build-only is selected')
    for key in ('hours', 'save_hours', 'eta_sigma', 'state_scale', 'domain_sigma', 'pos_tol_km', 'vel_tol_kms', 'ut_alpha'):
        if not math.isfinite(getattr(args, key)) or getattr(args, key) <= 0: raise ValueError(f'{key} must be positive and finite')
    if args.count < 2 or args.jobs < 1 or args.probes < 1 or args.order < 1 or args.max_order < args.order:
        raise ValueError('invalid count/jobs/probes/order limits')
    if any((args.srp_model, args.attitude, args.srp_error)) and len(args.cases) != 1:
        raise ValueError('physical overrides require exactly one --cases label for unambiguous metadata')
    args.output = args.output.resolve()
    args.opm_root, args.mode_root = args.opm_root.resolve(), args.mode_root.resolve()
    if args.output.exists(): raise FileExistsError('output directory already exists; use a new path')
    if not args.kernel_dir.is_dir() and not args.prepare_only: raise ValueError('missing --kernel-dir / POD_KERNEL_DIR')
    selected = [(label, *parse_method(label, args.order)) for label in args.methods]
    if len(set(args.methods)) != len(args.methods): raise ValueError('duplicate propagation methods')
    if len(set(args.cases)) != len(args.cases) or len(set(args.orbits)) != len(args.orbits):
        raise ValueError('duplicate cases/orbits')
    if args.opm and len(args.orbits) != 1: raise ValueError('--opm supports one orbit at a time')
    args.output.mkdir(parents=True, exist_ok=False)
    design, physical, ids = load_samples(args.samples, args.count, args.seed, args.domain_sigma,
                                         args.eta_sigma, args.eta_prior)
    np.savez_compressed(args.output/'input_design.npz', design=design, sample_id=ids)
    env = environment()
    jobs = []
    for orbit in args.orbits:
        opm = args.opm.resolve() if args.opm else args.opm_root/orbit/(orbit+'_init.opm.json')
        if not opm.is_file(): raise FileNotFoundError(opm)
        for case_name in args.cases:
            base = args.output/orbit/case_name
            pop, cfg, points, provenance = prepare_case(args, orbit, case_name, base, design, physical, ids, opm)
            for label, method, order, adaptive in selected:
                jobs.append((orbit, case_name, label, method, order, adaptive, pop, cfg, points, provenance, opm, base))
    invocation = [sys.executable, str(Path(__file__).resolve()), *(argv if argv is not None else sys.argv[1:])]
    (args.output/'commands.sh').write_text('#!/usr/bin/env bash\nset -euo pipefail\n'+shlex.join(invocation)+'\n', encoding='utf-8')
    if args.prepare_only:
        print(f'Prepared {len(jobs)} jobs in {args.output}', flush=True)
        return 0
    exe, build_seconds = build_executable(args, env)
    campaign = dict(schema='srp_unified_campaign_v1', source_base='boxwing-modes-p1 3f41851 plus portable driver',
                    executable_sha256=sha(exe), build_seconds=build_seconds, input_design_sha256=array_sha(design),
                    seed=args.seed, count=args.count, eta_prior=args.eta_prior,
                    eta_support=args.domain_sigma*args.eta_sigma, command=invocation, jobs=[])
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(run_one, job, args, exe, env) for job in jobs]
        for future in concurrent.futures.as_completed(futures):
            try: campaign['jobs'].append(future.result())
            except Exception as exc: campaign['jobs'].append(dict(status='driver_failed', reason=str(exc)))
            write_json(args.output/'campaign.json', campaign)
    fields = ('orbit', 'case', 'method', 'status', 'last_continuously_valid_hour', 'all_hours_numerically_valid', 'elapsed_seconds')
    with (args.output/'campaign_summary.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fields, extrasaction='ignore'); writer.writeheader(); writer.writerows(campaign['jobs'])
    return int(any(j['status'] != 'completed' for j in campaign['jobs']))


if __name__ == '__main__':
    try: sys.exit(main())
    except Exception as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        sys.exit(2)
