#!/usr/bin/env python3
"""Stream and independently audit layered 8D ADS probability output."""
import csv
import hashlib
import importlib.util
import json
import math
from pathlib import Path

import h5py
import numpy as np

_spec = importlib.util.spec_from_file_location('modes_probability_samples',
                                             Path(__file__).with_name('generate_srp_modes_probability_samples.py'))
SAMPLER = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(SAMPLER)
POINT_FIELDS = ('hour,epoch_et,keyframe_kind,sample_id,patch_id,x_km,y_km,z_km,vx_kms,vy_kms,vz_kms,'
                'dx_km,dy_km,dz_km,dvx_kms,dvy_kms,dvz_kms').split(',')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def read_statistics(path):
    with Path(path).open(newline='') as stream:
        reader = csv.DictReader(stream)
        require(reader.fieldnames and len(reader.fieldnames) == len(set(reader.fieldnames)), 'bad statistics header')
        result = []
        for row in reader:
            require(None not in row and all(value is not None for value in row.values()), 'bad statistics fields')
            converted = {}
            for key, value in row.items():
                if not value.strip():
                    require(key in ('first_failure_hour', 'first_unconfirmed_hour') or
                            key.startswith('skewness') or key.startswith('excess_kurtosis'), 'missing statistic '+key)
                    converted[key] = None
                else:
                    converted[key] = float(value)
                    require(math.isfinite(converted[key]), 'nonfinite statistic '+key)
            result.append(converted)
    require(result, 'missing statistics')
    return result


def vector(row, prefix, size=6):
    return np.array([row[f'{prefix}{i}'] for i in range(1, size+1)], dtype=float)


def covariance(row):
    return np.array([[row[f'cov_{i}{j}'] for j in range(1, 7)] for i in range(1, 7)])


def audit_statistics_row(row):
    """Algebraic consistency at every hour, without inferring N moments from saved V points."""
    n = int(row['statistics_sample_count']); cov = covariance(row)
    mean = vector(row, 'mean_error'); low = vector(row, 'min'); high = vector(row, 'max')
    std = vector(row, 'std'); reference = vector(row, 'reference')
    require(np.all(np.diag(cov) >= 0) and np.all(std >= 0), 'negative variance/standard deviation')
    _near(std, np.sqrt(np.diag(cov)), 'hourly standard deviation/covariance', rtol=1e-12, atol=1e-28)
    _near(vector(row, 'mean_state'), reference+mean, 'hourly mean translation', rtol=0,
          atol=np.array([1e-9]*3+[2e-14]*3))
    scale = np.maximum(1, np.maximum(np.maximum(abs(low), abs(high)), abs(mean)))
    tolerance = scale*2e-14
    require(np.all(low <= mean+tolerance) and np.all(mean <= high+tolerance), 'hourly min/mean/max ordering')
    sigma = std*np.sqrt((n-1)/n)
    variance_valid = sigma > 100*np.finfo(float).eps*scale
    for i in range(6):
        for name, flag, count_ok in (('skewness', 'skewness_valid', n > 2),
                                     ('excess_kurtosis', 'kurtosis_valid', n > 3)):
            valid = bool(variance_valid[i] and count_ok)
            require(row[f'{flag}{i+1}'] == int(valid), 'hourly standardized moment validity mismatch')
            require((row[f'{name}{i+1}'] is not None) == valid, 'hourly undefined moment must be blank')
    require(row['patches'] >= 1 and row['patches'].is_integer() and row['depth_limited'].is_integer()
            and 0 <= row['depth_limited'] <= row['patches'], 'invalid patch/depth count')
    for flag, first in (('ads_accuracy_valid', 'first_failure_hour'), ('precision_confirmed', 'first_unconfirmed_hour')):
        require(row[flag] in (0, 1), 'invalid cumulative validity bit')
        require((row[first] is None) == bool(row[flag]), 'first failure/precision flag mismatch')
        if row[first] is not None:
            require(0 <= row[first] <= row['hour'], 'invalid first failure/precision time')
    for name, block in (('position', slice(0, 3)), ('velocity', slice(3, 6))):
        quantiles = [row[name+'_'+q] for q in ('p50', 'p95', 'p99', 'p999')]
        require(0 <= quantiles[0] and np.all(np.diff(quantiles) >= 0) and quantiles[-1] <= row[name+'_max']
                and 0 <= row[name+'_norm_rms'] <= row[name+'_max']+1e-14, 'invalid hourly norm statistics')
        require(row[name+'_worst_id'].is_integer() and 1 <= row[name+'_worst_id'] <= n, 'invalid worst sample ID')
        _near(vector(row, name+'_axis_rms', 3),
              np.sqrt(np.maximum(0, np.linalg.eigvalsh(cov[block, block])[::-1])), 'hourly covariance RMS axes', atol=1e-9)


def convergence_report(rows, count):
    if count != 50000:
        return {'status': 'not_applicable_smoke' if count < 50000 else 'not_evaluated_nonstandard_count',
                'confirmed': False, 'metrics': {}, 'p999_is_rare_event_guarantee': False}
    grouped = {(row['hour'], int(row['statistics_sample_count'])): row for row in rows}
    hours = sorted({row['hour'] for row in rows})
    require(hours and all((h, n) in grouped for h in hours for n in (10000, 20000, 50000)),
            'missing convergence prefixes')
    accum = {}

    def add(name, delta, denominator, mandatory=True):
        old = accum.setdefault(name, [0., 0., mandatory])
        old[0] = max(old[0], float(delta)); old[1] = max(old[1], float(denominator))

    for hour in hours:
        small, large = grouped[hour, 20000], grouped[hour, 50000]
        a, b = covariance(small), covariance(large)
        for name, block in (('position', slice(0, 3)), ('velocity', slice(3, 6))):
            add(name+'_mean', np.linalg.norm((vector(small, 'mean_error')-vector(large, 'mean_error'))[block]),
                np.sqrt(max(0., np.trace(b[block, block]))))
            add(name+'_covariance', np.linalg.norm(a[block, block]-b[block, block]), np.linalg.norm(b[block, block]))
            for quantile in ('p50', 'p95', 'p99', 'p999'):
                key = name+'_'+quantile
                add(key, abs(small[key]-large[key]), abs(large[key]), quantile != 'p999')
        add('position_velocity_covariance', np.linalg.norm(a[:3, 3:]-b[:3, 3:]),
            np.sqrt(np.linalg.norm(b[:3, :3])*np.linalg.norm(b[3:, 3:])))
    metrics = {}
    for name, (delta, scale, mandatory) in accum.items():
        ratio = delta/scale if scale > 0 else (0. if delta == 0 else None)
        metrics[name] = {'maximum_absolute_change': delta, 'normalization_scale': scale,
                         'relative_change': ratio, 'degenerate': scale == delta == 0,
                         'reason': 'zero_scale_nonzero_change' if ratio is None else None,
                         'passed': ratio is not None and ratio <= .05, 'mandatory': mandatory}
    confirmed = all(item['passed'] for item in metrics.values() if item['mandatory'])
    return {'status': 'confirmed' if confirmed else 'unconfirmed', 'confirmed': confirmed,
            'comparison': '20000_to_50000', 'threshold': .05, 'metrics': metrics,
            'p999_effective_tail_sample_count': 50, 'p999_is_rare_event_guarantee': False,
            'interpretation': 'nested-prefix empirical stability, not independent confidence intervals'}


def iter_frames(path):
    """Read at most one frame into memory; strict sequential identities are checked by caller."""
    with Path(path).open(newline='') as stream:
        require(stream.readline().strip() == 'schema=srp_modes_probability_points_v1', 'invalid points schema')
        reader = csv.DictReader(stream)
        require(reader.fieldnames == POINT_FIELDS, 'invalid point fields/units')
        current, records = None, []
        for row in reader:
            require(None not in row and all(value is not None and value.strip() for value in row.values()),
                    'empty/extra point fields')
            number = [float(row[name]) for name in POINT_FIELDS if name != 'keyframe_kind']
            require(all(math.isfinite(x) for x in number), 'nonfinite point')
            identity = (number[0], number[1], row['keyframe_kind'])
            if current is not None and identity != current:
                yield current, np.asarray(records)
                records = []
            current = identity; records.append(number[2:])
        if current is not None:
            yield current, np.asarray(records)


def _near(actual, expected, label, rtol=2e-9, atol=1e-13):
    require(np.allclose(actual, expected, rtol=rtol, atol=atol), label+' mismatch')


def audit_keyframe(errors, reference, rows):
    """Independent NumPy implementation; native implementation uses compensated sums."""
    for row in rows:
        n = int(row['statistics_sample_count']); values = errors[:n]
        mean, cov = values.mean(0), np.cov(values, rowvar=False, ddof=1)
        scales = np.sqrt(np.maximum(np.diag(cov), 0))
        _near(vector(row, 'mean_error'), mean, 'full-population mean', atol=1e-13)
        _near(vector(row, 'mean_state'), reference+mean, 'absolute mean', atol=1e-13)
        _near(covariance(row), cov, 'full covariance', atol=np.outer(scales, scales)*2e-12+1e-28)
        _near(vector(row, 'std'), scales, 'standard deviation')
        _near(vector(row, 'min'), values.min(0), 'component min')
        _near(vector(row, 'max'), values.max(0), 'component max')
        centered = values-mean; m2 = (centered**2).mean(0)
        # Match numerical-zero rule of compute_group_stats, not an arbitrary Gaussian fit.
        scale = np.maximum(1, np.maximum(np.abs(values).max(0), np.abs(mean)))
        valid_variance = np.sqrt(m2) > 100*np.finfo(float).eps*scale
        for i in range(6):
            sv, kv = valid_variance[i] and n > 2, valid_variance[i] and n > 3
            require(row[f'skewness_valid{i+1}'] == int(sv) and row[f'kurtosis_valid{i+1}'] == int(kv),
                    'standardized moment validity mismatch')
            for name, valid in (('skewness', sv), ('excess_kurtosis', kv)):
                if not valid:
                    require(row[f'{name}{i+1}'] is None, 'invalid moment must be empty')
            if sv:
                skew = np.sqrt(n*(n-1))/(n-2)*(centered[:, i]**3).mean()/m2[i]**1.5
                _near(row[f'skewness{i+1}'], skew, 'skewness', atol=1e-10)
            if kv:
                g2 = (centered[:, i]**4).mean()/m2[i]**2-3
                kurt = (n-1)/((n-2)*(n-3))*((n+1)*g2+6)
                _near(row[f'excess_kurtosis{i+1}'], kurt, 'kurtosis', atol=1e-10)
        for name, block in (('position', slice(0, 3)), ('velocity', slice(3, 6))):
            norms = np.linalg.norm(values[:, block], axis=1)
            _near([row[f'{name}_{q}'] for q in ('p50', 'p95', 'p99', 'p999')],
                  np.quantile(norms, [.5, .95, .99, .999], method='linear'), 'norm quantiles')
            _near(row[f'{name}_norm_rms'], np.sqrt((norms**2).mean()), 'norm RMS')
            _near(row[f'{name}_max'], norms.max(), 'norm maximum')
            require(row[f'{name}_worst_id'] == int(np.argmax(norms))+1, 'worst sample identity mismatch')
            _near(vector(row, name+'_axis_rms', 3),
                  np.sqrt(np.maximum(0, np.linalg.eigvalsh(cov[block, block])[::-1])), 'RMS axes', atol=1e-9)


def pack_probability(prefix, samples_path, manifest_path, *, metadata=None):
    prefix = Path(prefix); metadata = dict(metadata or {})
    population, manifest = SAMPLER.load_sample_bundle(samples_path, manifest_path)
    n = len(population)
    stats_path = Path(str(prefix)+'_probability_hourly_stats.csv')
    visual_path = Path(str(prefix)+'_probability_visualization.csv')
    key_path = Path(str(prefix)+'_probability_keyframes.csv')
    stats = read_statistics(stats_path)
    hours = sorted({row['hour'] for row in stats}); counts = SAMPLER.prefixes(n)
    index = {}; expected = {(h, c) for h in hours for c in counts}
    for row in stats:
        identity = (row['hour'], row['statistics_sample_count'])
        require(identity in expected and identity not in index, 'statistics population/time identity mismatch')
        index[identity] = row
    require(set(index) == expected, 'missing statistics prefixes')
    v = int(stats[0]['saved_visualization_count']); ids = SAMPLER.choose_ids(n, v)
    require(len(ids) == v, 'invalid visualization count')
    references = np.array([vector(index[h, n], 'reference') for h in hours])
    epochs = np.array([index[h, n]['epoch_et'] for h in hours])
    for ti, h in enumerate(hours):
        for c in counts:
            row = index[h, c]
            require(row['saved_visualization_count'] == v, 'statistics/visualization counts changed')
            _near(vector(row, 'reference'), references[ti], 'prefix reference', rtol=0, atol=0)
            require(row['epoch_et'] == epochs[ti], 'prefix epoch mismatch')
            cov = covariance(row)
            audit_statistics_row(row)
            _near(cov, cov.T, 'covariance symmetry', rtol=1e-12, atol=1e-28)
            scale = np.asarray(metadata.get('covariance_scale', np.sqrt(np.maximum(np.diag(cov), 0))))
            require(scale.shape == (6,) and np.all(np.isfinite(scale)) and np.all(scale >= 0), 'invalid covariance units/scale')
            scale = np.where(scale > 0, scale, 1.)
            eigenvalues = np.linalg.eigvalsh(cov/np.outer(scale, scale))
            require(eigenvalues[0] >= -1e-10*max(1, eigenvalues[-1]), 'non-PSD scaled covariance')
            require(row['ads_accuracy_valid'] in (0, 1) and row['precision_confirmed'] in (0, 1), 'bad validity flags')
    require(len(hours) > 1 and hours[0] == 0 and np.all(np.diff(hours) > 0), 'bad probability time grid')
    _near(np.diff(hours), np.diff(hours)[0], 'hourly grid', rtol=0, atol=1e-10)
    _near(epochs, epochs[0]+np.array(hours)*3600, 'absolute epochs', rtol=0, atol=2e-6)
    input_hashes = {str(path): SAMPLER.digest(path) for path in
                    (Path(samples_path), Path(manifest_path), stats_path, visual_path, key_path)}
    published_hashes = {Path(path).name: value for path, value in input_hashes.items()}
    require(len(published_hashes) == len(input_hashes), 'input filenames must have unique identities')

    def checked_frames(path, role):
        seen = set()
        for identity, data in iter_frames(path):
            h, epoch, kind = identity
            require(h in hours and h not in seen and (not seen or h > max(seen)), 'duplicate/off-grid/nonmonotonic point frame')
            ti = hours.index(h); row = index[h, n]
            frame_ids = ids if role == 'visualization' else np.arange(1, n+1)
            require(data.shape == (len(frame_ids), 14), 'point frame count mismatch')
            require(np.array_equal(data[:, 0], frame_ids), 'duplicate/missing/nonsequential sample IDs')
            require(np.all(data[:, 1] == np.floor(data[:, 1])) and np.all(data[:, 1] >= 1)
                    and np.all(data[:, 1] <= row['patches']), 'invalid owning patch ID')
            require(epoch == epochs[ti], 'point epoch mismatch')
            require(kind == 'hourly_visualization' if role == 'visualization' else
                    bool(kind) and all(part in {'start','end','fixed','first_split','before_first_split',
                        'first_depth_limit','before_first_depth_limit','first_validation_failure',
                        'before_first_validation_failure'} for part in kind.split('|')), 'invalid frame role/label')
            _near(data[:, 2:8]-data[:, 8:14], references[ti], 'state-error reference', rtol=0,
                  atol=np.array([2e-9]*3+[2e-13]*3))
            if role == 'keyframes':
                audit_keyframe(data[:, 8:14], references[ti], [index[h, c] for c in counts])
            seen.add(h)
            yield ti, kind, data
        if role == 'visualization':
            require(seen == set(hours), 'missing hourly visualization frames')
        else:
            required = {hours[0], hours[-1]} | (set(hours) & {0, 24, 72, 168, 240, 360})
            require(required <= seen, 'missing fixed full keyframes')

    # First pass validates inputs before creating HDF5. Second pass remains frame-bounded.
    for _ in checked_frames(visual_path, 'visualization'): pass
    key_info = [(ti, kind) for ti, kind, _ in checked_frames(key_path, 'keyframes')]
    require(key_info, 'missing keyframes')
    target = Path(str(prefix)+'_probability_points.h5')
    require(not target.exists(), 'probability HDF5 already exists')
    temp = target.with_name(target.name+'.tmp')
    require(not temp.exists(), 'HDF5 staging already exists')
    convergence = convergence_report(stats, n)
    valid = [bool(index[h, n]['ads_accuracy_valid']) for h in hours]
    confirmed = [bool(index[h, n]['precision_confirmed']) for h in hours]
    with h5py.File(temp, 'x') as output:
        output.attrs.update(schema='srp_modes_probability_hdf5_v1', statistics_sample_count=n,
                            visualization_sample_count=v, keyframe_sample_count=n,
                            state_units='km,km,km,km/s,km/s,km/s', covariance_units='rr:km^2;rv:km^2/s;vv:km^2/s^2',
                            covariance_denominator='n-1', quantiles='linear_(n-1)*p',
                            prior=json.dumps(manifest['prior'], sort_keys=True),
                            reference='same_run_cannonball_real_nominal',
                            metadata_json=json.dumps(metadata, sort_keys=True, allow_nan=False),
                            input_sha256=json.dumps(published_hashes, sort_keys=True),
                            numpy_version=np.__version__, h5py_version=h5py.__version__,
                            keyframe_audit='independent_numpy_full_population',
                            regular_hour_audit='native_full_population_stats;CSV_to_HDF5;not_recomputed_from_visualization')
        output['hour'] = hours; output['epoch_et'] = epochs; output['reference_state'] = references
        for name, data in (('sample_id', population[:, 0].astype(np.int64)), ('z_state', population[:, 1:7]),
                           ('unit_coordinates', population[:, 7:15]), ('mode_parameters', population[:, 15:17]),
                           ('physical_weights', population[:, 17:19]), ('weight', np.full(n, 1/n))):
            output['population/'+name] = data
        output['population'].attrs['point_role'] = 'statistics_population_initial_definition'
        for role, size, frames, frame_path in (('visualization', v, len(hours), visual_path),
                                               ('keyframes', n, len(key_info), key_path)):
            group = output.create_group(role)
            group.attrs['point_role'] = 'fixed_hourly_subset' if role == 'visualization' else 'full_population_selected_hours'
            group['sample_id'] = ids if role == 'visualization' else np.arange(1, n+1)
            group['weight'] = np.full(size, 1/size)
            if role == 'keyframes':
                group['hour_index'] = np.array([ti for ti, _ in key_info], dtype=np.int64)
                group.create_dataset('kind', data=[kind for _, kind in key_info], dtype=h5py.string_dtype())
            datasets = {name: group.create_dataset(name, shape=(frames, size, 6), dtype='f8',
                            chunks=(1, min(size, 2048), 6), compression='gzip', compression_opts=4, shuffle=True)
                        for name in ('state', 'error')}
            patches = group.create_dataset('patch_id', shape=(frames, size), dtype='i4', chunks=(1, min(size, 2048)),
                                            compression='gzip', compression_opts=4, shuffle=True)
            for fi, (ti, _, data) in enumerate(checked_frames(frame_path, role)):
                slot = ti if role == 'visualization' else fi
                datasets['state'][slot] = data[:, 2:8]; datasets['error'][slot] = data[:, 8:14]
                patches[slot] = data[:, 1].astype(np.int32)
                require(np.array_equal(datasets['state'][slot], data[:, 2:8]) and
                        np.array_equal(datasets['error'][slot], data[:, 8:14]) and
                        np.array_equal(patches[slot], data[:, 1]), 'HDF5/CSV write verification failed')
        output['moments/sample_count'] = counts
        output['moments'].attrs['point_role'] = 'full_population_prefix_statistics_not_visualization'
        for name, getter in (('mean_state', lambda row: vector(row, 'mean_state')),
                             ('mean_error', lambda row: vector(row, 'mean_error')), ('covariance', covariance)):
            output['moments/'+name] = np.array([[getter(index[h, c]) for c in counts] for h in hours])
        output['validity/ads_accuracy_valid'] = np.asarray(valid, dtype=bool)
        output['validity/precision_confirmed'] = np.asarray(confirmed, dtype=bool)
        for fi, (ti, _) in enumerate(key_info):
            for name in ('state', 'error', 'patch_id'):
                require(np.array_equal(output['visualization/'+name][ti],
                                       output['keyframes/'+name][fi, ids-1]), 'subset/keyframe states differ')
    require({path: SAMPLER.digest(path) for path in input_hashes} == input_hashes, 'inputs changed during HDF5 packing')
    with h5py.File(temp, 'r') as output:
        require(output.attrs['schema'] == 'srp_modes_probability_hdf5_v1' and
                output['moments/covariance'].shape == (len(hours), len(counts), 6, 6), 'HDF5 reopen audit failed')
    temp.rename(target)
    numerical = all(valid) and all(confirmed)
    convergence.update(numerical_validation_passed=numerical,
                       probability_output_accepted=numerical and convergence['confirmed'])
    convergence_path = Path(str(prefix)+'_probability_convergence.json')
    with convergence_path.open('x') as stream:
        json.dump(convergence, stream, indent=2, allow_nan=False); stream.write('\n')
    return {'statistics_sample_count': n, 'visualization_sample_count': v, 'keyframe_sample_count': n,
            'hour_count': len(hours), 'keyframe_count': len(key_info), 'statistics_row_count': len(stats),
            'numerical_validation_passed': numerical, 'sample_convergence_confirmed': convergence['confirmed'],
            'sample_convergence_status': convergence['status'],
            'probability_output_accepted': convergence['probability_output_accepted'],
            'points_hdf5': target.name, 'points_hdf5_sha256': SAMPLER.digest(target),
            'convergence_file': convergence_path.name, 'input_sha256': published_hashes,
            'visualization_ids_sha256': SAMPLER.ids_digest(ids)}


def audit_probability_run(staging, report, population, requested_probes, reconstructor, history, read_rows):
    """Bind extra Real checks and all hourly saved polynomials to the native run."""
    staging = Path(staging)
    n = len(population); ids = SAMPLER.choose_ids(n, requested_probes)
    raw = report.get('probability_outputs', {})
    require(raw.get('statistics_sample_count') == n and raw.get('validation_count') == len(ids),
            'native probability population/probe count mismatch')
    require(raw.get('validation_failed') is False, 'native probability validation failed')
    cadence = report['save_hours']
    probes = read_rows(staging/'c2_probability_probes.csv')
    require(len(probes) == len(history)*len(ids), 'probability probe/hour count mismatch')
    lookup = {}; maxima = np.zeros((len(history), 2)); allowed_ids = set(ids.tolist())
    for row in probes:
        sample_id = int(row['sample_id']); hour = float(row['hour']); step = round(hour/cadence)
        identity = (step, sample_id)
        require(float(row['sample_id']) == sample_id and sample_id in allowed_ids and
                0 <= step < len(history) and abs(hour-step*cadence) < 1e-10 and identity not in lookup,
                'invalid/duplicate probability probe identity')
        _near([float(row['s']), float(row['u'])], population[sample_id-1, 15:17], 'fixed probability probe s/u', rtol=0, atol=1e-14)
        ads = np.array([float(row['ads_'+axis]) for axis in ('x','y','z','vx','vy','vz')])
        truth = np.array([float(row['real_'+axis]) for axis in ('x','y','z','vx','vy','vz')])
        errors = [np.linalg.norm(ads[:3]-truth[:3]), np.linalg.norm(ads[3:]-truth[3:])]
        _near([float(row['ads_numerical_error_pos_km']), float(row['ads_numerical_error_vel_kms'])], errors,
              'probability ADS/Real error', atol=1e-15)
        require(errors[0] <= report['position_tolerance_km'] and errors[1] <= report['velocity_tolerance_kms'],
                'probability independent Real gate failed')
        maxima[step] = np.maximum(maxima[step], errors)
        lookup[identity] = ads
    _near(raw['max_position_error_km'], maxima[:, 0].max(), 'probability report position gate')
    _near(raw['max_velocity_error_kms'], maxima[:, 1].max(), 'probability report velocity gate', atol=1e-15)
    stats = read_statistics(staging/'c2_probability_hourly_stats.csv')
    stats_index = {}; invalid = False; unconfirmed = False
    for row in stats:
        step = round(row['hour']/cadence)
        require(0 <= step < len(history) and abs(row['hour']-step*cadence) < 1e-10,
                'probability statistics off native grid')
        original = history[step]
        reference = [float(original['nominal_'+axis]) for axis in ('x_km','y_km','z_km','vx_kms','vy_kms','vz_kms')]
        _near(vector(row, 'reference'), reference, 'same-run Cannonball reference', rtol=0, atol=0)
        require(row['epoch_et'] == float(original['epoch_et']) and row['patches'] == float(original['patches'])
                and row['depth_limited'] == float(original['depth_limited']), 'native probability history mismatch')
        require(row['saved_visualization_count'] == raw['visualization_sample_count'], 'native visualization count mismatch')
        _near([row['original_probe_pos_km'], row['original_probe_vel_kms']],
              [float(original['ads_numerical_error_pos_km']), float(original['ads_numerical_error_vel_kms'])], 'original probe stats')
        _near([row['probability_probe_pos_km'], row['probability_probe_vel_kms']], maxima[step], 'probability probe stats')
        stats_index[step, int(row['statistics_sample_count'])] = row
    for step, original in enumerate(history):
        invalid = invalid or maxima[step, 0] > report['position_tolerance_km'] or maxima[step, 1] > report['velocity_tolerance_kms']
        invalid = invalid or float(original['ads_numerical_error_pos_km']) > report['position_tolerance_km'] or \
            float(original['ads_numerical_error_vel_kms']) > report['velocity_tolerance_kms'] or \
            float(original['zero_mode_pos_km']) > 1e-6 or float(original['zero_mode_vel_kms']) > 1e-9
        unconfirmed = unconfirmed or float(original['depth_limited']) > 0
        for count in SAMPLER.prefixes(n):
            require((step, count) in stats_index, 'probability statistics prefix missing')
            row = stats_index[step, count]
            require(row['ads_accuracy_valid'] == int(not invalid) and row['precision_confirmed'] == int(not unconfirmed),
                    'probability cumulative validity mismatch')
    seen = set(); checks = 0
    for prefix in report['snapshot_prefixes']:
        path = Path(prefix)
        if not path.is_absolute(): path = staging/path
        snapshot = reconstructor.load_snapshot(path)
        step = round(snapshot.time_seconds/(3600*cadence))
        require(step not in seen and 0 <= step < len(history), 'duplicate/off-grid probability snapshot')
        for sample_id in ids:
            patch_id, state = snapshot.evaluate_global(population[sample_id-1, 7:15])
            _near(state, lookup[step, sample_id], 'probability live/snapshot state', rtol=0,
                  atol=np.array([1e-8]*3+[1e-12]*3))
            require(patch_id >= 1, 'snapshot probability point is outside patches')
            checks += 1
        seen.add(step)
    require(seen == set(range(len(history))), 'probability hourly full snapshots missing')
    # Event/full frame contract is independently checked against original split/depth/gate history.
    required = {0: {'start', 'fixed'}, len(history)-1: {'end'}}
    for hour in (0,24,72,168,240,360):
        step = round(hour/cadence)
        if step < len(history) and abs(step*cadence-hour) < 1e-10:
            required.setdefault(step, set()).add('fixed')
    for label, predicate in (
        ('first_split', lambda i: any(float(history[i][f'split_dim{k}']) > 0 for k in range(1,9))),
        ('first_depth_limit', lambda i: float(history[i]['depth_limited']) > 0),
        ('first_validation_failure', lambda i: stats_index[i,n]['ads_accuracy_valid'] == 0)):
        first = next((i for i in range(len(history)) if predicate(i)), None)
        if first is not None:
            required.setdefault(first, set()).add(label)
            if first > 0: required.setdefault(first-1, set()).add('before_'+label)
    actual = {}
    for (hour, _, kind), _ in iter_frames(staging/'c2_probability_keyframes.csv'):
        step = round(hour/cadence)
        require(step not in actual, 'duplicate full keyframe')
        actual[step] = set(kind.split('|'))
    require(actual == required, 'fixed/event full keyframe labels or counts mismatch')
    return {'validation_count': len(ids), 'snapshot_probe_checks': checks,
            'max_position_error_km': float(maxima[:,0].max()), 'max_velocity_error_kms': float(maxima[:,1].max()),
            'numerical_validation_passed': not invalid and not unconfirmed,
            'hourly_full_snapshot_count': len(seen), 'validation_ids': ids.tolist()}
