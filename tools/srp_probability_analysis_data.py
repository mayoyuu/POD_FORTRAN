"""Read saved SRP probability products without changing or propagating inputs."""
import csv
import json
import math
from pathlib import Path
import re

import h5py
import numpy as np

CAMPAIGNS = ('boxwing_modes_DRO-2_init_c2_probability_15d_20261016',
             'boxwing_modes_six_orbit_probability_20261016',
             'boxwing_modes_remaining10_probability_20261017')
ORBITS = tuple(f'{family}-{i}' for family in ('DRO','L1Halo','L1PL','L2Halo','L2PL') for i in (1,2,3))
UNITS = 'km,km,km,km/s,km/s,km/s'
COVARIANCE_UNITS = 'rr:km^2;rv:km^2/s;vv:km^2/s^2'
LIMITS = {'original_probe_position': .1, 'original_probe_velocity': 1e-6,
          'probability_probe_position': .1, 'probability_probe_velocity': 1e-6,
          'zero_mode_position': 1e-6, 'zero_mode_velocity': 1e-9}
ERROR_FIELDS = {'original_probe_position':'original_probe_pos_km',
                'original_probe_velocity':'original_probe_vel_kms',
                'probability_probe_position':'probability_probe_pos_km',
                'probability_probe_velocity':'probability_probe_vel_kms',
                'zero_mode_position':'zero_mode_pos_km','zero_mode_velocity':'zero_mode_vel_kms'}


def require(ok, message):
    if not ok: raise ValueError(message)


def discover_cases(root: Path) -> dict[str, Path]:
    root = Path(root); found = {'DRO-2': root/CAMPAIGNS[0]}
    require(found['DRO-2'].is_dir(), f'missing campaign {found["DRO-2"]}')
    for campaign in CAMPAIGNS[1:]:
        results = root/campaign/'results'
        require(results.is_dir(), f'missing campaign results {results}')
        for folder in sorted(results.iterdir()):
            if not folder.is_dir(): continue
            match = re.fullmatch(r'\.(.+)\.staging-.+', folder.name)
            name = match.group(1) if match else folder.name
            if name not in ORBITS: continue
            require(name not in found, f'duplicate orbit {name}: {found.get(name)} and {folder}')
            found[name] = folder
    require(set(found) == set(ORBITS), f'missing orbits: {sorted(set(ORBITS)-set(found))}')
    return {name: found[name] for name in ORBITS}


def _read_numeric(path, max_hour, warnings):
    """Validate all complete rows; only an identifiable out-of-window tail may warn."""
    result = []
    with path.open(newline='', encoding='utf-8-sig') as stream:
        reader = csv.DictReader(stream)
        require(reader.fieldnames and len(set(reader.fieldnames)) == len(reader.fieldnames), f'{path}: invalid header')
        records = iter(reader)
        for row in records:
            line = reader.line_num
            hour = None
            try:
                hour = float(row.get('hour',''))
                require(math.isfinite(hour), 'invalid hour')
                require(None not in row and all(v is not None for v in row.values()), 'truncated/extra fields')
                converted = {}
                for key, value in row.items():
                    if not value.strip():
                        require(key in ('first_failure_hour','first_unconfirmed_hour') or
                                key.startswith(('skewness','excess_kurtosis')), f'blank {key}')
                        converted[key] = None
                    else:
                        converted[key] = float(value)
                        require(math.isfinite(converted[key]), f'nonfinite {key}')
                result.append(converted)
            except (ValueError, TypeError) as error:
                # A damaged record is tolerated only if it is the physical final record
                # and its explicit finite hour is strictly after the complete window.
                tail = list(records)
                if not tail and hour is not None and math.isfinite(hour) and hour > max_hour:
                    warnings.append(f'{path.name}:{line}: ignored malformed trailing record at hour {hour}: {error}')
                else:
                    raise ValueError(f'{path}:{line}: {error}') from error
    require(result, f'{path}: empty data')
    return result


def _vector(row, prefix):
    return np.array([row[f'{prefix}{i}'] for i in range(1,7)])


def load_case(folder: Path, max_hour: float = 300) -> dict:
    folder = Path(folder); warnings = []
    require(math.isfinite(max_hour) and max_hour >= 0 and max_hour == int(max_hour), 'max-hour must be a nonnegative integer')
    stats_path = folder/'c2_probability_hourly_stats.csv'
    try:
        stats = _read_numeric(stats_path,max_hour,warnings)
        grouped = {}
        for row in stats:
            key = (row['hour'], row['statistics_sample_count'])
            require(key not in grouped, f'duplicate hour/sample-count {key}')
            require(row['hour'] >= 0 and row['hour'].is_integer(), f'off-grid hour {row["hour"]}')
            require(row['statistics_sample_count'] in (10000,20000,50000), f'unsupported sample count {key}')
            require(row['saved_visualization_count'] == 2000, 'expected 2000 saved visualization samples')
            require(row['patches'].is_integer() and row['patches'] >= 1 and row['depth_limited'].is_integer()
                    and 0 <= row['depth_limited'] <= row['patches'], 'invalid depth/patch counts')
            require(all(row[field] >= 0 for field in ERROR_FIELDS.values() if field.startswith(('original','probability'))),
                    'negative numerical error')
            grouped[key] = row
            cov = np.array([[row[f'cov_{i}{j}'] for j in range(1,7)] for i in range(1,7)])
            std = _vector(row,'std')
            require(np.all(np.diag(cov) >= 0) and np.all(std >= 0), 'negative variance/std')
            require(np.allclose(std,np.sqrt(np.diag(cov)),rtol=1e-10,atol=1e-28), 'std/covariance mismatch')
            require(np.allclose(cov,cov.T,rtol=1e-12,atol=1e-28), 'asymmetric covariance')
            # Normalize component units before eigendecomposition: never mix km and km/s.
            scale = np.where(std > 0,std,1)
            eigen = np.linalg.eigvalsh(cov/np.outer(scale,scale))
            require(eigen.min() >= -1e-9*max(1,eigen.max()), 'non-PSD scaled covariance')
            for flag in ('ads_accuracy_valid','precision_confirmed'):
                require(row[flag] in (0,1), f'invalid {flag}')
        hours = np.arange(int(max_hour)+1,dtype=float)
        for h in hours:
            prefixes = [grouped.get((h,n)) for n in (10000,20000,50000)]
            require(all(r is not None for r in prefixes), f'missing hourly sample prefix at {h}')
            for key in ('ads_accuracy_valid','precision_confirmed','depth_limited',
                        'first_failure_hour','first_unconfirmed_hour'):
                require(len({r[key] for r in prefixes}) == 1, f'prefix {key} mismatch at {h}')
        primary = [grouped[h,50000] for h in hours]
    except (ValueError,KeyError) as error:
        raise ValueError(f'{stats_path}: {error}') from error
    history_path = folder/'c2_history.csv'
    hist = _read_numeric(history_path,max_hour,warnings); history = {}
    for row in hist:
        require(row['hour'] not in history,f'{history_path}: duplicate hour {row["hour"]}')
        history[row['hour']] = row
    require(all(h in history for h in hours),f'{history_path}: missing required hourly data')
    diagnostics = []; passing = []; alive = True; first_invalid = None; first_causes = []
    for row in primary:
        h = row['hour']; hi = history[h]
        diagnostic = dict(row)
        try:
            diagnostic.update(zero_mode_pos_km=hi['zero_mode_pos_km'],zero_mode_vel_kms=hi['zero_mode_vel_kms'],
                              history_original_probe_pos_km=hi['ads_numerical_error_pos_km'],
                              history_original_probe_vel_kms=hi['ads_numerical_error_vel_kms'])
            history_error_fields = ('zero_mode_pos_km','zero_mode_vel_kms',
                                    'ads_numerical_error_pos_km','ads_numerical_error_vel_kms')
            require(all(hi[field] >= 0 for field in history_error_fields),
                    f'{history_path}: negative history error at hour {h}')
            causes = [name for name, field in ERROR_FIELDS.items() if diagnostic[field] > LIMITS[name]]
            # Both original-probe records are retained; disagreement cannot hide a failure.
            for kind,field in [('position','ads_numerical_error_pos_km'),('velocity','ads_numerical_error_vel_kms')]:
                if hi[field] > LIMITS['original_probe_'+kind] and 'original_probe_'+kind not in causes:
                    causes.append('original_probe_'+kind)
            if row['depth_limited'] > 0: causes.append('depth_limit')
            bad_flag = not row['ads_accuracy_valid'] or not row['precision_confirmed']
            if bad_flag and not causes: causes.append('unknown_recorded_failure')
        except KeyError as error:
            raise ValueError(f'{folder}: missing diagnostic field {error}') from error
        invalid = bool(causes)
        if invalid and first_invalid is None: first_invalid, first_causes = h, causes.copy()
        alive = alive and not invalid
        passing.append(alive)
        diagnostic.update(recorded_prefix=int(alive), direct_failure_causes=';'.join(causes),
                          first_invalid=int(h == first_invalid), failure_observed=int(invalid))
        diagnostics.append(diagnostic)
    report_path = folder/'c2_report.json'
    report = json.loads(report_path.read_text(encoding='utf-8')) if report_path.exists() else None
    if report is None:
        state = 'interrupted_unfinalized'; warnings.append('native final report missing; run is interrupted/unfinalized')
    elif report.get('formal_result') is True: state = 'completed_accepted'
    elif report.get('validation_failed') or report.get('precision_unconfirmed') or report.get('probability_outputs',{}).get('validation_failed'):
        state = 'completed_failed'
    else: state = 'completed_unaccepted'
    mask = np.array(passing,dtype=bool)
    return dict(folder=folder, hours=hours, rows=primary, prefix_rows=stats, diagnostics=diagnostics,
                mean_state=np.array([_vector(r,'mean_state') for r in primary]),
                mean_error=np.array([_vector(r,'mean_error') for r in primary]),
                covariance=np.array([[[r[f'cov_{i}{j}'] for j in range(1,7)] for i in range(1,7)] for r in primary]),
                std=np.array([_vector(r,'std') for r in primary]),recorded_prefix=mask,
                last_passing_hour=float(hours[mask][-1]) if mask.any() else None,
                first_invalid_hour=first_invalid,first_failure_causes=first_causes,
                observed_saved_last_hour=max(r['hour'] for r in stats), requested_window_endpoint=float(max_hour),
                final_report_state=state,formal_result=report.get('formal_result') if report else None,
                source_revision=report.get('source_revision') if report else None, warnings=warnings, final_report=report)


def load_cloud(folder: Path, hour: float) -> tuple[np.ndarray,np.ndarray]:
    folder = Path(folder); h5path = folder/'c2_probability_points.h5'
    if h5path.exists():
        with h5py.File(h5path,'r') as stream:
            require(stream.attrs.get('schema') == 'srp_modes_probability_hdf5_v1',f'{h5path}: unsupported schema')
            indices = np.flatnonzero(stream['hour'][:] == hour)
            require(len(indices) == 1,f'{h5path}: missing/duplicate hour {hour}')
            ids = stream['visualization/sample_id'][:]
            values = stream['visualization/error'][int(indices[0]),:,:]
    else:
        path = folder/'c2_probability_visualization.csv'; ids = []; values = []
        with path.open(newline='',encoding='utf-8-sig') as stream:
            require(stream.readline().strip() == 'schema=srp_modes_probability_points_v1',f'{path}: unsupported schema')
            reader = csv.DictReader(stream)
            for row in reader:
                try:
                    h = float(row['hour'])
                    if h != hour: continue
                    ids.append(float(row['sample_id']))
                    values.append([float(row[name]) for name in ('dx_km','dy_km','dz_km','dvx_kms','dvy_kms','dvz_kms')])
                except (ValueError,TypeError,KeyError) as error:
                    raise ValueError(f'{path}:{reader.line_num}: malformed requested cloud: {error}') from error
        ids, values = np.asarray(ids), np.asarray(values)
    require(values.shape == (2000,6),f'{folder}: hour {hour} requires exactly 2000 saved vectors; got {values.shape}')
    require(ids.shape == (2000,) and np.all(np.isfinite(ids)) and np.all(ids == np.floor(ids))
            and np.all((ids >= 1)&(ids <= 50000)) and len(set(ids)) == 2000,f'{folder}: missing/duplicate/invalid sample IDs')
    require(np.all(np.isfinite(values)),f'{folder}: nonfinite cloud at hour {hour}')
    return ids.astype(np.int64),values
