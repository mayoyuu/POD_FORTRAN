#!/usr/bin/env python3
"""Fixed 8D prior: six componentwise truncated normals and independent uniform s/u."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

SCHEMA = 'srp_modes_probability_samples_v1'
HEADER = (['sample_id'] + [f'z{i}' for i in range(1, 7)] +
          [f'q{i}' for i in range(1, 9)] + ['s', 'u', 'xi_E', 'xi_M'])
PRIOR = {
    'state': 'independent_standard_normal_conditioned_componentwise_abs_z_le_3',
    'modes': 'independent_uniform_s_u_on_0_1_independent_of_state',
    'trajectory_parameters': 'fixed', 'sigma_limit': 3,
    'mapping': 'q1:6=z/3;q7=2s-1;q8=2u-1;xi_E=s;xi_M=(1-s)*u',
    'state_map': 'x0=mu0+L*z;P0=L*L^T;lower_cholesky',
    'triangle_uniform': False,
}


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def choose_ids(count, requested):
    if not isinstance(count, int) or not isinstance(requested, int) or min(count, requested) < 1:
        raise ValueError('sample and visualization counts must be positive integers')
    size = min(count, requested)
    if size == 1:
        return np.array([1], dtype=np.int64)
    return 1 + np.floor(np.arange(size)*(count-1)/(size-1)+.5).astype(np.int64)


def ids_digest(ids):
    return hashlib.sha256(('\n'.join(map(str, ids))+'\n').encode('ascii')).hexdigest()


def prefixes(count):
    return [n for n in (10000, 20000, 50000) if n < count] + [count]


def generate_samples(count, seed=20261015):
    if not isinstance(count, int) or count < 2:
        raise ValueError('at least two probability samples are required')
    if not isinstance(seed, int) or seed < 0:
        raise ValueError('seed must be a nonnegative integer')
    children = np.random.SeedSequence(seed).spawn(2)
    state_rng, mode_rng = [np.random.Generator(np.random.PCG64DXSM(child)) for child in children]
    accepted, remaining, rejected = [], count, 0
    # Fixed-sized candidate chunks make accepted prefixes independent of requested N.
    while remaining:
        candidates = state_rng.standard_normal((8192, 6))
        indices = np.flatnonzero(np.all(np.abs(candidates) <= 3, axis=1))
        take = min(remaining, len(indices))
        accepted.append(candidates[indices[:take]])
        consumed = int(indices[take-1])+1 if take == remaining else len(candidates)
        rejected += consumed - take
        remaining -= take
    z = np.vstack(accepted)
    modes = mode_rng.random((count, 2))
    values = np.column_stack((np.arange(1, count+1), z, z/3, 2*modes-1,
                              modes, modes[:, 0], (1-modes[:, 0])*modes[:, 1]))
    return values, {'seed': seed, 'generator': 'numpy.random.PCG64DXSM',
                    'numpy_version': np.__version__, 'state_spawn_key': list(children[0].spawn_key),
                    'modes_spawn_key': list(children[1].spawn_key), 'rejected_state_rows': rejected}


def write_sample_bundle(csv_path, manifest_path, *, count=50000, seed=20261015, visualization_count=2000):
    csv_path, manifest_path = Path(csv_path), Path(manifest_path)
    if csv_path.resolve() == manifest_path.resolve():
        raise ValueError('CSV and manifest paths must differ')
    if csv_path.exists() or manifest_path.exists():
        raise FileExistsError('probability sample bundle already exists')
    ids = choose_ids(count, visualization_count)
    values, metadata = generate_samples(count, seed)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open('x', newline='', encoding='utf-8') as stream:
        stream.write('schema='+SCHEMA+'\n')
        writer = csv.writer(stream, lineterminator='\n')
        writer.writerow(HEADER)
        for row in values:
            writer.writerow([int(row[0])] + [format(value, '.17g') for value in row[1:]])
    metadata.update(schema=SCHEMA, sample_count=count, prior=PRIOR,
                    prefix_counts=prefixes(count), population_weight=1/count,
                    visualization_count=len(ids), visualization_ids=ids.tolist(),
                    visualization_ids_sha256=ids_digest(ids), csv_sha256=digest(csv_path))
    with manifest_path.open('x', encoding='utf-8') as stream:
        json.dump(metadata, stream, indent=2, allow_nan=False)
        stream.write('\n')
    return metadata


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError('duplicate manifest key: '+key)
        result[key] = value
    return result


def load_sample_bundle(csv_path, manifest_path, expected_count=None):
    def invalid(value):
        raise ValueError('nonfinite manifest value '+value)
    with Path(manifest_path).open(encoding='utf-8') as stream:
        meta = json.load(stream, object_pairs_hook=_unique_object, parse_constant=invalid)
    if meta.get('schema') != SCHEMA or meta.get('prior') != PRIOR:
        raise ValueError('unsupported sample schema or prior')
    if meta.get('csv_sha256') != digest(csv_path):
        raise ValueError('probability sample CSV hash mismatch')
    count = meta.get('sample_count')
    if type(count) is not int or count < 2 or (expected_count is not None and count != expected_count):
        raise ValueError('probability sample count mismatch')
    if (meta.get('generator') != 'numpy.random.PCG64DXSM' or meta.get('state_spawn_key') != [0]
            or meta.get('modes_spawn_key') != [1] or meta.get('prefix_counts') != prefixes(count)
            or meta.get('population_weight') != 1/count):
        raise ValueError('invalid population identity')
    ids = choose_ids(count, meta.get('visualization_count'))
    if meta.get('visualization_ids') != ids.tolist() or meta.get('visualization_ids_sha256') != ids_digest(ids):
        raise ValueError('visualization ID identity mismatch')
    data = np.empty((count, len(HEADER)))
    with Path(csv_path).open(newline='', encoding='utf-8') as stream:
        if stream.readline().strip() != 'schema='+SCHEMA:
            raise ValueError('invalid sample schema')
        reader = csv.reader(stream)
        if next(reader, None) != HEADER:
            raise ValueError('invalid sample fields')
        n = 0
        for n, row in enumerate(reader, 1):
            if n > count or len(row) != len(HEADER):
                raise ValueError('invalid sample record/count')
            data[n-1] = [float(x) for x in row]
        if n != count:
            raise ValueError('missing sample records')
    if not np.all(np.isfinite(data)) or not np.array_equal(data[:, 0], np.arange(1, count+1)):
        raise ValueError('nonfinite sample or duplicate/nonsequential IDs')
    z, q, modes = data[:, 1:7], data[:, 7:15], data[:, 15:17]
    if np.any(np.abs(z) > 3) or np.any(np.abs(q) > 1) or np.any(modes < 0) or np.any(modes > 1):
        raise ValueError('sample outside prior domain')
    expected = np.column_stack((z/3, 2*modes-1, modes[:, 0], (1-modes[:, 0])*modes[:, 1]))
    actual = np.column_stack((q, data[:, 17:19]))
    if not np.allclose(actual, expected, rtol=0, atol=2e-14):
        raise ValueError('invalid probability parameter mapping')
    return data, meta


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--csv', type=Path, required=True)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--count', type=int, default=50000)
    parser.add_argument('--seed', type=int, default=20261015)
    parser.add_argument('--visualization-count', type=int, default=2000)
    args = parser.parse_args(argv)
    try:
        metadata = write_sample_bundle(args.csv, args.manifest, count=args.count, seed=args.seed,
                                       visualization_count=args.visualization_count)
    except (OSError, ValueError) as exc:
        parser.exit(2, str(exc)+'\n')
    print(json.dumps(metadata, indent=2))


if __name__ == '__main__': main()
