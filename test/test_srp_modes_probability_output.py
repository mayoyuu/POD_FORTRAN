import csv
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

import h5py
import numpy as np
from test_srp_modes_probability_samples import load_generator

ROOT = Path(__file__).resolve().parents[1]


def load_packer():
    path = ROOT/'tools/pack_srp_modes_probability_hdf5.py'
    assert path.is_file(), 'probability HDF5 packer is not implemented'
    spec = importlib.util.spec_from_file_location('probability_packer', path)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def fixture(folder, scenario=''):
    prefix = Path(folder)/'c2'
    paths = list((ROOT/'build').glob('gfortran_*/test/test_srp_modes_probability_io'))
    assert len(paths) == 1, 'compile native probability IO test first'
    env = os.environ.copy()
    env['LD_LIBRARY_PATH'] = ':'.join(map(str, [ROOT/'lib', Path.home()/'local_dace/lib',
                                              Path.home()/'local_openblas/lib']))+':'+env.get('LD_LIBRARY_PATH', '')
    result = subprocess.run([str(paths[0]), str(prefix), scenario], cwd=ROOT, env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout+result.stderr
    samples, manifest = Path(folder)/'population.csv', Path(folder)/'population.json'
    load_generator().write_sample_bundle(samples, manifest, count=8, visualization_count=2)
    return prefix, samples, manifest


class ProbabilityOutput(unittest.TestCase):
    def test_regular_hour_row_invariants_are_checked_without_full_point_cloud(self):
        packer = load_packer()
        for column in ('mean_state1', 'std1', 'skewness_valid1'):
            with self.subTest(column=column), tempfile.TemporaryDirectory() as folder:
                prefix, samples, manifest = fixture(folder)
                path = Path(str(prefix)+'_probability_hourly_stats.csv')
                with path.open() as stream:
                    reader = csv.DictReader(stream); fields = reader.fieldnames; rows = list(reader)
                rows[1][column] = '999' if column != 'skewness_valid1' else '0'
                with path.open('w', newline='') as stream:
                    writer = csv.DictWriter(stream, fieldnames=fields); writer.writeheader(); writer.writerows(rows)
                with self.assertRaises(ValueError): packer.pack_probability(prefix, samples, manifest)

    def test_published_input_identity_has_no_expired_staging_paths(self):
        packer = load_packer()
        with tempfile.TemporaryDirectory() as folder:
            prefix, samples, manifest = fixture(folder)
            result = packer.pack_probability(prefix, samples, manifest)
            self.assertTrue(all(Path(name).name == name for name in result['input_sha256']),
                            'published identities must survive staging directory rename')
            with h5py.File(str(prefix)+'_probability_points.h5') as data:
                self.assertTrue(all(Path(name).name == name for name in json.loads(data.attrs['input_sha256'])))

    def test_split_population_uses_global_modes_and_deterministic_patch_ids(self):
        packer = load_packer()
        with tempfile.TemporaryDirectory() as folder:
            prefix, samples, manifest = fixture(folder, 'split')
            frames = list(packer.iter_frames(str(prefix)+'_probability_keyframes.csv'))
            self.assertEqual(set(frames[0][1][:,1]), {1,2,3,4})

    def test_numerically_zero_variance_is_not_reported_as_valid_gaussian_moments(self):
        packer = load_packer()
        with tempfile.TemporaryDirectory() as folder:
            prefix, samples, manifest = fixture(folder, 'nearzero')
            row = packer.read_statistics(str(prefix)+'_probability_hourly_stats.csv')[0]
            self.assertGreater(row['cov_66'], 0)
            self.assertEqual(row['skewness_valid6'], 0, 'must match compute_group_stats numerical-zero threshold')
            self.assertIsNone(row['skewness6'])
            packer.pack_probability(prefix, samples, manifest)

    def test_first_events_keep_previous_full_frame_once_with_combined_labels(self):
        packer = load_packer()
        with tempfile.TemporaryDirectory() as folder:
            prefix, samples, manifest = fixture(folder, 'events')
            frames = list(packer.iter_frames(str(prefix)+'_probability_keyframes.csv'))
            self.assertEqual([identity[0] for identity, _ in frames], [0, 1, 2])
            self.assertEqual(set(frames[0][0][2].split('|')), {'start','fixed','before_first_split'})
            self.assertEqual(set(frames[1][0][2].split('|')), {'first_split','before_first_depth_limit',
                                                           'before_first_validation_failure'})
            result = packer.pack_probability(prefix, samples, manifest)
            self.assertEqual(result['keyframe_count'], 3)
            self.assertFalse(result['numerical_validation_passed'])

    def test_native_moments_and_layered_hdf5(self):
        packer = load_packer()
        with tempfile.TemporaryDirectory() as folder:
            prefix, samples, manifest = fixture(folder)
            result = packer.pack_probability(prefix, samples, manifest, metadata={'frame':'J2000', 'center':'EARTH'})
            self.assertEqual(result['hour_count'], 3)
            self.assertEqual(result['keyframe_count'], 2)
            with h5py.File(str(prefix)+'_probability_points.h5') as data:
                self.assertEqual(data['visualization/state'].shape, (3, 2, 6))
                self.assertEqual(data['keyframes/state'].shape, (2, 8, 6))
                self.assertEqual(data['moments/covariance'].shape, (3, 1, 6, 6))
                np.testing.assert_allclose(data['population/weight'][:].sum(), 1)
                np.testing.assert_allclose(data['visualization/weight'][:].sum(), 1)
                self.assertAlmostEqual(data['moments/covariance'][0, 0, 0, 3], .06, places=12)
                np.testing.assert_allclose(data['moments/covariance'][0, 0],
                                           np.cov(data['keyframes/error'][0], rowvar=False, ddof=1), atol=1e-12)
            stats = packer.read_statistics(str(prefix)+'_probability_hourly_stats.csv')
            self.assertIsNone(stats[0]['skewness6'])
            self.assertEqual(stats[0]['skewness_valid6'], 0)
            self.assertEqual(result['sample_convergence_status'], 'not_applicable_smoke')
            self.assertFalse(result['probability_output_accepted'])

    def test_rejects_duplicate_missing_nonfinite_and_wrong_reference(self):
        packer = load_packer()
        for change in ('duplicate', 'missing', 'nonfinite', 'reference'):
            with self.subTest(change=change), tempfile.TemporaryDirectory() as folder:
                prefix, samples, manifest = fixture(folder)
                path = Path(str(prefix)+'_probability_visualization.csv')
                lines = path.read_text().splitlines()
                if change == 'duplicate': lines[3] = lines[2]
                elif change == 'missing': lines.pop()
                else:
                    fields = lines[2].split(',')
                    fields[5] = 'nan' if change == 'nonfinite' else '999999'
                    lines[2] = ','.join(fields)
                path.write_text('\n'.join(lines)+'\n')
                with self.assertRaises(ValueError): packer.pack_probability(prefix, samples, manifest)
                self.assertFalse(Path(str(prefix)+'_probability_points.h5').exists())

    def test_rejects_wrong_full_covariance_or_population_count(self):
        packer = load_packer()
        for column in ('cov_14', 'statistics_sample_count'):
            with self.subTest(column=column), tempfile.TemporaryDirectory() as folder:
                prefix, samples, manifest = fixture(folder)
                path = Path(str(prefix)+'_probability_hourly_stats.csv')
                with path.open() as stream:
                    reader = csv.DictReader(stream); fields = reader.fieldnames; rows = list(reader)
                rows[0][column] = '999'
                with path.open('w', newline='') as stream:
                    writer = csv.DictWriter(stream, fieldnames=fields); writer.writeheader(); writer.writerows(rows)
                with self.assertRaises(ValueError): packer.pack_probability(prefix, samples, manifest)

    def test_convergence_does_not_certify_smoke(self):
        packer = load_packer()
        self.assertEqual(packer.convergence_report([], 8)['status'], 'not_applicable_smoke')
        # Zero dispersion is a supported degenerate result, not division by zero.
        with tempfile.TemporaryDirectory() as folder:
            prefix, _, _ = fixture(folder)
            rows = packer.read_statistics(str(prefix)+'_probability_hourly_stats.csv')[:1]
            toy = []
            for n in (10000, 20000, 50000):
                row = dict(rows[0]); row['statistics_sample_count'] = n
                for key in row:
                    if key.startswith(('cov_', 'mean_error', 'position_p', 'velocity_p')): row[key] = 0.
                toy.append(row)
            result = packer.convergence_report(toy, 50000)
            self.assertTrue(result['confirmed'])
            self.assertTrue(all(item['degenerate'] for item in result['metrics'].values()))
            toy[-2]['mean_error1'] = 1.
            result = packer.convergence_report(toy, 50000)
            self.assertFalse(result['confirmed'])
            self.assertIsNone(result['metrics']['position_mean']['relative_change'])
            json.dumps(result, allow_nan=False)


if __name__ == '__main__': unittest.main()
