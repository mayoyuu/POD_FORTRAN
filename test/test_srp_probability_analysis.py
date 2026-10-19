import csv
import importlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import h5py
import numpy as np

TOOLS = Path(__file__).resolve().parents[1] / 'tools'
sys.path.insert(0, str(TOOLS))


def fixture(folder, failure=None):
    folder.mkdir(parents=True, exist_ok=True)
    rows, history = [], []
    for h in range(3):
        for n in (10000, 20000, 50000):
            r = dict(hour=h, epoch_et=h*3600, statistics_sample_count=n,
                     saved_visualization_count=2000, patches=1, depth_limited=0,
                     original_probe_pos_km=0, original_probe_vel_kms=0,
                     probability_probe_pos_km=0, probability_probe_vel_kms=0,
                     ads_accuracy_valid=1, precision_confirmed=1,
                     first_failure_hour='', first_unconfirmed_hour='')
            for i in range(1, 7):
                r.update({f'reference{i}': 10*i, f'mean_state{i}': 10*i+n/50000,
                          f'mean_error{i}': n/50000, f'std{i}': i})
                for j in range(1, 7):
                    r[f'cov_{i}{j}'] = i*i if i == j else .01*i*j
            for kind in ('position', 'velocity'):
                for k, q in enumerate(('p50', 'p95', 'p99', 'p999', 'max')):
                    r[f'{kind}_{q}'] = k+1
            if h == 2 and failure:
                r['ads_accuracy_valid'] = 0
                r['first_failure_hour'] = 2
                if failure == 'probability': r['probability_probe_pos_km'] = .2
                if failure == 'depth': r.update(depth_limited=1, precision_confirmed=0, first_unconfirmed_hour=2)
            rows.append(r)
        history.append(dict(hour=h, ads_numerical_error_pos_km=0, ads_numerical_error_vel_kms=0,
                            zero_mode_pos_km=2e-6 if h == 2 and failure == 'zero' else 0,
                            zero_mode_vel_kms=0))
    write_csv(folder/'c2_probability_hourly_stats.csv', rows)
    write_csv(folder/'c2_history.csv', history)
    return rows, history


def write_csv(path, rows):
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def cloud(folder, h5=False, hour=1, bad=False):
    ids = np.arange(1, 2001)
    if bad: ids[-1] = ids[0]
    if h5:
        with h5py.File(folder/'c2_probability_points.h5', 'w') as f:
            f.attrs['schema'] = 'srp_modes_probability_hdf5_v1'
            f['hour'] = [0, 1, 2]; f['visualization/sample_id'] = ids
            f['visualization/error'] = np.ones((3, 2000, 6))
    else:
        with (folder/'c2_probability_visualization.csv').open('w', newline='') as f:
            f.write('schema=srp_modes_probability_points_v1\n')
            fields = ['hour','sample_id','dx_km','dy_km','dz_km','dvx_kms','dvy_kms','dvz_kms']
            w = csv.writer(f); w.writerow(fields)
            for h in range(3):
                for i in ids: w.writerow([h, i]+[h]*6)


class AnalysisTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(); self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.data = importlib.import_module('srp_probability_analysis_data')

    def test_primary_covariance_and_missing_report(self):
        fixture(self.root)
        c = self.data.load_case(self.root, 2)
        np.testing.assert_equal(c['mean_error'], np.ones((3,6)))
        self.assertEqual(c['covariance'].shape, (3,6,6))
        self.assertEqual(c['covariance'][0,1,4], .1)
        self.assertEqual(c['last_passing_hour'], 2)
        self.assertEqual(c['final_report_state'], 'interrupted_unfinalized')
        self.assertIsNone(c['formal_result'])

    def test_direct_and_unknown_failure_causes(self):
        for failure, cause in [('zero','zero_mode_position'), ('probability','probability_probe_position'),
                               ('depth','depth_limit'), ('unknown','unknown_recorded_failure')]:
            with self.subTest(failure=failure):
                fixture(self.root, failure)
                c = self.data.load_case(self.root, 2)
                self.assertEqual(c['last_passing_hour'], 1)
                self.assertEqual(c['first_invalid_hour'], 2)
                self.assertIn(cause, c['first_failure_causes'])

    def test_strict_bad_statistics(self):
        for problem in ('gap','duplicate','nonfinite','std','prefix','truncated'):
            with self.subTest(problem=problem):
                rows, _ = fixture(self.root)
                if problem == 'gap': rows = [r for r in rows if r['hour'] != 1]
                if problem == 'duplicate': rows.append(rows[-1])
                if problem == 'nonfinite': rows[0]['mean_error1'] = 'nan'
                if problem == 'std': rows[-1]['std1'] = 4
                if problem == 'prefix': rows[0]['ads_accuracy_valid'] = 0
                write_csv(self.root/'c2_probability_hourly_stats.csv', rows)
                if problem == 'truncated':
                    with (self.root/'c2_probability_hourly_stats.csv').open('a') as f: f.write('2,')
                with self.assertRaisesRegex(ValueError, 'c2_probability_hourly_stats.csv'):
                    self.data.load_case(self.root,2)

    def test_malformed_trailing_warns_only_beyond_window(self):
        fixture(self.root)
        with (self.root/'c2_probability_hourly_stats.csv').open('a') as f: f.write('3,')
        self.assertTrue(self.data.load_case(self.root,2)['warnings'])

    def test_unidentified_trailing_record_is_never_ignored(self):
        rows, _ = fixture(self.root)
        rows.append(dict(rows[-1],hour=3))
        write_csv(self.root/'c2_probability_hourly_stats.csv',rows)
        with (self.root/'c2_probability_hourly_stats.csv').open('a') as f: f.write('broken,')
        with self.assertRaises(ValueError): self.data.load_case(self.root,2)

    def test_no_recovery_after_recorded_failure(self):
        rows, _ = fixture(self.root)
        for r in rows:
            if r['hour'] == 1:
                r.update(ads_accuracy_valid=0,first_failure_hour=1,probability_probe_vel_kms=2e-6)
        write_csv(self.root/'c2_probability_hourly_stats.csv',rows)
        c = self.data.load_case(self.root,2)
        np.testing.assert_equal(c['recorded_prefix'],[True,False,False])

    def test_original_probe_failure_and_native_report_preserved(self):
        rows, _ = fixture(self.root)
        for r in rows:
            if r['hour'] == 2: r['original_probe_pos_km'] = .2
        write_csv(self.root/'c2_probability_hourly_stats.csv',rows)
        (self.root/'c2_report.json').write_text(json.dumps({'formal_result':False,'validation_failed':True}))
        c = self.data.load_case(self.root,2)
        self.assertEqual(c['final_report_state'],'completed_failed')
        self.assertIs(c['formal_result'],False)
        self.assertIn('original_probe_position',c['first_failure_causes'])

    def test_nonfinite_cloud_rejected(self):
        cloud(self.root,True)
        with h5py.File(self.root/'c2_probability_points.h5','a') as f: f['visualization/error'][1,0,0] = np.nan
        with self.assertRaises(ValueError): self.data.load_cloud(self.root,1)

    def test_invalid_grid_counts_and_negative_errors_rejected(self):
        for field,value in [('hour',.5),('depth_limited',-.5),('saved_visualization_count',1999),
                            ('original_probe_pos_km',-.1)]:
            with self.subTest(field=field):
                rows,_ = fixture(self.root)
                rows[0][field] = value
                write_csv(self.root/'c2_probability_hourly_stats.csv',rows)
                with self.assertRaises(ValueError): self.data.load_case(self.root,2)

    def test_csv_and_hdf5_cloud_exact_hour_and_identity(self):
        for use_h5 in (False, True):
            with self.subTest(h5=use_h5):
                cloud(self.root, use_h5)
                ids, values = self.data.load_cloud(self.root, 1)
                self.assertEqual(values.shape, (2000,6)); self.assertEqual(ids[-1],2000)
                with self.assertRaises(ValueError): self.data.load_cloud(self.root, 3)
                cloud(self.root,use_h5,bad=True)
                with self.assertRaises(ValueError): self.data.load_cloud(self.root,1)
                (self.root/('c2_probability_points.h5' if use_h5 else 'c2_probability_visualization.csv')).unlink()

    def test_discovery_and_duplicates(self):
        build_campaign(self.root)
        self.assertEqual(len(self.data.discover_cases(self.root)),15)
        extra = self.root/self.data.CAMPAIGNS[1]/'results'/'DRO-1'
        extra.mkdir()
        with self.assertRaisesRegex(ValueError, 'duplicate'): self.data.discover_cases(self.root)

    def test_cli_artifacts_and_npz(self):
        build_campaign(self.root)
        run = subprocess.run([sys.executable, str(TOOLS/'analyze_boxwing_probability.py'), '--root',
                              str(self.root), '--max-hour','2'], capture_output=True,text=True)
        self.assertEqual(run.returncode,0,run.stderr)
        out = self.root/'analysis'
        self.assertTrue((out/'campaign_summary.csv').exists())
        self.assertTrue((out/'provenance.json').exists())
        with np.load(out/'DRO-2'/'moments.npz') as archive:
            self.assertEqual(archive['covariance'].shape,(3,6,6))
        self.assertEqual(len(list(out.glob('*/*.png'))),60)


def build_campaign(root):
    names = ['boxwing_modes_DRO-2_init_c2_probability_15d_20261016',
             'boxwing_modes_six_orbit_probability_20261016', 'boxwing_modes_remaining10_probability_20261017']
    orbit_names = [f'{family}-{i}' for family in ('DRO','L1Halo','L1PL','L2Halo','L2PL') for i in (1,2,3)]
    for orbit in orbit_names:
        path = root/names[0] if orbit == 'DRO-2' else root/names[2]/'results'/f'.{orbit}.staging-test'
        fixture(path); cloud(path,h5=True)
    (root/names[1]/'results').mkdir(parents=True)


if __name__ == '__main__': unittest.main()
