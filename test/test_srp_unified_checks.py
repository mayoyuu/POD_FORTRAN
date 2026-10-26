import contextlib
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'tools'))
import check_t1_eta_prior as eta_check
import check_t1_c0_coverage as cov_check

POINTS = 'sample_id,x,y,z,vx,vy,vz,eta,s,u,wm,wc'
MEANS = [f'mean{i}' for i in range(1, 7)]
COVARIANCE = [f'cov{i}{j}' for i in range(1, 7) for j in range(1, 7)]


def write_case(root, orbit, case, eta, prior='uniform', support=1.0, state='state-sha'):
    directory = Path(root)/orbit/case
    directory.mkdir(parents=True)
    n = len(eta)
    data = np.column_stack((np.arange(1, n+1), np.zeros((n, 6)), eta, np.full(n, .5),
                            np.full(n, .5), np.full(n, 1/n), np.full(n, 1/n)))
    np.savetxt(directory/'initial_points.csv', data, delimiter=',', header=POINTS, comments='',
               fmt=['%d']+['%.17e']*11)
    (directory/'population.json').write_text(json.dumps(dict(
        eta_prior=prior, eta_support=support, input_state_sha256=state,
        joint_population_sha256=f'joint-{case}-{np.sum(eta):.17e}')))


def campaign(root, t1_eta, prior='uniform', support=1.0, t1_state='state-sha'):
    for case in ('C0', 'C1', 'C2', 'C3'):
        write_case(root, 'O1', case, np.zeros(len(t1_eta)), support=support)
    write_case(root, 'O1', 'T1', t1_eta, prior=prior, support=support, state=t1_state)
    return Path(root)


def write_moments(path, means, covs, hours=None, valid=1.):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    hours = np.arange(len(means), dtype=float) if hours is None else np.asarray(hours, dtype=float)
    rows = [[h, *mean, *np.asarray(cov).reshape(-1), valid, valid]
            for h, mean, cov in zip(hours, means, covs)]
    header = ','.join(['hour', *MEANS, *COVARIANCE, 'accuracy_valid', 'precision_confirmed'])
    np.savetxt(path, np.array(rows), delimiter=',', header=header, comments='', fmt='%.12e')


def coverage_campaign(root, shift, extra, shrink=0., hours=None):
    base = np.diag([100., 100., 100., 1e-4, 1e-4, 1e-4])
    n = 3 if hours is None else len(hours)
    means0 = [np.array([700000., 0., 0., 0., .6, 0.]) for _ in range(n)]
    covs0 = [base.copy() for _ in range(n)]
    means1 = [mean+np.asarray(shift, dtype=float)*i for i, mean in enumerate(means0)]
    covs1 = [base*(1.+shrink*i)+np.asarray(extra, dtype=float)*i for i in range(n)]
    write_moments(Path(root)/'O1'/'C0'/'UT'/'hourly_moments.csv', means0, covs0, hours)
    write_moments(Path(root)/'O1'/'T1'/'UT'/'hourly_moments.csv', means1, covs1, hours)
    return Path(root)


def audit(root, expect='uniform', alpha=.001):
    return eta_check.check_orbit(Path(root), 'O1', expect, None, 3., alpha)


def failed_names(report):
    return [c['name'] for c in report['checks'] if not c['ok']]


class EtaPriorAuditTests(unittest.TestCase):
    def test_uniform_prior_passes_every_check(self):
        rng = np.random.default_rng(7)
        with tempfile.TemporaryDirectory() as tmp:
            report = audit(campaign(tmp, rng.uniform(-1., 1., 4000)))
        self.assertTrue(report['passed'], failed_names(report))
        self.assertAlmostEqual(report['statistics']['eta_std'], 1/np.sqrt(3.), delta=.02)
        self.assertGreater(report['statistics']['ks_p_value'], .001)

    def test_truncated_gaussian_fails_the_uniform_audit(self):
        rng = np.random.default_rng(11)
        eta = np.clip(rng.normal(0., 1/3., 4000), -1., 1.)
        with tempfile.TemporaryDirectory() as tmp:
            report = audit(campaign(tmp, eta))
        self.assertFalse(report['passed'])
        self.assertIn('T1 eta is uniform on the ADS domain', failed_names(report))
        self.assertIn('T1 eta std matches uniform variance', failed_names(report))

    def test_truncated_gaussian_tail_is_detected_with_enough_samples(self):
        rng = np.random.default_rng(41)
        eta = np.clip(rng.normal(0., 1/3., 20000), -1., 1.)
        with tempfile.TemporaryDirectory() as tmp:
            report = audit(campaign(tmp, eta))
        self.assertIn('T1 populates the SRP~0 tail at the uniform rate', failed_names(report))

    def test_eta_outside_the_physical_bound_fails(self):
        rng = np.random.default_rng(13)
        eta = rng.uniform(-1.05, 1., 4000)
        with tempfile.TemporaryDirectory() as tmp:
            report = audit(campaign(tmp, eta))
        self.assertIn('T1 eta inside [-1, 1]', failed_names(report))

    def test_missing_prior_metadata_fails(self):
        rng = np.random.default_rng(17)
        with tempfile.TemporaryDirectory() as tmp:
            report = audit(campaign(tmp, rng.uniform(-1., 1., 2000), prior=None))
        self.assertIn('T1 records eta_prior', failed_names(report))

    def test_prior_name_mismatch_fails(self):
        rng = np.random.default_rng(19)
        with tempfile.TemporaryDirectory() as tmp:
            report = audit(campaign(tmp, rng.uniform(-1., 1., 2000), prior='normal'))
        self.assertIn('T1 prior name matches request', failed_names(report))

    def test_nonzero_eta_in_c0_fails_the_shared_design_check(self):
        rng = np.random.default_rng(23)
        with tempfile.TemporaryDirectory() as tmp:
            write_case(tmp, 'O1', 'C0', np.full(2000, .01))
            for case in ('C1', 'C2', 'C3'):
                write_case(tmp, 'O1', case, np.zeros(2000))
            write_case(tmp, 'O1', 'T1', rng.uniform(-1., 1., 2000))
            report = audit(tmp)
        self.assertIn('C0 keeps eta exactly zero', failed_names(report))

    def test_different_state_hash_fails_the_shared_design_check(self):
        rng = np.random.default_rng(29)
        with tempfile.TemporaryDirectory() as tmp:
            report = audit(campaign(tmp, rng.uniform(-1., 1., 2000), t1_state='other-sha'))
        self.assertIn('C0..C3 and T1 share one bit-identical 6D design', failed_names(report))

    def test_legacy_normal_prior_passes_when_requested(self):
        rng = np.random.default_rng(31)
        eta = np.clip(rng.normal(0., .3, 4000), -.9, .9)
        with tempfile.TemporaryDirectory() as tmp:
            report = audit(campaign(tmp, eta, prior='normal', support=.9), expect='normal')
        self.assertTrue(report['passed'], failed_names(report))
        self.assertIsNone(report['statistics']['expected_std'])

    def test_uniform_ks_detects_a_shifted_interval(self):
        rng = np.random.default_rng(37)
        eta = rng.uniform(0., 1., 4000)
        with tempfile.TemporaryDirectory() as tmp:
            report = audit(campaign(tmp, eta))
        self.assertIn('T1 eta is uniform on the ADS domain', failed_names(report))


class CoverageCheckTests(unittest.TestCase):
    def args(self, **overrides):
        values = dict(max_m6=6., lambda_tol=-.05, volume_tol=.9, relative_floor=1e-12)
        values.update(overrides)
        return type('Args', (), values)()

    def test_covariance_growth_passes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = coverage_campaign(tmp, [2., 0., 0., 0., 0., 0.],
                                     np.diag([4., 1., 1., 0., 0., 0.]))
            report = cov_check.check_pair(root, 'O1', 'UT', self.args())
        self.assertTrue(report['passed'], report['checks'])
        self.assertEqual(report['hour_count'], 3)
        self.assertAlmostEqual(report['worst']['m6']['m6'], .4, places=6)
        self.assertGreaterEqual(report['worst']['lambda_relative']['lambda_relative'], 0.)

    def test_shrinking_covariance_and_large_shift_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = coverage_campaign(tmp, [40., 0., 0., 0., 0., 0.],
                                     np.diag([-20., 0., 0., 0., 0., 0.]))
            report = cov_check.check_pair(root, 'O1', 'UT', self.args())
        self.assertFalse(report['passed'])
        names = [c['name'] for c in report['checks'] if not c['ok']]
        self.assertIn('mean shift stays inside the C0 cloud', names)
        self.assertIn('Sigma_T1 - Sigma_C0 is positive semidefinite', names)
        self.assertIn('T1 covariance volume is not smaller than C0', names)

    def test_mismatched_hour_grid_is_reported(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = coverage_campaign(tmp, np.zeros(6), np.zeros((6, 6)))
            base = np.diag([100., 100., 100., 1e-4, 1e-4, 1e-4])
            write_moments(root/'O1'/'T1'/'UT'/'hourly_moments.csv',
                          [np.zeros(6)]*2, [base]*2, hours=[0., 2.])
            report = cov_check.check_pair(root, 'O1', 'UT', self.args())
        self.assertFalse(report['passed'])
        self.assertIn('hour grids differ', report['reason'])

    def test_invalid_hours_are_skipped(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = coverage_campaign(tmp, [2., 0., 0., 0., 0., 0.],
                                     np.diag([4., 1., 1., 0., 0., 0.]))
            write_moments(root/'O1'/'T1'/'UT'/'hourly_moments.csv',
                          [np.zeros(6)]*3, [np.eye(6)]*3, valid=0.)
            report = cov_check.check_pair(root, 'O1', 'UT', self.args())
        self.assertFalse(report['passed'])
        self.assertIn('numerically valid', report['reason'])

    def test_missing_case_is_reported(self):
        with tempfile.TemporaryDirectory() as tmp:
            report = cov_check.check_pair(Path(tmp), 'O1', 'UT', self.args())
        self.assertFalse(report['passed'])
        self.assertIn('hourly_moments.csv', report['reason'])

    def test_cli_round_trip_writes_json(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = coverage_campaign(tmp, [2., 0., 0., 0., 0., 0.],
                                     np.diag([4., 1., 1., 0., 0., 0.]))
            out = Path(tmp)/'report.json'
            with contextlib.redirect_stdout(io.StringIO()):
                code = cov_check.main(['--root', str(root), '--method', 'UT', '--json', str(out)])
            payload = json.loads(out.read_text())
        self.assertEqual(code, 0)
        self.assertEqual(payload['schema'], 't1_c0_coverage_v1')
        self.assertTrue(payload['pairs'][0]['passed'])


if __name__ == '__main__':
    unittest.main()
