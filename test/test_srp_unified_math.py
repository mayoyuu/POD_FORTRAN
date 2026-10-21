import sys
import unittest
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'tools'))
from srp_unified_math import shared_design, joint_sigma_points, weighted_moments, henze_zirkler


class JointPropagationMathTests(unittest.TestCase):
    def test_shared_design_is_reproducible_bounded_and_non_degenerate(self):
        a = shared_design(2000, 20261021, 3.0)
        b = shared_design(2000, 20261021, 3.0)
        np.testing.assert_array_equal(a, b)
        self.assertEqual(a.shape, (2000, 9))
        self.assertTrue(np.all(np.abs(a[:, :7]) <= 3))
        self.assertTrue(np.all((a[:, 7:] >= 0) & (a[:, 7:] <= 1)))
        self.assertTrue(np.all(np.std(a, axis=0) > 0))

    def test_joint_ut_recovers_correlated_state_and_parameter_moments(self):
        rng = np.random.default_rng(71)
        x = rng.normal(size=(3000, 8))
        x[:, 6] = .2*x[:, 0] + .3*x[:, 6]
        x[:, 7] = .1*x[:, 2] + .5*x[:, 7]
        points, wm, wc = joint_sigma_points(x)
        mean, cov = weighted_moments(points, wm, wc)
        self.assertEqual(points.shape, (17, 8))
        np.testing.assert_allclose(mean, x.mean(axis=0), atol=1e-13)
        np.testing.assert_allclose(cov, np.cov(x, rowvar=False), atol=1e-13)
        self.assertGreater(abs(cov[0, 6]), .1)

    def test_ut_keeps_bounded_uniform_modes_inside_domain(self):
        d = shared_design(5000, 12)
        x = np.column_stack((d[:, :6], d[:, 7:9]))
        points, _, _ = joint_sigma_points(x)
        self.assertTrue(np.all((points[:, 6:] >= 0) & (points[:, 6:] <= 1)))

    def test_ut_handles_singular_covariance_without_jitter(self):
        x = np.column_stack((np.linspace(-1, 1, 100), np.ones(100), np.zeros(100)))
        points, wm, wc = joint_sigma_points(x)
        mean, cov = weighted_moments(points, wm, wc)
        np.testing.assert_allclose(mean, x.mean(axis=0), atol=1e-13)
        np.testing.assert_allclose(cov, np.cov(x, rowvar=False), atol=1e-13)
        self.assertTrue(np.all(points[:, 1] == 1))

    def test_hz_is_affine_invariant_and_detects_non_gaussian_mixture(self):
        rng = np.random.default_rng(173)
        x = rng.normal(size=(600, 3))
        a = np.array([[2., .2, .4], [.1, 3., -.3], [.3, .6, 4.]])
        r1 = henze_zirkler(x)
        r2 = henze_zirkler(x @ a + [10., -30., 500.])
        self.assertEqual(r1['status'], 'ok')
        self.assertAlmostEqual(r1['statistic'], r2['statistic'], places=10)
        self.assertAlmostEqual(r1['p_value'], r2['p_value'], places=10)
        x[:300, 0] -= 4
        x[300:, 0] += 4
        self.assertLess(henze_zirkler(x)['p_value'], .001)

    def test_hz_reports_singular_input_without_fake_pvalue(self):
        r = henze_zirkler(np.ones((100, 6)))
        self.assertEqual(r['status'], 'singular_covariance')
        self.assertIsNone(r['p_value'])

    def test_validation_rejects_invalid_inputs(self):
        with self.assertRaises(ValueError): shared_design(1, 5)
        with self.assertRaises(ValueError): shared_design(5, 5, -1)
        with self.assertRaises(ValueError): joint_sigma_points(np.full((3, 2), np.nan))


if __name__ == '__main__':
    unittest.main()
