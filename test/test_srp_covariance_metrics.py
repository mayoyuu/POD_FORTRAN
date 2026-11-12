"""Tests for endpoint position-covariance metrics used by SRP plots."""
from pathlib import Path
import math
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from srp_covariance_metrics import position_covariance_metrics


LU_KM = 384400.0


def covariance_state(position_block):
    cov = [0.0] * 36
    for i in range(3):
        for j in range(3):
            cov[6 * i + j] = position_block[i][j]
    return cov


class PositionCovarianceMetricsTests(unittest.TestCase):
    def test_dominant_z_direction_is_reported_as_longest_axis(self):
        metrics = position_covariance_metrics(covariance_state([
            [1.0, 0.0, 0.0],
            [0.0, 4.0, 0.0],
            [0.0, 0.0, 9.0],
        ]))

        self.assertEqual(metrics.axes_km, (3.0 * LU_KM, 2.0 * LU_KM, LU_KM))
        self.assertAlmostEqual(metrics.total_rms_km, math.sqrt(14.0) * LU_KM)
        self.assertEqual(metrics.longest_direction, (0.0, 0.0, 1.0))

    def test_position_block_is_symmetrized_before_decomposition(self):
        metrics = position_covariance_metrics(covariance_state([
            [2.0, 1.0 + 1.0e-12, 0.0],
            [1.0 - 1.0e-12, 2.0, 0.0],
            [0.0, 0.0, 1.0],
        ]), lu_km=1.0)

        self.assertAlmostEqual(metrics.axes_km[0], math.sqrt(3.0), places=13)
        self.assertAlmostEqual(metrics.axes_km[1], 1.0, places=13)
        self.assertAlmostEqual(metrics.axes_km[2], 1.0, places=13)
        self.assertAlmostEqual(sum(v * v for v in metrics.longest_direction),
                               1.0, places=13)
        self.assertGreaterEqual(max(metrics.longest_direction,
                                    key=lambda value: abs(value)), 0.0)

    def test_decomposition_is_scale_invariant_for_tiny_covariance(self):
        scale = 1.0e-24
        metrics = position_covariance_metrics(covariance_state([
            [2.0 * scale, scale, 0.0],
            [scale, 2.0 * scale, 0.0],
            [0.0, 0.0, scale],
        ]), lu_km=1.0)

        self.assertAlmostEqual(metrics.axes_km[0], math.sqrt(3.0 * scale),
                               delta=1.0e-27)
        self.assertAlmostEqual(metrics.axes_km[1], math.sqrt(scale),
                               delta=1.0e-27)
        self.assertAlmostEqual(metrics.axes_km[2], math.sqrt(scale),
                               delta=1.0e-27)

    def test_roundoff_negative_eigenvalue_is_clamped(self):
        metrics = position_covariance_metrics(covariance_state([
            [2.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, -1.0e-12],
        ]), lu_km=1.0, psd_rtol=1.0e-10)

        self.assertEqual(metrics.axes_km[2], 0.0)

    def test_materially_indefinite_covariance_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "positive semidefinite"):
            position_covariance_metrics(covariance_state([
                [2.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, -0.1],
            ]), lu_km=1.0)


if __name__ == "__main__":
    unittest.main()
