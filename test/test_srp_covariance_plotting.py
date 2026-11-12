"""Integration tests for SRP covariance plotting tools."""
from pathlib import Path
import json
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
import analyze_srp_beta_scan_summary as beta_summary
import analyze_srp_cov_axis as cov_axis
import analyze_srp_split_vs_axis as split_axis


LU_KM = 384400.0


def covariance_state():
    cov = [0.0] * 36
    cov[0] = 1.0
    cov[7] = 4.0
    cov[14] = 9.0
    return cov


class CovariancePlottingTests(unittest.TestCase):
    def test_all_plotters_use_true_largest_eigenvalue(self):
        for module in (cov_axis, beta_summary, split_axis):
            with self.subTest(module=module.__name__):
                self.assertEqual(module.axis_km(covariance_state()), 3.0 * LU_KM)

    def test_covariance_scan_exposes_three_axes_total_rms_and_direction(self):
        with tempfile.TemporaryDirectory() as tmp:
            case_dir = Path(tmp) / "Q2" / "H1N_A005000"
            case_dir.mkdir(parents=True)
            (case_dir / "p00__a0.json").write_text(json.dumps({
                "status": 0,
                "t_end_days": 15.0,
                "cov_state": covariance_state(),
            }))
            jobs = cov_axis.scan_results(Path(tmp))

        self.assertEqual(len(jobs), 1)
        job = jobs[0]
        self.assertEqual(job["axis_km"], 3.0 * LU_KM)
        self.assertEqual(job["axes_km"], (3.0 * LU_KM, 2.0 * LU_KM, LU_KM))
        self.assertAlmostEqual(job["total_rms_km"], (14.0 ** 0.5) * LU_KM)
        self.assertEqual(job["longest_direction"], (0.0, 0.0, 1.0))

    def test_covariance_plot_accepts_metric_selector(self):
        for metric in ("longest", "second", "third", "total-rms"):
            args = cov_axis.parser().parse_args(["--metric", metric])
            self.assertEqual(args.metric, metric)


if __name__ == "__main__":
    unittest.main()
