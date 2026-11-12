"""Unit tests for fixed-checkpoint ADS batch CLI wiring."""
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
import run_crtbp_srp_study_batch as batch


class FixedCheckpointAdsCliTests(unittest.TestCase):
    def test_defaults_and_validation(self):
        args = batch.parser().parse_args([])
        self.assertEqual(args.n_segments, 48)
        self.assertEqual(args.tol_div, 1.0e-3)
        batch.validate_ads_args(args)

        args.tol_div = 0.0
        with self.assertRaisesRegex(ValueError, "tol-div"):
            batch.validate_ads_args(args)

    def test_build_command_appends_explicit_block(self):
        args = batch.parser().parse_args([
            "--n-segments", "24",
            "--tol-div", "1e-2",
        ])
        cfg = dict(batch.CONFIGS["Q1_L1_main"])
        orbit = {"mu": 0.012153614091891635, "period_nd": 2.0}
        phase = {"state_nd": [1.0, 0.0, 0.0, 0.0, 0.1, 0.0]}
        cmd = batch.build_command(
            Path("/tmp/run_crtbp_srp_study"), orbit, phase, 0, 0.0,
            "Q1_L1_main", cfg, args, Path("/tmp/result.json"))
        self.assertEqual(cmd[-4:], ["n_segments", "24", "tol_div", "0.01"])

    def test_theta_phi_places_block_after_angle_block(self):
        args = batch.parser().parse_args([
            "--attitude", "theta_phi", "--phi-deg", "-45", "45",
        ])
        cfg = dict(batch.CONFIGS["Q2_L2_t10"])
        cfg.update(attitude_coordinates="theta_phi",
                   theta_deg=(0.0, 10.0), phi_deg=(-45.0, 45.0))
        orbit = {"mu": 0.012153614091891635, "period_nd": 2.0}
        phase = {"state_nd": [1.0, 0.0, 0.0, 0.0, 0.1, 0.0]}
        cmd = batch.build_command(
            Path("/tmp/run_crtbp_srp_study"), orbit, phase, 0, 0.0,
            "Q2_L2_t10", cfg, args, Path("/tmp/result.json"))
        self.assertEqual(cmd[-9:-4], ["theta_phi", "0", "10", "-45", "45"])
        self.assertEqual(cmd[-4:], ["n_segments", "48", "tol_div", "0.001"])

    def test_skip_existing_requires_matching_ads_controller(self):
        args = batch.parser().parse_args([])
        cfg = dict(batch.CONFIGS["Q1_L1_main"])
        result = {
            "attitude_coordinates": "uv",
            "n_segments": 48,
            "tol_div": 1.0e-3,
        }
        self.assertTrue(batch._matching_request(result, cfg, args))
        result["n_segments"] = 24
        self.assertFalse(batch._matching_request(result, cfg, args))
        result["n_segments"] = 48
        result["tol_div"] = 1.0e-2
        self.assertFalse(batch._matching_request(result, cfg, args))


if __name__ == "__main__":
    unittest.main()
