"""Unit tests for adaptive ADS time-step batch CLI wiring."""
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
import run_crtbp_srp_study_batch as batch


class AdaptiveAdsCliTests(unittest.TestCase):
    def test_defaults_and_validation(self):
        args = batch.parser().parse_args([])
        self.assertEqual(args.ads_dt_initial_seconds, 21600.0)
        self.assertEqual(args.ads_dt_min_seconds, 3600.0)
        self.assertEqual(args.ads_dt_max_seconds, 43200.0)
        batch.validate_ads_time_args(args)

        args.ads_dt_min_seconds = 30000.0
        with self.assertRaisesRegex(ValueError, "ADS time steps"):
            batch.validate_ads_time_args(args)

    def test_build_command_appends_explicit_ads_time_block(self):
        args = batch.parser().parse_args([
            "--ads-dt-initial-seconds", "7200",
            "--ads-dt-min-seconds", "30",
            "--ads-dt-max-seconds", "14400",
        ])
        cfg = dict(batch.CONFIGS["Q1_L1_main"])
        orbit = {"mu": 0.012153614091891635, "period_nd": 2.0}
        phase = {"state_nd": [1.0, 0.0, 0.0, 0.0, 0.1, 0.0]}
        cmd = batch.build_command(
            Path("/tmp/run_crtbp_srp_study"), orbit, phase, 0, 0.0,
            "Q1_L1_main", cfg, args, Path("/tmp/result.json"))
        self.assertEqual(cmd[-4:], ["ads_time", "7200", "30", "14400"])

    def test_theta_phi_places_ads_block_after_angle_block(self):
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
        self.assertEqual(cmd[-4:], ["ads_time", "21600", "3600", "43200"])

    def test_skip_existing_requires_matching_ads_controller(self):
        args = batch.parser().parse_args([])
        cfg = dict(batch.CONFIGS["Q1_L1_main"])
        result = {
            "attitude_coordinates": "uv",
            "ads_dt_initial_seconds": 21600.0,
            "ads_dt_min_seconds": 3600.0,
            "ads_dt_max_seconds": 43200.0,
        }
        self.assertTrue(batch._matching_request(result, cfg, args))
        result["ads_dt_min_seconds"] = 60.0
        self.assertFalse(batch._matching_request(result, cfg, args))
        result.pop("ads_dt_min_seconds")
        self.assertFalse(batch._matching_request(result, cfg, args))


if __name__ == "__main__":
    unittest.main()
