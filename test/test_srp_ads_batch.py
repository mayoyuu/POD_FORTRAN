import importlib.util
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace


sys.dont_write_bytecode = True

SCRIPT = Path(__file__).parents[1] / "tools" / "run_srp_ads_batch.py"
SPEC = importlib.util.spec_from_file_location("run_srp_ads_batch", SCRIPT)
batch = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(batch)


def write_opm(path: Path, pos_rss_km: float, vel_rss_mps: float) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pos_var = pos_rss_km**2 / 3.0
    vel_var = (vel_rss_mps / 1000.0) ** 2 / 3.0
    data = {
        "CX_X": pos_var,
        "CY_Y": pos_var,
        "CZ_Z": pos_var,
        "CX_DOT_X_DOT": vel_var,
        "CY_DOT_Y_DOT": vel_var,
        "CZ_DOT_Z_DOT": vel_var,
    }
    path.write_text(json.dumps(data), encoding="utf-8")


def campaign_args(eta_sigma: float = 0.3) -> SimpleNamespace:
    return SimpleNamespace(
        days=15.0,
        save_hours=1.0,
        eta_sigma=eta_sigma,
        pos_tol_km=0.1,
        vel_tol_kms=1e-6,
        da_order=4,
        max_depth=8,
    )


def compatible_report(prefix: Path, args: SimpleNamespace, opm: Path, config: Path) -> None:
    report = {
        "output_schema_version": 2,
        "n_variables": 6 if args.eta_sigma == 0.0 else 7,
        "nominal_srp_enabled": True,
        "srp_uncertainty_enabled": args.eta_sigma > 0.0,
        "duration_hours": args.days * 24.0,
        "save_hours": args.save_hours,
        "eta_sigma": args.eta_sigma,
        "position_tolerance_km": args.pos_tol_km,
        "velocity_tolerance_kms": args.vel_tol_kms,
        "da_order": args.da_order,
        "max_split_depth": args.max_depth,
        "opm_sha256": batch.file_sha256(opm),
        "config_sha256": batch.file_sha256(config),
    }
    Path(str(prefix) + "_report.json").write_text(json.dumps(report), encoding="utf-8")
    Path(str(prefix) + "_split_events.csv").write_text(
        "interval_start_hour,interval_end_hour,parent_depth,component,direction,direction_score7\n", encoding="utf-8"
    )


class BatchDiscoveryTest(unittest.TestCase):
    def test_discovers_mirrors_and_filters_expected_initial_error(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            opm = root / "OPM"
            output = root / "SRP" / "SRP_UCERTAINTY_PROP_WITHINIT"
            good = opm / "DRO-1" / "DRO-1_init.opm.json"
            bad = opm / "supp_batch" / "nested" / "DRO-1_corrected_init.opm.json"
            ignored = opm / "DRO-2" / "DRO-2_final.opm.json"
            write_opm(good, 10.0, 0.03)
            write_opm(bad, 1.0, 0.01)
            write_opm(ignored, 10.0, 0.03)

            found = batch.find_initial_opms(opm)
            self.assertEqual(found, [good])
            pos_rss, vel_rss = batch.covariance_rss(good)
            self.assertTrue(math.isclose(pos_rss, 10.0, rel_tol=1e-12))
            self.assertTrue(math.isclose(vel_rss, 0.03, rel_tol=1e-12))
            self.assertTrue(batch.matches_expected_error(good, 10.0, 0.03))
            self.assertFalse(batch.matches_expected_error(bad, 10.0, 0.03))
            self.assertEqual(
                batch.output_prefix_for(good, opm, output),
                output / "DRO-1" / "DRO-1_init",
            )


class CompatibilityAuditTest(unittest.TestCase):
    def test_requires_schema_dimension_split_log_and_current_hashes(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            prefix = root / "orbit"
            opm = root / "init.opm.json"
            config = root / "config.txt"
            opm.write_text("opm", encoding="utf-8")
            config.write_text("config", encoding="utf-8")
            args = campaign_args(0.0)

            status, reasons = batch.audit_report(prefix, args, opm, config)
            self.assertEqual(status, "missing_result")
            self.assertIn("report_missing_or_invalid", reasons)

            compatible_report(prefix, args, opm, config)
            status, reasons = batch.audit_report(prefix, args, opm, config)
            self.assertEqual((status, reasons), ("compatible", []))

            report_path = Path(str(prefix) + "_report.json")
            report = json.loads(report_path.read_text(encoding="utf-8"))
            del report["output_schema_version"]
            report_path.write_text(json.dumps(report), encoding="utf-8")
            status, reasons = batch.audit_report(prefix, args, opm, config)
            self.assertEqual(status, "incompatible")
            self.assertIn("output_schema_version", reasons)

    def test_successful_temp_prefix_promotion_replaces_old_set(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            prefix = root / "orbit"
            temp_prefix = root / "orbit.__new__"
            Path(str(prefix) + "_report.json").write_text("old", encoding="utf-8")
            stale = Path(str(prefix) + "_obsolete.csv")
            stale.write_text("stale", encoding="utf-8")
            Path(str(temp_prefix) + "_report.json").write_text("new", encoding="utf-8")
            Path(str(temp_prefix) + "_history.csv").write_text("history", encoding="utf-8")

            promoted = batch.promote_prefix(temp_prefix, prefix)

            self.assertEqual(
                Path(str(prefix) + "_report.json").read_text(encoding="utf-8"), "new"
            )
            self.assertTrue(Path(str(prefix) + "_history.csv").is_file())
            self.assertFalse(stale.exists())
            self.assertEqual(promoted, 2)

    def test_eta_zero_dry_run_is_six_dimensional(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            opm_root = root / "OPM"
            output_root = root / "ETA0"
            config = root / "config.txt"
            write_opm(opm_root / "DRO-1" / "DRO-1_init.opm.json", 10.0, 0.03)
            config.write_text("test config", encoding="utf-8")

            code = batch.main(
                [
                    "--opm-root", str(opm_root),
                    "--output-root", str(output_root),
                    "--config", str(config),
                    "--eta-sigma", "0",
                    "--dry-run",
                ]
            )

            self.assertEqual(code, 0)
            summary = json.loads((output_root / "batch_summary.json").read_text())
            self.assertEqual(summary["eta_ads_domain"], [0.0, 0.0])
            self.assertEqual(summary["n_variables"], 6)


if __name__ == "__main__":
    unittest.main()
