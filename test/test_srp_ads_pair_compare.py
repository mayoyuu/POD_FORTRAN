import csv
import importlib.util
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path


sys.dont_write_bytecode = True
SCRIPT = Path(__file__).parents[1] / "tools" / "compare_srp_ads_pairs.py"


def load_module():
    spec = importlib.util.spec_from_file_location("compare_srp_ads_pairs", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def history_row(hour: float, treatment: bool) -> dict:
    nominal = {
        "nominal_x_km": 1000.0 + hour,
        "nominal_y_km": 2000.0,
        "nominal_z_km": 3000.0,
        "nominal_vx_kms": 1.0,
        "nominal_vy_kms": 2.0,
        "nominal_vz_kms": 3.0,
    }
    effect = hour if treatment else 0.0
    row = {
        "hour": hour,
        "epoch_et": 100.0 + 3600.0 * hour,
        **nominal,
    }
    components = ("dx_km", "dy_km", "dz_km", "dvx_kms", "dvy_kms", "dvz_kms")
    for index, component in enumerate(components):
        low = effect if index == 0 else (0.001 * effect if index == 3 else 0.0)
        high = 3.0 * effect if index == 0 else (0.001 * effect if index == 3 else 0.0)
        row[f"lower_{component}"] = low
        row[f"upper_{component}"] = high
    row.update(
        {
            "pos_max_km": 3.0 * effect,
            "vel_max_kms": 0.001 * effect,
            "pos_axis1_km": 3.0 * effect,
            "pos_axis2_km": 2.0 * effect,
            "pos_axis3_km": effect,
            "vel_axis1_kms": 0.003 * effect,
            "vel_axis2_kms": 0.002 * effect,
            "vel_axis3_kms": 0.001 * effect,
            "patches": 2 if treatment else 1,
            "depth_limited": 0,
            "split_dim1": 1 if treatment else 0,
            "split_dim2": 0,
            "split_dim3": 0,
            "split_dim4": 0,
            "split_dim5": 0,
            "split_dim6": 0,
            "split_dim7": 1 if treatment else 0,
            "bfs_iterations": 2 if treatment else 1,
            "max_queue_size": 3 if treatment else 1,
            "step_elapsed_seconds": 0.2 if treatment else 0.1,
            "max_validation_pos_km": 1e-4 if treatment else 5e-5,
            "max_validation_vel_kms": 1e-8 if treatment else 5e-9,
        }
    )
    return row


def create_pair(root: Path) -> tuple[Path, Path, Path]:
    treatment = root / "treatment"
    control = root / "control"
    output = root / "effect"
    for prefix in (treatment, control):
        prefix.parent.mkdir(parents=True, exist_ok=True)

    common = {
        "output_schema_version": 2,
        "nominal_srp_enabled": True,
        "duration_hours": 1.0,
        "save_hours": 1.0,
        "position_tolerance_km": 0.1,
        "velocity_tolerance_kms": 1e-6,
        "da_order": 4,
        "max_split_depth": 8,
        "shape_point_count": 512,
        "validation_point_count": 207,
        "opm_sha256": "same-opm",
        "config_sha256": "same-config",
        "source_revision": "abc123",
        "first_split_hour": None,
        "first_split_interval_start_hour": None,
        "first_unconfirmed_hour": None,
        "first_failed_hour": None,
        "previous_passing_hour": None,
    }
    treatment_report = {
        **common,
        "n_variables": 7,
        "srp_uncertainty_enabled": True,
        "eta_sigma": 0.3,
        "first_eta_split_hour": None,
        "first_eta_split_interval_start_hour": None,
    }
    control_report = {
        **common,
        "n_variables": 6,
        "srp_uncertainty_enabled": False,
        "eta_sigma": 0.0,
        "first_eta_split_hour": None,
        "first_eta_split_interval_start_hour": None,
    }
    Path(str(treatment) + "_report.json").write_text(json.dumps(treatment_report))
    Path(str(control) + "_report.json").write_text(json.dumps(control_report))

    initial_t, initial_c = [], []
    for point_id in range(1, 513):
        u7 = -1.0 if point_id <= 64 else (1.0 if point_id <= 128 else 0.0)
        common_point = {
            "point_id": point_id,
            **{f"u{k}": 0.0 for k in range(1, 7)},
            "u7": u7,
            "x0_km": 1000.0,
            "y0_km": 2000.0,
            "z0_km": 3000.0,
            "vx0_kms": 1.0,
            "vy0_kms": 2.0,
            "vz0_kms": 3.0,
            "validation": 1 if point_id <= 207 else 0,
        }
        initial_t.append({**common_point, "eta": 0.3 * u7})
        initial_c.append({**common_point, "eta": 0.0})
    write_csv(Path(str(treatment) + "_initial_points.csv"), initial_t)
    write_csv(Path(str(control) + "_initial_points.csv"), initial_c)

    write_csv(
        Path(str(treatment) + "_history.csv"),
        [history_row(0.0, True), history_row(1.0, True)],
    )
    write_csv(
        Path(str(control) + "_history.csv"),
        [history_row(0.0, False), history_row(1.0, False)],
    )

    shape_t, shape_c = [], []
    for hour in (0.0, 1.0):
        for point_id in range(1, 513):
            if hour == 0.0:
                dx = dvx = 0.0
            else:
                dx = 1.0 if point_id <= 64 else (3.0 if point_id <= 128 else 2.0)
                dvx = 0.001
            shape_t.append(
                {
                    "hour": hour,
                    "point_id": point_id,
                    "dx_km": dx,
                    "dy_km": 0.0,
                    "dz_km": 0.0,
                    "dvx_kms": dvx,
                    "dvy_kms": 0.0,
                    "dvz_kms": 0.0,
                }
            )
            shape_c.append(
                {
                    "hour": hour,
                    "point_id": point_id,
                    "dx_km": 0.0,
                    "dy_km": 0.0,
                    "dz_km": 0.0,
                    "dvx_kms": 0.0,
                    "dvy_kms": 0.0,
                    "dvz_kms": 0.0,
                }
            )
    write_csv(Path(str(treatment) + "_shape.csv"), shape_t)
    write_csv(Path(str(control) + "_shape.csv"), shape_c)
    return treatment, control, output


class PairComparisonTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.compare = load_module()

    def test_streamed_metrics_quantile_and_odd_even_decomposition(self):
        with tempfile.TemporaryDirectory() as folder:
            treatment, control, output = create_pair(Path(folder))
            report = self.compare.compare_orbit(treatment, control, output)

            with Path(str(output) + "_paired_history.csv").open(newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 2)
            final = rows[1]
            self.assertAlmostEqual(float(final["delta_r_max_km"]), 3.0)
            self.assertAlmostEqual(float(final["delta_r_rms_km"]), math.sqrt(4.25))
            self.assertAlmostEqual(float(final["delta_r_design_p95_km"]), 3.0)
            self.assertEqual(int(final["delta_r_worst_point_id"]), 65)
            self.assertAlmostEqual(float(final["odd_dx_km_max"]), 1.0)
            self.assertAlmostEqual(float(final["even_dx_km_max"]), 2.0)
            self.assertAlmostEqual(float(final["treatment_pos_volume_proxy_km3"]), 6.0)
            self.assertTrue(report["pairing_integrity"]["complete"])
            self.assertEqual(report["checkpoint_count"], 2)
            self.assertAlmostEqual(
                report["metrics"]["delta_r_max_km"]["endpoint"], 3.0
            )
            self.assertFalse(Path(str(output) + "_paired_shape.csv").exists())

    def test_rejects_hash_mismatch(self):
        with tempfile.TemporaryDirectory() as folder:
            treatment, control, output = create_pair(Path(folder))
            path = Path(str(control) + "_report.json")
            report = json.loads(path.read_text())
            report["opm_sha256"] = "different"
            path.write_text(json.dumps(report))
            with self.assertRaisesRegex(self.compare.PairingError, "opm_sha256"):
                self.compare.compare_orbit(treatment, control, output)

    def test_rejects_point_id_mismatch(self):
        with tempfile.TemporaryDirectory() as folder:
            treatment, control, output = create_pair(Path(folder))
            path = Path(str(control) + "_shape.csv")
            with path.open(newline="") as handle:
                rows = list(csv.DictReader(handle))
            rows[512]["point_id"] = "999"
            write_csv(path, rows)
            with self.assertRaisesRegex(self.compare.PairingError, "point_id"):
                self.compare.compare_orbit(treatment, control, output)


if __name__ == "__main__":
    unittest.main()
