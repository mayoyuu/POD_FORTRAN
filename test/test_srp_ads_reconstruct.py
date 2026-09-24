import csv
import importlib.util
import json
import math
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


sys.dont_write_bytecode = True
ROOT = Path(__file__).parents[1]
SCRIPT = ROOT / "tools" / "reconstruct_srp_ads_snapshot.py"


def load_module():
    spec = importlib.util.spec_from_file_location("reconstruct_srp_ads_snapshot", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_csv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def coefficient_row(
    patch_id: int,
    component: int,
    coefficient: float,
    powers: tuple[int, ...] = (0, 0, 0, 0, 0, 0, 0),
) -> dict:
    return {
        "patch_id": patch_id,
        "component": component,
        **{f"power{k}": powers[k - 1] for k in range(1, 8)},
        "coefficient": coefficient,
    }


def write_map(prefix: Path, eta_sigma: float) -> None:
    rows = [
        {"key": "state_unit", "i": 0, "j": 0, "value": "km_and_km_per_s"},
        {
            "key": "polynomial_variables",
            "i": 0,
            "j": 0,
            "value": "local_patch_coordinates_minus_one_to_one",
        },
        {"key": "coordinate_mode", "i": 0, "j": 0, "value": "whitened"},
        {"key": "da_order", "i": 0, "j": 0, "value": 4},
        {"key": "domain_sigma", "i": 0, "j": 0, "value": 3.0},
        {"key": "eta_sigma", "i": 0, "j": 0, "value": eta_sigma},
    ]
    write_csv(Path(str(prefix) + "_map.csv"), ["key", "i", "j", "value"], rows)


def write_unsplit_snapshot(prefix: Path, n_variables: int = 7) -> None:
    eta_sigma = 0.3 if n_variables == 7 else 0.0
    write_map(prefix, eta_sigma)
    center = [0.0] * 7
    width = [2.0] * n_variables + [0.0] * (7 - n_variables)
    meta = {
        "patch_id": 1,
        "time_seconds": 3600.0,
        **{f"center{k}": center[k - 1] for k in range(1, 8)},
        **{f"width{k}": width[k - 1] for k in range(1, 8)},
        "split_count": 0,
    }
    write_csv(
        Path(str(prefix) + "_meta.csv"),
        [
            "patch_id",
            "time_seconds",
            *[f"center{k}" for k in range(1, 8)],
            *[f"width{k}" for k in range(1, 8)],
            "split_count",
        ],
        [meta],
    )
    write_csv(
        Path(str(prefix) + "_splits.csv"),
        ["patch_id", "sequence", "split_direction"],
        [],
    )
    coeffs = [
        coefficient_row(1, 1, 10.0),
        coefficient_row(1, 1, 2.0, (1, 0, 0, 0, 0, 0, 0)),
        coefficient_row(1, 2, 20.0),
        coefficient_row(1, 3, 30.0),
        coefficient_row(1, 4, 1.0),
        coefficient_row(1, 5, 2.0),
        coefficient_row(1, 6, 3.0),
    ]
    if n_variables == 7:
        coeffs.append(
            coefficient_row(
                1, 1, -3.0, (0, 1, 0, 0, 0, 0, 1)
            )
        )
    write_csv(
        Path(str(prefix) + "_coeff.csv"),
        [
            "patch_id",
            "component",
            *[f"power{k}" for k in range(1, 8)],
            "coefficient",
        ],
        coeffs,
    )


def write_split_snapshot(prefix: Path) -> None:
    write_map(prefix, 0.3)
    rows = []
    for patch_id, center, direction in ((1, -0.5, -1), (2, 0.5, 1)):
        rows.append(
            {
                "patch_id": patch_id,
                "time_seconds": 7200.0,
                **{
                    f"center{k}": center if k == 1 else 0.0
                    for k in range(1, 8)
                },
                **{
                    f"width{k}": 1.0 if k == 1 else 2.0
                    for k in range(1, 8)
                },
                "split_count": 1,
            }
        )
    write_csv(
        Path(str(prefix) + "_meta.csv"),
        [
            "patch_id",
            "time_seconds",
            *[f"center{k}" for k in range(1, 8)],
            *[f"width{k}" for k in range(1, 8)],
            "split_count",
        ],
        rows,
    )
    write_csv(
        Path(str(prefix) + "_splits.csv"),
        ["patch_id", "sequence", "split_direction"],
        [
            {"patch_id": 1, "sequence": 1, "split_direction": -1},
            {"patch_id": 2, "sequence": 1, "split_direction": 1},
        ],
    )
    coeffs = []
    for patch_id, center in ((1, -0.5), (2, 0.5)):
        coeffs.extend(
            [
                coefficient_row(patch_id, 1, center),
                coefficient_row(
                    patch_id, 1, 0.5, (1, 0, 0, 0, 0, 0, 0)
                ),
                coefficient_row(patch_id, 2, 20.0),
                coefficient_row(patch_id, 3, 30.0),
                coefficient_row(patch_id, 4, 1.0),
                coefficient_row(patch_id, 5, 2.0),
                coefficient_row(patch_id, 6, 3.0),
            ]
        )
    write_csv(
        Path(str(prefix) + "_coeff.csv"),
        [
            "patch_id",
            "component",
            *[f"power{k}" for k in range(1, 8)],
            "coefficient",
        ],
        coeffs,
    )


def write_parent_files(run_prefix: Path) -> None:
    initial_rows = []
    shape_rows = []
    points = [
        (1, [0.2, 0.3, 0.0, 0.0, 0.0, 0.0, 0.5], [-0.05, 0, 0, 0, 0, 0]),
        (2, [0.0] * 7, [0.0] * 6),
    ]
    for point_id, unit, deviation in points:
        initial_rows.append(
            {
                "point_id": point_id,
                **{f"u{k}": unit[k - 1] for k in range(1, 8)},
            }
        )
        shape_rows.append(
            {
                "hour": 1.0,
                "point_id": point_id,
                "dx_km": deviation[0],
                "dy_km": deviation[1],
                "dz_km": deviation[2],
                "dvx_kms": deviation[3],
                "dvy_kms": deviation[4],
                "dvz_kms": deviation[5],
            }
        )
    write_csv(
        Path(str(run_prefix) + "_initial_points.csv"),
        ["point_id", *[f"u{k}" for k in range(1, 8)]],
        initial_rows,
    )
    write_csv(
        Path(str(run_prefix) + "_shape.csv"),
        [
            "hour",
            "point_id",
            "dx_km",
            "dy_km",
            "dz_km",
            "dvx_kms",
            "dvy_kms",
            "dvz_kms",
        ],
        shape_rows,
    )
    write_csv(
        Path(str(run_prefix) + "_history.csv"),
        [
            "hour",
            "nominal_x_km",
            "nominal_y_km",
            "nominal_z_km",
            "nominal_vx_kms",
            "nominal_vy_kms",
            "nominal_vz_kms",
        ],
        [
            {
                "hour": 1.0,
                "nominal_x_km": 10.0,
                "nominal_y_km": 20.0,
                "nominal_z_km": 30.0,
                "nominal_vx_kms": 1.0,
                "nominal_vy_kms": 2.0,
                "nominal_vz_kms": 3.0,
            }
        ],
    )


class SnapshotReconstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.reconstruct = load_module()

    def test_loads_and_evaluates_known_seven_dimensional_polynomial(self):
        with tempfile.TemporaryDirectory() as folder:
            prefix = Path(folder) / "orbit_patch_s0003600"
            write_unsplit_snapshot(prefix, 7)
            snapshot = self.reconstruct.load_snapshot(prefix)
            patch_id, state = snapshot.evaluate_global(
                (0.2, 0.3, 0.0, 0.0, 0.0, 0.0, 0.5)
            )
            self.assertEqual(snapshot.n_variables, 7)
            self.assertEqual(snapshot.da_order, 4)
            self.assertEqual(patch_id, 1)
            self.assertAlmostEqual(state[0], 9.95)
            self.assertEqual(state[1:], (20.0, 30.0, 1.0, 2.0, 3.0))
            self.assertAlmostEqual(snapshot.physical_eta(0.5), 0.45)

    def test_replays_split_history_and_uses_first_patch_on_boundary(self):
        with tempfile.TemporaryDirectory() as folder:
            prefix = Path(folder) / "orbit_patch_s0007200"
            write_split_snapshot(prefix)
            snapshot = self.reconstruct.load_snapshot(prefix)
            left_center, left_width = self.reconstruct.replay_split_history([-1], 7)
            right_center, right_width = self.reconstruct.replay_split_history([1], 7)
            self.assertEqual(left_center[0], -0.5)
            self.assertEqual(right_center[0], 0.5)
            self.assertEqual(left_width[0], 1.0)
            self.assertEqual(right_width[0], 1.0)
            patch_id, state = snapshot.evaluate_global((0.0,) * 7)
            self.assertEqual(patch_id, 1)
            self.assertAlmostEqual(state[0], 0.0)

    def test_generates_deterministic_mixed_samples_and_pads_six_dimensions(self):
        first = self.reconstruct.generate_local_samples(7, 150, 0.5)
        second = self.reconstruct.generate_local_samples(7, 150, 0.5)
        self.assertEqual(first, second)
        self.assertEqual(len(first), 150)
        self.assertEqual(sum(sample.kind == "corner" for sample in first), 128)
        self.assertEqual(sum(sample.kind == "axis" for sample in first), 14)
        self.assertEqual(sum(sample.kind == "center" for sample in first), 1)
        self.assertGreater(sum(sample.kind == "face" for sample in first), 0)
        self.assertGreater(sum(sample.kind == "interior" for sample in first), 0)
        six = self.reconstruct.generate_local_samples(6, 82, 0.5)
        self.assertTrue(all(len(sample.coordinates) == 7 for sample in six))
        self.assertTrue(all(sample.coordinates[6] == 0.0 for sample in six))

    def test_loads_padded_six_dimensional_control_and_rejects_power7(self):
        with tempfile.TemporaryDirectory() as folder:
            prefix = Path(folder) / "control_patch_s0003600"
            write_unsplit_snapshot(prefix, 6)
            snapshot = self.reconstruct.load_snapshot(prefix)
            self.assertEqual(snapshot.n_variables, 6)
            self.assertEqual(snapshot.physical_eta(0.8), 0.0)
            _, state = snapshot.evaluate_global((0.2, 0.3, 0, 0, 0, 0))
            self.assertAlmostEqual(state[0], 10.4)

            coeff_path = Path(str(prefix) + "_coeff.csv")
            with coeff_path.open(newline="") as handle:
                rows = list(csv.DictReader(handle))
            rows[0]["power7"] = "1"
            write_csv(coeff_path, list(rows[0]), rows)
            with self.assertRaisesRegex(
                self.reconstruct.SnapshotError, "power7"
            ):
                self.reconstruct.load_snapshot(prefix)

    def test_rejects_metadata_split_history_disagreement(self):
        with tempfile.TemporaryDirectory() as folder:
            prefix = Path(folder) / "bad_patch_s0007200"
            write_split_snapshot(prefix)
            meta_path = Path(str(prefix) + "_meta.csv")
            with meta_path.open(newline="") as handle:
                rows = list(csv.DictReader(handle))
            rows[0]["center1"] = "-0.4"
            write_csv(meta_path, list(rows[0]), rows)
            with self.assertRaisesRegex(
                self.reconstruct.SnapshotError, "split history"
            ):
                self.reconstruct.load_snapshot(prefix)

    def test_reconstructs_csv_and_verifies_saved_shape_cloud(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            run_prefix = root / "orbit"
            prefix = root / "orbit_patch_s0003600"
            write_unsplit_snapshot(prefix, 7)
            write_parent_files(run_prefix)
            output = root / "dense.csv"
            report_path = root / "dense.json"
            report = self.reconstruct.reconstruct_snapshot(
                prefix=prefix,
                output=output,
                report_path=report_path,
                points_per_patch=150,
                boundary_fraction=0.5,
                history_path=None,
                absolute_only=False,
                verify=True,
            )
            with output.open(newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 150)
            self.assertIn("local7", rows[0])
            self.assertIn("u7", rows[0])
            self.assertIn("dx_km", rows[0])
            self.assertEqual(report["output_schema_version"], 1)
            self.assertEqual(report["n_variables"], 7)
            self.assertEqual(report["verification"]["status"], "passed")
            self.assertLess(report["verification"]["max_position_error_km"], 1e-12)
            self.assertTrue(report_path.is_file())

    def test_absolute_only_standalone_snapshot_needs_no_parent_files(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            prefix = root / "standalone"
            output = root / "absolute.csv"
            write_unsplit_snapshot(prefix, 7)
            report = self.reconstruct.reconstruct_snapshot(
                prefix=prefix,
                output=output,
                points_per_patch=150,
                absolute_only=True,
                verify=False,
            )
            with output.open(newline="") as handle:
                header = next(csv.reader(handle))
            self.assertNotIn("dx_km", header)
            self.assertEqual(report["verification"]["status"], "skipped")
            self.assertNotIn("history_file", report)

    def test_cli_writes_default_report(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            run_prefix = root / "orbit"
            prefix = root / "orbit_patch_s0003600"
            write_unsplit_snapshot(prefix, 7)
            write_parent_files(run_prefix)
            output = root / "cli.csv"
            completed = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--prefix",
                    str(prefix),
                    "--points-per-patch",
                    "150",
                    "--output",
                    str(output),
                ],
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stdout)
            self.assertTrue(output.is_file())
            report_path = Path(str(prefix) + "_reconstruction.json")
            report = json.loads(report_path.read_text())
            self.assertEqual(report["total_samples"], 150)


if __name__ == "__main__":
    unittest.main()
