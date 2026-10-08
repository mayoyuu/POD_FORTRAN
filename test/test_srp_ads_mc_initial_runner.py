#!/usr/bin/env python3
import csv
import json
import os
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PYTHON = Path("/home/songyu/venvs/astro/bin/python")


class InitialOnlyRunnerTest(unittest.TestCase):
    def run_fpm(self, arguments):
        env = os.environ.copy()
        lib_paths = [
            str(ROOT / "lib"),
            str(Path.home() / "local_dace" / "lib"),
            str(Path.home() / "local_openblas" / "lib"),
        ]
        env["LIBRARY_PATH"] = ":".join(lib_paths + [env.get("LIBRARY_PATH", "")])
        env["LD_LIBRARY_PATH"] = ":".join(lib_paths + [env.get("LD_LIBRARY_PATH", "")])
        return subprocess.run(
            ["fpm", "run", "--target", "run_srp_ads_mc_initial_only", "--", *arguments],
            cwd=ROOT,
            env=env,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
        )

    def test_cli_rejects_invalid_contracts(self):
        help_result = self.run_fpm(["--help"])
        self.assertEqual(help_result.returncode, 0, help_result.stdout[-2000:])
        self.assertIn("run_srp_ads_mc_initial_only", help_result.stdout)

        missing = self.run_fpm([])
        self.assertNotEqual(missing.returncode, 0)
        self.assertIn("required", missing.stdout)

        positive_eta = self.run_fpm(["--eta-sigma", "0.3"])
        self.assertNotEqual(positive_eta.returncode, 0)
        self.assertIn("unknown option", positive_eta.stdout)

    def test_hour_zero_earth_output_is_six_dimensional(self):
        opm_candidates = (
            ROOT / "OPM" / "L1Halo-1" / "L1Halo-1_init.opm.json",
            ROOT.parent.parent / "OPM" / "L1Halo-1" / "L1Halo-1_init.opm.json",
        )
        opm = next((path for path in opm_candidates if path.is_file()), opm_candidates[0])
        self.assertTrue(opm.is_file())

        with tempfile.TemporaryDirectory() as folder:
            folder = Path(folder)
            samples = folder / "samples.csv"
            manifest = folder / "samples.json"
            subprocess.run(
                [
                    str(PYTHON),
                    str(ROOT / "tools" / "generate_srp_mc_samples.py"),
                    "--csv", str(samples),
                    "--manifest", str(manifest),
                    "--count", "32",
                ],
                cwd=ROOT,
                check=True,
                stdout=subprocess.DEVNULL,
            )
            prefix = folder / "earth"
            result = self.run_fpm([
                "-opm", str(opm),
                "-o", str(prefix),
                "--mc-samples", str(samples),
                "--mc-count", "32",
                "--visualization-count", "16",
                "--hours", "0",
                "-cfg", str(ROOT / "config" / "config_box_wing_earth.txt"),
            ])
            self.assertEqual(result.returncode, 0, result.stdout[-4000:])

            report = json.loads(Path(str(prefix) + "_mc_report.json").read_text())
            self.assertEqual(report["output_schema_version"], "srp_ads_mc_initial_only_v1")
            self.assertEqual(report["active_dimensions"], 6)
            self.assertEqual(report["dace_storage_dimensions"], 6)
            self.assertEqual(report["srp_eta"], 0.0)
            self.assertEqual(report["srp_model"], "box_wing")
            self.assertEqual(report["srp_attitude_mode"], "earth")
            self.assertEqual(report["sample_count"], 32)

            with Path(str(prefix) + "_mc_hourly_stats.csv").open(newline="") as stream:
                hourly = list(csv.DictReader(stream))
            self.assertEqual(len(hourly), 3)
            self.assertEqual({row["group"] for row in hourly}, {"earth"})
            self.assertEqual({int(row["active_dimensions"]) for row in hourly}, {6})
            self.assertEqual({int(row["split_dim7"]) for row in hourly}, {0})

            point_lines = Path(str(prefix) + "_mc_keyframe_points.csv").read_text().splitlines()
            self.assertEqual(point_lines[0], "schema=srp_mc_initial_points_v1")
            rows = list(csv.DictReader(point_lines[1:]))
            self.assertEqual(len(rows), 32)
            self.assertEqual([int(row["sample_id"]) for row in rows], list(range(1, 33)))
            for row in rows:
                self.assertAlmostEqual(
                    float(row["x_km"]) - float(row["dx_km"]),
                    float(hourly[0]["center_x_km"]),
                    places=10,
                )

            with Path(str(prefix) + "_split_events.csv").open(newline="") as stream:
                event_header = next(csv.reader(stream))
            self.assertIn("direction_score7", event_header)

    def test_forced_split_records_full_first_event_and_limits(self):
        opm = ROOT.parent.parent / "OPM" / "L1Halo-1" / "L1Halo-1_init.opm.json"
        with tempfile.TemporaryDirectory() as folder:
            folder = Path(folder)
            samples = folder / "samples.csv"
            manifest = folder / "samples.json"
            subprocess.run([
                str(PYTHON), str(ROOT / "tools" / "generate_srp_mc_samples.py"),
                "--csv", str(samples), "--manifest", str(manifest), "--count", "32",
            ], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
            prefix = folder / "forced"
            result = self.run_fpm([
                "-opm", str(opm), "-o", str(prefix),
                "--mc-samples", str(samples), "--mc-count", "32",
                "--visualization-count", "16", "--hours", "0.01",
                "--save-hours", "0.01", "--max-depth", "1",
                "--pos-tol-km", "1e-20", "--vel-tol-kms", "1e-25",
                "-cfg", str(ROOT / "config" / "config_box_wing_earth.txt"),
            ])
            self.assertEqual(result.returncode, 0, result.stdout[-4000:])
            with Path(str(prefix) + "_split_events.csv").open(newline="") as stream:
                events = list(csv.DictReader(stream))
            self.assertGreaterEqual(len(events), 1)
            self.assertIn(int(events[0]["direction"]), range(1, 7))
            report = json.loads(Path(str(prefix) + "_mc_report.json").read_text())
            first = report["events"]["first_split"]
            self.assertTrue(first["observed"])
            self.assertIsInstance(first["full_criterion"], dict)
            self.assertEqual(first["full_criterion"]["direction"],
                             int(events[0]["direction"]))
            self.assertTrue(report["events"]["first_unconfirmed"]["observed"])
            self.assertTrue(
                report["events"]["first_validation_failure"]["observed"]
            )

    def test_rejects_cannonball_config(self):
        opm_candidates = (
            ROOT / "OPM" / "L1Halo-1" / "L1Halo-1_init.opm.json",
            ROOT.parent.parent / "OPM" / "L1Halo-1" / "L1Halo-1_init.opm.json",
        )
        opm = next((path for path in opm_candidates if path.is_file()), opm_candidates[0])
        with tempfile.TemporaryDirectory() as folder:
            folder = Path(folder)
            samples = folder / "samples.csv"
            manifest = folder / "samples.json"
            subprocess.run([
                str(PYTHON), str(ROOT / "tools" / "generate_srp_mc_samples.py"),
                "--csv", str(samples), "--manifest", str(manifest), "--count", "2",
            ], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
            result = self.run_fpm([
                "-opm", str(opm), "-o", str(folder / "bad"),
                "--mc-samples", str(samples), "--mc-count", "2", "--hours", "0",
                "-cfg", str(ROOT / "config" / "config.txt"),
            ])
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("requires srp_model=box_wing", result.stdout)


if __name__ == "__main__":
    unittest.main()
