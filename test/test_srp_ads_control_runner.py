import csv
import json
import os
import subprocess
import tempfile
import unittest
from pathlib import Path


class SrpAdsControlRunnerTest(unittest.TestCase):
    def test_eta_zero_runs_as_six_dimensional_paired_control(self):
        executable = Path(os.environ["POD_SRP_ADS_EXE"]).resolve()
        project = Path(os.environ["POD_PROJECT_ROOT"]).resolve()
        with tempfile.TemporaryDirectory() as folder:
            prefix = Path(folder) / "eta0"
            command = [
                str(executable),
                "-opm", str(project / "input" / "DROb_20251210_9.opm"),
                "-o", str(prefix),
                "--eta-sigma", "0",
                "--hours", "0.0002777777777777778",
                "--save-hours", "0.0002777777777777778",
                "--da-order", "4",
                "--max-depth", "1",
                "--pos-tol-km", "0.1",
                "--vel-tol-kms", "1e-6",
                "-cfg", str(project / "config" / "config.txt"),
            ]
            completed = subprocess.run(
                command,
                cwd=project,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stdout)
            report = json.loads(Path(str(prefix) + "_report.json").read_text())
            self.assertEqual(report["output_schema_version"], 2)
            self.assertEqual(report["n_variables"], 6)
            self.assertFalse(report["srp_uncertainty_enabled"])
            self.assertTrue(report["nominal_srp_enabled"])
            self.assertEqual(report["eta_sigma"], 0.0)

            with Path(str(prefix) + "_initial_points.csv").open(newline="") as stream:
                initial = list(csv.DictReader(stream))
            self.assertEqual(len(initial), 512)
            self.assertEqual([int(row["point_id"]) for row in initial], list(range(1, 513)))
            self.assertTrue(all(float(row["eta"]) == 0.0 for row in initial))
            self.assertEqual(float(initial[0]["u7"]), -1.0)
            self.assertEqual(float(initial[127]["u7"]), 1.0)

            with Path(str(prefix) + "_history.csv").open(newline="") as stream:
                history = list(csv.DictReader(stream))
            self.assertEqual(len(history), 2)
            self.assertTrue(all(int(row["split_dim7"]) == 0 for row in history))

            with Path(str(prefix) + "_shape.csv").open(newline="") as stream:
                self.assertEqual(sum(1 for _ in stream), 1 + 2 * 512)

            with Path(str(prefix) + "_split_events.csv").open(newline="") as stream:
                event_header = next(csv.reader(stream))
            self.assertIn("interval_start_hour", event_header)
            self.assertIn("direction_score7", event_header)


if __name__ == "__main__":
    unittest.main()
