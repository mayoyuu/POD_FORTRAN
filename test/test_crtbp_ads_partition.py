"""Offline postprocessing checks; no Fortran execution or orbit propagation."""
import csv
import json
import math
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
import plot_srp_ads_partition as boxes
import analyze_srp_study_progress as progress


def result():
    # beta splits at day 5; its left child splits along theta at day 14.
    # Events are deliberately out of storage order; i specifies replay order.
    return {
        "attitude_coordinates": "theta_phi",
        "parameter_names": ["beta", "chi", "theta", "phi"],
        "parameter_units": ["m2/kg", "1", "rad", "rad"],
        "parameter_lower": [.08, .6, 0, -math.pi / 4],
        "parameter_upper": [.12, .6, math.pi / 18, math.pi / 4],
        "t_end_days": 20, "t_end_nd": 20 / progress.ND_TO_DAYS,
        "n_splits_total": 2, "n_patches": 3, "n_diverged": 0,
        "first_divergence_time": -1, "first_split_time": 5 / progress.ND_TO_DAYS,
        "status": 0, "alpha0_rad": 0, "split_count_dir": [1, 0, 1, 0],
        "split_events": [
            {"i": 2, "td": 14, "ax": 3, "d": 1,
             "call": [-.5, 0, 0, 0], "wall": [1, 0, 2, 2]},
            {"i": 1, "td": 5, "ax": 1, "d": 0,
             "call": [0, 0, 0, 0], "wall": [2, 0, 2, 2]},
        ],
    }


class PartitionTests(unittest.TestCase):
    def test_time_replay_conserves_root_volume(self):
        for day, count in [(0, 1), (13.9, 2), (14, 3), (20, 3)]:
            snapshot = boxes.reconstruct(result(), day)
            self.assertEqual(len(snapshot["boxes"]), count)
            self.assertEqual(snapshot["partition_kind"], "active")
            volume = sum(math.prod(b["width"][j] / 2 for j in [0, 2, 3])
                         for b in snapshot["boxes"])
            self.assertAlmostEqual(volume, 1)

    def test_slice_degrees_and_half_open_hidden_boundary(self):
        snap = boxes.reconstruct(result(), 14)
        # beta=.1 belongs only to the right child, not both closed boxes.
        selected, fixed = boxes.select_slice(snap, ["theta", "phi"], {"beta": .1})
        self.assertEqual(len(selected), 1)
        selected, _ = boxes.select_slice(snap, ["theta", "phi"], {"beta": .09})
        self.assertEqual(len(selected), 2)
        selected, fixed = boxes.select_slice(snap, ["beta", "theta"], {"phi": 20})
        self.assertAlmostEqual(fixed["phi"], math.radians(20))
        self.assertEqual(len(selected), 3)
        with self.assertRaises(ValueError):
            boxes.select_slice(snap, ["theta", "phi"], {"beta": .2})

    def test_reject_missing_or_inconsistent_log(self):
        for mutate in (
            lambda r: r["split_events"][0].pop("wall"),
            lambda r: r.update(n_splits_total=3),
            lambda r: r["split_events"][0].update(call=[.5, 0, .5, 0]),
            lambda r: r.update(n_patches=4),
        ):
            r = result()
            mutate(r)
            with self.assertRaises(ValueError):
                boxes.reconstruct(r, 20)
        with self.assertRaises(ValueError):
            boxes.reconstruct(result(), 21)

    def test_divergence_is_not_an_active_partition(self):
        r = result()
        r.update(n_diverged=1, n_patches=2, first_divergence_time=10 / progress.ND_TO_DAYS)
        self.assertEqual(boxes.reconstruct(r, 9)["partition_kind"], "active")
        self.assertEqual(boxes.reconstruct(r, 14)["partition_kind"], "historical")
        r["first_divergence_time"] = -1
        self.assertEqual(boxes.reconstruct(r, 0)["partition_kind"], "historical")

    def test_q4_hidden_chi_and_legacy_uv(self):
        r = result()
        r["parameter_lower"][1], r["parameter_upper"][1] = 0, 1
        for event in r["split_events"]:
            event["wall"][1] = 2
        snap = boxes.reconstruct(r, 14)
        chosen, fixed = boxes.select_slice(snap, ["theta", "phi"], {"beta": .09, "chi": .6})
        self.assertEqual(len(chosen), 2)
        self.assertEqual(fixed["chi"], .6)
        r = result()
        for name in ("parameter_lower", "parameter_upper", "parameter_names", "parameter_units"):
            r.pop(name)
        r.update(attitude_coordinates="uv", beta_lower=.08, beta_upper=.12, tilt_max_deg=10)
        snap = boxes.reconstruct(r, 14)
        self.assertEqual(snap["parameter_names"], ["beta", "chi", "u", "v"])
        self.assertEqual(len(boxes.select_slice(snap, ["u", "v"], {"beta": .09})[0]), 2)

    def test_csv_cli_without_matplotlib(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "cfg/H1S_A005000/p00__a0.json"
            source.parent.mkdir(parents=True)
            source.write_text(json.dumps(result()))
            out = root / "analysis"
            self.assertEqual(boxes.main(["--json", str(source), "--time-days", "14",
                "--axes", "theta", "phi", "--slice", "beta=0.09",
                "--outdir", str(out), "--no-png"]), 0)
            csv_path, = out.rglob("*boxes.csv")
            with csv_path.open() as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(len(rows), 3)
            self.assertEqual(sum(int(r["in_slice"]) for r in rows), 2)
            self.assertIn("theta_lower_rad", rows[0])


class SummaryTests(unittest.TestCase):
    def test_reject_uv_output_under_angle_config(self):
        import contextlib
        import io
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            d = root / "Q2_L2_t10__theta_phi_t0_10_p-45_45/H1S_A005000"
            d.mkdir(parents=True)
            r = result()
            r.update(attitude_coordinates="uv", parameter_names=["beta", "chi", "u", "v"])
            (d / "p00__a0.json").write_text(json.dumps(r))
            with contextlib.redirect_stderr(io.StringIO()) as warnings:
                jobs = progress.scan_results(root)
            self.assertEqual(jobs, [])
            self.assertIn("mismatch", warnings.getvalue())

    def test_multiple_solar_phases_have_separate_summaries(self):
        import contextlib
        import io
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            d = root / "results/cfg/H1S_A005000"
            d.mkdir(parents=True)
            for solar in [0, 1]:
                r = result()
                r["alpha0_rad"] = solar
                (d / f"p00__a{solar}.json").write_text(json.dumps(r))
            with contextlib.redirect_stdout(io.StringIO()):
                status = progress.main(["--results-root", str(root), "--orbit-prefix", "H1S",
                                        "--outdir", str(root / "out"), "--no-png"])
            self.assertEqual(status, 0)
            with (root / "out/direction_summary_H1S.csv").open() as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(len(rows), 2)
            self.assertEqual({row["solar_phase_rad"] for row in rows}, {"0", "1"})

    def test_angle_counts_and_solar_phase_filter(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            d = root / "cfg/H1S_A005000"
            d.mkdir(parents=True)
            for solar in [0, 1]:
                r = result()
                r["alpha0_rad"] = solar
                (d / f"p00__a{solar}.json").write_text(json.dumps(r))
            jobs = progress.scan_results(root)
            selected = progress.select_jobs(jobs, "H1S", [], solar_phase=0)
            self.assertEqual(len(selected), 1)
            dest = root / "stats.csv"
            progress.write_stats_csv(selected, dest)
            with dest.open() as stream:
                row, = list(csv.DictReader(stream))
            self.assertAlmostEqual(float(row["first_split_time_days"]), 5)
            self.assertEqual(row["split_theta"], "1")
            self.assertEqual(row["split_u"], "")
            self.assertEqual(row["solar_phase_rad"], "0")


if __name__ == "__main__":
    unittest.main()
