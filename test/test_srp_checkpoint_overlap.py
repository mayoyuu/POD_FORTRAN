"""Regression checks at the checkpoint analysis CLI."""
import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
import analyze_srp_cov_containment as analysis


class CheckpointOverlap(unittest.TestCase):
    def test_overlap_denominators_and_time_matching(self):
        def cp(t, axes, x=0):
            return dict(t_nd=t, mean=[x, 0, 0, 0, 0, 0],
                        cov=np.diag([a*a for a in axes] + [1]*3).ravel().tolist())

        a = [cp(t, [1, 1, 1]) for t in range(1, 6)]
        b = [cp(1, [1, 1, 1], 1), cp(2, [0.5]*3), cp(3, [2]*3),
             cp(4, [0, 1, 1]), cp(6, [1]*3)]
        del b[2]["mean"]
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for name, series in (("q1", a), ("q2", b)):
                case = root / name / "H1N_A005000" / "p00__a0.json"
                case.parent.mkdir(parents=True)
                case.write_text(json.dumps(dict(checkpoints=series, status=0)))
            out = root / "analysis"
            self.assertEqual(analysis.main([
                "--inner-results", str(root / "q1"), "--outer-results", str(root / "q2"),
                "--orbit-prefix", "H1N", "--outdir", str(out), "--no-png",
                "--overlap-samples", "16384"]), 0)
            with (out / "cov_overlap_checkpoints_H1N.csv").open() as stream:
                rows = list(csv.DictReader(stream))
            with (out / "cov_containment_H1N.csv").open() as stream:
                summary, = csv.DictReader(stream)
        self.assertEqual((summary["n_volume"], summary["n_actual_means"]), ("3", "2"))
        self.assertAlmostEqual(float(rows[0]["centered_overlap_q1"]), 1)
        self.assertAlmostEqual(float(rows[0]["actual_means_overlap_q1"]), 5/16, delta=.02)
        self.assertAlmostEqual(float(rows[1]["volume_ratio_q2_q1"]), 1/8)
        self.assertAlmostEqual(float(rows[1]["centered_overlap_q1"]), 1/8)
        self.assertEqual(rows[2]["actual_means_overlap_q1"], "")
        self.assertEqual(rows[2]["centered_overlap_q1"], "1.0")
        self.assertIn("invalid_q2_cov", rows[3]["status"])
        self.assertEqual([r["status"] for r in rows[4:]], ["unmatched_time"]*2)


if __name__ == "__main__":
    unittest.main()
