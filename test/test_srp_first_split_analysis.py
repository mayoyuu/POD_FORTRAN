"""Tests for selecting event-detected first-split times in analysis tools."""
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))

import analyze_srp_beta_scan_summary as beta_summary
import analyze_srp_split_vs_axis as split_axis
import analyze_srp_study_progress as progress
from srp_first_split_time import select_first_split_time


ND_TO_DAYS = progress.ND_TO_DAYS


def covariance_state():
    cov = [0.0] * 36
    cov[0] = cov[7] = cov[14] = 1.0
    return cov


class FirstSplitSelectionTests(unittest.TestCase):
    def test_event_time_has_priority_over_checkpoint_time(self):
        value, source = select_first_split_time({
            "first_split_time": 2.0,
            "first_split_time_accurate": 1.25,
        })
        self.assertEqual(value, 1.25)
        self.assertEqual(source, "event")

    def test_invalid_event_time_falls_back_to_checkpoint(self):
        value, source = select_first_split_time({
            "first_split_time": 2.0,
            "first_split_time_accurate": -1.0,
        })
        self.assertEqual(value, 2.0)
        self.assertEqual(source, "checkpoint")

    def test_no_split_and_missing_values_are_distinguished(self):
        self.assertEqual(select_first_split_time({
            "first_split_time": -1.0,
            "first_split_time_accurate": -1.0,
        }), (None, "never"))
        self.assertEqual(select_first_split_time({}), (None, "missing"))

    def test_progress_uses_event_time(self):
        job = {"first_split_time": 2.0, "first_split_time_accurate": 1.25}
        self.assertEqual(progress.first_split_days(job), 1.25 * ND_TO_DAYS)

    def test_progress_grid_preserves_never_split_sentinel(self):
        job = {
            "orbit": "H1N_A005000",
            "phase_index": 0,
            "first_split_time": -1.0,
            "first_split_time_accurate": -1.0,
        }
        self.assertEqual(progress.build_first_split_grid([job])[5000][0], -1.0)

    def test_beta_summary_uses_event_time(self):
        row = beta_summary.summarize_group([{
            "status": 0,
            "cov_state": covariance_state(),
            "amplitude": 5000,
            "phase_index": 0,
            "first_split_time": 2.0,
            "first_split_time_accurate": 1.25,
        }])
        self.assertEqual(row["split_min"], 1.25 * ND_TO_DAYS)

    def test_split_axis_uses_event_time_and_omits_never_split(self):
        base = {
            "status": 0,
            "cov_state": covariance_state(),
            "config": "cfg",
            "family": "H1N",
            "orbit": "H1N_A005000",
            "amplitude": 5000,
            "phase_index": 0,
        }
        event = dict(base, first_split_time=2.0,
                     first_split_time_accurate=1.25)
        never = dict(base, phase_index=1, first_split_time=-1.0,
                     first_split_time_accurate=-1.0)
        points, _, _ = split_axis.valid_points([event, never])
        self.assertEqual(len(points), 1)
        self.assertEqual(points[0][0], 1.25 * ND_TO_DAYS)


if __name__ == "__main__":
    unittest.main()
