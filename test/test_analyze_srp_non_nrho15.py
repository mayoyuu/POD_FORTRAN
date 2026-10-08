import importlib.util
import unittest
from pathlib import Path

import pandas as pd

SCRIPT = Path(__file__).resolve().parents[1] / "tools" / "analyze_srp_non_nrho15.py"
SPEC = importlib.util.spec_from_file_location("analyze_srp_non_nrho15", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class PairedRelativeMaximumTests(unittest.TestCase):
    def test_paired_delta_is_shared_numerator_and_group_axis_is_denominator(self):
        treatment = pd.DataFrame({
            "hour": [0.0, 1.0],
            "pos_longest_rms_axis_km": [2.0, 4.0],
            "vel_longest_rms_axis_km_s": [0.5, 1.0],
        })
        control = pd.DataFrame({
            "hour": [0.0, 1.0],
            "pos_longest_rms_axis_km": [1.0, 8.0],
            "vel_longest_rms_axis_km_s": [0.25, 2.0],
        })
        delta = pd.DataFrame({
            "hour": [0.0, 1.0],
            "delta_r_max_km": [10.0, 20.0],
            "delta_v_max_km_s": [1.0, 4.0],
        })

        treatment_result, control_result = MODULE.attach_paired_relative_max(
            treatment, control, delta
        )

        self.assertEqual(treatment_result["rel_max_pos_num_km"].tolist(), [10.0, 20.0])
        self.assertEqual(control_result["rel_max_pos_num_km"].tolist(), [10.0, 20.0])
        self.assertEqual(treatment_result["rel_max_pos"].tolist(), [5.0, 5.0])
        self.assertEqual(control_result["rel_max_pos"].tolist(), [10.0, 2.5])
        self.assertEqual(treatment_result["rel_max_vel"].tolist(), [2.0, 4.0])
        self.assertEqual(control_result["rel_max_vel"].tolist(), [4.0, 2.0])


if __name__ == "__main__":
    unittest.main()
