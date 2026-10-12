#!/usr/bin/env python3
import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
import pack_srp_mc_initial_points_hdf5 as packer


HEADER = [
    "hour", "keyframe_kind", "sample_id",
    "x_km", "y_km", "z_km", "vx_kms", "vy_kms", "vz_kms",
    "dx_km", "dy_km", "dz_km", "dvx_kms", "dvy_kms", "dvz_kms",
]


class InitialPointHdf5PackerTest(unittest.TestCase):
    def make_csv(self, path, *, duplicate=False, missing=False, nonfinite=False):
        with Path(path).open("w", newline="") as stream:
            stream.write("schema=srp_mc_initial_points_v1\n")
            writer = csv.DictWriter(stream, fieldnames=HEADER)
            writer.writeheader()
            for hour in (0.0, 1.0):
                for sample_id in (2, 5, 9):
                    if missing and hour == 1.0 and sample_id == 9:
                        continue
                    state = np.arange(6, dtype=float) + 10.0 * hour + sample_id
                    error = state - np.arange(6, dtype=float) - 100.0
                    if nonfinite and hour == 1.0 and sample_id == 5:
                        error[3] = float("nan")
                    row = {
                        "hour": hour,
                        "keyframe_kind": "hourly_visualization",
                        "sample_id": sample_id,
                    }
                    for name, value in zip(HEADER[3:9], state):
                        row[name] = value
                    for name, value in zip(HEADER[9:15], error):
                        row[name] = value
                    writer.writerow(row)
                    if duplicate and hour == 0.0 and sample_id == 2:
                        writer.writerow(row)

    def test_pack_preserves_state_error_and_metadata(self):
        with tempfile.TemporaryDirectory() as folder:
            folder = Path(folder)
            source = folder / "points.csv"
            output = folder / "points.h5"
            report = folder / "report.json"
            self.make_csv(source)
            report.write_text(json.dumps({
                "output_schema_version": "srp_ads_mc_initial_only_v1",
                "srp_model": "box_wing",
                "srp_attitude_mode": "earth",
            }))
            result = packer.pack_initial_points(
                source, output, 3, "visualization", report
            )
            self.assertEqual(result["frame_count"], 2)
            self.assertEqual(result["point_count"], 3)
            with h5py.File(output, "r") as h5:
                self.assertEqual(
                    h5.attrs["output_schema_version"],
                    "srp_mc_initial_points_hdf5_v1",
                )
                self.assertEqual(h5.attrs["source_schema_version"],
                                 "srp_mc_initial_points_v1")
                self.assertEqual(h5.attrs["point_role"], "visualization")
                self.assertEqual(h5.attrs["srp_model"], "box_wing")
                self.assertEqual(h5.attrs["srp_attitude_mode"], "earth")
                self.assertIn("source_csv_sha256", h5.attrs)
                self.assertIn("report_sha256", h5.attrs)
                np.testing.assert_array_equal(h5["sample_id"][0], [2, 5, 9])
                self.assertEqual(h5["state"].shape, (2, 3, 6))
                self.assertEqual(h5["error"].shape, (2, 3, 6))
                np.testing.assert_allclose(h5["state"][1, 1], np.arange(6) + 15)
                np.testing.assert_allclose(h5["error"][1, 1], np.arange(6) + 15 - np.arange(6) - 100)
                np.testing.assert_allclose(h5["hour"][:], [0.0, 1.0])

    def assert_rejected_atomically(self, **kwargs):
        with tempfile.TemporaryDirectory() as folder:
            source = Path(folder) / "points.csv"
            output = Path(folder) / "points.h5"
            self.make_csv(source, **kwargs)
            with self.assertRaises(packer.PointFormatError):
                packer.pack_initial_points(source, output, 3, "visualization")
            self.assertFalse(output.exists())

    def test_missing_point_is_rejected_atomically(self):
        self.assert_rejected_atomically(missing=True)

    def test_duplicate_point_is_rejected_atomically(self):
        self.assert_rejected_atomically(duplicate=True)

    def test_nonfinite_value_is_rejected_atomically(self):
        self.assert_rejected_atomically(nonfinite=True)

    def test_changed_id_order_is_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            source = Path(folder) / "points.csv"
            output = Path(folder) / "points.h5"
            self.make_csv(source)
            lines = source.read_text().splitlines()
            lines[-2], lines[-1] = lines[-1], lines[-2]
            source.write_text("\n".join(lines) + "\n")
            with self.assertRaises(packer.PointFormatError):
                packer.pack_initial_points(source, output, 3, "visualization")
            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
