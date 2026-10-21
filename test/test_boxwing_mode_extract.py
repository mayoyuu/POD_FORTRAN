import csv
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

from test_audit_boxwing_mode_table import load_auditor


ROOT = Path(__file__).resolve().parents[1]


def read_rows(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


class ExtractorTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)
        self.prefix = self.folder / "candidate"

    def run_extract(self, *extra):
        env = os.environ.copy()
        libraries = [ROOT / "lib", Path.home() / "local_dace/lib", Path.home() / "local_openblas/lib"]
        for key in ("LIBRARY_PATH", "LD_LIBRARY_PATH"):
            env[key] = ":".join(map(str, libraries)) + ":" + env.get(key, "")
        return subprocess.run(
            ["fpm", "run", "--target", "run_boxwing_mode_extract", "--", *map(str, extra)],
            cwd=ROOT, env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        )

    def arguments(self, hours=1, cadence=0.5, cfg="config_box_wing_earth.txt"):
        return ["-opm", ROOT / "OPM/DRO-2/DRO-2_init.opm.json", "-cfg", ROOT / "config" / cfg,
                "-o", self.prefix, "--hours", hours, "--sample-hours", cadence]

    def test_help_and_invalid_cadence(self):
        result = self.run_extract("--help")
        self.assertEqual(result.returncode, 0, result.stdout[-3000:])
        self.assertIn("sample-hours", result.stdout)
        for value in (0, -1, "NaN", 0.3):
            result = self.run_extract(*self.arguments(cadence=value))
            self.assertNotEqual(result.returncode, 0)
            self.assertFalse(self.prefix.with_suffix(".csv").exists())

    def test_spice_epoch_inspector(self):
        result = self.run_extract("--utc-to-tdb", "2027-01-01T00:01:00")
        self.assertEqual(result.returncode, 0, result.stdout[-3000:])
        value = next(line.split("=", 1)[1] for line in result.stdout.splitlines()
                     if line.startswith("P1_EPOCH_TDB="))
        self.assertAlmostEqual(float(value), 852033729.1839125, delta=1e-6)

    def test_common_reference_smoke_and_no_overwrite(self):
        result = self.run_extract(*self.arguments())
        self.assertEqual(result.returncode, 0, result.stdout[-4000:])
        csv_file = self.prefix.with_suffix(".csv")
        meta_file = Path(str(self.prefix) + "_meta.json")
        self.assertTrue(csv_file.is_file() and meta_file.is_file(), result.stdout)
        metadata = json.loads(meta_file.read_text())
        self.assertFalse(metadata["audit_complete"])
        self.assertEqual(metadata["reference_model"], "cannonball")
        rows = read_rows(csv_file)
        self.assertEqual([float(row["time_hours"]) for row in rows], [0, 0.5, 1])
        opm = json.loads((ROOT / "OPM/DRO-2/DRO-2_init.opm.json").read_text())
        for name, key in zip(("ref_x_km", "ref_y_km", "ref_z_km", "ref_vx_km_s", "ref_vy_km_s", "ref_vz_km_s"),
                             ("X", "Y", "Z", "X_DOT", "Y_DOT", "Z_DOT")):
            self.assertAlmostEqual(float(rows[0][name]), opm[key], delta=1e-10)
        audit = load_auditor().audit_table(csv_file, meta_file, expected_hours=1)
        self.assertTrue(audit["gate_passed"])
        self.assertAlmostEqual(float(rows[0]["eE_norm"]), 0.248269417459696, delta=1e-12)
        self.assertAlmostEqual(float(rows[0]["eM_norm"]), 0.161447268471512, delta=1e-12)
        # Independent pre-P1 Cannonball reference at this OPM epoch, not a
        # second execution of the extractor. Detects accidental BW integration.
        reference = [-299090.27439841500, -171951.86271499429, -108274.48865325656,
                     0.30935104398452318, -1.0383452259538237, -0.51963322762264963]
        for name, expected in zip(("ref_x_km", "ref_y_km", "ref_z_km",
                                   "ref_vx_km_s", "ref_vy_km_s", "ref_vz_km_s"), reference):
            self.assertAlmostEqual(float(rows[-1][name]), expected, delta=1e-7)
        for name, expected in zip(("a_cb_x_km_s2", "a_cb_y_km_s2", "a_cb_z_km_s2"),
                                   (-2.5236441738494006e-11, 1.3032706830717093e-10, 5.6461693709338575e-11)):
            # The historical experiment used 1367/c; the complete production
            # config stores a slightly different rounded pressure. SRP is linear in P.
            pressure_ratio = metadata["physical_parameters"]["srp_pressure_1au_n_m2"] / (1367/299792458)
            self.assertAlmostEqual(float(rows[0][name]), expected*pressure_ratio, delta=1e-22)
        original = csv_file.read_bytes()
        repeated = self.run_extract(*self.arguments())
        self.assertNotEqual(repeated.returncode, 0)
        self.assertEqual(csv_file.read_bytes(), original)

    def test_sampling_cadence_and_input_model_do_not_change_reference(self):
        a = self.run_extract(*self.arguments(cadence=1))
        self.assertEqual(a.returncode, 0, a.stdout[-3000:])
        self.assertTrue(self.prefix.with_suffix(".csv").is_file(), a.stdout)
        baseline = read_rows(self.prefix.with_suffix(".csv"))
        self.prefix = self.folder / "earth_input"
        b = self.run_extract(*self.arguments(cadence=0.5, cfg="config_box_wing_moon.txt"))
        self.assertEqual(b.returncode, 0, b.stdout[-3000:])
        other = read_rows(self.prefix.with_suffix(".csv"))
        for key in ("ref_x_km", "ref_y_km", "ref_z_km", "ref_vx_km_s", "ref_vy_km_s", "ref_vz_km_s"):
            self.assertAlmostEqual(float(baseline[-1][key]), float(other[-1][key]), delta=1e-7)
        self.assertEqual(json.loads(Path(str(self.prefix) + "_meta.json").read_text())["sample_hours"], 0.5)


if __name__ == "__main__":
    unittest.main()
