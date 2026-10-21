"""Algebraic fixtures exercise the P1 data contract, not a force-model mock."""
import csv
import importlib.util
import json
import math
from pathlib import Path
import shutil
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
EPOCH0 = 852033729.1839125
SCRIPT = ROOT / "tools" / "audit_boxwing_mode_table.py"
HEADER = [
    "time_hours", "epoch_tdb_s",
    "ref_x_km", "ref_y_km", "ref_z_km", "ref_vx_km_s", "ref_vy_km_s", "ref_vz_km_s",
    "a_cb_x_km_s2", "a_cb_y_km_s2", "a_cb_z_km_s2", "a_cb_norm_km_s2",
    "a_bwe_x_km_s2", "a_bwe_y_km_s2", "a_bwe_z_km_s2",
    "a_bwm_x_km_s2", "a_bwm_y_km_s2", "a_bwm_z_km_s2",
    "eE_x", "eE_y", "eE_z", "eE_norm", "eM_x", "eM_y", "eM_z", "eM_norm",
]
PHYSICAL = {
    "srp_cannonball_cr": 1.2213, "srp_cannonball_area_m2": 30.0, "srp_mass_kg": 1200.0,
    "srp_box_dimensions_m": [2.0, 3.0, 4.0], "srp_box_optical": [0.3, 0.4, 0.3],
    "srp_array_total_area_m2": 24.0, "srp_array_tracking_mode": "single_axis",
    "srp_array_hinge_axis_body": [0.0, 1.0, 0.0],
    "srp_array_reference_normal_body": [1.0, 0.0, 0.0],
    "srp_array_front_optical": [0.85, 0.08, 0.07], "srp_array_back_optical": [0.6, 0.2, 0.2],
    "srp_primary_axis_body": [0.0, 0.0, 1.0], "srp_secondary_axis_body": [0.0, 1.0, 0.0],
    "srp_roll_reference": "orbit_normal", "srp_pressure_1au_n_m2": 1367.0 / 299792458.0,
    "srp_geometry_tolerance": 1e-12,
}


def load_auditor():
    spec = importlib.util.spec_from_file_location("audit_boxwing_modes", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_fixture(folder, cadence=1.0, earth=0.2, hours=360.0):
    opm = folder / "input.opm.json"
    cfg = folder / "config.txt"
    opm.write_text(json.dumps({"EPOCH": "2027-01-01T00:01:00", "X": 300000.0,
                              "Y": 1000.0, "Z": 20000.0, "X_DOT": 0.1,
                              "Y_DOT": 0.5, "Z_DOT": 0.2}) + "\n")
    cfg.write_text("\n".join(
        f"{key} = {' '.join(map(str, value)) if isinstance(value, list) else value}"
        for key, value in PHYSICAL.items()
    ) + "\n")
    cb = [1e-10, -2e-10, 3e-10]
    norm = math.sqrt(sum(x * x for x in cb))
    rows = []
    for i in range(round(hours / cadence) + 1):
        hour = i * cadence
        values = [hour, EPOCH0 + hour * 3600.0, 300000.0, 1000.0, 20000.0, 0.1, 0.5, 0.2]
        values += cb + [norm]
        values += [cb[0] + earth * norm, cb[1], cb[2]]
        values += [cb[0], cb[1] + 0.1 * norm, cb[2]]
        values += [earth, 0.0, 0.0, earth, 0.0, 0.1, 0.0, 0.1]
        rows.append(dict(zip(HEADER, values)))
    path = folder / "candidate.csv"
    write_rows(path, rows)
    meta = {
        "schema_version": "boxwing_srp_modes_v1", "audit_complete": False,
        "frame": "J2000", "center": "EARTH", "timescale": "TDB", "illumination": "full_sun",
        "reference_model": "cannonball", "evaluation_contract": "common_reference_state",
        "units": {"position": "km", "velocity": "km/s", "acceleration": "km/s2", "modes": "dimensionless"},
        "opm_file": str(opm), "config_file": str(cfg), "epoch0_tdb_s": EPOCH0,
        "executable_file": str(next(ROOT.glob("build/*/app/run_boxwing_mode_extract"))),
        "duration_hours": hours, "sample_hours": cadence, "row_count": len(rows),
        "physical_parameters": PHYSICAL,
        "integrator": {"method": "RKF78", "rel_tol": 1e-12, "abs_tol": 1e-12,
                       "dt_min_s": 1e-6, "dt_max_s": 3600.0},
        "max_eE_norm": earth, "max_eE_time_hours": 0.0,
        "max_eM_norm": 0.1, "max_eM_time_hours": 0.0,
    }
    metadata = folder / "candidate_meta.json"
    metadata.write_text(json.dumps(meta))
    return path, metadata, rows, meta


def write_rows(path, rows, header=HEADER):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=header)
        writer.writeheader()
        writer.writerows(rows)


class AuditorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if SCRIPT.is_file():
            cls.audit = load_auditor()

    def setUp(self):
        self.assertTrue(SCRIPT.is_file(), "P1 auditor has not been implemented")
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)
        self.csv, self.meta, self.rows, self.metadata = write_fixture(self.folder)

    def test_recomputes_gate_and_norm_preserving_rtn(self):
        result = self.audit.audit_table(self.csv, self.meta)
        self.assertTrue(result["gate_passed"])
        self.assertEqual(result["row_count"], 361)
        self.assertAlmostEqual(result["max_eE_norm"], 0.2)
        self.assertLessEqual(result["max_rtn_relative_error"], 1e-12)

    def test_half_hour_cadence_has_721_rows(self):
        self.csv, self.meta, _, _ = write_fixture(self.folder, cadence=0.5)
        self.assertEqual(self.audit.audit_table(self.csv, self.meta)["row_count"], 721)

    def test_bad_csv_values_are_rejected(self):
        mutations = {"epoch_tdb_s": 0.0, "time_hours": 0.5, "eE_x": 0.3,
                     "a_cb_norm_km_s2": 1e-10, "eM_norm": 0.2, "ref_x_km": float("nan")}
        for key, value in mutations.items():
            with self.subTest(key=key):
                rows = [row.copy() for row in self.rows]
                rows[1][key] = value
                write_rows(self.csv, rows)
                with self.assertRaises(self.audit.AuditError):
                    self.audit.audit_table(self.csv, self.meta)

    def test_bad_header_and_row_count_are_rejected(self):
        self.csv.write_text(self.csv.read_text().replace("eE_x", "wrong", 1))
        with self.assertRaises(self.audit.AuditError):
            self.audit.audit_table(self.csv, self.meta)
        write_rows(self.csv, self.rows[:-1])
        with self.assertRaises(self.audit.AuditError):
            self.audit.audit_table(self.csv, self.meta)

    def test_bad_metadata_identity_and_maximum_are_rejected(self):
        for key, value in {"frame": "RTN", "timescale": "UTC", "max_eE_norm": 0.25,
                           "evaluation_contract": "independent_orbits", "sample_hours": 0.5}.items():
            with self.subTest(key=key):
                meta = dict(self.metadata, **{key: value})
                self.meta.write_text(json.dumps(meta))
                with self.assertRaises(self.audit.AuditError):
                    self.audit.audit_table(self.csv, self.meta)

    def test_physical_config_mismatch_is_rejected(self):
        (self.folder / "config.txt").write_text("srp_mass_kg = 1000\n")
        with self.assertRaises(self.audit.AuditError):
            self.audit.audit_table(self.csv, self.meta)

    def test_singular_rtn_state_is_rejected(self):
        for name in ("ref_vx_km_s", "ref_vy_km_s", "ref_vz_km_s"):
            self.rows[0][name] = 0.0
        write_rows(self.csv, self.rows)
        with self.assertRaises(self.audit.AuditError):
            self.audit.audit_table(self.csv, self.meta)

    def test_reference_initial_state_must_match_the_source_opm(self):
        self.rows[0]["ref_x_km"] += 1.0
        write_rows(self.csv, self.rows)
        with self.assertRaises(self.audit.AuditError):
            self.audit.audit_table(self.csv, self.meta)

    def test_uniform_epoch_shift_is_rejected_against_opm_epoch(self):
        for row in self.rows:
            row["epoch_tdb_s"] += 86400.0
        self.metadata["epoch0_tdb_s"] += 86400.0
        write_rows(self.csv, self.rows)
        self.meta.write_text(json.dumps(self.metadata))
        with self.assertRaises(self.audit.AuditError):
            self.audit.audit_table(self.csv, self.meta)

    def test_failed_physical_gate_is_diagnostic_and_cannot_publish(self):
        self.csv, self.meta, _, _ = write_fixture(self.folder, earth=0.31)
        result = self.audit.audit_table(self.csv, self.meta)
        self.assertFalse(result["gate_passed"])
        output = self.folder / "published"
        with self.assertRaises(self.audit.AuditError):
            self.audit.finalize_table(self.csv, self.meta, output)
        self.assertFalse(output.exists())
        self.assertFalse(json.loads(self.meta.read_text())["audit_complete"])

    def test_finalize_binds_hashes_and_never_overwrites(self):
        output = self.folder / "published"
        self.audit.finalize_table(self.csv, self.meta, output)
        csv_file = output / "boxwing_modes.csv"
        meta_file = output / "boxwing_modes_meta.json"
        finalized = json.loads(meta_file.read_text())
        self.assertTrue(finalized["audit_complete"])
        self.assertTrue(finalized["gate_passed"])
        self.assertEqual(finalized["csv_sha256"], self.audit.file_sha256(csv_file))
        self.assertTrue(self.audit.audit_table(csv_file, meta_file)["gate_passed"])
        with self.assertRaises(self.audit.AuditError):
            self.audit.finalize_table(self.csv, self.meta, output)

    def test_finalized_input_and_csv_tampering_are_rejected(self):
        output = self.folder / "published"
        self.audit.finalize_table(self.csv, self.meta, output)
        csv_file = output / "boxwing_modes.csv"
        meta_file = output / "boxwing_modes_meta.json"
        original = csv_file.read_bytes()
        csv_file.write_bytes(original + b"\n")
        with self.assertRaises(self.audit.AuditError):
            self.audit.audit_table(csv_file, meta_file)
        csv_file.write_bytes(original)
        (self.folder / "input.opm.json").write_text('{"EPOCH":"2028-01-01"}')
        with self.assertRaises(self.audit.AuditError):
            self.audit.audit_table(csv_file, meta_file)

    def test_runtime_executable_hash_is_bound_and_checked(self):
        executable = self.folder / "extractor"
        shutil.copyfile(self.metadata["executable_file"], executable)
        executable.chmod(0o755)
        self.metadata["executable_file"] = str(executable)
        self.meta.write_text(json.dumps(self.metadata))
        output = self.folder / "published"
        report = self.audit.finalize_table(self.csv, self.meta, output)
        self.assertEqual(report.get("executable_sha256"), self.audit.file_sha256(executable))
        self.assertEqual(report["executable_file"], "extractor")
        self.assertTrue((output / "extractor").is_file())
        # Subsequent builds must not invalidate an already published mode table.
        executable.write_bytes(b"changed executable")
        self.assertTrue(self.audit.audit_table(output / "boxwing_modes.csv", output / "boxwing_modes_meta.json")["gate_passed"])
        (output / "extractor").write_bytes(b"changed frozen executable")
        with self.assertRaises(self.audit.AuditError):
            self.audit.audit_table(output / "boxwing_modes.csv", output / "boxwing_modes_meta.json")


if __name__ == "__main__":
    unittest.main()
