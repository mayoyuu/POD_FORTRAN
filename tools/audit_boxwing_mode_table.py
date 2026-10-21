#!/usr/bin/env python3
"""Validate common-state P1 modes; publish only data passing the 30% gate."""
import argparse
import csv
from datetime import datetime, timezone
from functools import lru_cache
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "boxwing_srp_modes_v1"
HEADER = [
    "time_hours", "epoch_tdb_s",
    "ref_x_km", "ref_y_km", "ref_z_km", "ref_vx_km_s", "ref_vy_km_s", "ref_vz_km_s",
    "a_cb_x_km_s2", "a_cb_y_km_s2", "a_cb_z_km_s2", "a_cb_norm_km_s2",
    "a_bwe_x_km_s2", "a_bwe_y_km_s2", "a_bwe_z_km_s2",
    "a_bwm_x_km_s2", "a_bwm_y_km_s2", "a_bwm_z_km_s2",
    "eE_x", "eE_y", "eE_z", "eE_norm", "eM_x", "eM_y", "eM_z", "eM_norm",
]
PHYSICAL_KEYS = (
    "srp_cannonball_cr", "srp_cannonball_area_m2", "srp_mass_kg", "srp_box_dimensions_m",
    "srp_box_optical", "srp_array_total_area_m2", "srp_array_tracking_mode",
    "srp_array_hinge_axis_body", "srp_array_reference_normal_body",
    "srp_array_front_optical", "srp_array_back_optical", "srp_primary_axis_body",
    "srp_secondary_axis_body", "srp_roll_reference", "srp_pressure_1au_n_m2", "srp_geometry_tolerance",
)


class AuditError(ValueError):
    """A P1 dataset violates its numerical, identity or publication contract."""


def require(condition, message):
    if not condition:
        raise AuditError(message)


def file_sha256(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def cross(a, b):
    return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]


def norm(values):
    return math.hypot(*values)


def close(actual, expected, label, absolute=0.0):
    require(math.isfinite(actual) and math.isfinite(expected)
            and math.isclose(actual, expected, rel_tol=1e-12, abs_tol=absolute),
            f"{label}: got {actual!r}, expected {expected!r}")


def rtn_rows(position, velocity):
    radius = norm(position)
    angular = cross(position, velocity)
    hnorm = norm(angular)
    require(radius > 0 and hnorm > 1e-15*radius*norm(velocity), "reference state has a singular RTN frame")
    radial = [x/radius for x in position]
    normal = [x/hnorm for x in angular]
    return radial, cross(normal, radial), normal


def validate_physical_config(meta, config_file):
    config = {}
    for line in config_file.read_text(encoding="utf-8").splitlines():
        line = line.split("#", 1)[0].strip()
        if "=" in line:
            key, value = line.split("=", 1)
            config[key.strip()] = value.strip()
    physical = meta.get("physical_parameters", {})
    require(set(physical) == set(PHYSICAL_KEYS), "physical parameter keys do not match schema")
    for key in PHYSICAL_KEYS:
        require(key in config, f"source config is missing {key}")
        value = physical[key]
        try:
            if isinstance(value, str):
                require(value == config[key], f"source config disagrees with {key}")
            elif isinstance(value, list):
                fields = [float(x) for x in config[key].replace(",", " ").split()]
                require(len(fields) == len(value), f"invalid source config vector {key}")
                for actual, expected in zip(value, fields):
                    close(actual, expected, key)
            else:
                close(value, float(config[key]), key)
        except (TypeError, ValueError) as exc:
            raise AuditError(f"invalid physical parameter {key}: {exc}") from exc


@lru_cache(maxsize=16)
def spice_epoch_tdb(executable, executable_sha256, utc_epoch):
    """Use the actual project's SPICE conversion, without approximating UTC/TDB."""
    require(isinstance(utc_epoch, str) and bool(utc_epoch.strip()), "source OPM EPOCH must be a string")
    env = os.environ.copy()
    libraries = (ROOT/"lib", Path.home()/"local_dace/lib", Path.home()/"local_openblas/lib")
    env["LD_LIBRARY_PATH"] = ":".join(map(str, libraries)) + ":" + env.get("LD_LIBRARY_PATH", "")
    try:
        result = subprocess.run([executable, "--utc-to-tdb", utc_epoch], cwd=ROOT, env=env,
                                capture_output=True, text=True, timeout=30, check=True)
        values = [line.split("=", 1)[1] for line in result.stdout.splitlines()
                  if line.startswith("P1_EPOCH_TDB=")]
        require(len(values) == 1, "extractor did not return the SPICE OPM epoch")
        value = float(values[0])
        require(math.isfinite(value), "SPICE OPM epoch is nonfinite")
        return value
    except (OSError, subprocess.SubprocessError, ValueError) as exc:
        if isinstance(exc, AuditError):
            raise
        raise AuditError(f"cannot verify OPM epoch using SPICE: {exc}") from exc


def audit_table(csv_path, meta_path, expected_hours=360.0):
    """Return a recomputed audit; a valid dataset may still fail the physical gate."""
    csv_path, meta_path = Path(csv_path).resolve(), Path(meta_path).resolve()
    require(math.isfinite(expected_hours) and expected_hours > 0, "expected_hours must be positive")
    try:
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        expected = {
            "schema_version": SCHEMA, "frame": "J2000", "center": "EARTH", "timescale": "TDB",
            "illumination": "full_sun", "reference_model": "cannonball",
            "evaluation_contract": "common_reference_state",
            "units": {"position": "km", "velocity": "km/s", "acceleration": "km/s2", "modes": "dimensionless"},
        }
        for key, value in expected.items():
            require(meta.get(key) == value, f"invalid metadata {key}")
        require(type(meta.get("audit_complete")) is bool, "audit_complete must be boolean")
        epoch0, cadence = float(meta["epoch0_tdb_s"]), float(meta["sample_hours"])
        require(math.isfinite(epoch0) and math.isfinite(cadence) and cadence > 0, "invalid epoch or cadence")
        close(float(meta["duration_hours"]), expected_hours, "duration_hours", 1e-12)
        segments = round(expected_hours/cadence)
        close(segments*cadence, expected_hours, "duration must contain whole sampling intervals", 1e-12)
        count = segments + 1
        require(type(meta["row_count"]) is int and meta["row_count"] == count, "metadata row count is wrong")
        integrator = meta["integrator"]
        require(integrator.get("method") == "RKF78", "requires RKF78")
        for key, value in {"rel_tol": 1e-12, "abs_tol": 1e-12, "dt_min_s": 1e-6, "dt_max_s": 3600.0}.items():
            close(float(integrator[key]), value, f"integrator {key}")
        hashes = {}
        source_kinds = ["opm", "config", "executable"]
        for kind in source_kinds:
            path = Path(meta[f"{kind}_file"])
            if not path.is_absolute():
                path = meta_path.parent/path
            path = path.resolve()
            require(path.is_file(), f"missing source {kind}: {path}")
            hashes[f"{kind}_file"] = str(path)
            hashes[f"{kind}_sha256"] = file_sha256(path)
            if f"{kind}_sha256" in meta:
                require(hashes[f"{kind}_sha256"] == meta[f"{kind}_sha256"], f"{kind} SHA-256 mismatch")
        validate_physical_config(meta, Path(hashes["config_file"]))
        opm = json.loads(Path(hashes["opm_file"]).read_text(encoding="utf-8"))
        require("EPOCH" in opm, "source OPM must contain EPOCH")
        opm_epoch = spice_epoch_tdb(hashes["executable_file"], hashes["executable_sha256"], opm["EPOCH"])
        require(abs(epoch0-opm_epoch) <= 1e-6, "metadata TDB epoch differs from source OPM EPOCH")
        hashes["csv_sha256"] = file_sha256(csv_path)
        if "csv_sha256" in meta:
            require(hashes["csv_sha256"] == meta["csv_sha256"], "CSV SHA-256 mismatch")
        if meta["audit_complete"]:
            require(all(key in meta for key in ("csv_sha256", "opm_sha256", "config_sha256")),
                    "finalized metadata is missing SHA-256 identities")
        with csv_path.open(newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream)
            require(reader.fieldnames == HEADER, "CSV header does not match P1 schema")
            rows = list(reader)
        require(len(rows) == count, f"expected {count} CSV rows, found {len(rows)}")
        maxima = {"eE": (-1.0, 0.0), "eM": (-1.0, 0.0)}
        max_rtn_error = 0.0
        previous = -math.inf
        for index, raw in enumerate(rows):
            require(set(raw) == set(HEADER) and None not in raw.values(), f"row {index}: wrong field count")
            row = {key: float(value) for key, value in raw.items()}
            require(all(math.isfinite(value) for value in row.values()), f"row {index}: nonfinite value")
            hour = row["time_hours"]
            require(hour > previous, f"row {index}: time is not strictly increasing")
            close(hour, index*cadence, f"row {index} time_hours", 1e-12)
            require(abs(row["epoch_tdb_s"] - (epoch0+3600*hour)) <= 1e-6, f"row {index}: incorrect TDB epoch")
            previous = hour
            position = [row[f"ref_{axis}_km"] for axis in "xyz"]
            velocity = [row[f"ref_v{axis}_km_s"] for axis in "xyz"]
            if index == 0:
                for key, actual in zip(("X", "Y", "Z", "X_DOT", "Y_DOT", "Z_DOT"), position+velocity):
                    expected_initial = float(opm[key])
                    tolerance = 1e-9 if "DOT" not in key else 1e-12
                    require(math.isfinite(expected_initial) and abs(actual-expected_initial) <= tolerance,
                            f"reference initial {key} differs from source OPM")
            rotation = rtn_rows(position, velocity)
            cb = [row[f"a_cb_{axis}_km_s2"] for axis in "xyz"]
            cb_norm = norm(cb)
            require(cb_norm > 0, f"row {index}: zero Cannonball acceleration")
            close(row["a_cb_norm_km_s2"], cb_norm, f"row {index}: Cannonball norm")
            for mode, force in (("eE", "a_bwe"), ("eM", "a_bwm")):
                vector = [row[f"{mode}_{axis}"] for axis in "xyz"]
                for axis, baseline, value in zip("xyz", cb, vector):
                    close(value, (row[f"{force}_{axis}_km_s2"]-baseline)/cb_norm,
                          f"row {index}: {mode}_{axis}", 1e-14)
                length = norm(vector)
                close(row[f"{mode}_norm"], length, f"row {index}: {mode} norm", 1e-14)
                rotated = [math.fsum(a*b for a, b in zip(direction, vector)) for direction in rotation]
                error = abs(norm(rotated)-length)/length if length else 0.0
                require(error <= 1e-12, f"row {index}: RTN rotation changed {mode} norm")
                max_rtn_error = max(max_rtn_error, error)
                if length > maxima[mode][0]:
                    maxima[mode] = length, hour
        for mode, (value, hour) in maxima.items():
            close(float(meta[f"max_{mode}_norm"]), value, f"metadata {mode} maximum", 1e-14)
            close(float(meta[f"max_{mode}_time_hours"]), hour, f"metadata {mode} maximum time", 1e-9)
        passed = all(value <= 0.3 for value, _ in maxima.values())
        if meta["audit_complete"]:
            require(passed and meta.get("gate_passed") is True, "finalized dataset failed the physical gate")
        result = {
            **meta, **hashes, "row_count": count, "start_epoch_tdb_s": epoch0,
            "end_epoch_tdb_s": epoch0+3600*expected_hours, "opm_epoch_verified": True,
            "gate_limit": 0.3, "gate_passed": passed,
            "max_rtn_relative_error": max_rtn_error,
            "audit_time_utc": datetime.now(timezone.utc).isoformat(),
        }
        for mode, (value, hour) in maxima.items():
            result[f"max_{mode}_norm"] = value
            result[f"max_{mode}_time_hours"] = hour
        return result
    except (OSError, KeyError, TypeError, ValueError, OverflowError) as exc:
        if isinstance(exc, AuditError):
            raise
        raise AuditError(str(exc)) from exc


def finalize_table(csv_path, meta_path, output_dir, expected_hours=360.0):
    """Publish CSV and finalized metadata together, without replacing a prior run."""
    report = audit_table(csv_path, meta_path, expected_hours)
    require(report["gate_passed"], "30% physical gate failed; preserve staging files and do not start P2")
    output = Path(output_dir).absolute()
    require(not output.exists(), f"output directory already exists: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    lock = output.parent/f".{output.name}.publish.lock"
    try:
        lock_fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError as exc:
        raise AuditError(f"publication is locked: {lock}") from exc
    staging = None
    try:
        require(not output.exists(), f"output directory already exists: {output}")
        staging = Path(tempfile.mkdtemp(prefix=f".{output.name}.staging-", dir=output.parent))
        shutil.copyfile(csv_path, staging/"boxwing_modes.csv")
        require(file_sha256(staging/"boxwing_modes.csv") == report["csv_sha256"], "CSV changed during publication")
        source_executable = report["executable_file"]
        shutil.copy2(source_executable, staging/"extractor")
        require(file_sha256(staging/"extractor") == report["executable_sha256"],
                "extractor changed during publication")
        report.update(source_executable_file=source_executable, executable_file="extractor")
        revision = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"],
                                  capture_output=True, text=True, check=True).stdout.strip()
        report.update(audit_complete=True, source_revision=revision,
                      source_files_sha256={str(path.relative_to(ROOT)): file_sha256(path)
                                           for path in (ROOT/"app/run_boxwing_mode_extract.f90", Path(__file__))
                                           if path.is_file()})
        (staging/"boxwing_modes_meta.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n", encoding="utf-8")
        audit_table(staging/"boxwing_modes.csv", staging/"boxwing_modes_meta.json", expected_hours)
        require(not output.exists(), f"output directory already exists: {output}")
        staging.rename(output)
        staging = None
    finally:
        os.close(lock_fd)
        lock.unlink()
        if staging is not None:
            shutil.rmtree(staging)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", required=True, type=Path)
    parser.add_argument("--meta", required=True, type=Path)
    parser.add_argument("--expected-hours", type=float, default=360.0)
    parser.add_argument("--finalize", action="store_true")
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    if args.finalize and args.output_dir is None:
        parser.error("--finalize requires --output-dir")
    try:
        report = (finalize_table(args.csv, args.meta, args.output_dir, args.expected_hours) if args.finalize
                  else audit_table(args.csv, args.meta, args.expected_hours))
    except (AuditError, OSError, subprocess.SubprocessError) as exc:
        parser.exit(2, f"P1 audit failed: {exc}\n")
    print(json.dumps(report, indent=2, allow_nan=False))
    return 0 if report["gate_passed"] else 3


if __name__ == "__main__":
    raise SystemExit(main())
