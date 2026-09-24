#!/usr/bin/env python3
"""Run paired SRP-ADS histories for matching initial OPM files in parallel."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
import time
import uuid
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

OUTPUT_SCHEMA_VERSION = 2
POS_KEYS = ("CX_X", "CY_Y", "CZ_Z")
VEL_KEYS = ("CX_DOT_X_DOT", "CY_DOT_Y_DOT", "CZ_DOT_Z_DOT")
SPLIT_EVENT_COLUMNS = {
    "interval_start_hour",
    "interval_end_hour",
    "parent_depth",
    "component",
    "direction",
    "direction_score7",
}


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def find_initial_opms(root: Path) -> list[Path]:
    """Only accept init OPMs directly below one orbit directory."""
    paths = set(root.glob("*/*_init.opm"))
    paths.update(root.glob("*/*_init.opm.json"))
    return sorted(path for path in paths if path.is_file())


def covariance_rss(path: Path) -> tuple[float, float]:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"cannot parse JSON OPM: {exc}") from exc
    values = []
    for key in POS_KEYS + VEL_KEYS:
        if key not in data:
            raise ValueError(f"missing top-level covariance element {key}")
        try:
            value = float(data[key])
        except (TypeError, ValueError) as exc:
            raise ValueError(f"non-numeric covariance element {key}") from exc
        if not math.isfinite(value) or value < 0.0:
            raise ValueError(f"invalid covariance variance {key}={value}")
        values.append(value)
    return math.sqrt(sum(values[:3])), 1000.0 * math.sqrt(sum(values[3:]))


def matches_expected_error(
    path: Path, pos_km: float, vel_mps: float, relative_tolerance: float = 1e-9
) -> bool:
    try:
        actual_pos, actual_vel = covariance_rss(path)
    except ValueError:
        return False
    return math.isclose(actual_pos, pos_km, rel_tol=relative_tolerance, abs_tol=1e-12) and (
        math.isclose(actual_vel, vel_mps, rel_tol=relative_tolerance, abs_tol=1e-12)
    )


def output_prefix_for(path: Path, opm_root: Path, output_root: Path) -> Path:
    relative = path.relative_to(opm_root)
    name = relative.name
    for suffix in (".opm.json", ".opm"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    return output_root / relative.parent / name


def mirror_directories(opm_root: Path, output_root: Path) -> None:
    output_root.mkdir(parents=True, exist_ok=True)
    for opm_path in find_initial_opms(opm_root):
        output_prefix_for(opm_path, opm_root, output_root).parent.mkdir(
            parents=True, exist_ok=True
        )


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    with temp.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    temp.replace(path)


def write_json(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    temp.replace(path)


def preflight(
    opm_root: Path,
    output_root: Path,
    expected_pos: float,
    expected_vel: float,
    relative_tolerance: float,
) -> tuple[list[dict], list[dict]]:
    rows, eligible = [], []
    for path in find_initial_opms(opm_root):
        prefix = output_prefix_for(path, opm_root, output_root)
        row = {
            "opm": str(path),
            "relative_opm": str(path.relative_to(opm_root)),
            "output_prefix": str(prefix),
            "position_rss_km": "",
            "velocity_rss_mps": "",
            "status": "",
            "reason": "",
        }
        try:
            pos, vel = covariance_rss(path)
            row["position_rss_km"] = f"{pos:.16g}"
            row["velocity_rss_mps"] = f"{vel:.16g}"
            if math.isclose(pos, expected_pos, rel_tol=relative_tolerance, abs_tol=1e-12) and (
                math.isclose(vel, expected_vel, rel_tol=relative_tolerance, abs_tol=1e-12)
            ):
                row["status"] = "eligible"
                eligible.append(row)
            else:
                row["status"] = "initial_error_mismatch"
                row["reason"] = f"expected RSS {expected_pos:g} km and {expected_vel:g} m/s"
        except ValueError as exc:
            row["status"] = "invalid_opm"
            row["reason"] = str(exc)
        rows.append(row)
    return rows, eligible


def find_executable(project_root: Path) -> Path:
    candidates = [
        path
        for path in (project_root / "build").glob("*/app/run_srp_ads_history")
        if path.is_file() and os.access(path, os.X_OK)
    ]
    if not candidates:
        raise RuntimeError("cannot locate build/*/app/run_srp_ads_history")
    return max(candidates, key=lambda path: path.stat().st_mtime)


def _close_number(actual: object, expected: float) -> bool:
    try:
        return math.isclose(float(actual), expected, rel_tol=1e-12, abs_tol=1e-15)
    except (TypeError, ValueError):
        return False


def audit_report(
    prefix: Path,
    args: argparse.Namespace,
    opm_path: Path,
    config_path: Path,
) -> tuple[str, list[str]]:
    report_path = Path(str(prefix) + "_report.json")
    try:
        report = json.loads(report_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return "missing_result", ["report_missing_or_invalid"]

    reasons: list[str] = []
    expected_ints = {
        "output_schema_version": OUTPUT_SCHEMA_VERSION,
        "n_variables": 6 if args.eta_sigma == 0.0 else 7,
        "da_order": args.da_order,
        "max_split_depth": args.max_depth,
    }
    for key, value in expected_ints.items():
        try:
            matches = int(report[key]) == value
        except (KeyError, TypeError, ValueError):
            matches = False
        if not matches:
            reasons.append(key)

    expected_numbers = {
        "duration_hours": args.days * 24.0,
        "save_hours": args.save_hours,
        "eta_sigma": args.eta_sigma,
        "position_tolerance_km": args.pos_tol_km,
        "velocity_tolerance_kms": args.vel_tol_kms,
    }
    for key, value in expected_numbers.items():
        if key not in report or not _close_number(report[key], value):
            reasons.append(key)

    expected_bools = {
        "nominal_srp_enabled": True,
        "srp_uncertainty_enabled": args.eta_sigma > 0.0,
    }
    for key, value in expected_bools.items():
        if report.get(key) is not value:
            reasons.append(key)

    expected_hashes = {
        "opm_sha256": file_sha256(opm_path),
        "config_sha256": file_sha256(config_path),
    }
    for key, value in expected_hashes.items():
        if report.get(key) != value:
            reasons.append(key)

    split_path = Path(str(prefix) + "_split_events.csv")
    try:
        with split_path.open(newline="", encoding="utf-8") as handle:
            header = set(next(csv.reader(handle)))
        if not SPLIT_EVENT_COLUMNS.issubset(header):
            reasons.append("split_events_columns")
    except (OSError, StopIteration):
        reasons.append("split_events_file")

    return ("compatible", []) if not reasons else ("incompatible", reasons)


def report_matches(
    prefix: Path,
    args: argparse.Namespace,
    opm_path: Path | None = None,
    config_path: Path | None = None,
) -> bool:
    if opm_path is None or config_path is None:
        return False
    return audit_report(prefix, args, opm_path, config_path)[0] == "compatible"


def inject_provenance(
    prefix: Path,
    opm_path: Path,
    config_path: Path,
    executable: Path,
    source_revision: str,
) -> None:
    report_path = Path(str(prefix) + "_report.json")
    report = json.loads(report_path.read_text(encoding="utf-8"))
    report.update(
        {
            "output_schema_version": OUTPUT_SCHEMA_VERSION,
            "source_opm": str(opm_path.resolve()),
            "opm_sha256": file_sha256(opm_path),
            "source_config": str(config_path.resolve()),
            "config_sha256": file_sha256(config_path),
            "executable_sha256": file_sha256(executable),
            "source_revision": source_revision,
        }
    )
    write_json(report_path, report)


def promote_prefix(temp_prefix: Path, final_prefix: Path) -> int:
    temp_files = sorted(
        path
        for path in temp_prefix.parent.glob(temp_prefix.name + "_*")
        if path.is_file()
    )
    if not temp_files:
        raise RuntimeError(f"temporary run produced no artifacts: {temp_prefix}")
    old_files = {
        path
        for path in final_prefix.parent.glob(final_prefix.name + "_*")
        if path.is_file()
    }
    targets: set[Path] = set()
    for source in temp_files:
        suffix = source.name[len(temp_prefix.name) :]
        target = final_prefix.parent / (final_prefix.name + suffix)
        source.replace(target)
        targets.add(target)
    for stale in sorted(old_files - targets):
        stale.unlink()
    return len(targets)


def stage_executable(executable: Path, output_root: Path) -> Path:
    digest = file_sha256(executable)
    runtime = output_root / ".batch_runtime"
    runtime.mkdir(parents=True, exist_ok=True)
    target = runtime / f"run_srp_ads_history_{digest[:16]}"
    if not target.exists() or file_sha256(target) != digest:
        temp = target.with_suffix(".tmp")
        shutil.copy2(executable, temp)
        temp.chmod(temp.stat().st_mode | 0o111)
        temp.replace(target)
    return target


def source_revision(project_root: Path) -> str:
    completed = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=project_root,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        check=False,
    )
    return completed.stdout.strip() if completed.returncode == 0 else "unknown"


def run_one(
    row: dict,
    executable: Path,
    project_root: Path,
    config: Path,
    args: argparse.Namespace,
    environment: dict[str, str],
    revision: str = "unknown",
) -> dict:
    prefix = Path(row["output_prefix"])
    opm_path = Path(row["opm"])
    prefix.parent.mkdir(parents=True, exist_ok=True)
    audit_status, audit_reasons = audit_report(prefix, args, opm_path, config)
    final_log = Path(str(prefix) + "_run.log")

    if audit_status == "compatible" and not args.force:
        return {
            **row,
            "audit_status": audit_status,
            "audit_reasons": "",
            "run_status": "skipped_complete",
            "return_code": 0,
            "elapsed_seconds": 0.0,
            "log": str(final_log),
        }
    if audit_status == "incompatible" and not (args.refresh_incompatible or args.force):
        return {
            **row,
            "audit_status": audit_status,
            "audit_reasons": ";".join(audit_reasons),
            "run_status": "conflict_existing_report",
            "return_code": 2,
            "elapsed_seconds": 0.0,
            "log": str(final_log),
        }

    temp_prefix = prefix.with_name(prefix.name + f".__new__{uuid.uuid4().hex}")
    temp_report = Path(str(temp_prefix) + "_report.json")
    temp_log = Path(str(temp_prefix) + "_run.log")
    command = [
        str(executable),
        "-opm", str(opm_path.resolve()),
        "-o", str(temp_prefix.resolve()),
        "--eta-sigma", str(args.eta_sigma),
        "--days", str(args.days),
        "--save-hours", str(args.save_hours),
        "--da-order", str(args.da_order),
        "--max-depth", str(args.max_depth),
        "--pos-tol-km", str(args.pos_tol_km),
        "--vel-tol-kms", str(args.vel_tol_kms),
        "-cfg", str(config),
    ]
    started = time.monotonic()
    with temp_log.open("w", encoding="utf-8") as handle:
        handle.write("COMMAND: " + " ".join(command) + "\n")
        handle.flush()
        completed = subprocess.run(
            command,
            cwd=project_root,
            env=environment,
            stdout=handle,
            stderr=subprocess.STDOUT,
            check=False,
        )
    elapsed = round(time.monotonic() - started, 3)
    run_status = "failed"
    result_log = temp_log

    if completed.returncode == 0 and temp_report.exists():
        try:
            inject_provenance(
                temp_prefix, opm_path, config, executable, revision
            )
            verified, reasons = audit_report(
                temp_prefix, args, opm_path, config
            )
            if verified != "compatible":
                raise RuntimeError(
                    "new result failed compatibility audit: " + ";".join(reasons)
                )
            promote_prefix(temp_prefix, prefix)
            run_status = "refreshed" if audit_status == "incompatible" else "succeeded"
            result_log = final_log
        except Exception as exc:
            with temp_log.open("a", encoding="utf-8") as handle:
                handle.write(f"PROMOTION ERROR: {exc}\n")

    return {
        **row,
        "audit_status": audit_status,
        "audit_reasons": ";".join(audit_reasons),
        "run_status": run_status,
        "return_code": completed.returncode if run_status != "failed" else max(1, completed.returncode),
        "elapsed_seconds": elapsed,
        "log": str(result_log),
    }


def make_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--opm-root", type=Path, default=root / "OPM")
    parser.add_argument(
        "--output-root",
        type=Path,
        default=root / "SRP" / "SRP_UCERTAINTY_PROP_WITHINIT",
    )
    parser.add_argument("--jobs", type=int, default=min(4, max(1, os.cpu_count() or 1)))
    parser.add_argument("--eta-sigma", type=float, default=0.3)
    parser.add_argument("--days", type=float, default=15.0)
    parser.add_argument("--save-hours", type=float, default=1.0)
    parser.add_argument("--da-order", type=int, default=4)
    parser.add_argument("--max-depth", type=int, default=8)
    parser.add_argument("--pos-tol-km", type=float, default=0.1)
    parser.add_argument("--vel-tol-kms", type=float, default=1e-6)
    parser.add_argument("--expected-position-rss-km", type=float, default=10.0)
    parser.add_argument("--expected-velocity-rss-mps", type=float, default=0.03)
    parser.add_argument("--error-relative-tolerance", type=float, default=1e-9)
    parser.add_argument("--config", type=Path, default=root / "config" / "config.txt")
    parser.add_argument("--executable", type=Path)
    parser.add_argument("--no-build", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--audit-only",
        action="store_true",
        help="write compatibility audit without building or propagating",
    )
    parser.add_argument(
        "--refresh-incompatible",
        action="store_true",
        help="rerun incompatible or legacy results and replace only after success",
    )
    parser.add_argument("--force", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = make_parser().parse_args(argv)
    project_root = Path(__file__).resolve().parents[1]
    args.opm_root = args.opm_root.resolve()
    args.output_root = args.output_root.resolve()
    args.config = args.config.resolve()
    if args.jobs < 1:
        raise SystemExit("--jobs must be at least 1")
    if not math.isfinite(args.eta_sigma) or args.eta_sigma < 0.0:
        raise SystemExit("--eta-sigma must be finite and nonnegative")
    positive = (
        args.days,
        args.save_hours,
        args.pos_tol_km,
        args.vel_tol_kms,
        args.expected_position_rss_km,
        args.expected_velocity_rss_mps,
    )
    if not all(math.isfinite(value) and value > 0.0 for value in positive):
        raise SystemExit("numeric campaign parameters must be finite and positive")
    if args.da_order < 1 or args.max_depth < 0:
        raise SystemExit("invalid DA order or maximum split depth")
    if not args.opm_root.is_dir():
        raise SystemExit(f"OPM root does not exist: {args.opm_root}")
    if not args.config.is_file():
        raise SystemExit(f"configuration file does not exist: {args.config}")

    mirror_directories(args.opm_root, args.output_root)
    rows, eligible = preflight(
        args.opm_root,
        args.output_root,
        args.expected_position_rss_km,
        args.expected_velocity_rss_mps,
        args.error_relative_tolerance,
    )
    preflight_path = args.output_root / "batch_preflight.csv"
    write_csv(
        preflight_path,
        rows,
        [
            "opm", "relative_opm", "output_prefix", "position_rss_km",
            "velocity_rss_mps", "status", "reason",
        ],
    )
    rejected = len(rows) - len(eligible)
    print(f"Preflight: discovered={len(rows)} eligible={len(eligible)} rejected={rejected}")
    print(f"Preflight CSV: {preflight_path}")
    if not eligible:
        return 2

    audit_rows = []
    for row in eligible:
        status, reasons = audit_report(
            Path(row["output_prefix"]), args, Path(row["opm"]), args.config
        )
        audit_rows.append(
            {
                **row,
                "audit_status": status,
                "audit_reasons": ";".join(reasons),
            }
        )
    audit_path = args.output_root / "batch_compatibility_audit.csv"
    write_csv(
        audit_path,
        audit_rows,
        [
            "opm", "relative_opm", "output_prefix", "position_rss_km",
            "velocity_rss_mps", "status", "reason", "audit_status", "audit_reasons",
        ],
    )

    eta_domain = [0.0, 0.0] if args.eta_sigma == 0.0 else [
        -3.0 * args.eta_sigma,
        3.0 * args.eta_sigma,
    ]
    common = {
        "output_schema_version": OUTPUT_SCHEMA_VERSION,
        "discovered": len(rows),
        "eligible": len(eligible),
        "rejected_preflight": rejected,
        "jobs": args.jobs,
        "eta_sigma": args.eta_sigma,
        "eta_ads_domain": eta_domain,
        "n_variables": 6 if args.eta_sigma == 0.0 else 7,
        "nominal_srp_enabled": True,
        "srp_uncertainty_enabled": args.eta_sigma > 0.0,
        "days": args.days,
        "save_hours": args.save_hours,
        "da_order": args.da_order,
        "max_split_depth": args.max_depth,
        "position_tolerance_km": args.pos_tol_km,
        "velocity_tolerance_kms": args.vel_tol_kms,
        "expected_position_rss_km": args.expected_position_rss_km,
        "expected_velocity_rss_mps": args.expected_velocity_rss_mps,
        "config_sha256": file_sha256(args.config),
        "compatibility_audit_csv": str(audit_path),
    }
    if args.audit_only:
        write_json(
            args.output_root / "batch_summary.json",
            {"status": "audit_complete", **common},
        )
        return 0
    if args.dry_run:
        write_json(args.output_root / "batch_summary.json", {"status": "dry_run", **common})
        for row in audit_rows:
            print(
                f"DRY RUN [{row['audit_status']}]: "
                f"{row['relative_opm']} -> {row['output_prefix']}"
            )
        return 0

    build_log = args.output_root / "batch_build.log"
    if args.executable:
        executable = args.executable.resolve()
    else:
        if not args.no_build:
            with build_log.open("w", encoding="utf-8") as handle:
                built = subprocess.run(
                    ["fpm", "build"],
                    cwd=project_root,
                    stdout=handle,
                    stderr=subprocess.STDOUT,
                    check=False,
                )
            if built.returncode:
                print(f"Build failed; see {build_log}", file=sys.stderr)
                return 2
        try:
            executable = find_executable(project_root)
        except RuntimeError as exc:
            print(exc, file=sys.stderr)
            return 2
    if not executable.is_file():
        print(f"Executable does not exist: {executable}", file=sys.stderr)
        return 2

    executable = stage_executable(executable, args.output_root)
    revision = source_revision(project_root)
    common.update(
        {
            "executable": str(executable),
            "executable_sha256": file_sha256(executable),
            "source_revision": revision,
        }
    )
    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = "1"
    environment["OPENBLAS_NUM_THREADS"] = "1"
    print(f"Staged executable: {executable}")
    print(f"Launching {len(eligible)} orbit(s) with {args.jobs} parallel worker(s)")
    started = time.monotonic()
    results = []
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = {
            pool.submit(
                run_one,
                row,
                executable,
                project_root,
                args.config,
                args,
                environment,
                revision,
            ): row
            for row in eligible
        }
        for index, future in enumerate(as_completed(futures), 1):
            row = futures[future]
            try:
                result = future.result()
            except Exception as exc:
                result = {
                    **row,
                    "audit_status": "unknown",
                    "audit_reasons": "",
                    "run_status": "worker_exception",
                    "return_code": 3,
                    "elapsed_seconds": 0.0,
                    "log": "",
                    "reason": str(exc),
                }
            results.append(result)
            print(f"[{index}/{len(eligible)}] {result['run_status']}: {row['relative_opm']}")

    results.sort(key=lambda row: row["relative_opm"])
    results_path = args.output_root / "batch_results.csv"
    write_csv(
        results_path,
        results,
        [
            "opm", "relative_opm", "output_prefix", "position_rss_km",
            "velocity_rss_mps", "status", "reason", "audit_status",
            "audit_reasons", "run_status", "return_code", "elapsed_seconds", "log",
        ],
    )
    successful = ("succeeded", "refreshed", "skipped_complete")
    failures = [row for row in results if row["run_status"] not in successful]
    summary = {
        "status": "complete" if not failures else "completed_with_failures",
        **common,
        "succeeded": sum(row["run_status"] == "succeeded" for row in results),
        "refreshed": sum(row["run_status"] == "refreshed" for row in results),
        "skipped_complete": sum(row["run_status"] == "skipped_complete" for row in results),
        "failed": len(failures),
        "elapsed_seconds": round(time.monotonic() - started, 3),
        "results_csv": str(results_path),
    }
    write_json(args.output_root / "batch_summary.json", summary)
    print(f"Results CSV: {results_path}")
    print(f"Summary: {args.output_root / 'batch_summary.json'}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
