#!/usr/bin/env python3
"""Run 15-day SRP-ADS histories for matching initial OPM files in parallel."""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

POS_KEYS = ("CX_X", "CY_Y", "CZ_Z")
VEL_KEYS = ("CX_DOT_X_DOT", "CY_DOT_Y_DOT", "CZ_DOT_Z_DOT")


def find_initial_opms(root: Path) -> list[Path]:
    paths = set(root.rglob("*_init.opm"))
    paths.update(root.rglob("*_init.opm.json"))
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
    for source in sorted(path for path in opm_root.rglob("*") if path.is_dir()):
        (output_root / source.relative_to(opm_root)).mkdir(parents=True, exist_ok=True)


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    with temp.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    temp.replace(path)


def write_json(path: Path, data: dict) -> None:
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


def report_matches(prefix: Path, args: argparse.Namespace) -> bool:
    try:
        report = json.loads(
            Path(str(prefix) + "_report.json").read_text(encoding="utf-8")
        )
    except (OSError, json.JSONDecodeError):
        return False
    expected = {
        "duration_hours": args.days * 24.0,
        "save_hours": args.save_hours,
        "eta_sigma": args.eta_sigma,
        "position_tolerance_km": args.pos_tol_km,
        "velocity_tolerance_kms": args.vel_tol_kms,
        "da_order": args.da_order,
        "max_split_depth": args.max_depth,
    }
    for key, value in expected.items():
        if key not in report:
            return False
        if isinstance(value, int):
            if int(report[key]) != value:
                return False
        elif not math.isclose(float(report[key]), value, rel_tol=1e-12, abs_tol=1e-15):
            return False
    return True


def run_one(
    row: dict,
    executable: Path,
    project_root: Path,
    config: Path,
    args: argparse.Namespace,
    environment: dict[str, str],
) -> dict:
    prefix = Path(row["output_prefix"])
    prefix.parent.mkdir(parents=True, exist_ok=True)
    report = Path(str(prefix) + "_report.json")
    log = Path(str(prefix) + "_run.log")
    if report.exists() and not args.force:
        status = "skipped_complete" if report_matches(prefix, args) else "conflict_existing_report"
        return {
            **row,
            "run_status": status,
            "return_code": 0 if status == "skipped_complete" else 2,
            "elapsed_seconds": 0.0,
            "log": str(log),
        }
    command = [
        str(executable),
        "-opm", str(Path(row["opm"]).resolve()),
        "-o", str(prefix.resolve()),
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
    with log.open("w", encoding="utf-8") as handle:
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
    status = "succeeded" if completed.returncode == 0 and report.exists() else "failed"
    return {
        **row,
        "run_status": status,
        "return_code": completed.returncode,
        "elapsed_seconds": elapsed,
        "log": str(log),
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
    positive = (
        args.eta_sigma, args.days, args.save_hours, args.pos_tol_km,
        args.vel_tol_kms, args.expected_position_rss_km, args.expected_velocity_rss_mps,
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
    common = {
        "discovered": len(rows),
        "eligible": len(eligible),
        "rejected_preflight": rejected,
        "jobs": args.jobs,
        "eta_sigma": args.eta_sigma,
        "eta_ads_domain": [-3.0 * args.eta_sigma, 3.0 * args.eta_sigma],
        "days": args.days,
        "save_hours": args.save_hours,
        "da_order": args.da_order,
        "max_split_depth": args.max_depth,
        "position_tolerance_km": args.pos_tol_km,
        "velocity_tolerance_kms": args.vel_tol_kms,
        "expected_position_rss_km": args.expected_position_rss_km,
        "expected_velocity_rss_mps": args.expected_velocity_rss_mps,
    }
    if args.dry_run:
        write_json(args.output_root / "batch_summary.json", {"status": "dry_run", **common})
        for row in eligible:
            print(f"DRY RUN: {row['relative_opm']} -> {row['output_prefix']}")
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

    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = "1"
    environment["OPENBLAS_NUM_THREADS"] = "1"
    print(f"Executable: {executable}")
    print(f"Launching {len(eligible)} orbit(s) with {args.jobs} parallel worker(s)")
    started = time.monotonic()
    results = []
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = {
            pool.submit(
                run_one, row, executable, project_root, args.config, args, environment
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
            "velocity_rss_mps", "status", "reason", "run_status",
            "return_code", "elapsed_seconds", "log",
        ],
    )
    failures = [
        row for row in results
        if row["run_status"] not in ("succeeded", "skipped_complete")
    ]
    summary = {
        "status": "complete" if not failures else "completed_with_failures",
        **common,
        "succeeded": sum(row["run_status"] == "succeeded" for row in results),
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
