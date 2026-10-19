#!/usr/bin/env python3
"""Run no-solar-radiation-pressure initial-error ADS Monte Carlo jobs."""
from __future__ import annotations

import argparse
import csv
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

import h5py

from pack_srp_mc_initial_points_hdf5 import pack_initial_points
from run_srp_mc_campaign import (
    ORBITS as _BASE_ORBITS,
    CampaignError,
    atomic_csv,
    atomic_json,
    audit_opm,
    file_sha256,
    publish_directory,
)


ORBITS = _BASE_ORBITS
NO_SRP_GROUP = "no_srp"
OUTPUT_SCHEMA = "no_srp_initial_mc_campaign_v1"


def discover_opms(opm_root: str | Path) -> dict[str, Path]:
    root = Path(opm_root)
    result = {}
    problems = []
    for orbit in ORBITS:
        candidates = [
            root / orbit / f"{orbit}_init.opm",
            root / orbit / f"{orbit}_init.opm.json",
        ]
        present = [path for path in candidates if path.is_file()]
        if len(present) != 1:
            problems.append(
                f"{orbit}: expected one direct init OPM, found {len(present)}"
            )
        else:
            result[orbit] = present[0]
    if problems:
        raise CampaignError("; ".join(problems))
    return result


def find_executable(project_root: Path) -> Path:
    candidates = [
        path
        for path in (project_root / "build").glob(
            "*/app/run_srp_ads_mc_initial_only"
        )
        if path.is_file() and os.access(path, os.X_OK)
    ]
    if not candidates:
        raise CampaignError(
            "cannot locate build/*/app/run_srp_ads_mc_initial_only"
        )
    return max(candidates, key=lambda path: path.stat().st_mtime)


def stage_executable(executable: Path, experiment_root: Path) -> Path:
    digest = file_sha256(executable)
    runtime = experiment_root / ".runtime"
    runtime.mkdir(parents=True, exist_ok=True)
    target = runtime / f"run_srp_ads_mc_initial_only_{digest[:16]}"
    if not target.exists() or file_sha256(target) != digest:
        temporary = target.with_suffix(".tmp")
        shutil.copy2(executable, temporary)
        temporary.chmod(temporary.stat().st_mode | 0o111)
        temporary.replace(target)
    return target


def make_work_items(opms, config):
    return [
        (NO_SRP_GROUP, orbit, opms[orbit], config)
        for orbit in ORBITS
        if orbit in opms
    ]


def build_command(
    executable: Path,
    opm: Path,
    prefix: Path,
    samples: Path,
    config: Path,
    args,
) -> list[str]:
    return [
        str(executable),
        "-opm", str(opm.resolve()),
        "-o", str(prefix.resolve()),
        "--mc-samples", str(samples.resolve()),
        "--mc-count", str(args.mc_count),
        "--visualization-count", str(args.visualization_count),
        "--days", str(args.days),
        "--save-hours", str(args.save_hours),
        "--da-order", str(args.da_order),
        "--max-depth", str(args.max_depth),
        "--pos-tol-km", str(args.pos_tol_km),
        "--vel-tol-kms", str(args.vel_tol_kms),
        "-cfg", str(config.resolve()),
    ]


def inject_provenance(
    report_path: Path,
    attitude: str,
    orbit: str,
    opm: Path,
    config: Path,
    samples: Path,
    executable: Path,
    args,
) -> dict:
    report = json.loads(report_path.read_text(encoding="utf-8"))
    report.update({
        "orbit_name": orbit,
        "campaign_attitude": attitude,
        "campaign_force_model": NO_SRP_GROUP,
        "source_opm": str(opm.resolve()),
        "opm_sha256": file_sha256(opm),
        "source_config": str(config.resolve()),
        "config_sha256": file_sha256(config),
        "sample_file": str(samples.resolve()),
        "sample_sha256": file_sha256(samples),
        "sample_seed": 20261005,
        "sample_pairing_dimensions": 6,
        "executable_sha256": file_sha256(executable),
        "campaign_schema_version": OUTPUT_SCHEMA,
        "expected_hour_count": int(round(args.days * 24.0 / args.save_hours)) + 1,
    })
    atomic_json(report_path, report)
    return report


def _csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def _all_numeric_finite(rows, text_fields=()):
    for row in rows:
        for key, value in row.items():
            if key in text_fields:
                continue
            try:
                if not math.isfinite(float(value)):
                    return False
            except (TypeError, ValueError):
                return False
    return True


def organize_outputs(stage: Path, prefix: Path, attitude: str, orbit: str, args):
    visualization = Path(str(prefix) + "_mc_visualization_points.csv")
    keyframes = Path(str(prefix) + "_mc_keyframe_points.csv")
    report = Path(str(prefix) + "_mc_report.json")
    pack_initial_points(
        visualization,
        stage / f"{orbit}_visualization_points.h5",
        min(args.visualization_count, args.mc_count),
        "visualization",
        report,
    )
    pack_initial_points(
        keyframes,
        stage / f"{orbit}_keyframe_points.h5",
        args.mc_count,
        "keyframe",
        report,
    )
    mappings = {
        Path(str(prefix) + "_mc_hourly_stats.csv"):
            stage / f"{orbit}_hourly_mc_stats.csv",
        Path(str(prefix) + "_split_events.csv"):
            stage / f"{orbit}_split_events.csv",
        report: stage / f"{orbit}_events.json",
        Path(str(prefix) + "_run.log"): stage / f"{orbit}_run.log",
    }
    for source, target in mappings.items():
        if not source.is_file():
            raise CampaignError(f"missing staged artifact: {source.name}")
        source.replace(target)
    visualization.unlink()
    keyframes.unlink()

    snapshot_dir = stage / "patch_snapshots"
    snapshot_dir.mkdir()
    for source in sorted(stage.glob(f"{orbit}_patch_s*")):
        source.replace(snapshot_dir / source.name)
    if not any(snapshot_dir.iterdir()):
        raise CampaignError("no patch snapshots were produced")


def audit_result(
    directory: str | Path, attitude: str, orbit: str, args
) -> tuple[bool, list[str]]:
    directory = Path(directory)
    reasons = []
    expected_hours = int(round(args.days * 24.0 / args.save_hours)) + 1
    required = {
        "hourly": directory / f"{orbit}_hourly_mc_stats.csv",
        "events": directory / f"{orbit}_events.json",
        "split": directory / f"{orbit}_split_events.csv",
        "visual": directory / f"{orbit}_visualization_points.h5",
        "keyframe": directory / f"{orbit}_keyframe_points.h5",
        "log": directory / f"{orbit}_run.log",
        "snapshots": directory / "patch_snapshots",
    }
    for key, path in required.items():
        if key == "snapshots":
            if not path.is_dir() or not any(path.iterdir()):
                reasons.append("missing_snapshots")
        elif not path.is_file():
            reasons.append(f"missing_{key}")
    if reasons:
        return False, reasons

    try:
        hourly = _csv_rows(required["hourly"])
        if len(hourly) != expected_hours * 3:
            reasons.append("hourly_row_count")
        expected_grid = {
            float(index) * args.save_hours for index in range(expected_hours)
        }
        if {float(row["hour"]) for row in hourly} != expected_grid:
            reasons.append("hourly_time_grid")
        if {row["group"] for row in hourly} != {attitude}:
            reasons.append("hourly_group")
        if {int(row["sample_count"]) for row in hourly} != {
            10000, 20000, 50000
        }:
            reasons.append("hourly_prefixes")
        if {int(row["active_dimensions"]) for row in hourly} != {6}:
            reasons.append("hourly_dimensions")
        if any(int(row["split_dim7"]) != 0 for row in hourly):
            reasons.append("seventh_direction_split")
        if not _all_numeric_finite(hourly, ("group",)):
            reasons.append("hourly_nonfinite")
    except (OSError, ValueError, KeyError, TypeError):
        reasons.append("hourly_parse")
        hourly = []

    try:
        report = json.loads(required["events"].read_text(encoding="utf-8"))
        if report.get("output_schema_version") != "srp_ads_mc_initial_only_v1":
            reasons.append("report_schema")
        if report.get("nominal_srp_enabled") is not False:
            reasons.append("report_nominal_srp")
        if report.get("srp_parameter_uncertainty_enabled") is not False:
            reasons.append("report_srp_uncertainty")
        if report.get("srp_model") != "cannonball":
            reasons.append("report_model")
        if report.get("srp_attitude_mode") != attitude:
            reasons.append("report_attitude")
        if int(report.get("active_dimensions", -1)) != 6:
            reasons.append("report_dimensions")
        if float(report.get("srp_eta", math.nan)) != 0.0:
            reasons.append("report_eta")
        if int(report.get("sample_count", -1)) != args.mc_count:
            reasons.append("report_sample_count")
        if report.get("orbit_name") != orbit:
            reasons.append("report_orbit")
        diagnostics = report.get("geometry_fallback_diagnostics")
        if not isinstance(diagnostics, dict):
            reasons.append("report_geometry_fallback_diagnostics")
        else:
            for key in ("real_count", "da_count"):
                value = diagnostics.get(key)
                if not isinstance(value, int) or value < 0:
                    reasons.append(f"report_geometry_{key}")
            for key in ("first_real_message", "first_da_message"):
                if not isinstance(diagnostics.get(key), str):
                    reasons.append(f"report_geometry_{key}")
            if diagnostics.get("real_count", 0) > 0 and not diagnostics.get("first_real_message"):
                reasons.append("report_geometry_real_message")
            if diagnostics.get("da_count", 0) > 0 and not diagnostics.get("first_da_message"):
                reasons.append("report_geometry_da_message")
        for key in (
            "opm_sha256", "config_sha256", "sample_sha256",
            "executable_sha256",
        ):
            value = report.get(key)
            if (
                not isinstance(value, str)
                or len(value) != 64
                or any(character not in "0123456789abcdef"
                       for character in value.lower())
            ):
                reasons.append(f"report_{key}")
    except (OSError, json.JSONDecodeError, TypeError, ValueError):
        reasons.append("report_parse")
        report = {}

    try:
        event_rows = _csv_rows(required["split"])
        finest = [
            row for row in hourly if int(row["sample_count"]) == 50000
        ]
        expected_events = sum(
            sum(int(row[f"split_dim{dimension}"]) for dimension in range(1, 7))
            for row in finest
        )
        if len(event_rows) != expected_events:
            reasons.append("split_event_count")
        if event_rows and not _all_numeric_finite(event_rows):
            reasons.append("split_event_nonfinite")
        first = report.get("events", {}).get("first_split", {})
        if bool(first.get("observed")) != bool(event_rows):
            reasons.append("first_split_presence")
    except (OSError, KeyError, TypeError, ValueError):
        reasons.append("split_event_parse")

    try:
        with h5py.File(required["visual"], "r") as h5:
            if h5.attrs.get("output_schema_version") !=                     "srp_mc_initial_points_hdf5_v1":
                reasons.append("visual_schema")
            expected_shape = (
                expected_hours,
                min(args.visualization_count, args.mc_count),
                6,
            )
            if h5["state"].shape != expected_shape or                     h5["error"].shape != expected_shape:
                reasons.append("visual_shape")
            if len(set(map(int, h5["sample_id"][0]))) != expected_shape[1]:
                reasons.append("visual_ids")
        with h5py.File(required["keyframe"], "r") as h5:
            if h5["state"].shape[1:] != (args.mc_count, 6):
                reasons.append("keyframe_shape")
            hours = set(map(float, h5["hour"][:]))
            final_hour = args.days * 24.0
            if 0.0 not in hours or final_hour not in hours:
                reasons.append("keyframe_fixed_hours")
    except (OSError, KeyError, ValueError):
        reasons.append("hdf5_parse")
    return not reasons, reasons


def write_file_manifest(directory: Path, attitude: str, orbit: str):
    files = []
    for path in sorted(directory.rglob("*")):
        if path.is_file() and path.name != f"{orbit}_file_manifest.json":
            files.append({
                "relative_path": str(path.relative_to(directory)),
                "size_bytes": path.stat().st_size,
                "sha256": file_sha256(path),
            })
    atomic_json(directory / f"{orbit}_file_manifest.json", {
        "output_schema_version": "no_srp_initial_orbit_manifest_v1",
        "attitude": attitude,
        "orbit_name": orbit,
        "files": files,
    })

def run_one(
    attitude: str,
    orbit: str,
    opm: Path,
    executable: Path,
    project_root: Path,
    experiment_root: Path,
    samples: Path,
    config: Path,
    args,
    environment: dict[str, str],
) -> dict:
    results_root = experiment_root / "results_no_srp_initial_only"
    final = results_root / orbit
    if final.is_dir() and not args.force:
        valid, reasons = audit_result(final, attitude, orbit, args)
        if valid:
            return {
                "attitude": attitude,
                "orbit": orbit,
                "status": "skipped_complete",
                "elapsed_seconds": 0.0,
                "reason": "",
                "result_directory": str(final),
            }
        if not args.refresh_incompatible:
            return {
                "attitude": attitude,
                "orbit": orbit,
                "status": "incompatible_existing",
                "elapsed_seconds": 0.0,
                "reason": ";".join(reasons),
                "result_directory": str(final),
            }

    staging_root = results_root / ".staging"
    staging_root.mkdir(parents=True, exist_ok=True)
    stage = staging_root / f"{orbit}.{uuid.uuid4().hex}"
    stage.mkdir()
    prefix = stage / orbit
    log = Path(str(prefix) + "_run.log")
    command = build_command(
        executable, opm, prefix, samples, config, args
    )
    started = time.monotonic()
    with log.open("w", encoding="utf-8") as stream:
        stream.write("COMMAND: " + " ".join(command) + "\n")
        stream.flush()
        completed = subprocess.run(
            command,
            cwd=project_root,
            env=environment,
            stdout=stream,
            stderr=subprocess.STDOUT,
            check=False,
        )
    elapsed = round(time.monotonic() - started, 3)
    try:
        report_path = Path(str(prefix) + "_mc_report.json")
        if completed.returncode != 0 or not report_path.is_file():
            raise CampaignError(
                f"propagator failed with return code {completed.returncode}"
            )
        inject_provenance(
            report_path, attitude, orbit, opm, config, samples,
            executable, args
        )
        organize_outputs(stage, prefix, attitude, orbit, args)
        valid, reasons = audit_result(stage, attitude, orbit, args)
        if not valid:
            raise CampaignError(";".join(reasons))
        write_file_manifest(stage, attitude, orbit)
        publish_directory(
            stage,
            final,
            lambda path: audit_result(path, attitude, orbit, args),
        )
        return {
            "attitude": attitude,
            "orbit": orbit,
            "status": "succeeded",
            "elapsed_seconds": elapsed,
            "reason": "",
            "result_directory": str(final),
        }
    except Exception as exc:
        with log.open("a", encoding="utf-8") as stream:
            stream.write(f"CAMPAIGN ERROR: {exc}\n")
        return {
            "attitude": attitude,
            "orbit": orbit,
            "status": "failed",
            "elapsed_seconds": elapsed,
            "reason": str(exc),
            "result_directory": str(stage),
        }


def run_jobs(
    work_items,
    executable: Path,
    project_root: Path,
    experiment_root: Path,
    samples: Path,
    args,
    environment: dict[str, str],
    jobs: int,
):
    results = []
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        futures = {
            pool.submit(
                run_one,
                attitude,
                orbit,
                opm,
                executable,
                project_root,
                experiment_root,
                samples,
                config,
                args,
                environment,
            ): (attitude, orbit)
            for attitude, orbit, opm, config in work_items
        }
        for index, future in enumerate(as_completed(futures), 1):
            attitude, orbit = futures[future]
            try:
                row = future.result()
            except Exception as exc:
                row = {
                    "attitude": attitude,
                    "orbit": orbit,
                    "status": "worker_exception",
                    "elapsed_seconds": 0.0,
                    "reason": str(exc),
                    "result_directory": "",
                }
            results.append(row)
            print(
                f"[{index}/{len(work_items)}] {attitude}/{orbit}: "
                f"{row['status']}",
                flush=True,
            )
    order = {
        (attitude, orbit): index
        for index, (attitude, orbit, _, _) in enumerate(work_items)
    }
    results.sort(key=lambda row: order[(row["attitude"], row["orbit"])])
    return results


def audit_all(experiment_root: Path, args):
    rows = []
    for orbit in ORBITS:
        directory = (
            experiment_root / "results_no_srp_initial_only" / orbit
        )
        valid, reasons = (
            audit_result(directory, NO_SRP_GROUP, orbit, args)
            if directory.is_dir()
            else (False, ["missing_result_directory"])
        )
        rows.append({
            "attitude": NO_SRP_GROUP,
            "orbit": orbit,
            "complete": valid,
            "reasons": ";".join(reasons),
        })
    atomic_csv(
        experiment_root / "audit" / "no_srp_output_completeness.csv",
        rows,
        ["attitude", "orbit", "complete", "reasons"],
    )
    return rows


def make_parser() -> argparse.ArgumentParser:
    project_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=project_root)
    parser.add_argument("--opm-root", type=Path, default=project_root / "OPM")
    parser.add_argument(
        "--experiment-root", type=Path,
        default=project_root / "20261005SRP",
    )
    parser.add_argument("--samples", type=Path)
    parser.add_argument("--executable", type=Path)
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--mc-count", type=int, default=50000)
    parser.add_argument("--visualization-count", type=int, default=2048)
    parser.add_argument("--days", type=float, default=15.0)
    parser.add_argument("--save-hours", type=float, default=1.0)
    parser.add_argument("--da-order", type=int, default=4)
    parser.add_argument("--max-depth", type=int, default=8)
    parser.add_argument("--pos-tol-km", type=float, default=0.1)
    parser.add_argument("--vel-tol-kms", type=float, default=1e-6)
    parser.add_argument("--expected-position-rss-km", type=float, default=10.0)
    parser.add_argument("--expected-velocity-rss-mps", type=float, default=0.03)
    parser.add_argument("--error-relative-tolerance", type=float, default=1e-9)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--audit-only", action="store_true")
    parser.add_argument("--refresh-incompatible", action="store_true")
    parser.add_argument("--force", action="store_true")
    return parser


def main(argv=None) -> int:
    args = make_parser().parse_args(argv)
    args.project_root = args.project_root.resolve()
    args.opm_root = args.opm_root.resolve()
    args.experiment_root = args.experiment_root.resolve()
    if args.jobs < 1 or args.mc_count < 2 or args.visualization_count < 1:
        raise SystemExit("jobs/counts must be positive and MC count at least two")
    if (
        args.days <= 0.0
        or args.save_hours <= 0.0
        or args.da_order < 1
        or args.max_depth < 0
    ):
        raise SystemExit("invalid propagation duration, cadence, order, or depth")
    intervals = args.days * 24.0 / args.save_hours
    if not math.isclose(intervals, round(intervals), abs_tol=1e-12):
        raise SystemExit(
            "duration must contain an integral number of save intervals"
        )

    config = (args.project_root / "config" / "config_no_srp.txt").resolve()
    if not config.is_file():
        raise SystemExit(f"missing no-SRP configuration: {config}")
    try:
        opms = discover_opms(args.opm_root)
    except CampaignError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    for name in ("audit", "config", "logs", "population", "samples"):
        (args.experiment_root / name).mkdir(parents=True, exist_ok=True)
    audits = []
    for orbit, opm in opms.items():
        try:
            row = audit_opm(
                opm,
                args.expected_position_rss_km,
                args.expected_velocity_rss_mps,
                args.error_relative_tolerance,
            )
        except CampaignError as exc:
            row = {
                "opm": str(opm),
                "opm_sha256": "",
                "position_rss_km": "",
                "velocity_rss_mps": "",
                "positive_definite": False,
                "status": "invalid_opm",
                "reason": str(exc),
            }
        row["orbit"] = orbit
        audits.append(row)
    atomic_csv(
        args.experiment_root / "audit" / "no_srp_input_audit.csv",
        audits,
        [
            "orbit", "opm", "opm_sha256", "position_rss_km",
            "velocity_rss_mps", "positive_definite", "status", "reason",
        ],
    )
    rejected = [row for row in audits if row["status"] != "eligible"]
    if rejected:
        print(
            f"ERROR: {len(rejected)} whitelist OPM inputs failed preflight",
            file=sys.stderr,
        )
        return 2

    experiment_config = {
        "output_schema_version": OUTPUT_SCHEMA,
        "orbits": list(ORBITS),
        "groups": [NO_SRP_GROUP],
        "nominal_srp_enabled": False,
        "srp_model": "disabled",
        "srp_eta": 0.0,
        "active_dimensions": 6,
        "mc_count": args.mc_count,
        "visualization_count": args.visualization_count,
        "sample_seed": 20261005,
        "days": args.days,
        "save_hours": args.save_hours,
        "da_order": args.da_order,
        "max_split_depth": args.max_depth,
        "position_tolerance_km": args.pos_tol_km,
        "velocity_tolerance_kms": args.vel_tol_kms,
        "expected_position_rss_km": args.expected_position_rss_km,
        "expected_velocity_rss_mps": args.expected_velocity_rss_mps,
        "config": {
            "path": str(config),
            "sha256": file_sha256(config),
        },
    }
    atomic_json(
        args.experiment_root / "config" / "no_srp_experiment_config.json",
        experiment_config,
    )
    if args.dry_run:
        print(json.dumps({
            "status": "dry_run",
            "eligible_orbits": len(audits),
            "work_item_count": len(audits),
        }))
        return 0

    if args.audit_only:
        rows = audit_all(args.experiment_root, args)
        return 0 if all(row["complete"] for row in rows) else 1

    if args.samples is None:
        samples = (
            args.experiment_root / "samples" /
            f"mc_standard_normal_{args.mc_count}.csv"
        )
        sample_manifest = (
            args.experiment_root / "samples" / "sample_manifest.json"
        )
        if not samples.is_file():
            generated = subprocess.run([
                sys.executable,
                str(args.project_root / "tools" /
                    "generate_srp_mc_samples.py"),
                "--csv", str(samples),
                "--manifest", str(sample_manifest),
                "--count", str(args.mc_count),
                "--seed", "20261005",
                "--eta-sigma", "0.3",
            ], cwd=args.project_root, check=False)
            if generated.returncode != 0:
                return 2
    else:
        samples = args.samples.resolve()
        if not samples.is_file():
            raise SystemExit(f"sample file does not exist: {samples}")

    try:
        executable = (
            args.executable.resolve()
            if args.executable is not None
            else find_executable(args.project_root)
        )
        if not executable.is_file():
            raise CampaignError(f"executable does not exist: {executable}")
        executable = stage_executable(executable, args.experiment_root)
    except CampaignError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = "1"
    environment["OPENBLAS_NUM_THREADS"] = "1"
    library_paths = [
        str(args.project_root / "lib"),
        str(Path.home() / "local_dace" / "lib"),
        str(Path.home() / "local_openblas" / "lib"),
    ]
    environment["LD_LIBRARY_PATH"] = ":".join(
        library_paths + [environment.get("LD_LIBRARY_PATH", "")]
    )
    environment["LIBRARY_PATH"] = ":".join(
        library_paths + [environment.get("LIBRARY_PATH", "")]
    )

    work_items = make_work_items(opms, config)
    started = time.monotonic()
    results = run_jobs(
        work_items, executable, args.project_root, args.experiment_root,
        samples, args, environment, args.jobs
    )
    atomic_csv(
        args.experiment_root / "audit" / "no_srp_campaign_results.csv",
        results,
        [
            "attitude", "orbit", "status", "elapsed_seconds",
            "reason", "result_directory",
        ],
    )
    audits_out = audit_all(args.experiment_root, args)
    completed = sum(bool(row["complete"]) for row in audits_out)
    failed = [
        row for row in results
        if row["status"] not in ("succeeded", "skipped_complete")
    ]
    expected = len(work_items)
    batch = {
        "output_schema_version": OUTPUT_SCHEMA,
        "status": (
            "complete"
            if not failed and completed == expected
            else "completed_with_failures"
        ),
        "work_items_requested": expected,
        "completed_count": completed,
        "failed_count": len(failed),
        "elapsed_seconds": round(time.monotonic() - started, 3),
        "executable": str(executable),
        "executable_sha256": file_sha256(executable),
        "sample_file": str(samples),
        "sample_sha256": file_sha256(samples),
        "campaign_results_csv": str(
            args.experiment_root / "audit" / "no_srp_campaign_results.csv"
        ),
    }
    atomic_json(
        args.experiment_root / "population" / "no_srp_batch_summary.json", batch
    )
    return 1 if failed or completed != expected else 0


if __name__ == "__main__":
    raise SystemExit(main())
