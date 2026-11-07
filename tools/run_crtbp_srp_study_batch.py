#!/usr/bin/env python3
"""Batch dispatcher for the single-case SRP sensitivity study.

Reads CRTBP_SRP/crtbp_initial_conditions.json, enumerates
(orbit x orbit-phase x solar-phase x experiment-config) cases, and invokes the
Fortran single-case executable ``run_crtbp_srp_study`` once per case in
parallel. JSON reading and parallel dispatch stay here in Python; the Fortran
``app/`` handles only one case per invocation (see app/run_crtbp_srp_study.f90).

Each case runs one time-advancing ADS propagation of a rectangular SRP
parameter box, and returns endpoint mean/covariance, split diagnostics, and
(where beta=0 is an edge of the box) the d(x)/d(beta) sensitivity at beta=0.

Base 19-argument contract (order is significant):
    1-6  x0(6) = x y z vx vy vz  (nondimensional)
    7    mu                     (CRTBP mass ratio)
    8    period_nd              (orbit period, nondimensional)
    9    alpha0_rad             (solar initial phase)
    10   model                  ("cannonball" | "flatplate")
    11   beta_lower  (m2/kg)
    12   beta_upper  (m2/kg)
    13   tilt_max_deg            (0 => u,v inactive)
    14   chi_lower
    15   chi_upper               (==chi_lower => chi inactive)
    16   da_order                (e.g. 6)
    17   err_toll                (truncation-error threshold, ND)
    18   n_split_max             (max split depth)
    19   output_json_path
Angle mode appends: theta_phi theta_lower_deg theta_upper_deg phi_lower_deg phi_upper_deg.
Use --attitude theta_phi --phi-deg LOWER UPPER for Q2/Q4. Theta defaults to
[0, config tilt_deg], or is explicitly replaced by --theta-deg LOWER UPPER.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import json
import math
import os
from pathlib import Path
import shlex
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_JSON = ROOT / "CRTBP_SRP" / "crtbp_initial_conditions.json"

# ---------------------------------------------------------------------------
# Experiment configurations: one entry = one ADS task over one parameter box.
# chi=(0.6, 0.6) fixes chi at the baseline, so Q2/Q3 are 3D (beta,u,v); an
# open chi interval makes Q4 4D (beta,chi,u,v). tilt_deg=0 leaves u,v inactive.
# ---------------------------------------------------------------------------
CONFIGS = {
    # Q1 — Layer 1: 1D cannonball strength sensitivity
    "Q1_L1_small": dict(layer="L1", model="cannonball", beta=(0.00, 0.02),
                        tilt_deg=0.0, chi=(0.6, 0.6)),
    "Q1_L1_main":  dict(layer="L1", model="cannonball", beta=(0.00, 0.10),
                        tilt_deg=0.0, chi=(0.6, 0.6)),
    "Q1_L1_ext":   dict(layer="L1", model="cannonball", beta=(0.00, 0.20),
                        tilt_deg=0.0, chi=(0.6, 0.6)),
    # Q2 — Layer 2: 3D flatplate, strength vs attitude (chi fixed 0.6)
    "Q2_L2_t5":        dict(layer="L2", model="flatplate", beta=(0.08, 0.12),
                            tilt_deg=5.0, chi=(0.6, 0.6)),
    "Q2_L2_t10":       dict(layer="L2", model="flatplate", beta=(0.08, 0.12),
                            tilt_deg=10.0, chi=(0.6, 0.6)),
    "Q2_L2_t20":       dict(layer="L2", model="flatplate", beta=(0.08, 0.12),
                            tilt_deg=20.0, chi=(0.6, 0.6)),
    "Q2_L2_ext_main":  dict(layer="L2", model="flatplate", beta=(0.00, 0.10),
                            tilt_deg=10.0, chi=(0.6, 0.6)),
    "Q2_L2_ext_ext":   dict(layer="L2", model="flatplate", beta=(0.00, 0.20),
                            tilt_deg=10.0, chi=(0.6, 0.6)),
    # Q4 — Layer 3: 4D flatplate, reflection-attitude coupling (C_eq fixed)
    "Q4_L3_chi01":     dict(layer="L3", model="flatplate", beta=(0.08, 0.12),
                            tilt_deg=10.0, chi=(0.0, 1.0)),
    "Q4_L3_chi_narrow": dict(layer="L3", model="flatplate", beta=(0.08, 0.12),
                             tilt_deg=10.0, chi=(0.4, 0.8)),
    "Q4_L3_ext":       dict(layer="L3", model="flatplate", beta=(0.00, 0.10),
                            tilt_deg=10.0, chi=(0.0, 1.0)),
}


def _angle_bounds(values, flag: str, *, theta: bool = False) -> tuple:
    lo, hi = values
    if not (math.isfinite(lo) and math.isfinite(hi)) or lo > hi:
        raise ValueError(f"{flag}: need finite lower <= upper")
    if theta and not (0.0 <= lo <= hi < 90.0):
        raise ValueError(f"{flag}: need 0 <= lower <= upper < 90 degrees")
    return lo, hi


def select_configs(args) -> dict:
    """Resolve angle boxes without changing the legacy CONFIGS definitions.

    Angles have separate IDs (including bounds) so skip-existing cannot reuse
    a UV result or a different angle domain. Phi stays unwrapped throughout.
    """
    if args.attitude == "uv":
        if args.theta_deg is not None or args.phi_deg is not None:
            raise ValueError("--theta-deg/--phi-deg require --attitude theta_phi")
        return {key: dict(CONFIGS[key]) for key in args.configs}
    if args.beta is not None:
        raise ValueError("--beta selects cannonball Q1; use --configs Q2... or Q4... for theta_phi")
    if args.phi_deg is None:
        raise ValueError("--attitude theta_phi requires --phi-deg LOWER UPPER (degrees)")
    phi = _angle_bounds(args.phi_deg, "--phi-deg")
    if args.theta_deg is not None:
        _angle_bounds(args.theta_deg, "--theta-deg", theta=True)
    selected = {}
    angle_count = 0
    for key in args.configs:
        cfg = dict(CONFIGS[key])
        if cfg["model"] == "flatplate":
            theta = _angle_bounds(args.theta_deg if args.theta_deg is not None
                                  else (0.0, cfg["tilt_deg"]), "--theta-deg", theta=True)
            cfg.update(attitude_coordinates="theta_phi", theta_deg=theta, phi_deg=phi)
            key += f"__theta_phi_t{fmt(theta[0])}_{fmt(theta[1])}_p{fmt(phi[0])}_{fmt(phi[1])}"
            angle_count += 1
        selected[key] = cfg
    if not angle_count:
        raise ValueError("--attitude theta_phi needs at least one flatplate Q2/Q4 configuration")
    return selected


def write_json(path: Path, value) -> None:
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
                    encoding="utf-8")


def fmt(v: float) -> str:
    """Free-format double that Fortran list-directed read accepts."""
    if not math.isfinite(v):
        raise ValueError(f"non-finite value {v!r} cannot be passed to Fortran")
    return f"{v:.17g}"


def environment() -> dict:
    """Compiler/runtime library paths for DACE/SPICE/OpenBLAS, mirroring
    setup_env.sh (CPLUS_INCLUDE_PATH for the DACE headers, LIBRARY_PATH for
    linking, LD_LIBRARY_PATH for the runtime loader)."""
    env = os.environ.copy()
    home = Path.home()
    dace = Path(env.get("DACE_ROOT", home / "local_dace"))
    spice = Path(env.get("SPICE_ROOT", home / "spice_fortran" / "toolkit"))
    openblas = home / "local_openblas"
    libdirs = [str(ROOT / "lib"), str(dace / "lib"), str(spice / "lib"), str(openblas / "lib")]
    for key, entries in (("CPLUS_INCLUDE_PATH", [str(dace / "include")]),
                         ("LIBRARY_PATH", libdirs),
                         ("LD_LIBRARY_PATH", libdirs)):
        cur = env.get(key, "")
        env[key] = ":".join(entries) + (":" + cur if cur else "")
    return env


def resolve_executable(args, env: dict) -> Path:
    if args.executable:
        path = args.executable.resolve(strict=True)
        if not os.access(path, os.X_OK):
            raise ValueError(f"executable is not executable: {path}")
        return path
    marker = "__POD_SRP_STUDY_TARGET__ "
    log_path = ROOT / "build" / "srp-study-build.log"
    log_path.parent.mkdir(exist_ok=True)
    with log_path.open("w") as log:
        result = subprocess.run(
            ["fpm", "run", "--target", "run_crtbp_srp_study",
             "--runner", "echo " + marker],
            cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError("Fortran build failed: " + log_path.read_text()[-4000:])
    targets = [line[len(marker):].strip() for line in log_path.read_text().splitlines()
               if line.strip().startswith(marker)]
    if len(targets) != 1:
        raise RuntimeError("fpm did not report exactly one run_crtbp_srp_study target")
    binary = Path(targets[0])
    if not binary.is_absolute():
        binary = ROOT / binary
    binary = binary.resolve(strict=True)
    if not binary.is_file() or not os.access(binary, os.X_OK):
        raise RuntimeError("fpm target is not an executable file: " + str(binary))
    return binary


def load_orbits(path: Path) -> list:
    data = json.loads(path.read_text(encoding="utf-8-sig"))
    orbits = data.get("orbits")
    if not isinstance(orbits, list) or not orbits:
        raise ValueError(f"no orbits in {path}")
    return orbits


def git_commit() -> str:
    try:
        out = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT,
                             capture_output=True, text=True, check=True)
        return out.stdout.strip()
    except Exception:
        return "unknown"


def enumerate_cases(args, orbits: list, configs=None):
    """Yield (case_key, orbit, phase, phase_index, solar_phase, cfg_id, cfg)."""
    configs = CONFIGS if configs is None else configs
    orbit_filter = None if not args.orbits else set(args.orbits)
    phase_filter = None if not args.phases else set(args.phases)
    prefix = args.orbit_prefix
    count = 0
    for orbit in orbits:
        oid = orbit.get("orbit_id")
        if orbit_filter is not None and oid not in orbit_filter:
            continue
        if prefix and not oid.startswith(prefix):
            continue
        mu = orbit.get("mu")
        period_nd = orbit.get("period_nd")
        phases = orbit.get("phases")
        if not isinstance(phases, list) or not phases:
            raise ValueError(f"orbit {oid}: missing phases")
        if mu is None or period_nd is None:
            raise ValueError(f"orbit {oid}: missing mu or period_nd")
        for phase in phases:
            pi = phase.get("phase_index")
            state = phase.get("state_nd")
            if pi is None or not isinstance(state, list) or len(state) != 6:
                raise ValueError(f"orbit {oid}: malformed phase")
            if phase_filter is not None and pi not in phase_filter:
                continue
            for solar_phase in args.solar_phases:
                for cfg_id in args.configs:
                    cfg = configs[cfg_id]
                    key = f"{cfg_id}__{oid}__p{int(pi):02d}__a{fmt(solar_phase)}"
                    yield key, orbit, phase, pi, solar_phase, cfg_id, cfg
                    count += 1
                    if args.limit and count >= args.limit:
                        return


def build_command(exe: Path, orbit, phase, pi, solar_phase, cfg_id, cfg, args, out_path: Path) -> list:
    beta_lo, beta_hi = cfg["beta"]
    chi_lo, chi_hi = cfg["chi"]
    cmd = [str(exe)]
    for comp in phase["state_nd"]:
        cmd.append(fmt(float(comp)))
    cmd += [
        fmt(float(orbit["mu"])),
        fmt(float(orbit["period_nd"])),
        fmt(solar_phase),
        cfg["model"],
        fmt(beta_lo),
        fmt(beta_hi),
        fmt(0.0 if cfg.get("attitude_coordinates") == "theta_phi" else cfg["tilt_deg"]),
        fmt(chi_lo),
        fmt(chi_hi),
        str(args.da_order),
        fmt(args.err_toll),
        str(args.n_split_max),
        str(out_path),
    ]
    if cfg.get("attitude_coordinates") == "theta_phi":
        # Fortran performs degree->radian conversion once, in its domain factory.
        cmd += ["theta_phi", *(fmt(v) for v in cfg["theta_deg"]), *(fmt(v) for v in cfg["phi_deg"])]
    return cmd


def read_app_result(path: Path) -> dict:
    try:
        return json.loads(path.read_text(encoding="utf-8-sig"))
    except Exception:
        return {}


def _matching_coordinates(res: dict, cfg: dict) -> bool:
    # A pre-angle binary may silently ignore trailing arguments. Do not count
    # its UV result as a completed theta/phi case or skip it on a later run.
    return res.get("attitude_coordinates", "uv") == cfg.get("attitude_coordinates", "uv")


def _fill_from_result(record: dict, res: dict) -> None:
    """Copy app-derived fields out of a per-case result JSON into a campaign record."""
    for field in ("parameter_names", "parameter_units", "parameter_lower",
                  "parameter_upper", "parameter_da_index", "parameter_prior"):
        if field in res:
            record[field] = res[field]
    record["app_status"] = res.get("status")
    record["wall_time_s"] = res.get("wall_time_s")
    record["n_patches"] = res.get("n_patches")
    record["n_splits_total"] = res.get("n_splits_total")
    record["n_max_active"] = res.get("n_max_active")
    record["first_split_time"] = res.get("first_split_time")
    record["n_diverged"] = res.get("n_diverged")
    record["first_divergence_time"] = res.get("first_divergence_time")
    record["split_count_dir"] = res.get("split_count_dir")
    record["t_end_days"] = res.get("t_end_days")
    record["has_beta_sens"] = res.get("beta_sensitivity") is not None


def run_one(job, args, exe: Path, env: dict) -> dict:
    key, orbit, phase, pi, solar_phase, cfg_id, cfg = job
    oid = orbit.get("orbit_id")
    out_path = args.outdir / "results" / cfg_id / oid / f"p{int(pi):02d}__a{fmt(solar_phase)}.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    record = {
        "config": cfg_id, "layer": cfg["layer"], "model": cfg["model"],
        "orbit": oid, "family": orbit.get("family"),
        "libration_point": orbit.get("libration_point"),
        "phase_index": pi, "phase_fraction": phase.get("phase_fraction"),
        "solar_phase_rad": solar_phase, "beta": list(cfg["beta"]),
        "tilt_deg": cfg["tilt_deg"], "chi": list(cfg["chi"]),
        "attitude_coordinates": cfg.get("attitude_coordinates", "uv"),
        "theta_deg": cfg.get("theta_deg"), "phi_deg": cfg.get("phi_deg"),
        "da_order": args.da_order, "err_toll": args.err_toll,
        "n_split_max": args.n_split_max, "result_path": str(out_path.relative_to(args.outdir)),
    }

    if args.skip_existing and out_path.is_file():
        res = read_app_result(out_path)
        # A completed case is any Fortran-written JSON carrying a "status" key,
        # regardless of its value: status<0 (-2 empty manifold, -4 sample in a
        # diverged region) is still a valid, finished result, not a rerun trigger.
        if res and "status" in res and _matching_coordinates(res, cfg):
            record["status"] = "skipped_existing"
            record["elapsed_s"] = 0.0
            _fill_from_result(record, res)
            print(f"{cfg_id} {oid} p{int(pi):02d} a={fmt(solar_phase)}: "
                  f"skipped_existing (already on disk)", flush=True)
            return record

    cmd = build_command(exe, orbit, phase, pi, solar_phase, cfg_id, cfg, args, out_path)
    (out_path.parent / (out_path.stem + ".sh")).write_text(
        "#!/usr/bin/env bash\nset -euo pipefail\n" + shlex.join(cmd) + "\n", encoding="utf-8")
    started = time.perf_counter()
    with (out_path.parent / (out_path.stem + ".log")).open("w") as log:
        proc = subprocess.run(cmd, cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
    elapsed = time.perf_counter() - started
    record["elapsed_s"] = elapsed
    if proc.returncode:
        record["status"] = "process_failed"
        record["returncode"] = proc.returncode
    else:
        res = read_app_result(out_path)
        record["status"] = "completed" if res else "no_output"
        if res and not _matching_coordinates(res, cfg):
            record["status"] = "configuration_mismatch"
            record["reason"] = "result attitude mode differs from request; rebuild the Fortran executable"
        _fill_from_result(record, res)
    print(f"{cfg_id} {oid} p{int(pi):02d} a={fmt(solar_phase)}: "
          f"{record['status']}; {elapsed:.2f}s", flush=True)
    return record


def parser():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--json", type=Path, default=DEFAULT_JSON)
    p.add_argument("--configs", nargs="+", default=list(CONFIGS), choices=list(CONFIGS),
                   help="experiment configs to run")
    p.add_argument("--beta", type=float, nargs=2, metavar=("LOWER", "UPPER"),
                   help="run a custom cannonball (Q1) config over this beta range; "
                        "overrides --configs")
    p.add_argument("--attitude", choices=("uv", "theta_phi"), default="uv",
                   help="attitude coordinates for selected flatplate Q2/Q4 configs (default: uv)")
    p.add_argument("--theta-deg", type=float, nargs=2, metavar=("LOWER", "UPPER"),
                   help="theta interval in degrees for theta_phi; default: [0, config tilt_deg]")
    p.add_argument("--phi-deg", type=float, nargs=2, metavar=("LOWER", "UPPER"),
                   help="required unwrapped phi interval in degrees for theta_phi, e.g. -45 45")
    p.add_argument("--orbits", nargs="+", help="subset of orbit_id values (default: all)")
    p.add_argument("--orbit-prefix", type=str, default=None,
                   help="only orbits whose orbit_id starts with this prefix (e.g. PL1)")
    p.add_argument("--phases", type=int, nargs="+",
                   help="subset of phase_index values (default: all)")
    p.add_argument("--solar-phases", type=float, nargs="+", default=[0.0],
                   help="solar initial phases alpha0 in radians (default: 0.0)")
    p.add_argument("--da-order", type=int, default=8)
    p.add_argument("--err-toll", type=float, default=1.0e-8)
    p.add_argument("--n-split-max", type=int, default=12)
    p.add_argument("--jobs", type=int, default=1)
    p.add_argument("--limit", type=int, default=0,
                   help="stop after N cases (0 = no limit); useful for smoke tests")
    p.add_argument("--skip-existing", action="store_true",
                   help="reuse finished per-case result JSONs with a matching attitude mode")
    p.add_argument("--outdir", type=Path, default=ROOT / "results" / "srp_sensitivity")
    p.add_argument("--executable", type=Path, help="path to run_crtbp_srp_study (skip fpm build)")
    p.add_argument("--dry-run", action="store_true", help="print commands without building, running, or writing output files")
    return p


def main(argv=None) -> int:
    args = parser().parse_args(argv)
    if args.da_order < 1:
        raise ValueError("--da-order must be positive")
    if not math.isfinite(args.err_toll) or args.err_toll <= 0:
        raise ValueError("--err-toll must be positive and finite")
    if args.n_split_max < 1:
        raise ValueError("--n-split-max must be positive")
    if args.jobs < 1:
        raise ValueError("--jobs must be positive")
    if any(not math.isfinite(a) for a in args.solar_phases):
        raise ValueError("--solar-phases must be finite")
    if args.beta is not None:
        beta_lo, beta_hi = args.beta
        if not (math.isfinite(beta_lo) and math.isfinite(beta_hi)):
            raise ValueError("--beta: values must be finite")
        if beta_lo < 0.0 or beta_hi <= beta_lo:
            raise ValueError("--beta: need 0 <= lower < upper")
        cfg_id = f"Q1_beta_{beta_lo:g}_{beta_hi:g}"
        CONFIGS[cfg_id] = dict(layer="L1", model="cannonball",
                               beta=(beta_lo, beta_hi), tilt_deg=0.0, chi=(0.6, 0.6))
        args.configs = [cfg_id]
    args.json = args.json.resolve()
    if not args.json.is_file():
        raise FileNotFoundError(args.json)

    orbits = load_orbits(args.json)
    selected_configs = select_configs(args)
    args.configs = list(selected_configs)
    cases = list(enumerate_cases(args, orbits, selected_configs))
    if not cases:
        raise ValueError("no cases matched the requested filters")

    if args.dry_run:
        # Command preview only: never build or launch Fortran during --dry-run.
        exe = args.executable.resolve() if args.executable else Path("run_crtbp_srp_study")
        for job in cases:
            _, orbit, phase, pi, solar_phase, cfg_id, cfg = job
            oid = orbit.get("orbit_id")
            out_path = args.outdir / "results" / cfg_id / oid / f"p{int(pi):02d}__a{fmt(solar_phase)}.json"
            print(shlex.join(build_command(exe, orbit, phase, pi, solar_phase, cfg_id, cfg, args, out_path)))
        return 0

    commit = git_commit()
    args.outdir.mkdir(parents=True, exist_ok=True)
    write_json(args.outdir / "manifest.json", {
        "schema": "crtbp_srp_study_batch_v1",
        "git_commit": commit,
        "json_source": str(args.json),
        "configs": list(args.configs),
        "config_definitions": selected_configs,
        "solar_phases_rad": args.solar_phases,
        "da_order": args.da_order, "err_toll": args.err_toll,
        "n_split_max": args.n_split_max,
        "orbit_count": len({c[1].get("orbit_id") for c in cases}),
        "case_count": len(cases),
    })

    env = environment()
    exe = resolve_executable(args, env)
    print(f"binary: {exe}", flush=True)
    print(f"cases: {len(cases)} (jobs={args.jobs})", flush=True)

    campaign = {"schema": "crtbp_srp_study_campaign_v1", "git_commit": commit,
                "config_definitions": selected_configs, "case_count": len(cases), "jobs": []}
    started = time.perf_counter()
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(run_one, job, args, exe, env) for job in cases]
        for future in concurrent.futures.as_completed(futures):
            try:
                campaign["jobs"].append(future.result())
            except Exception as exc:
                campaign["jobs"].append({"status": "driver_failed", "reason": str(exc)})
            write_json(args.outdir / "campaign.json", campaign)
    campaign["total_elapsed_s"] = time.perf_counter() - started
    write_json(args.outdir / "campaign.json", campaign)

    fields = ("config", "orbit", "family", "libration_point", "phase_index",
              "solar_phase_rad", "attitude_coordinates", "theta_deg", "phi_deg",
              "status", "elapsed_s", "wall_time_s", "app_status",
              "n_patches", "n_splits_total", "n_max_active", "first_split_time",
              "n_diverged", "first_divergence_time",
              "t_end_days", "has_beta_sens")
    with (args.outdir / "campaign_summary.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(campaign["jobs"])
    failed = [j for j in campaign["jobs"] if j.get("status") not in ("completed", "skipped_existing")]
    print(f"done: {len(campaign['jobs'])} cases, {len(failed)} not completed, "
          f"{campaign['total_elapsed_s']:.1f}s wall", flush=True)
    return 1 if failed else 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(2)
