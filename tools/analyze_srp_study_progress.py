#!/usr/bin/env python3
"""Progress analysis for the CRTBP + SRP sensitivity study (partial results).

The batch runner (run_crtbp_srp_study_batch.py) writes campaign.json *incrementally*
as each case finishes, so it is a reliable snapshot of partial progress. This script
reads that snapshot and, for one orbit family (default Halo North L1, orbit_id prefix
``H1N`` = 16 amplitudes x 24 phases = 384 cases per config), produces per-config:

  * heatmaps of ``first_split_time`` (days) and ``n_splits_total`` over amplitude x phase
  * a per-case table of split statistics (total splits + per-axis split counts)
  * a per-config split-direction summary (beta/chi/u/v or theta/phi totals)

If matplotlib + numpy are available the heatmap is saved as PNG; otherwise a plain
ASCII grid is printed instead, so the script always works on a bare server.

Usage:
    python3 tools/analyze_srp_study_progress.py \
        --campaign results/srp_sensitivity/campaign.json \
        --orbit-prefix H1N \
        --outdir results/srp_sensitivity/analysis

    # if campaign.json is missing, scan the raw result JSON files instead:
    python3 tools/analyze_srp_study_progress.py \
        --results-root results/srp_sensitivity/results \
        --orbit-prefix H1N
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTDIR = ROOT / "results" / "srp_sensitivity"

TU_S = 375237.080249162
SEC_PER_DAY = 86400.0
ND_TO_DAYS = TU_S / SEC_PER_DAY
N_PHASE = 24

try:
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAS_MPL = True
except Exception:  # pragma: no cover - optional dependency
    HAS_MPL = False


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
def parse_amplitude(orbit_id: str) -> int:
    """Extract the amplitude (km) encoded as the ``_A<digits>`` suffix of an orbit_id."""
    tag = orbit_id.split("_A")[-1]
    return int(tag)


def load_campaign(path: Path) -> list:
    if not path.is_file():
        return []
    data = json.loads(path.read_text(encoding="utf-8-sig"))
    jobs = data.get("jobs")
    return jobs if isinstance(jobs, list) else []


def scan_results(root: Path) -> list:
    """Fallback: read each per-case result JSON under root/<config>/<orbit>/p*.json."""
    jobs = []
    if not root.is_dir():
        return jobs
    for cfg_dir in sorted(root.iterdir()):
        if not cfg_dir.is_dir():
            continue
        for orbit_dir in sorted(cfg_dir.iterdir()):
            if not orbit_dir.is_dir():
                continue
            for f in sorted(orbit_dir.glob("p*.json")):
                try:
                    res = json.loads(f.read_text(encoding="utf-8-sig"))
                except (OSError, ValueError) as exc:
                    print(f"WARNING: skipping unreadable result {f}: {exc}", file=sys.stderr)
                    continue
                # A legacy executable can leave UV output in an angle-tagged
                # directory even though the batch marks configuration_mismatch.
                if "__theta_phi_" in cfg_dir.name and (
                    res.get("attitude_coordinates") != "theta_phi"
                    or res.get("parameter_names") != ["beta", "chi", "theta", "phi"]
                ):
                    print(f"WARNING: configuration mismatch; skipping {f}", file=sys.stderr)
                    continue
                jobs.append({
                    "config": cfg_dir.name,
                    "orbit": orbit_dir.name,
                    "phase_index": int(f.stem.split("__")[0][1:]),
                    "solar_phase_rad": res.get("alpha0_rad"),
                    "parameter_names": res.get("parameter_names"),
                    "attitude_coordinates": res.get("attitude_coordinates", "uv"),
                    "t_end_days": res.get("t_end_days"),
                    "result_path": str(f),
                    "first_split_time": res.get("first_split_time"),
                    "n_splits_total": res.get("n_splits_total"),
                    "n_patches": res.get("n_patches"),
                    "n_max_active": res.get("n_max_active"),
                    "n_diverged": res.get("n_diverged"),
                    "split_count_dir": res.get("split_count_dir"),
                    "wall_time_s": res.get("wall_time_s"),
                    "app_status": res.get("status"),
                    "status": "completed",
                })
    return jobs


def select_jobs(jobs: list, orbit_prefix: str, configs: list, solar_phase=None) -> list:
    out = []
    for job in jobs:
        if orbit_prefix and not str(job.get("orbit", "")).startswith(orbit_prefix):
            continue
        if configs and job.get("config") not in configs:
            continue
        if solar_phase is not None:
            phase = job.get("solar_phase_rad")
            if phase is None or not math.isclose(float(phase), solar_phase, rel_tol=0, abs_tol=1e-10):
                continue
        if job.get("status") not in ("completed", "skipped_existing"):
            continue
        try:
            parse_amplitude(str(job["orbit"]))
            int(job.get("phase_index"))
        except (KeyError, TypeError, ValueError):
            continue
        out.append(job)
    return out


# ---------------------------------------------------------------------------
# Grid assembly
# ---------------------------------------------------------------------------
def sorted_amplitudes(jobs: list) -> list:
    return sorted({parse_amplitude(str(j["orbit"])) for j in jobs})


def sorted_phases(jobs: list) -> list:
    return sorted({int(j["phase_index"]) for j in jobs})


def build_grid(jobs: list, key: str) -> dict:
    """Return {amp: {phase: value}} for the given metric key."""
    amps = sorted_amplitudes(jobs)
    phases = sorted_phases(jobs)
    grid = {a: {p: None for p in phases} for a in amps}
    for job in jobs:
        amp = parse_amplitude(str(job["orbit"]))
        ph = int(job["phase_index"])
        val = job.get(key)
        if val is None:
            continue
        grid[amp][ph] = val
    return grid


def build_dir_grid(jobs: list, k: int) -> dict:
    """Return {amp: {phase: split_count_dir[k]}} for direction k (0=beta,1=chi,2=u,3=v)."""
    amps = sorted_amplitudes(jobs)
    phases = sorted_phases(jobs)
    grid = {a: {p: None for p in phases} for a in amps}
    for job in jobs:
        amp = parse_amplitude(str(job["orbit"]))
        ph = int(job["phase_index"])
        scd = job.get("split_count_dir") or [0, 0, 0, 0]
        grid[amp][ph] = int(scd[k]) if len(scd) > k else 0
    return grid


def first_split_days(job: dict):
    v = job.get("first_split_time")
    if v is None:
        return None
    v = float(v)
    if v < 0.0:
        return None  # -1 => no split over the whole arc
    return v * ND_TO_DAYS


# ---------------------------------------------------------------------------
# Output: per-case statistics table
# ---------------------------------------------------------------------------
def parameter_names(job):
    """Physical slot names; packed DA indices are not the output column order."""
    names = job.get("parameter_names")
    if names:
        return list(names)
    return ["beta", "chi", "theta", "phi"] if job.get("attitude_coordinates") == "theta_phi" else ["beta", "chi", "u", "v"]


def direction_counts(job):
    counts = job.get("split_count_dir") or [0, 0, 0, 0]
    return dict(zip(parameter_names(job), counts))


def write_stats_csv(jobs: list, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    header = ("config", "orbit", "amplitude_km", "phase_index", "phase_fraction",
              "first_split_time_nd", "first_split_time_days",
              "n_splits_total", "n_patches", "n_max_active", "n_diverged",
              "split_beta", "split_chi", "split_u", "split_v", "split_theta", "split_phi",
              "solar_phase_rad", "t_end_days", "result_path",
              "app_status", "wall_time_s")
    rows = []
    for job in sorted(jobs, key=lambda j: (j["config"], parse_amplitude(str(j["orbit"])),
                                           int(j["phase_index"]))):
        rows.append({
            "config": job.get("config"),
            "orbit": job.get("orbit"),
            "amplitude_km": parse_amplitude(str(job["orbit"])),
            "phase_index": int(job["phase_index"]),
            "phase_fraction": job.get("phase_fraction", ""),
            "first_split_time_nd": _fmt_float(job.get("first_split_time")),
            "first_split_time_days": _fmt_float(first_split_days(job)),
            "n_splits_total": job.get("n_splits_total", ""),
            "n_patches": job.get("n_patches", ""),
            "n_max_active": job.get("n_max_active", ""),
            "n_diverged": job.get("n_diverged", ""),
            **{f"split_{name}": direction_counts(job).get(name, "")
               for name in ("beta", "chi", "u", "v", "theta", "phi")},
            "solar_phase_rad": job.get("solar_phase_rad", ""),
            "t_end_days": job.get("t_end_days", ""),
            "result_path": job.get("result_path", ""),
            "app_status": job.get("app_status", ""),
            "wall_time_s": _fmt_float(job.get("wall_time_s")),
        })
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=header, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def _fmt_float(v):
    if v is None or v == "":
        return ""
    try:
        return f"{float(v):.6g}"
    except (TypeError, ValueError):
        return v


# ---------------------------------------------------------------------------
# Output: heatmap
# ---------------------------------------------------------------------------
def png_heatmap(amps, phases, matrix, title, path, fmt="{:.1f}", colorbar_label="first_split_time (days)"):
    fig, ax = plt.subplots(figsize=(max(10, 0.45 * len(phases)), max(6, 0.35 * len(amps))))
    arr = np.array(matrix, dtype=float)
    masked = np.ma.masked_invalid(arr)
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad(color="0.85")
    im = ax.imshow(masked, aspect="auto", origin="lower", cmap=cmap)
    ax.set_xticks(range(len(phases)))
    ax.set_xticklabels([str(p) for p in phases], fontsize=7)
    ax.set_yticks(range(len(amps)))
    ax.set_yticklabels([f"{a/1000:.0f}k" for a in amps], fontsize=7)
    ax.set_xlabel("phase_index")
    ax.set_ylabel("amplitude (x1000 km)")
    ax.set_title(title)
    for i in range(len(amps)):
        for j in range(len(phases)):
            v = arr[i][j]
            if np.isfinite(v):
                ax.text(j, i, fmt.format(v), ha="center", va="center", fontsize=5.5)
    fig.colorbar(im, ax=ax, label=colorbar_label)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def ascii_heatmap(amps, phases, str_matrix, title, legend, stream):
    width = 7
    print(f"\n=== {title} ===", file=stream)
    header = "amp\\ph " + "".join(f"{p:>{width}}" for p in phases)
    print(header, file=stream)
    print("-" * len(header), file=stream)
    for i, a in enumerate(amps):
        row = f"{a/1000:>5.1f}k "
        for j in range(len(phases)):
            row += str_matrix[i][j]
        print(row, file=stream)
    print(legend, file=stream)


def _emit_heatmap(metric, title, amps, phases, float_matrix, str_matrix, cfg, outdir,
                  png, fmt, legend):
    if png and HAS_MPL:
        path = outdir / f"heatmap_{metric}_{cfg}.png"
        png_heatmap(amps, phases, float_matrix, title, path, fmt=fmt,
                    colorbar_label="first_split_time (days)" if metric == "first_split_time" else "split count")
        print(f"  {metric} heatmap -> {path}")
    else:
        ascii_heatmap(amps, phases, str_matrix, title, legend, sys.stdout)
        if png and not HAS_MPL:
            print("  (matplotlib not installed; install it to get PNG heatmaps)",
                  file=sys.stderr)


def render_direction_summary(jobs) -> dict:
    cfg = jobs[0]["config"]
    total = 0
    dirs = [0, 0, 0, 0]
    n_div = 0
    for j in jobs:
        scd = j.get("split_count_dir") or [0, 0, 0, 0]
        for k in range(4):
            dirs[k] += int(scd[k])
            total += int(scd[k])
        n_div += int(j.get("n_diverged") or 0)
    names = parameter_names(jobs[0])
    parts = "  ".join(f"{n}={dirs[k]}" for k, n in enumerate(names))
    print(f"  split directions (total={total} splits across {len(jobs)} cases): {parts}; "
          f"diverged={n_div}")
    return {"config": cfg, "n_cases": len(jobs), "total_splits": total,
            **{f"split_{name}": dict(zip(names, dirs)).get(name, "")
               for name in ("beta", "chi", "u", "v", "theta", "phi")},
            "solar_phase_rad": jobs[0].get("solar_phase_rad", ""), "n_diverged": n_div}


def render_config(jobs, outdir: Path, png: bool) -> dict:
    cfg = jobs[0]["config"]
    if any(parameter_names(j) != parameter_names(jobs[0]) for j in jobs):
        raise ValueError(f"mixed attitude parameterizations in config {cfg}")
    cells = [(j["orbit"], j["phase_index"]) for j in jobs]
    if len(set(cells)) != len(cells):
        raise ValueError("duplicate orbit/phase cells; select one --solar-phase")
    amps = sorted_amplitudes(jobs)
    phases = sorted_phases(jobs)
    n_total = len(amps) * len(phases)
    n_done = len(jobs)

    print(f"\nconfig {cfg}: {n_done}/{n_total} cases completed "
          f"(amplitudes {amps[0]}..{amps[-1]} km, phases {phases[0]}..{phases[-1]})")

    # ---- first_split_time (days) ----
    fs_grid = build_grid(jobs, "first_split_time")
    fs_float, fs_str = [], []
    n_no_split = 0
    for a in amps:
        frow, srow = [], []
        for p in phases:
            v = fs_grid[a][p]
            if v is None:
                frow.append(float("nan")); srow.append("  .   ")
            else:
                v = float(v)
                if v < 0.0:
                    frow.append(float("nan")); srow.append("  inf  "); n_no_split += 1
                else:
                    d = v * ND_TO_DAYS
                    frow.append(d); srow.append(f"{d:7.2f}")
        fs_float.append(frow); fs_str.append(srow)
    _emit_heatmap("first_split_time",
                  f"{cfg}: first_split_time (days) [{n_done}/{n_total} completed, "
                  f"{n_no_split} never split]",
                  amps, phases, fs_float, fs_str, cfg, outdir, png, "{:.1f}",
                  "(. = not completed; inf = completed but never split; else days)")

    # ---- n_splits_total ----
    ns_grid = build_grid(jobs, "n_splits_total")
    ns_float, ns_str = [], []
    for a in amps:
        frow, srow = [], []
        for p in phases:
            v = ns_grid[a][p]
            if v is None:
                frow.append(float("nan")); srow.append("  .   ")
            else:
                n = int(v)
                frow.append(float(n)); srow.append(f"{n:7d}")
        ns_float.append(frow); ns_str.append(srow)
    _emit_heatmap("n_splits_total",
                  f"{cfg}: n_splits_total [{n_done}/{n_total} completed]",
                  amps, phases, ns_float, ns_str, cfg, outdir, png, "{:.0f}",
                  "(. = not completed; else n_splits_total)")

    # ---- per-direction split counts (beta/chi/u/v) ----
    dir_names = parameter_names(jobs[0])
    for k, dname in enumerate(dir_names):
        dgrid = build_dir_grid(jobs, k)
        d_float, d_str = [], []
        for a in amps:
            frow, srow = [], []
            for p in phases:
                v = dgrid[a][p]
                if v is None:
                    frow.append(float("nan")); srow.append("  .   ")
                else:
                    n = int(v)
                    frow.append(float(n)); srow.append(f"{n:7d}")
            d_float.append(frow); d_str.append(srow)
        _emit_heatmap(f"split_{dname}",
                      f"{cfg}: splits along '{dname}' [{n_done}/{n_total} completed]",
                      amps, phases, d_float, d_str, cfg, outdir, png, "{:.0f}",
                      "(. = not completed; else splits along this axis)")

    return render_direction_summary(jobs)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def parser():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--campaign", type=Path,
                   default=None,
                   help="fallback campaign.json; default: next to the selected results directory")
    p.add_argument("--results-root", type=Path,
                   default=DEFAULT_OUTDIR / "results",
                   help="fallback: scan raw result JSONs under this root if campaign.json is empty")
    p.add_argument("--orbit-prefix", type=str, default="H1N",
                   help="orbit_id prefix to analyze (default: H1N = Halo North L1)")
    p.add_argument("--configs", nargs="+", default=[],
                   help="restrict to these config ids (default: all configs present)")
    p.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR / "analysis",
                   help="where to write CSVs and PNGs")
    p.add_argument("--solar-phase", type=float, default=None,
                   help="select a single solar initial phase in radians (not attitude phi)")
    p.add_argument("--no-png", action="store_true", help="force ASCII output even if matplotlib is available")
    return p


def main(argv=None) -> int:
    args = parser().parse_args(argv)

    # Prefer the raw result files: campaign.json is overwritten by every batch run,
    # so it only reflects the *last* run (e.g. a partial Q2), not all completed work.
    # Accept either the campaign directory or its results/ child.
    result_root = args.results_root
    if (result_root / "results").is_dir():
        result_root = result_root / "results"
    jobs = scan_results(result_root)
    if not jobs:
        campaign = args.campaign or result_root.parent / "campaign.json"
        print(f"no result files under {result_root}; falling back to "
              f"campaign.json ({campaign})", file=sys.stderr)
        jobs = load_campaign(campaign)
    if not jobs:
        print("no completed cases found", file=sys.stderr)
        return 1

    sel = select_jobs(jobs, args.orbit_prefix, args.configs, args.solar_phase)
    if not sel:
        print(f"no completed cases for orbit prefix '{args.orbit_prefix}'", file=sys.stderr)
        return 1

    configs = sorted({j["config"] for j in sel})
    print(f"orbit prefix: {args.orbit_prefix}  configs: {configs}  "
          f"completed cases: {len(sel)}")

    args.outdir.mkdir(parents=True, exist_ok=True)

    # Per-case statistics table (all configs combined)
    stats_path = args.outdir / f"split_stats_{args.orbit_prefix}.csv"
    write_stats_csv(sel, stats_path)
    print(f"split statistics table -> {stats_path}")

    # Per-config heatmap + direction summary
    dir_summaries = []
    for cfg in configs:
        cfg_jobs = [j for j in sel if j["config"] == cfg]
        # A solar phase is an independent experiment, never the same heatmap cell.
        groups = {}
        for job in cfg_jobs:
            key = (str(job["orbit"]).split("_A")[0], job.get("solar_phase_rad"))
            groups.setdefault(key, []).append(job)
        for (family, solar), group in groups.items():
            outdir = args.outdir
            if len(groups) > 1:
                outdir = outdir / family / f"alpha0_{solar}"
            outdir.mkdir(parents=True, exist_ok=True)
            summary = render_config(group, outdir, png=not args.no_png)
            summary["family"] = family
            dir_summaries.append(summary)

    dir_path = args.outdir / f"direction_summary_{args.orbit_prefix}.csv"
    with dir_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(dir_summaries[0].keys()))
        w.writeheader()
        w.writerows(dir_summaries)
    print(f"\nsplit direction summary -> {dir_path}")

    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(2)
