#!/usr/bin/env python3
"""Progress analysis for the CRTBP + SRP sensitivity study (partial results).

The batch runner (run_crtbp_srp_study_batch.py) writes campaign.json *incrementally*
as each case finishes, so it is a reliable snapshot of partial progress. This script
reads that snapshot and, for one orbit family (default Halo North L1, orbit_id prefix
``H1N`` = 16 amplitudes x 24 phases = 384 cases per config), produces per-config:

  * a heatmap of ``first_split_time`` (days) over amplitude x phase
  * a per-case table of split statistics (total splits + per-axis split counts)

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
                except Exception:
                    continue
                jobs.append({
                    "config": cfg_dir.name,
                    "orbit": orbit_dir.name,
                    "phase_index": int(f.stem.split("__")[0][1:]),
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


def select_jobs(jobs: list, orbit_prefix: str, configs: list) -> list:
    out = []
    for job in jobs:
        if orbit_prefix and not str(job.get("orbit", "")).startswith(orbit_prefix):
            continue
        if configs and job.get("config") not in configs:
            continue
        if job.get("status") != "completed":
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
def write_stats_csv(jobs: list, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    header = ("config", "orbit", "amplitude_km", "phase_index", "phase_fraction",
              "first_split_time_nd", "first_split_time_days",
              "n_splits_total", "n_patches", "n_max_active", "n_diverged",
              "split_beta", "split_chi", "split_u", "split_v",
              "app_status", "wall_time_s")
    rows = []
    for job in sorted(jobs, key=lambda j: (j["config"], parse_amplitude(str(j["orbit"])),
                                           int(j["phase_index"]))):
        scd = job.get("split_count_dir") or [0, 0, 0, 0]
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
            "split_beta": scd[0], "split_chi": scd[1], "split_u": scd[2], "split_v": scd[3],
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
def png_heatmap(amps, phases, matrix, title, path, fmt="{:.1f}"):
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
    fig.colorbar(im, ax=ax, label="first_split_time (days)")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def ascii_heatmap(amps, phases, raw_grid, title, stream):
    width = 7
    print(f"\n=== {title} ===", file=stream)
    header = "amp\\ph " + "".join(f"{p:>{width}}" for p in phases)
    print(header, file=stream)
    print("-" * len(header), file=stream)
    for a in amps:
        row = f"{a/1000:>5.1f}k "
        for p in phases:
            v = raw_grid[a][p]
            if v is None:
                cell = "  .   "  # missing / not yet completed
            elif float(v) < 0.0:
                cell = "  inf  "  # completed, never split
            else:
                cell = f"{float(v)*ND_TO_DAYS:{width}.2f}"
            row += cell
        print(row, file=stream)
    print("(. = not completed; inf = completed but never split; else first_split_time in days)",
          file=stream)


def render_config(jobs, outdir: Path, png: bool) -> None:
    cfg = jobs[0]["config"]
    amps = sorted_amplitudes(jobs)
    phases = sorted_phases(jobs)
    n_total = len(amps) * len(phases)

    fs_grid = build_grid(jobs, "first_split_time")
    matrix = []
    n_no_split = 0
    for a in amps:
        row = []
        for p in phases:
            v = fs_grid[a][p]
            if v is None:
                row.append(float("nan"))
                continue
            v = float(v)
            if v < 0.0:
                row.append(float("nan"))
                n_no_split += 1
            else:
                row.append(v * ND_TO_DAYS)
        matrix.append(row)

    n_done = sum(1 for j in jobs)
    title = f"{cfg}: first_split_time (days)  [{n_done}/{n_total} completed]"
    print(f"\nconfig {cfg}: {n_done}/{n_total} cases completed "
          f"(amplitudes {amps[0]}..{amps[-1]} km, phases {phases[0]}..{phases[-1]}, "
          f"{n_no_split} never split)")

    if png and HAS_MPL:
        png_path = outdir / f"heatmap_first_split_{cfg}.png"
        png_heatmap(amps, phases, matrix, title, png_path)
        print(f"  heatmap -> {png_path}")
    else:
        ascii_heatmap(amps, phases, fs_grid, title, sys.stdout)
        if png and not HAS_MPL:
            print("  (matplotlib not installed; install it to get PNG heatmaps)",
                  file=sys.stderr)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def parser():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--campaign", type=Path,
                   default=DEFAULT_OUTDIR / "campaign.json",
                   help="path to the incremental campaign.json (default: results/srp_sensitivity/campaign.json)")
    p.add_argument("--results-root", type=Path,
                   default=DEFAULT_OUTDIR / "results",
                   help="fallback: scan raw result JSONs under this root if campaign.json is empty")
    p.add_argument("--orbit-prefix", type=str, default="H1N",
                   help="orbit_id prefix to analyze (default: H1N = Halo North L1)")
    p.add_argument("--configs", nargs="+", default=[],
                   help="restrict to these config ids (default: all configs present)")
    p.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR / "analysis",
                   help="where to write CSVs and PNGs")
    p.add_argument("--no-png", action="store_true", help="force ASCII output even if matplotlib is available")
    return p


def main(argv=None) -> int:
    args = parser().parse_args(argv)

    jobs = load_campaign(args.campaign)
    if not jobs:
        print(f"campaign.json empty/missing ({args.campaign}); scanning result files under "
              f"{args.results_root}", file=sys.stderr)
        jobs = scan_results(args.results_root)
    if not jobs:
        print("no completed cases found", file=sys.stderr)
        return 1

    sel = select_jobs(jobs, args.orbit_prefix, args.configs)
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

    # Per-config heatmap
    for cfg in configs:
        render_config([j for j in sel if j["config"] == cfg],
                      args.outdir, png=not args.no_png)

    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(2)
