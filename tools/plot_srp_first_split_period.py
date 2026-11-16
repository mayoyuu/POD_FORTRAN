#!/usr/bin/env python3
"""First-split-time heatmap in units of orbital period for the SRP study.

Companion to ``analyze_srp_study_progress.py``: the same first-split statistic
(``first_split_time_accurate`` preferred, then the checkpoint ``first_split_time``),
but expressed in *orbital periods* (``first_split_time_nd / period_nd``) instead
of days. A value of ``1.3`` means the first split occurred ~1.3 orbital periods
into the propagation.

Heatmap layout: x = phase_index, y = amplitude. Completed cases that never
split are shown as ``inf`` (ASCII) or masked (PNG); missing/incomplete cases as
``.`` / masked.

Usage:
    python3 tools/plot_srp_first_split_period.py \
        --results-root results/srp_sensitivity/results \
        --orbit-prefix H1N \
        --outdir results/srp_sensitivity/analysis
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

from srp_first_split_time import select_first_split_time

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RESULTS = ROOT / "results" / "srp_sensitivity" / "results"
DEFAULT_OUTDIR = ROOT / "results" / "srp_sensitivity" / "analysis"

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
    return int(orbit_id.split("_A")[-1])


def _positive_float(v):
    if v is None or isinstance(v, bool):
        return None
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    return x if (math.isfinite(x) and x > 0.0) else None


def scan_results(root: Path) -> list:
    """Walk root/<config>/<orbit>/p*.json and return first-split fields."""
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
                    "solar_phase_rad": res.get("alpha0_rad"),
                    "period_nd": res.get("period_nd"),
                    "first_split_time": res.get("first_split_time"),
                    "first_split_time_accurate": res.get("first_split_time_accurate"),
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
            sp = job.get("solar_phase_rad")
            if sp is None or not math.isclose(float(sp), solar_phase,
                                             rel_tol=0, abs_tol=1e-10):
                continue
        try:
            parse_amplitude(str(job["orbit"]))
            int(job.get("phase_index"))
        except (KeyError, TypeError, ValueError):
            continue
        out.append(job)
    return out


def sorted_amplitudes(jobs: list) -> list:
    return sorted({parse_amplitude(str(j["orbit"])) for j in jobs})


def sorted_phases(jobs: list) -> list:
    return sorted({int(j["phase_index"]) for j in jobs})


def build_period_grid(jobs: list) -> dict:
    """Return {amplitude: {phase: first_split_periods}}; -1.0 = never split."""
    amps = sorted_amplitudes(jobs)
    phases = sorted_phases(jobs)
    grid = {a: {p: None for p in phases} for a in amps}
    for job in jobs:
        amp = parse_amplitude(str(job["orbit"]))
        phase = int(job["phase_index"])
        first_split, source = select_first_split_time(job)
        if source == "never":
            grid[amp][phase] = -1.0
        elif source == "missing" or first_split is None:
            grid[amp][phase] = None
        else:
            period = _positive_float(job.get("period_nd"))
            grid[amp][phase] = None if period is None else first_split / period
    return grid


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------
def png_heatmap(amps, phases, matrix, title, path, colorbar_label):
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
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=5.5)
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


def render_config(jobs, outdir: Path, png: bool) -> None:
    cfg = jobs[0]["config"]
    cells = [(j["orbit"], j["phase_index"]) for j in jobs]
    if len(set(cells)) != len(cells):
        raise ValueError("duplicate orbit/phase cells; select one --solar-phase")
    amps = sorted_amplitudes(jobs)
    phases = sorted_phases(jobs)
    n_total = len(amps) * len(phases)
    n_done = len(jobs)

    grid = build_period_grid(jobs)
    fmat, smat = [], []
    n_no_split = 0
    for a in amps:
        frow, srow = [], []
        for p in phases:
            v = grid[a][p]
            if v is None:
                frow.append(float("nan")); srow.append("  .   ")
            elif v < 0.0:
                frow.append(float("nan")); srow.append("  inf  "); n_no_split += 1
            else:
                frow.append(v); srow.append(f"{v:7.2f}")
        fmat.append(frow); smat.append(srow)

    title = (f"{cfg}: first split time, event preferred (periods) "
             f"[{n_done}/{n_total} completed, {n_no_split} never split]")
    print(f"\nconfig {cfg}: {n_done}/{n_total} cases completed "
          f"(amplitudes {amps[0]}..{amps[-1]} km)")

    if png and HAS_MPL:
        path = outdir / f"heatmap_first_split_period_{cfg}.png"
        png_heatmap(amps, phases, fmat, title, path, "first split time (periods)")
        print(f"  first split period heatmap -> {path}")
    else:
        ascii_heatmap(amps, phases, smat, title,
                      "(. = not completed; inf = completed but never split; else periods)",
                      sys.stdout)
        if png and not HAS_MPL:
            print("  (matplotlib not installed; install for PNG heatmaps)",
                  file=sys.stderr)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def parser():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--results-root", type=Path, default=DEFAULT_RESULTS,
                   help="root of per-case result JSONs (default: results/srp_sensitivity/results)")
    p.add_argument("--orbit-prefix", type=str, default="H1N",
                   help="orbit_id prefix to analyze (default: H1N)")
    p.add_argument("--configs", nargs="+", default=[],
                   help="restrict to these config ids (default: all present)")
    p.add_argument("--solar-phase", type=float, default=None,
                   help="select a single solar initial phase in radians (default: all)")
    p.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR,
                   help="where to write PNG heatmaps (default: results/srp_sensitivity/analysis)")
    p.add_argument("--no-png", action="store_true",
                   help="force ASCII output even if matplotlib is available")
    return p


def main(argv=None) -> int:
    args = parser().parse_args(argv)
    jobs = scan_results(args.results_root)
    if not jobs:
        print(f"no result JSONs found under {args.results_root}", file=sys.stderr)
        return 1

    sel = select_jobs(jobs, args.orbit_prefix, args.configs, args.solar_phase)
    if not sel:
        print(f"no cases for orbit prefix '{args.orbit_prefix}'", file=sys.stderr)
        return 1

    configs = sorted({j["config"] for j in sel})
    print(f"orbit prefix: {args.orbit_prefix}  configs: {configs}  cases: {len(sel)}")

    args.outdir.mkdir(parents=True, exist_ok=True)
    for cfg in configs:
        render_config([j for j in sel if j["config"] == cfg], args.outdir,
                      png=not args.no_png)

    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(2)
