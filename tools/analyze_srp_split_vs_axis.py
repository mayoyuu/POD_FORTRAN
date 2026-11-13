#!/usr/bin/env python3
"""Scatter plot of first-split time vs. endpoint position longest-axis.

For each completed case with a valid covariance (status == 0), reads
the event-detected ``first_split_time_accurate`` when available (otherwise the
checkpoint ``first_split_time``; nondimensional; -1 = never split), together
with ``cov_state``
(6x6 row-major), computes the endpoint position longest axis in km, and plots
one point per case: x = first split time (days), y = longest axis (km).

Diverged cases (status == -4) have no valid covariance and are excluded from
the scatter (their counts are still reported).

Usage:
    python3 tools/analyze_srp_split_vs_axis.py \
        --results-root results/srp_sensitivity/results \
        --configs Q1_L1_small Q1_L1_main Q1_L1_ext \
        --color-by config \
        --outdir results/srp_sensitivity/analysis
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

from srp_covariance_metrics import position_covariance_metrics
from srp_first_split_time import select_first_split_time

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RESULTS = ROOT / "results" / "srp_sensitivity" / "results"
TU_S = 375237.080249162
SEC_PER_DAY = 86400.0
ND_TO_DAYS = TU_S / SEC_PER_DAY

try:
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAS_MPL = True
except Exception:  # pragma: no cover - optional dependency
    HAS_MPL = False


def parse_amplitude(orbit_id: str) -> int:
    return int(orbit_id.split("_A")[-1])


def family_of(orbit_id: str) -> str:
    return orbit_id.split("_A")[0]


def axis_km(cov_state: list) -> float:
    """Compatibility alias for the true longest one-sigma position axis."""
    return position_covariance_metrics(cov_state).longest_axis_km


def scan_results(root: Path) -> list:
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
                oid = orbit_dir.name
                jobs.append({
                    "config": cfg_dir.name,
                    "orbit": oid,
                    "family": family_of(oid),
                    "amplitude": parse_amplitude(oid),
                    "phase_index": int(f.stem.split("__")[0][1:]),
                    "status": res.get("status"),
                    "first_split_time": res.get("first_split_time"),
                    "first_split_time_accurate": res.get("first_split_time_accurate"),
                    "cov_state": res.get("cov_state"),
                })
    return jobs


def valid_points(jobs):
    """(first_split_days, axis_km, color, row) for status==0 cases with cov."""
    pts = []
    n_div = 0
    n_nocov = 0
    for j in jobs:
        if j["status"] is not None and j["status"] < 0:
            n_div += 1
            continue
        cov = j["cov_state"]
        if not isinstance(cov, list) or len(cov) != 36:
            n_nocov += 1
            continue
        fst, _ = select_first_split_time(j)
        if fst is None:
            continue
        pts.append((fst * ND_TO_DAYS, axis_km(cov), j))
    return pts, n_div, n_nocov


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
def _color_key(j, color_by):
    return j["config"] if color_by == "config" else j["family"]


def png_scatter(points, color_by, outdir, tag):
    groups = {}
    for x, y, j in points:
        groups.setdefault(_color_key(j, color_by), []).append((x, y))
    fig, ax = plt.subplots(figsize=(11, 7))
    keys = sorted(groups)
    cmap = plt.get_cmap("tab20")
    for i, k in enumerate(keys):
        xs = [p[0] for p in groups[k]]
        ys = [p[1] for p in groups[k]]
        ax.scatter(xs, ys, s=14, alpha=0.55, color=cmap(i % 20), label=k)
    ax.set_xlabel("selected first split time (event preferred, days)")
    ax.set_ylabel("endpoint position longest axis (km)")
    ax.set_title(f"selected first split time vs. longest axis (color = {color_by})")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, ncol=max(1, len(keys) // 6), loc="best")
    fig.tight_layout()
    path = outdir / f"split_vs_axis_{tag}_{color_by}.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  scatter plot -> {path}")


def ascii_report(points, n_div, n_nocov, stream):
    print(f"\n=== first split time vs. longest axis ===", file=stream)
    print(f"status==0 cases plotted: {len(points)}; "
          f"diverged (excluded): {n_div}; missing cov: {n_nocov}", file=stream)
    if not points:
        return
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    print(f"first split time (days): min={min(xs):.2f} max={max(xs):.2f} "
          f"mean={sum(xs)/len(xs):.2f}", file=stream)
    print(f"longest axis (km): min={min(ys):.3g} max={max(ys):.3g} "
          f"mean={sum(ys)/len(ys):.3g}", file=stream)
    # crude correlation
    mx = sum(xs) / len(xs); my = sum(ys) / len(ys)
    num = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    dx = sum((a - mx) ** 2 for a in xs)
    dy = sum((b - my) ** 2 for b in ys)
    if dx > 0 and dy > 0:
        print(f"Pearson r = {num / math.sqrt(dx * dy):.3f}", file=stream)


def write_csv(points, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    header = ("config", "family", "orbit", "amplitude_km", "phase_index",
              "first_split_days", "first_split_source",
              "first_split_checkpoint_nd", "first_split_accurate_nd", "axis_km")
    rows = []
    for x, y, j in sorted(points, key=lambda p: (p[2]["config"], p[2]["orbit"],
                                                 p[2]["phase_index"])):
        _, source = select_first_split_time(j)
        rows.append({"config": j["config"], "family": j["family"], "orbit": j["orbit"],
                     "amplitude_km": j["amplitude"], "phase_index": j["phase_index"],
                     "first_split_days": f"{x:.6g}", "first_split_source": source,
                     "first_split_checkpoint_nd": j.get("first_split_time"),
                     "first_split_accurate_nd": j.get("first_split_time_accurate"),
                     "axis_km": f"{y:.6g}"})
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=header, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def parser():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--results-root", type=Path, default=DEFAULT_RESULTS)
    p.add_argument("--orbit-prefix", type=str, default=None)
    p.add_argument("--configs", nargs="+", default=[])
    p.add_argument("--color-by", choices=["config", "family"], default="config")
    p.add_argument("--outdir", type=Path,
                   default=ROOT / "results" / "srp_sensitivity" / "analysis")
    p.add_argument("--no-png", action="store_true")
    return p


def main(argv=None) -> int:
    args = parser().parse_args(argv)
    jobs = scan_results(args.results_root)
    if not jobs:
        print(f"no result JSONs under {args.results_root}", file=sys.stderr)
        return 1

    sel = [j for j in jobs
           if (not args.orbit_prefix or j["orbit"].startswith(args.orbit_prefix))
           and (not args.configs or j["config"] in args.configs)]
    if not sel:
        print("no matching cases", file=sys.stderr)
        return 1

    points, n_div, n_nocov = valid_points(sel)
    args.outdir.mkdir(parents=True, exist_ok=True)
    tag = args.orbit_prefix or "all"

    csv_path = args.outdir / f"split_vs_axis_{tag}.csv"
    write_csv(points, csv_path)
    print(f"per-case table -> {csv_path}")

    if not args.no_png and HAS_MPL:
        png_scatter(points, args.color_by, args.outdir, tag)
    else:
        ascii_report(points, n_div, n_nocov, sys.stdout)
        if not args.no_png and not HAS_MPL:
            print("(matplotlib not installed; install for PNG scatter)",
                  file=sys.stderr)

    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(2)
