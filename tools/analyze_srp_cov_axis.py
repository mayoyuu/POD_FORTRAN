#!/usr/bin/env python3
"""Endpoint position-covariance principal-axis analysis for the SRP study.

The per-case result JSON saves only the endpoint covariance (``cov_state``), not
an intermediate-time covariance series. This tool computes all three one-sigma
position principal axes, total RMS spread, and the longest-axis direction.

Outputs, for a chosen orbit family (default H1N):
  * cov_axis_<prefix>.csv          -- per-case axes, total RMS, and direction
  * cov_axis_summary_<prefix>.csv  -- longest-axis summary over phases
  * selected metric heatmap over amplitude x phase per config (PNG or ASCII)

The default ``--metric longest`` preserves the original heatmap name and behavior.

Usage:
    python3 tools/analyze_srp_cov_axis.py \
        --results-root results/srp_sensitivity/results \
        --orbit-prefix H1N \
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

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RESULTS = ROOT / "results" / "srp_sensitivity" / "results"

try:
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAS_MPL = True
except Exception:  # pragma: no cover - optional dependency
    HAS_MPL = False


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def parse_amplitude(orbit_id: str) -> int:
    return int(orbit_id.split("_A")[-1])


def axis_km(cov_state: list) -> float:
    """Compatibility alias for the true longest one-sigma position axis."""
    return position_covariance_metrics(cov_state).longest_axis_km


def scan_results(root: Path) -> list:
    """Walk root/<config>/<orbit>/p*.json and return completed cases with cov_state."""
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
                cov = res.get("cov_state")
                if not isinstance(cov, list) or len(cov) != 36:
                    continue
                if res.get("status") != 0:
                    continue  # statistics failed; cov is uninitialized
                metrics = position_covariance_metrics(cov)
                jobs.append({
                    "config": cfg_dir.name,
                    "orbit": orbit_dir.name,
                    "phase_index": int(f.stem.split("__")[0][1:]),
                    "t_end_days": res.get("t_end_days"),
                    "axis_km": metrics.longest_axis_km,
                    "axes_km": metrics.axes_km,
                    "total_rms_km": metrics.total_rms_km,
                    "longest_direction": metrics.longest_direction,
                })
    return jobs


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------
def write_cases_csv(jobs: list, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    header = (
        "config", "orbit", "amplitude_km", "phase_index", "t_end_days",
        "axis_km", "axis1_km", "axis2_km", "axis3_km", "total_rms_km",
        "axis1_dir_x", "axis1_dir_y", "axis1_dir_z",
    )
    rows = []
    for j in sorted(jobs, key=lambda j: (j["config"], parse_amplitude(j["orbit"]),
                                         j["phase_index"])):
        axes = j["axes_km"]
        direction = j["longest_direction"]
        rows.append({
            "config": j["config"], "orbit": j["orbit"],
            "amplitude_km": parse_amplitude(j["orbit"]),
            "phase_index": j["phase_index"],
            "t_end_days": _fmt(j.get("t_end_days")),
            "axis_km": f"{j['axis_km']:.6g}",
            "axis1_km": f"{axes[0]:.6g}",
            "axis2_km": f"{axes[1]:.6g}",
            "axis3_km": f"{axes[2]:.6g}",
            "total_rms_km": f"{j['total_rms_km']:.6g}",
            "axis1_dir_x": f"{direction[0]:.12g}",
            "axis1_dir_y": f"{direction[1]:.12g}",
            "axis1_dir_z": f"{direction[2]:.12g}",
        })
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=header, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def _fmt(v):
    if v is None or v == "":
        return ""
    try:
        return f"{float(v):.6g}"
    except (TypeError, ValueError):
        return v


def summary_rows(jobs: list) -> list:
    """Per (config, amplitude): mean/min/max axis_km over phases."""
    groups = {}
    for j in jobs:
        key = (j["config"], parse_amplitude(j["orbit"]))
        groups.setdefault(key, []).append(j["axis_km"])
    rows = []
    for (cfg, amp), vals in sorted(groups.items()):
        rows.append({"config": cfg, "amplitude_km": amp, "n_phases": len(vals),
                     "axis_km_mean": f"{sum(vals)/len(vals):.6g}",
                     "axis_km_min": f"{min(vals):.6g}",
                     "axis_km_max": f"{max(vals):.6g}"})
    return rows


METRIC_INFO = {
    "longest": ("longest RMS axis (km)", "axis"),
    "second": ("second RMS axis (km)", "axis2"),
    "third": ("third RMS axis (km)", "axis3"),
    "total-rms": ("total RMS spread (km)", "total_rms"),
}


def metric_value(job, metric):
    if metric == "total-rms":
        return job["total_rms_km"]
    return job["axes_km"][{"longest": 0, "second": 1, "third": 2}[metric]]


def png_heatmap(amps, phases, matrix, title, path, colorbar_label):
    fig, ax = plt.subplots(figsize=(max(10, 0.45*len(phases)), max(6, 0.35*len(amps))))
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
            value = arr[i][j]
            if np.isfinite(value):
                ax.text(j, i, f"{value:.0f}", ha="center", va="center", fontsize=5.5)
    fig.colorbar(im, ax=ax, label=colorbar_label)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def ascii_heatmap(amps, phases, matrix, title, stream, metric_label):
    width = 11
    miss = f"{'.':^{width}}"
    print(f"\n=== {title} ===", file=stream)
    header = "amp\\ph " + "".join(f"{p:>{width}}" for p in phases)
    print(header, file=stream)
    print("-" * len(header), file=stream)
    for i, a in enumerate(amps):
        row = f"{a/1000:>5.1f}k "
        for j in range(len(phases)):
            value = matrix[i][j]
            row += miss if (value is None or
                             (isinstance(value, float) and math.isnan(value))) \
                else f"{value:{width}.1f}"
        print(row, file=stream)
    print(f"(. = not completed; else {metric_label})", file=stream)


def render_config(jobs, outdir: Path, png: bool, metric: str) -> None:
    cfg = jobs[0]["config"]
    amps = sorted({parse_amplitude(j["orbit"]) for j in jobs})
    phases = sorted({j["phase_index"] for j in jobs})
    n_total = len(amps) * len(phases)
    metric_label, filename_tag = METRIC_INFO[metric]

    grid = {(a, p): None for a in amps for p in phases}
    for job in jobs:
        key = (parse_amplitude(job["orbit"]), job["phase_index"])
        grid[key] = metric_value(job, metric)

    matrix = [[grid[(a, p)] for p in phases] for a in amps]
    title = f"{cfg}: position-covariance {metric_label} [{len(jobs)}/{n_total} completed]"
    print(f"\nconfig {cfg}: {len(jobs)}/{n_total} cases completed "
          f"(amplitudes {amps[0]}..{amps[-1]} km)")

    if png and HAS_MPL:
        path = outdir / f"heatmap_cov_{filename_tag}_{cfg}.png"
        png_heatmap(amps, phases, matrix, title, path, metric_label)
        print(f"  heatmap -> {path}")
    else:
        ascii_heatmap(amps, phases, matrix, title, sys.stdout, metric_label)
        if png and not HAS_MPL:
            print("  (matplotlib not installed; install it for PNG heatmaps)",
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
                   help="orbit_id prefix to analyze (default: H1N = Halo North L1)")
    p.add_argument("--configs", nargs="+", default=[],
                   help="restrict to these config ids (default: all present)")
    p.add_argument("--metric", choices=tuple(METRIC_INFO), default="longest",
                   help="heatmap metric (default: longest)")
    p.add_argument("--outdir", type=Path, default=ROOT / "results" / "srp_sensitivity" / "analysis",
                   help="where to write CSVs and PNGs")
    p.add_argument("--no-png", action="store_true", help="force ASCII even if matplotlib is available")
    return p


def main(argv=None) -> int:
    args = parser().parse_args(argv)
    jobs = scan_results(args.results_root)
    if not jobs:
        print(f"no result JSONs found under {args.results_root}", file=sys.stderr)
        return 1

    sel = [j for j in jobs
           if str(j["orbit"]).startswith(args.orbit_prefix)
           and (not args.configs or j["config"] in args.configs)]
    if not sel:
        print(f"no cases for orbit prefix '{args.orbit_prefix}'", file=sys.stderr)
        return 1

    configs = sorted({j["config"] for j in sel})
    print(f"orbit prefix: {args.orbit_prefix}  configs: {configs}  cases: {len(sel)}")

    args.outdir.mkdir(parents=True, exist_ok=True)

    cases_path = args.outdir / f"cov_axis_{args.orbit_prefix}.csv"
    write_cases_csv(sel, cases_path)
    print(f"per-case axis table -> {cases_path}")

    rows = summary_rows(sel)
    summary_path = args.outdir / f"cov_axis_summary_{args.orbit_prefix}.csv"
    with summary_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["config", "amplitude_km", "n_phases",
                                          "axis_km_mean", "axis_km_min", "axis_km_max"])
        w.writeheader()
        w.writerows(rows)
    print(f"per-amplitude summary -> {summary_path}")

    for cfg in configs:
        render_config([j for j in sel if j["config"] == cfg],
                      args.outdir, png=not args.no_png, metric=args.metric)

    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(2)
