#!/usr/bin/env python3
"""Endpoint position-covariance longest-axis analysis for the SRP study.

The per-case result JSON saves only the ENDPOINT mean/covariance (``cov_state``,
a 6x6 row-major matrix), not intermediate-time covariance. This script reads the
``cov_state`` of each completed case, takes the upper-left 3x3 position block
(x,y,z), and computes its longest principal axis = sqrt(largest eigenvalue),
converted to km (x LU_km = 384400). This is one value per case at t_end, not a
time series.

Outputs, for a chosen orbit family (default H1N):
  * cov_axis_<prefix>.csv          -- per case: config/amplitude/phase/t_end/axis_km
  * cov_axis_summary_<prefix>.csv  -- per config x amplitude: mean/min/max over phases
  * heatmap of axis_km over amplitude x phase per config (PNG or ASCII)

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

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RESULTS = ROOT / "results" / "srp_sensitivity" / "results"
LU_KM = 384400.0

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


def position_cov(cov_state: list) -> list:
    """Upper-left 3x3 (x,y,z) block of a row-major 6x6 covariance."""
    c = cov_state
    return [[c[0], c[1], c[2]],
            [c[6], c[7], c[8]],
            [c[12], c[13], c[14]]]


def largest_eigval(A: list) -> float:
    """Largest eigenvalue of a 3x3 symmetric PSD matrix (power iteration)."""
    v = [1.0, 0.0, 0.0]
    lam = 0.0
    for _ in range(200):
        w = [A[0][0]*v[0] + A[0][1]*v[1] + A[0][2]*v[2],
             A[1][0]*v[0] + A[1][1]*v[1] + A[1][2]*v[2],
             A[2][0]*v[0] + A[2][1]*v[1] + A[2][2]*v[2]]
        lam = v[0]*w[0] + v[1]*w[1] + v[2]*w[2]
        n = math.sqrt(w[0]*w[0] + w[1]*w[1] + w[2]*w[2])
        if n < 1e-30:
            return 0.0
        v = [w[0]/n, w[1]/n, w[2]/n]
    return lam


def axis_km(cov_state: list) -> float:
    lam = largest_eigval(position_cov(cov_state))
    if lam < 0.0:
        lam = 0.0
    return math.sqrt(lam) * LU_KM


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
                jobs.append({
                    "config": cfg_dir.name,
                    "orbit": orbit_dir.name,
                    "phase_index": int(f.stem.split("__")[0][1:]),
                    "t_end_days": res.get("t_end_days"),
                    "axis_km": axis_km(cov),
                })
    return jobs


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------
def write_cases_csv(jobs: list, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    header = ("config", "orbit", "amplitude_km", "phase_index", "t_end_days", "axis_km")
    rows = []
    for j in sorted(jobs, key=lambda j: (j["config"], parse_amplitude(j["orbit"]),
                                         j["phase_index"])):
        rows.append({"config": j["config"], "orbit": j["orbit"],
                     "amplitude_km": parse_amplitude(j["orbit"]),
                     "phase_index": j["phase_index"],
                     "t_end_days": _fmt(j.get("t_end_days")),
                     "axis_km": f"{j['axis_km']:.6g}"})
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


def png_heatmap(amps, phases, matrix, title, path):
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
            v = arr[i][j]
            if np.isfinite(v):
                ax.text(j, i, f"{v:.0f}", ha="center", va="center", fontsize=5.5)
    fig.colorbar(im, ax=ax, label="longest axis (km)")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def ascii_heatmap(amps, phases, matrix, title, stream):
    width = 11
    miss = f"{'.':^{width}}"
    print(f"\n=== {title} ===", file=stream)
    header = "amp\\ph " + "".join(f"{p:>{width}}" for p in phases)
    print(header, file=stream)
    print("-" * len(header), file=stream)
    for i, a in enumerate(amps):
        row = f"{a/1000:>5.1f}k "
        for j in range(len(phases)):
            v = matrix[i][j]
            row += miss if (v is None or (isinstance(v, float) and math.isnan(v))) \
                else f"{v:{width}.1f}"
        print(row, file=stream)
    print("(. = not completed; else longest axis in km)", file=stream)


def render_config(jobs, outdir: Path, png: bool) -> None:
    cfg = jobs[0]["config"]
    amps = sorted({parse_amplitude(j["orbit"]) for j in jobs})
    phases = sorted({j["phase_index"] for j in jobs})
    n_total = len(amps) * len(phases)

    grid = {(a, p): None for a in amps for p in phases}
    for j in jobs:
        grid[(parse_amplitude(j["orbit"]), j["phase_index"])] = j["axis_km"]

    matrix = [[grid[(a, p)] for p in phases] for a in amps]
    title = f"{cfg}: position-covariance longest axis (km) [{len(jobs)}/{n_total} completed]"
    print(f"\nconfig {cfg}: {len(jobs)}/{n_total} cases completed "
          f"(amplitudes {amps[0]}..{amps[-1]} km)")

    if png and HAS_MPL:
        path = outdir / f"heatmap_cov_axis_{cfg}.png"
        png_heatmap(amps, phases, matrix, title, path)
        print(f"  heatmap -> {path}")
    else:
        ascii_heatmap(amps, phases, matrix, title, sys.stdout)
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
                      args.outdir, png=not args.no_png)

    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(2)
