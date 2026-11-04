#!/usr/bin/env python3
"""Divergence distribution for the CRTBP + SRP sensitivity study.

``status == 0`` means the whole parameter box survived (valid covariance).
``status == -4`` means some fixed-seed samples landed in a diverged (dropped)
patch region, so endpoint statistics/covariance are unavailable.

This script reads per-case result JSONs and reports, per config and per orbit
family, how many cases diverged and how severely (``n_diverged``, first
divergence time), plus a status heatmap over amplitude x phase.

Usage:
    python3 tools/analyze_srp_divergence.py \
        --results-root results/srp_sensitivity/results \
        --configs Q1_L1_small Q1_L1_main Q1_L1_ext \
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
                    "n_diverged": res.get("n_diverged"),
                    "first_divergence_time": res.get("first_divergence_time"),
                    "n_patches": res.get("n_patches"),
                })
    return jobs


def select(jobs, orbit_prefix, configs):
    out = []
    for j in jobs:
        if orbit_prefix and not j["orbit"].startswith(orbit_prefix):
            continue
        if configs and j["config"] not in configs:
            continue
        out.append(j)
    return out


# ---------------------------------------------------------------------------
# Family summary
# ---------------------------------------------------------------------------
def family_summary(jobs) -> list:
    """Per (config, family): n_cases, n_diverged, rate, mean n_diverged among
    diverged cases, and mean/median first divergence time (days)."""
    groups = {}
    for j in jobs:
        groups.setdefault((j["config"], j["family"]), []).append(j)
    rows = []
    for (cfg, fam), g in sorted(groups.items()):
        n = len(g)
        div = [j for j in g if j["status"] is not None and j["status"] < 0]
        n_div = len(div)
        nd = [int(j["n_diverged"]) for j in div if j["n_diverged"] is not None]
        ft = [_days(j["first_divergence_time"]) for j in div
              if j["first_divergence_time"] is not None
              and float(j["first_divergence_time"]) >= 0.0]
        rows.append({
            "config": cfg, "family": fam, "n_cases": n, "n_diverged": n_div,
            "divergence_rate": n_div / n if n else 0.0,
            "mean_n_diverged": f"{sum(nd) / len(nd):.3g}" if nd else "",
            "mean_first_div_days": f"{sum(ft) / len(ft):.3g}" if ft else "",
            "median_first_div_days": f"{_median(ft):.3g}" if ft else "",
        })
    return rows


def _days(v):
    try:
        v = float(v)
    except (TypeError, ValueError):
        return None
    if v < 0.0:
        return None
    return v * ND_TO_DAYS


def _median(xs):
    xs = sorted(xs)
    n = len(xs)
    if n == 0:
        return ""
    if n % 2:
        return xs[n // 2]
    return 0.5 * (xs[n // 2 - 1] + xs[n // 2])


# ---------------------------------------------------------------------------
# Status heatmap (amplitude x phase): 0=valid, -4=diverged, missing=nan
# ---------------------------------------------------------------------------
def status_grid(jobs):
    amps = sorted({j["amplitude"] for j in jobs})
    phases = sorted({j["phase_index"] for j in jobs})
    grid = {(a, p): None for a in amps for p in phases}
    for j in jobs:
        key = (j["amplitude"], j["phase_index"])
        if j["status"] == 0:
            grid[key] = 0.0
        elif j["status"] is not None and j["status"] < 0:
            grid[key] = 1.0
        else:
            grid[key] = None
    return amps, phases, grid


def render_status_heatmap(jobs, cfg, outdir, png):
    amps, phases, grid = status_grid(jobs)
    matrix = [[grid[(a, p)] for p in phases] for a in amps]
    title = f"{cfg}: status (green=valid, red=diverged) [{len(jobs)} cases]"
    print(f"\nconfig {cfg}: {len(jobs)} cases, "
          f"{sum(1 for j in jobs if j['status'] is not None and j['status'] < 0)} diverged")
    if png and HAS_MPL:
        path = outdir / f"heatmap_divergence_{cfg}.png"
        _png_status(amps, phases, matrix, title, path)
        print(f"  status heatmap -> {path}")
    else:
        _ascii_status(amps, phases, matrix, title, sys.stdout)
        if png and not HAS_MPL:
            print("  (matplotlib not installed; install for PNG heatmaps)",
                  file=sys.stderr)


def _png_status(amps, phases, matrix, title, path):
    arr = np.array(matrix, dtype=float)
    fig, ax = plt.subplots(figsize=(max(10, 0.45 * len(phases)), max(6, 0.35 * len(amps))))
    cmap = plt.get_cmap("RdYlGn").copy()  # red=diverged, green=valid
    cmap.set_bad(color="0.85")
    im = ax.imshow(arr, aspect="auto", origin="lower", cmap=cmap, vmin=0.0, vmax=1.0)
    ax.set_xticks(range(len(phases)))
    ax.set_xticklabels([str(p) for p in phases], fontsize=7)
    ax.set_yticks(range(len(amps)))
    ax.set_yticklabels([f"{a/1000:.0f}k" for a in amps], fontsize=7)
    ax.set_xlabel("phase_index")
    ax.set_ylabel("amplitude (x1000 km)")
    ax.set_title(title)
    cbar = fig.colorbar(im, ax=ax, ticks=[0.0, 1.0])
    cbar.ax.set_yticklabels(["valid", "diverged"])
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def _ascii_status(amps, phases, matrix, title, stream):
    print(f"\n=== {title} ===", file=stream)
    header = "amp\\ph " + "".join(f"{p:>3}" for p in phases)
    print(header, file=stream)
    print("-" * len(header), file=stream)
    for i, a in enumerate(amps):
        row = f"{a/1000:>5.1f}k "
        for j in range(len(phases)):
            v = matrix[i][j]
            if v is None:
                row += "  ."
            elif v == 0.0:
                row += "  o"
            else:
                row += "  X"
        print(row, file=stream)
    print("(o = valid, X = diverged, . = missing)", file=stream)


# ---------------------------------------------------------------------------
# Diverged-case detail table
# ---------------------------------------------------------------------------
def write_diverged_csv(jobs, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    header = ("config", "family", "orbit", "amplitude_km", "phase_index",
              "status", "n_diverged", "first_divergence_time_nd",
              "first_divergence_days", "n_patches")
    rows = []
    for j in sorted(jobs, key=lambda j: (j["config"], j["family"], j["amplitude"],
                                         j["phase_index"])):
        if j["status"] is None or j["status"] >= 0:
            continue
        rows.append({
            "config": j["config"], "family": j["family"], "orbit": j["orbit"],
            "amplitude_km": j["amplitude"], "phase_index": j["phase_index"],
            "status": j["status"], "n_diverged": j["n_diverged"],
            "first_divergence_time_nd": _fmt(j["first_divergence_time"]),
            "first_divergence_days": _fmt(_days(j["first_divergence_time"])),
            "n_patches": j["n_patches"],
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


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def parser():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--results-root", type=Path, default=DEFAULT_RESULTS)
    p.add_argument("--orbit-prefix", type=str, default=None,
                   help="restrict to this orbit_id prefix (default: all families)")
    p.add_argument("--configs", nargs="+", default=[],
                   help="restrict to these config ids (default: all present)")
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

    sel = select(jobs, args.orbit_prefix, args.configs)
    if not sel:
        print("no matching cases", file=sys.stderr)
        return 1

    args.outdir.mkdir(parents=True, exist_ok=True)
    tag = args.orbit_prefix or "all"

    # family summary
    rows = family_summary(sel)
    summary_path = args.outdir / f"divergence_summary_{tag}.csv"
    with summary_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=[
            "config", "family", "n_cases", "n_diverged", "divergence_rate",
            "mean_n_diverged", "mean_first_div_days", "median_first_div_days"])
        w.writeheader()
        w.writerows(rows)
    print(f"family divergence summary -> {summary_path}")

    # diverged-case detail
    detail_path = args.outdir / f"divergence_cases_{tag}.csv"
    write_diverged_csv(sel, detail_path)
    print(f"diverged-case detail -> {detail_path}")

    # per-config status heatmap (amplitude is only unique within one family)
    if args.orbit_prefix:
        configs = sorted({j["config"] for j in sel})
        for cfg in configs:
            render_status_heatmap([j for j in sel if j["config"] == cfg],
                                  cfg, args.outdir, png=not args.no_png)
    else:
        print("\n(heatmap skipped: amplitudes collide across families; pass "
              "--orbit-prefix to render an amplitude x phase status heatmap)")

    # console family summary
    print("\n=== divergence rate by family (per config) ===")
    for r in rows:
        print(f"{r['config']:<16} {r['family']:<6} "
              f"{r['n_diverged']:>4}/{r['n_cases']:<4} "
              f"({100*r['divergence_rate']:.0f}%) "
              f"mean_n_div={r['mean_n_diverged']} "
              f"first_div_d={r['median_first_div_days']}d")

    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(2)
