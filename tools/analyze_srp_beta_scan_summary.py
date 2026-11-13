#!/usr/bin/env python3
"""Cross-config beta-scan summary of the SRP sensitivity study.

For every ``Q1_beta_<lo>_<hi>`` config and every orbit family, reads the
per-case result JSONs and reports the min/max of two heatmap metrics:

  * cov longest axis (km):  sqrt(largest eigenvalue of the 3x3 endpoint
    position covariance) x 384400, over status==0 cases only (diverged cases
    have no valid covariance and are counted separately).
  * first split time (days): event-detected first_split_time_accurate when
    available, otherwise checkpoint first_split_time, multiplied by TU_s/86400;
    "never split" (-1) cases are counted
    separately.

For each metric the exact (amplitude km, phase index) of the argmin/argmax is
also recorded, so you can see where in the amplitude x phase heatmap the
extrema concentrate. The config's beta center = (lo+hi)/2 and range width =
hi-lo are parsed from the directory name, so the summary shows how the extrema
move as the beta center and range change.

Outputs:
  * beta_scan_summary.csv          -- one row per (family, config)
  * trend_cov_axis.png             -- max/min axis vs beta range width, 7 subplots
  * trend_split_time.png           -- max/min split time vs beta range width
  * grid_axis_<family>.png         -- per-family (center x width) location grid
  * grid_split_<family>.png

Usage:
    python3 tools/analyze_srp_beta_scan_summary.py \
        --results-root results/srp_sensitivity/results \
        --outdir results/srp_sensitivity/analysis

    # restrict to an explicit set of configs (else all Q1_beta_* are auto-detected):
    python3 tools/analyze_srp_beta_scan_summary.py \
        --results-root results/srp_sensitivity/results \
        --configs Q1_beta_0_0.1 Q1_beta_0.01_0.09 Q1_beta_0.08_0.12
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

from srp_covariance_metrics import position_covariance_metrics
from srp_first_split_time import select_first_split_time

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RESULTS = ROOT / "results" / "srp_sensitivity" / "results"
TU_S = 375237.080249162
SEC_PER_DAY = 86400.0
ND_TO_DAYS = TU_S / SEC_PER_DAY

FAMILY_ORDER = ["H1N", "H1S", "H2N", "H2S", "PL1", "PL2", "DRO"]

try:
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    HAS_MPL = True
except Exception:  # pragma: no cover - optional dependency
    HAS_MPL = False


def family_of(orbit_id: str) -> str:
    return orbit_id.split("_A")[0]


def parse_amplitude(orbit_id: str) -> int:
    return int(orbit_id.split("_A")[-1])


def family_key(fam: str) -> int:
    return FAMILY_ORDER.index(fam) if fam in FAMILY_ORDER else len(FAMILY_ORDER)


# ---------------------------------------------------------------------------
# cov longest axis (identical to analyze_srp_cov_axis.py)
# ---------------------------------------------------------------------------
def axis_km(cov_state: list) -> float:
    """Compatibility alias for the true longest one-sigma position axis."""
    return position_covariance_metrics(cov_state).longest_axis_km


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
def scan_results(root: Path) -> list:
    """Walk root/<config>/<orbit>/p*.json and return every completed case."""
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


def parse_beta_config(config: str):
    """Return (lo, hi, center, width) for a Q1_beta_<lo>_<hi> name, else None."""
    if not config.startswith("Q1_beta_"):
        return None
    body = config[len("Q1_beta_"):]
    lo_str, _, hi_str = body.rpartition("_")
    if not lo_str or not hi_str:
        return None
    try:
        lo = float(lo_str)
        hi = float(hi_str)
    except ValueError:
        return None
    return lo, hi, 0.5 * (lo + hi), hi - lo


# ---------------------------------------------------------------------------
# Per-(family, config) summary
# ---------------------------------------------------------------------------
def _extremes(vals):
    """vals: list of (value, amp, phase). Return dict with min/max + coords."""
    if not vals:
        return {"min": None, "min_amp": None, "min_phase": None,
                "max": None, "max_amp": None, "max_phase": None}
    mn = min(vals, key=lambda t: t[0])
    mx = max(vals, key=lambda t: t[0])
    return {"min": mn[0], "min_amp": mn[1], "min_phase": mn[2],
            "max": mx[0], "max_amp": mx[1], "max_phase": mx[2]}


def summarize_group(jobs) -> dict:
    axis_vals = []   # (axis_km, amp, phase)
    n_diverged = 0
    for j in jobs:
        if j["status"] is not None and j["status"] < 0:
            n_diverged += 1
            continue
        cov = j["cov_state"]
        if not isinstance(cov, list) or len(cov) != 36:
            continue
        axis_vals.append((axis_km(cov), j["amplitude"], j["phase_index"]))

    split_vals = []  # (days, amp, phase)
    n_never = 0
    for j in jobs:
        fst, source = select_first_split_time(j)
        if fst is None:
            if source == "never":
                n_never += 1
            continue
        split_vals.append((fst * ND_TO_DAYS, j["amplitude"], j["phase_index"]))

    a = _extremes(axis_vals)
    s = _extremes(split_vals)
    return {
        "axis_min": a["min"], "axis_max": a["max"],
        "axis_min_amp": a["min_amp"], "axis_min_phase": a["min_phase"],
        "axis_max_amp": a["max_amp"], "axis_max_phase": a["max_phase"],
        "axis_n_valid": len(axis_vals), "axis_n_diverged": n_diverged,
        "split_min": s["min"], "split_max": s["max"],
        "split_min_amp": s["min_amp"], "split_min_phase": s["min_phase"],
        "split_max_amp": s["max_amp"], "split_max_phase": s["max_phase"],
        "split_n_split": len(split_vals), "split_n_never": n_never,
    }


# ---------------------------------------------------------------------------
# CSV
# ---------------------------------------------------------------------------
CSV_HEADER = ("family", "config", "beta_lo", "beta_hi", "beta_center", "beta_width",
              "axis_min_km", "axis_max_km", "axis_min_amp", "axis_min_phase",
              "axis_max_amp", "axis_max_phase", "axis_n_valid", "axis_n_diverged",
              "split_min_days", "split_max_days", "split_min_amp", "split_min_phase",
              "split_max_amp", "split_max_phase", "split_n_split", "split_n_never")


def _fmt(v):
    if v is None or v == "":
        return ""
    try:
        return f"{float(v):.6g}"
    except (TypeError, ValueError):
        return v


def _fmt_int(v):
    return "" if v is None else str(int(v))


def write_csv(rows, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=CSV_HEADER, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({
                "family": r["family"], "config": r["config"],
                "beta_lo": _fmt(r["beta_lo"]), "beta_hi": _fmt(r["beta_hi"]),
                "beta_center": _fmt(r["beta_center"]), "beta_width": _fmt(r["beta_width"]),
                "axis_min_km": _fmt(r["axis_min"]), "axis_max_km": _fmt(r["axis_max"]),
                "axis_min_amp": _fmt_int(r["axis_min_amp"]), "axis_min_phase": _fmt_int(r["axis_min_phase"]),
                "axis_max_amp": _fmt_int(r["axis_max_amp"]), "axis_max_phase": _fmt_int(r["axis_max_phase"]),
                "axis_n_valid": _fmt_int(r["axis_n_valid"]), "axis_n_diverged": _fmt_int(r["axis_n_diverged"]),
                "split_min_days": _fmt(r["split_min"]), "split_max_days": _fmt(r["split_max"]),
                "split_min_amp": _fmt_int(r["split_min_amp"]), "split_min_phase": _fmt_int(r["split_min_phase"]),
                "split_max_amp": _fmt_int(r["split_max_amp"]), "split_max_phase": _fmt_int(r["split_max_phase"]),
                "split_n_split": _fmt_int(r["split_n_split"]), "split_n_never": _fmt_int(r["split_n_never"]),
            })


# ---------------------------------------------------------------------------
# Trend plots
# ---------------------------------------------------------------------------
def plot_trend(rows, metric, outdir):
    if metric == "axis":
        ylabel = "cov longest axis (km)"
        fname = "trend_cov_axis.png"
        max_key, min_key = "axis_max", "axis_min"
    else:
        ylabel = "first split time (days)"
        fname = "trend_split_time.png"
        max_key, min_key = "split_max", "split_min"

    families = [f for f in FAMILY_ORDER if any(r["family"] == f for r in rows)]
    centers = sorted({r["beta_center"] for r in rows})
    cmap = plt.get_cmap("tab10")
    center_color = {c: cmap(i % 10) for i, c in enumerate(centers)}

    fig, axes = plt.subplots(2, 4, figsize=(22, 10))
    axes = axes.ravel()

    for i, fam in enumerate(families):
        ax = axes[i]
        fam_rows = [r for r in rows if r["family"] == fam]
        for c in centers:
            cr = sorted([r for r in fam_rows if abs(r["beta_center"] - c) < 1e-9],
                        key=lambda r: r["beta_width"])
            if not cr:
                continue
            xs = [r["beta_width"] for r in cr]
            pts_max = [(x, r[max_key]) for x, r in zip(xs, cr) if r[max_key] is not None]
            pts_min = [(x, r[min_key]) for x, r in zip(xs, cr) if r[min_key] is not None]
            if pts_max:
                ax.plot([p[0] for p in pts_max], [p[1] for p in pts_max],
                        "o-", color=center_color[c], lw=1.5)
            if pts_min:
                ax.plot([p[0] for p in pts_min], [p[1] for p in pts_min],
                        "s--", color=center_color[c], lw=1.2, alpha=0.7)
        ax.set_title(fam, fontsize=11)
        ax.set_xlabel("beta range width")
        ax.set_ylabel(ylabel, fontsize=8)
        ax.grid(alpha=0.3)
        ax.tick_params(labelsize=8)

    legend_handles = []
    for c in centers:
        col = center_color[c]
        legend_handles.append(Line2D([0], [0], color=col, lw=1.5, marker="o",
                                     label=f"center {c:.3g} max"))
        legend_handles.append(Line2D([0], [0], color=col, lw=1.2, ls="--", marker="s",
                                     label=f"center {c:.3g} min"))
    axes[7].axis("off")
    axes[7].legend(handles=legend_handles, loc="center", fontsize=9, title=ylabel)
    for j in range(len(families), 7):
        axes[j].axis("off")

    fig.suptitle(ylabel, fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(outdir / fname, dpi=150)
    plt.close(fig)
    print(f"  trend plot -> {outdir / fname}")


# ---------------------------------------------------------------------------
# Location grid
# ---------------------------------------------------------------------------
def plot_location_grid(rows, metric, outdir):
    """Per-family grid of where each metric's extrema concentrate.

    x-axis = beta center, y-axis = beta width; each cell shows the
    (amplitude, phase) of the max and min separately ("max A.. phi.." /
    "min A.. phi.."), or is left blank when that (center, width) config is
    absent or the metric has no valid cases. One PNG per family.
    """
    if metric == "axis":
        ylabel = "cov longest axis (km)"
        fname_prefix = "grid_axis"
        max_key, min_key = "axis_max", "axis_min"
        max_amp, min_amp = "axis_max_amp", "axis_min_amp"
        max_phase, min_phase = "axis_max_phase", "axis_min_phase"
    else:
        ylabel = "first split time (days)"
        fname_prefix = "grid_split"
        max_key, min_key = "split_max", "split_min"
        max_amp, min_amp = "split_max_amp", "split_min_amp"
        max_phase, min_phase = "split_max_phase", "split_min_phase"

    families = [f for f in FAMILY_ORDER if any(r["family"] == f for r in rows)]
    centers = sorted({round(r["beta_center"], 6) for r in rows})
    widths = sorted({round(r["beta_width"], 6) for r in rows})

    for fam in families:
        fam_rows = [r for r in rows if r["family"] == fam]
        lookup = {(round(r["beta_center"], 6), round(r["beta_width"], 6)): r
                  for r in fam_rows}

        fig, ax = plt.subplots(figsize=(max(8, 3.2 * len(centers)),
                                        max(5, 2.0 * len(widths))))
        for yi, w in enumerate(widths):
            for xi, c in enumerate(centers):
                r = lookup.get((c, w))
                lines = []
                if r is not None:
                    if r[max_key] is not None:
                        lines.append(f"max A{r[max_amp]} φ{r[max_phase]}")
                    if r[min_key] is not None:
                        lines.append(f"min A{r[min_amp]} φ{r[min_phase]}")
                ax.text(xi, yi, "\n".join(lines), ha="center", va="center",
                        fontsize=10)

        ax.set_xlim(-0.5, len(centers) - 0.5)
        ax.set_ylim(len(widths) - 0.5, -0.5)
        ax.set_xticks(range(len(centers)))
        ax.set_xticklabels([f"{c:.3g}" for c in centers], fontsize=10)
        ax.set_yticks(range(len(widths)))
        ax.set_yticklabels([f"{w:.3g}" for w in widths], fontsize=10)
        ax.set_xlabel("beta center", fontsize=11)
        ax.set_ylabel("beta width", fontsize=11)
        ax.set_title(f"{fam}: {ylabel} extrema location", fontsize=12)
        ax.set_xticks([i - 0.5 for i in range(len(centers) + 1)], minor=True)
        ax.set_yticks([i - 0.5 for i in range(len(widths) + 1)], minor=True)
        ax.grid(which="minor", color="0.8", lw=1)
        ax.tick_params(which="major", length=0)
        ax.tick_params(which="minor", length=0)

        fname = f"{fname_prefix}_{fam}.png"
        fig.tight_layout()
        fig.savefig(outdir / fname, dpi=150)
        plt.close(fig)
        print(f"  location grid -> {outdir / fname}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def parser():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--results-root", type=Path, default=DEFAULT_RESULTS)
    p.add_argument("--configs", nargs="+", default=[],
                   help="explicit config ids (default: auto-detect all Q1_beta_*)")
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

    if args.configs:
        configs = args.configs
    else:
        configs = sorted({j["config"] for j in jobs
                          if j["config"].startswith("Q1_beta_")})

    parsed = {}
    for cfg in configs:
        p = parse_beta_config(cfg)
        if p is None:
            print(f"skip config '{cfg}': not a Q1_beta_<lo>_<hi> name", file=sys.stderr)
            continue
        parsed[cfg] = p

    if not parsed:
        print("no Q1_beta_* configs found (pass --configs to specify them)",
              file=sys.stderr)
        return 1

    rows = []
    for cfg in configs:
        if cfg not in parsed:
            continue
        lo, hi, center, width = parsed[cfg]
        cfg_jobs = [j for j in jobs if j["config"] == cfg]
        fam_groups = {}
        for j in cfg_jobs:
            fam_groups.setdefault(j["family"], []).append(j)
        for fam, g in fam_groups.items():
            s = summarize_group(g)
            rows.append({"family": fam, "config": cfg,
                         "beta_lo": lo, "beta_hi": hi,
                         "beta_center": center, "beta_width": width, **s})

    rows.sort(key=lambda r: (family_key(r["family"]), r["beta_center"], r["beta_width"]))

    args.outdir.mkdir(parents=True, exist_ok=True)

    csv_path = args.outdir / "beta_scan_summary.csv"
    write_csv(rows, csv_path)
    print(f"summary table -> {csv_path}")

    print(f"\n{len(parsed)} configs x up to {len(FAMILY_ORDER)} families "
          f"= {len(rows)} rows")

    if not args.no_png and HAS_MPL:
        plot_trend(rows, "axis", args.outdir)
        plot_trend(rows, "split", args.outdir)
        plot_location_grid(rows, "axis", args.outdir)
        plot_location_grid(rows, "split", args.outdir)
    elif not args.no_png and not HAS_MPL:
        print("(matplotlib not installed; install it for PNG trend plots)",
              file=sys.stderr)

    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(2)
