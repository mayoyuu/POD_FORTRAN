#!/usr/bin/env python3
"""Single-case ADS split-geometry visualization for the CRTBP + SRP study.

Reads ONE per-case result JSON and renders where/when the automatic domain
splitting (ADS) happened during the propagation:

  * timeline   -- each split as a point (time vs normalized xi, colored by the
                  physical axis beta/chi/u/v that was split)
  * density 1D -- per-axis histogram of the split centers in PHYSICAL units,
                  showing where along each parameter the splits cluster
  * density 2D -- per-axis time-vs-parameter heatmap, showing the split-density
                  hotspots in time-parameter space

Each ``split_events`` entry records: time (``t`` nd / ``td`` days), axis ``ax``
(1=beta, 2=chi, 3=u, 4=v), depth ``d``, active patches ``na``, and the patch
center ``c`` / half-width ``w`` in GLOBAL normalized xi ([-1,1]) along the split
axis. Physical value = domain_center[ax] + domain_halfwidth[ax] * xi, where the
domain is rebuilt from the JSON fields beta_lower/beta_upper/tilt_max_deg/
chi_lower/chi_upper.

Note: only the split axis's center is stored, so we cannot reconstruct each
patch's full 4-D box; but for every split we DO know when it happened and where
along the split parameter it landed.

Usage (single result file):
    python3 tools/analyze_srp_split_geometry.py \
        --json results/srp_sensitivity/results/Q2_L2_t10/H1N_A05000/p03__a0.5.json

Usage (locate a case by config/orbit/phase):
    python3 tools/analyze_srp_split_geometry.py \
        --results-root results/srp_sensitivity/results \
        --config Q2_L2_t10 --orbit H1N_A05000 --phase 3
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RESULTS = ROOT / "results" / "srp_sensitivity" / "results"
LU_KM = 384400.0
TU_S = 375237.080249162

AXIS_NAMES = {1: "beta", 2: "chi", 3: "u", 4: "v"}
AXIS_UNITS = {1: "m^2/kg", 2: "", 3: "", 4: ""}
AXIS_COLORS = {1: "#1f77b4", 2: "#d62728", 3: "#2ca02c", 4: "#9467bd"}

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
def rebuild_domain(res: dict):
    """Return (center[4], halfwidth[4], active[4]) from the JSON physical fields.

    Axes are 0-indexed here: 0=beta, 1=chi, 2=u, 3=v.
    """
    bl = float(res["beta_lower"]); bu = float(res["beta_upper"])
    cl = float(res.get("chi_lower", 0.6)); cu = float(res.get("chi_upper", 0.6))
    tilt = float(res.get("tilt_max_deg", 0.0))

    center = [0.0] * 4
    hw = [0.0] * 4
    center[0] = 0.5 * (bl + bu); hw[0] = 0.5 * (bu - bl)
    center[1] = 0.5 * (cl + cu); hw[1] = 0.5 * (cu - cl)
    if tilt > 0.0:
        q = math.tan(math.radians(tilt)) / math.sqrt(2.0)
        center[2] = center[3] = 0.0
        hw[2] = hw[3] = q
    else:
        center[2] = center[3] = 0.0
        hw[2] = hw[3] = 0.0
    active = [hw[k] > 0.0 for k in range(4)]
    return center, hw, active


def physical_value(center, hw, ax, c):
    k = ax - 1
    return center[k] + hw[k] * c


def load_events(res: dict, center, hw):
    evs = res.get("split_events") or []
    out = []
    for e in evs:
        ax = int(e.get("ax", 0))
        c = float(e.get("c", 0.0))
        out.append({
            "i": int(e.get("i", 0)),
            "t": float(e.get("t", 0.0)),
            "td": float(e.get("td", 0.0)),
            "ax": ax,
            "d": int(e.get("d", 0)),
            "na": int(e.get("na", 0)),
            "c": c,                              # global xi along split axis
            "w": float(e.get("w", 0.0)),
            "crit": float(e.get("crit", 0.0)),
            "x_phys": physical_value(center, hw, ax, c),
        })
    return out


def locate_case(root: Path, config: str, orbit: str, phase: int) -> Path | None:
    d = root / config / orbit
    if not d.is_dir():
        return None
    for f in sorted(d.glob("p*__*.json")):
        stem = f.stem.split("__")[0]
        try:
            idx = int(stem[1:])
        except ValueError:
            continue
        if idx == phase:
            return f
    return None


def describe_case(res: dict) -> str:
    model = res.get("model", "?")
    beta = f"[{res.get('beta_lower')}, {res.get('beta_upper')}]"
    tilt = res.get("tilt_max_deg", 0.0)
    chi = f"[{res.get('chi_lower')}, {res.get('chi_upper')}]"
    return f"model={model} beta={beta} chi={chi} tilt={tilt} deg"


# ---------------------------------------------------------------------------
# Text summary (always printed)
# ---------------------------------------------------------------------------
def print_summary(res, center, hw, active, events):
    ax_names = [k + 1 for k in range(4) if active[k]]
    print(f"case: {describe_case(res)}")
    print(f"active axes: {[AXIS_NAMES[a] for a in ax_names]}")
    print(f"t_end = {res.get('t_end_days', float('nan'))} days, "
          f"n_splits_total = {res.get('n_splits_total')}, "
          f"n_patches = {res.get('n_patches')}, "
          f"n_max_active = {res.get('n_max_active')}, "
          f"n_diverged = {res.get('n_diverged')}")
    print(f"split_count_dir (beta/chi/u/v) = {res.get('split_count_dir')}")
    if not events:
        print("no split events recorded")
        return

    print("\nper-axis split centers (physical):")
    for a in ax_names:
        pts = [e["x_phys"] for e in events if e["ax"] == a]
        ts = [e["td"] for e in events if e["ax"] == a]
        n = len(pts)
        if n == 0:
            print(f"  {AXIS_NAMES[a]:>4}: 0 splits")
            continue
        lo = center[a - 1] - hw[a - 1]
        hi = center[a - 1] + hw[a - 1]
        print(f"  {AXIS_NAMES[a]:>4}: {n} splits  "
              f"center in [{min(pts):.4g}, {max(pts):.4g}] "
              f"(domain [{lo:.4g}, {hi:.4g}] {AXIS_UNITS[a]})  "
              f"t in [{min(ts):.4g}, {max(ts):.4g}] days")


# ---------------------------------------------------------------------------
# ASCII fallback
# ---------------------------------------------------------------------------
def ascii_render(res, center, hw, active, events, stream=sys.stdout):
    print_summary(res, center, hw, active, events)
    if not events:
        return

    ax_names = [k + 1 for k in range(4) if active[k]]
    NB = 12

    print("\nper-axis center histogram (physical; count per bin):")
    for a in ax_names:
        pts = [e["x_phys"] for e in events if e["ax"] == a]
        if not pts:
            continue
        lo = center[a - 1] - hw[a - 1]
        hi = center[a - 1] + hw[a - 1]
        bins = [0] * NB
        for x in pts:
            k = int((x - lo) / (hi - lo) * NB)
            k = min(max(k, 0), NB - 1)
            bins[k] += 1
        mx = max(bins) or 1
        print(f"  {AXIS_NAMES[a]:>4} [{lo:+.4g} .. {hi:+.4g}] {AXIS_UNITS[a]}")
        for k in range(NB):
            bar = "#" * max(1, round(40 * bins[k] / mx)) if bins[k] else "."
            x0 = lo + (hi - lo) * k / NB
            x1 = lo + (hi - lo) * (k + 1) / NB
            print(f"    [{x0:+.4g},{x1:+.4g}) {bins[k]:4d} {bar}")

    print("\ntimeline (sorted by split index i):")
    print(f"  {'i':>4} {'axis':>5} {'t_days':>9} {'xi_center':>10} {'phys':>10} {'depth':>6} {'n_act':>6}")
    for e in events:
        print(f"  {e['i']:>4} {AXIS_NAMES.get(e['ax'], '?'):>5} {e['td']:9.3f} "
              f"{e['c']:10.3f} {e['x_phys']:10.4g} {e['d']:6d} {e['na']:6d}")


# ---------------------------------------------------------------------------
# Matplotlib rendering
# ---------------------------------------------------------------------------
def fig_timeline(res, events, outdir, stem):
    fig, ax = plt.subplots(figsize=(11, 5))
    ax_names = sorted({e["ax"] for e in events})
    for a in ax_names:
        pts = [e for e in events if e["ax"] == a]
        ax.scatter([e["td"] for e in pts], [e["c"] for e in pts],
                   s=28, c=AXIS_COLORS.get(a, "#888"),
                   label=AXIS_NAMES.get(a, f"ax{a}"), alpha=0.85, edgecolors="k",
                   linewidths=0.3)
    ax.axhline(0.0, color="0.7", lw=0.6)
    ax.set_xlabel("time (days)")
    ax.set_ylabel("normalized parameter coordinate (global xi)")
    ax.set_ylim(-1.02, 1.02)
    ax.set_title(f"{stem}: split timeline ({len(events)} splits) — "
                 f"{describe_case(res)}")
    ax.legend(loc="best", ncol=len(ax_names))
    ax.grid(alpha=0.3)
    fig.tight_layout()
    path = outdir / f"{stem}_timeline.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def fig_density_1d(res, center, hw, active, events, outdir, stem):
    ax_names = [k + 1 for k in range(4) if active[k]]
    n = len(ax_names)
    fig, axes = plt.subplots(1, n, figsize=(4.5 * n, 3.6), squeeze=False)
    for j, a in enumerate(ax_names):
        ax = axes[0][j]
        pts = [e["x_phys"] for e in events if e["ax"] == a]
        lo = center[a - 1] - hw[a - 1]
        hi = center[a - 1] + hw[a - 1]
        unit = AXIS_UNITS[a]
        ax.hist(pts, bins=16, range=(lo, hi), color=AXIS_COLORS.get(a), edgecolor="k",
                linewidth=0.4)
        ax.set_xlabel(f"{AXIS_NAMES[a]} ({unit})" if unit else AXIS_NAMES[a])
        ax.set_ylabel("split count")
        ax.set_title(f"axis {AXIS_NAMES[a]} (n={len(pts)})")
        ax.axvline(center[a - 1], color="0.5", lw=0.6)
    fig.suptitle(f"{stem}: split-center density along each parameter", y=1.02)
    fig.tight_layout()
    path = outdir / f"{stem}_density_1d.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


def fig_density_2d(res, center, hw, active, events, outdir, stem):
    ax_names = [k + 1 for k in range(4) if active[k]]
    n = len(ax_names)
    fig, axes = plt.subplots(1, n, figsize=(4.8 * n, 3.8), squeeze=False)
    for j, a in enumerate(ax_names):
        ax = axes[0][j]
        pts = [e for e in events if e["ax"] == a]
        if len(pts) >= 3:
            xs = [e["td"] for e in pts]
            ys = [e["x_phys"] for e in pts]
            lo = center[a - 1] - hw[a - 1]
            hi = center[a - 1] + hw[a - 1]
            h = ax.hist2d(xs, ys, bins=[24, 16],
                          range=[[0.0, float(res.get("t_end_days", 1.0))], [lo, hi]],
                          cmap="viridis")
            fig.colorbar(h[3], ax=ax)
            ax.scatter(xs, ys, s=6, c="w", alpha=0.6, edgecolors="k", linewidths=0.2)
        ax.set_xlabel("time (days)")
        unit = AXIS_UNITS[a]
        ax.set_ylabel(f"{AXIS_NAMES[a]} ({unit})" if unit else AXIS_NAMES[a])
        ax.set_title(f"axis {AXIS_NAMES[a]}")
    fig.suptitle(f"{stem}: split density in time x parameter", y=1.02)
    fig.tight_layout()
    path = outdir / f"{stem}_density_2d.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def parser():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--json", type=Path, default=None,
                   help="path to a single per-case result JSON")
    p.add_argument("--results-root", type=Path, default=DEFAULT_RESULTS,
                   help="root of per-case result JSONs (default: results/srp_sensitivity/results)")
    p.add_argument("--config", type=str, default=None)
    p.add_argument("--orbit", type=str, default=None)
    p.add_argument("--phase", type=int, default=None)
    p.add_argument("--outdir", type=Path, default=ROOT / "results" / "srp_sensitivity" / "analysis",
                   help="where to write PNGs")
    p.add_argument("--no-png", action="store_true",
                   help="force ASCII output even if matplotlib is available")
    return p


def main(argv=None) -> int:
    args = parser().parse_args(argv)

    path = args.json
    if path is None:
        if None in (args.config, args.orbit, args.phase):
            print("need --json, or (--config --orbit --phase) with --results-root",
                  file=sys.stderr)
            return 1
        path = locate_case(args.results_root, args.config, args.orbit, args.phase)
        if path is None:
            print(f"no case found for config={args.config} orbit={args.orbit} "
                  f"phase={args.phase} under {args.results_root}", file=sys.stderr)
            return 1

    if not path.is_file():
        print(f"result JSON not found: {path}", file=sys.stderr)
        return 1

    try:
        res = json.loads(path.read_text(encoding="utf-8-sig"))
    except Exception as exc:
        print(f"failed to read {path}: {exc}", file=sys.stderr)
        return 1

    center, hw, active = rebuild_domain(res)
    events = load_events(res, center, hw)

    stem = path.stem
    if args.orbit and args.config:
        stem = f"{args.config}__{args.orbit}__p{args.phase:02d}"

    if not HAS_MPL or args.no_png:
        ascii_render(res, center, hw, active, events)
        if args.no_png and HAS_MPL:
            print("\n(--no-png: PNG output suppressed)", file=sys.stderr)
        return 0

    print_summary(res, center, hw, active, events)
    if not events:
        print("no split events -> nothing to plot")
        return 0

    args.outdir.mkdir(parents=True, exist_ok=True)
    p1 = fig_timeline(res, events, args.outdir, stem)
    p2 = fig_density_1d(res, center, hw, active, events, args.outdir, stem)
    p3 = fig_density_2d(res, center, hw, active, events, args.outdir, stem)
    print(f"\ntimeline     -> {p1}")
    print(f"density 1D   -> {p2}")
    print(f"density 2D   -> {p3}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(2)
