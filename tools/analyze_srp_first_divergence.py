#!/usr/bin/env python3
"""First-divergence-time heatmap for the CRTBP + SRP sensitivity study.

For each completed case, reads ``status`` and ``first_divergence_time``
(nondimensional; -1 = never diverged). A case that fully survived
(``status == 0``) is shown as "never"; a case where some fixed-seed samples
landed in a diverged patch region (``status == -4``) is shown as its first
divergence time in days (``first_divergence_time * TU_s / 86400``). A missing
(not-yet-run) case is left blank.

One heatmap (amplitude x phase) is rendered per config. Amplitudes are only
unique within a single orbit family, so a heatmap needs ``--orbit-prefix``;
without it only the per-case CSV is written.

Usage:
    python3 tools/analyze_srp_first_divergence.py \
        --results-root results/srp_sensitivity/results \
        --orbit-prefix H1N \
        --configs Q1_beta_0.08_0.12 \
        --outdir results/srp_sensitivity/analysis/H1N

    # all families into their own folders:
    for prefix in H1N H1S H2N H2S PL1 PL2 DRO; do
        python3 tools/analyze_srp_first_divergence.py \
            --results-root results/srp_sensitivity/results \
            --configs Q1_beta_0.08_0.12 \
            --orbit-prefix "$prefix" \
            --outdir "results/srp_sensitivity/analysis/$prefix"
    done
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
                    "first_divergence_time": res.get("first_divergence_time"),
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


def classify(job):
    """Return (category, days) for one case.

    category: "diverged" (has a divergence time), "never" (status==0 survived),
    "failed" (status<0 but no usable first_divergence_time), or "missing".
    """
    st = job["status"]
    if st == 0:
        return "never", None
    if st is not None and st < 0:
        v = job["first_divergence_time"]
        try:
            v = float(v)
        except (TypeError, ValueError):
            v = None
        if v is not None and v >= 0.0:
            return "diverged", v * ND_TO_DAYS
        return "failed", None
    return "missing", None


# ---------------------------------------------------------------------------
# Heatmap
# ---------------------------------------------------------------------------
def render_config(jobs, outdir: Path, png: bool) -> None:
    cfg = jobs[0]["config"]
    amps = sorted({j["amplitude"] for j in jobs})
    phases = sorted({j["phase_index"] for j in jobs})
    n_total = len(amps) * len(phases)
    n_done = len(jobs)

    grid = {(j["amplitude"], j["phase_index"]): classify(j) for j in jobs}

    n_div = n_never = n_failed = n_missing = 0
    float_rows, str_rows = [], []
    for a in amps:
        frow, srow = [], []
        for p in phases:
            cat, days = grid.get((a, p), ("missing", None))
            if cat == "diverged":
                frow.append(days); srow.append(f"{days:7.1f}"); n_div += 1
            elif cat == "never":
                frow.append(-1.0); srow.append("  --  "); n_never += 1
            elif cat == "failed":
                frow.append(-1.0); srow.append("  xx  "); n_failed += 1
            else:
                frow.append(float("nan")); srow.append("  .   "); n_missing += 1
        float_rows.append(frow); str_rows.append(srow)

    title = (f"{cfg}: first divergence time (days) "
             f"[{n_done}/{n_total} done, {n_div} diverged, {n_never} survived]")
    print(f"\nconfig {cfg}: {n_done}/{n_total} cases completed "
          f"({n_div} diverged, {n_never} survived, {n_failed} failed, "
          f"{n_missing} missing)")

    if png and HAS_MPL:
        path = outdir / f"heatmap_first_divergence_{cfg}.png"
        _png_heatmap(amps, phases, float_rows, title, path)
        print(f"  first-divergence heatmap -> {path}")
    else:
        _ascii_heatmap(amps, phases, str_rows, title, sys.stdout)
        if png and not HAS_MPL:
            print("  (matplotlib not installed; install it for PNG heatmaps)",
                  file=sys.stderr)


def _png_heatmap(amps, phases, float_rows, title, path):
    arr = np.array(float_rows, dtype=float)
    fig, ax = plt.subplots(figsize=(max(10, 0.45 * len(phases)), max(6, 0.35 * len(amps))))
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad(color="0.45")     # missing -> dark gray
    cmap.set_under(color="0.85")   # never/failed -> light gray
    vmax = np.nanmax(arr)
    if not np.isfinite(vmax) or vmax <= 0.0:
        vmax = 1.0
    im = ax.imshow(arr, aspect="auto", origin="lower", cmap=cmap,
                   vmin=0.0, vmax=vmax)
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
            if np.isfinite(v) and v >= 0.0:
                ax.text(j, i, f"{v:.1f}", ha="center", va="center", fontsize=5.5)
    cbar = fig.colorbar(im, ax=ax, extend="min")
    cbar.set_label("first divergence time (days)")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def _ascii_heatmap(amps, phases, str_rows, title, stream):
    width = 7
    print(f"\n=== {title} ===", file=stream)
    header = "amp\\ph " + "".join(f"{p:>{width}}" for p in phases)
    print(header, file=stream)
    print("-" * len(header), file=stream)
    for i, a in enumerate(amps):
        row = f"{a/1000:>5.1f}k "
        row += "".join(str_rows[i])
        print(row, file=stream)
    print("(. = missing; -- = survived, never diverged; "
          "xx = failed, no divergence time; else days)", file=stream)


# ---------------------------------------------------------------------------
# Per-case CSV
# ---------------------------------------------------------------------------
def write_cases_csv(jobs, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    header = ("config", "family", "orbit", "amplitude_km", "phase_index",
              "status", "first_divergence_time_nd", "first_divergence_days")
    rows = []
    for j in sorted(jobs, key=lambda j: (j["config"], j["family"], j["amplitude"],
                                         j["phase_index"])):
        cat, days = classify(j)
        rows.append({
            "config": j["config"], "family": j["family"], "orbit": j["orbit"],
            "amplitude_km": j["amplitude"], "phase_index": j["phase_index"],
            "status": j["status"],
            "first_divergence_time_nd": _fmt(j["first_divergence_time"]),
            "first_divergence_days": _fmt(days) if cat == "diverged" else "",
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
                   help="restrict to this orbit_id prefix (required for heatmaps)")
    p.add_argument("--configs", nargs="+", default=[])
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

    csv_path = args.outdir / f"first_divergence_{tag}.csv"
    write_cases_csv(sel, csv_path)
    print(f"per-case table -> {csv_path}")

    if args.orbit_prefix:
        configs = sorted({j["config"] for j in sel})
        for cfg in configs:
            render_config([j for j in sel if j["config"] == cfg],
                          args.outdir, png=not args.no_png)
    else:
        print("\n(heatmap skipped: amplitudes collide across families; pass "
              "--orbit-prefix to render an amplitude x phase heatmap)")

    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(2)
