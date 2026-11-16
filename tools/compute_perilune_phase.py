#!/usr/bin/env python3
"""Compute the perilune (closest-to-Moon) phase for each catalog orbit.

The Moon sits at (1-mu, 0, 0) in the Earth-Moon synodic frame. Each orbit in
``CRTBP_SRP/crtbp_initial_conditions.json`` carries 24 evenly spaced phase
states along one period, so the perilune phase is the phase of minimum distance
to the Moon, refined by a three-point parabolic fit around the minimum.

Output: a console table plus a CSV (``orbit_id, family, libration_point,
amplitude_km, mu, period_nd, perilune_phase_frac, perilune_time_nd,
r_min_nd, r_max_nd``). ``perilune_phase_frac`` is in [0,1) periods; 0.0 means
phase 0 (reference state) is perilune, 0.5 means apolune. Near-circular orbits
(e.g. DRO) have r_min ~ r_max, so the perilune location is ill-conditioned;
those rows are still emitted but the phase is not meaningful.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_JSON = ROOT / "CRTBP_SRP" / "crtbp_initial_conditions.json"
DEFAULT_OUT = ROOT / "CRTBP_SRP" / "perilune_phase.csv"


def distance_to_moon(state_nd, mu):
    """Distance from a synodic-frame state to the Moon center (LU)."""
    moon_x = 1.0 - mu
    dx = state_nd[0] - moon_x
    dy = state_nd[1]
    dz = state_nd[2]
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def refine_minimum(i_min, d):
    """Parabolic-refined offset (in sample units) around the minimum index."""
    n = len(d)
    a = d[(i_min - 1) % n]
    b = d[i_min]
    c = d[(i_min + 1) % n]
    denom = a - 2.0 * b + c
    if denom == 0.0:
        return 0.0
    delta = 0.5 * (a - c) / denom
    return max(-0.5, min(0.5, delta))


def compute_orbit(orbit) -> dict:
    mu = float(orbit["mu"])
    period_nd = float(orbit["period_nd"])
    phases = orbit["phases"]
    d = [distance_to_moon(p["state_nd"], mu) for p in phases]
    i_min = min(range(len(d)), key=lambda i: d[i])
    delta = refine_minimum(i_min, d)
    n = len(d)
    phase_frac = ((i_min + delta) / n) % 1.0
    if phase_frac < 1e-6 or phase_frac > 1.0 - 1e-6:
        phase_frac = 0.0
    oid = orbit["orbit_id"]
    try:
        amplitude = int(oid.split("_A")[-1])
    except (IndexError, ValueError):
        amplitude = int(orbit.get("target_amplitude_km", 0))
    return {
        "orbit_id": oid,
        "family": orbit.get("family", ""),
        "libration_point": orbit.get("libration_point", ""),
        "amplitude_km": amplitude,
        "mu": mu,
        "period_nd": period_nd,
        "perilune_phase_frac": phase_frac,
        "perilune_time_nd": phase_frac * period_nd,
        "r_min_nd": min(d),
        "r_max_nd": max(d),
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--json", type=Path, default=DEFAULT_JSON)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args(argv)

    data = json.loads(args.json.read_text(encoding="utf-8-sig"))
    orbits = data.get("orbits")
    if not isinstance(orbits, list) or not orbits:
        print(f"no orbits in {args.json}", file=sys.stderr)
        return 1

    rows = [compute_orbit(o) for o in orbits]

    header = ("orbit_id", "family", "libration_point", "amplitude_km", "mu",
              "period_nd", "perilune_phase_frac", "perilune_time_nd",
              "r_min_nd", "r_max_nd")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=header)
        w.writeheader()
        w.writerows(rows)
    print(f"perilune phase table -> {args.out}")

    print(f"\n{'orbit':<14} {'fam':<5} {'amp_km':>7} {'phase_frac':>11} "
          f"{'r_min':>8} {'r_max':>8} {'r_min/r_max':>10}")
    for r in sorted(rows, key=lambda r: r["orbit_id"]):
        ratio = r["r_min_nd"] / r["r_max_nd"]
        print(f"{r['orbit_id']:<14} {r['family']:<5} {r['amplitude_km']:>7} "
              f"{r['perilune_phase_frac']:>11.4f} {r['r_min_nd']:>8.4f} "
              f"{r['r_max_nd']:>8.4f} {ratio:>10.5f}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(2)
