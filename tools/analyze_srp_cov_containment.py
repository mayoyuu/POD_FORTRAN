#!/usr/bin/env python3
"""Cross-config covariance containment / volume-expansion heatmaps.

Compares the propagation-time covariance (the per-checkpoint ``cov`` written in
each case result JSON) between two SRP study configurations:

  * ``inner`` -- the reference config whose uncertainty should be CONTAINED
    (e.g. Q1 cannonball, beta only);
  * ``outer`` -- the config that should CONTAIN it (e.g. Q2 flatplate with a
    theta/phi attitude box).

For every orbit+phase present in both configs, and for each checkpoint in the
propagation, the position 3x3 covariance blocks are compared:

  * Containment test: the minimum generalized eigenvalue lambda_min of
    (Sigma_outer, Sigma_inner), i.e. the smallest eigenvalue of
    Sigma_inner^{-1} Sigma_outer. lambda_min >= 1  <=>  Sigma_outer >=
    Sigma_inner  <=>  the outer covariance ellipsoid contains the inner one.
  * Volume expansion factor: only over checkpoints where containment holds,
    sqrt(det Sigma_outer / det Sigma_inner) -- how many times the outer
    ellipsoid volume is larger than the inner one.

A case is SKIPPED (left blank) if either side diverged (``case_diverged`` or
``status < 0``). Per case the containment rate (fraction of checkpoints where
the outer contains the inner) and the mean/median/min/max of the volume factor
are aggregated and written to CSV, and rendered as one heatmap each
(amplitude x phase) when ``--orbit-prefix`` selects a single family.

Outputs (prefix = --orbit-prefix):
  * cov_containment_<prefix>.csv
  * heatmap_containment_rate_<prefix>.png
  * heatmap_vol_mean_<prefix>.png / _median_ / _min_ / _max_

Usage:
    python3 tools/analyze_srp_cov_containment.py \
        --inner-results results/srp_sensitivity/results/Q1_beta_0.08_0.12 \
        --outer-results results/srp_sensitivity/results/q2_theta_phi/results/Q2_L2_t10__theta_phi_t0_45_p-45_45 \
        --orbit-prefix H1N \
        --outdir results/srp_sensitivity/analysis/containment
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTDIR = ROOT / "results" / "srp_sensitivity" / "analysis" / "containment"

FAMILY_ORDER = ["H1N", "H1S", "H2N", "H2S", "PL1", "PL2", "DRO"]

try:
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAS_MPL = True
except Exception:  # pragma: no cover - optional dependency
    HAS_MPL = False


def family_of(orbit_id: str) -> str:
    return orbit_id.split("_A")[0]


def parse_amplitude(orbit_id: str) -> int:
    return int(orbit_id.split("_A")[-1])


# ---------------------------------------------------------------------------
# 3x3 linear algebra (pure Python, no numpy dependency for the core metric)
# ---------------------------------------------------------------------------
def pos_cov_3x3(cov36):
    """Extract the symmetrized 3x3 position block from a flattened 6x6 cov."""
    try:
        c = [float(x) for x in cov36]
    except (TypeError, ValueError):
        return None
    if len(c) != 36 or not all(math.isfinite(x) for x in c):
        return None
    a00, a01, a02 = c[0], c[1], c[2]
    a10, a11, a12 = c[6], c[7], c[8]
    a20, a21, a22 = c[12], c[13], c[14]
    return [[a00, 0.5 * (a01 + a10), 0.5 * (a02 + a20)],
            [0.5 * (a10 + a01), a11, 0.5 * (a12 + a21)],
            [0.5 * (a20 + a02), 0.5 * (a21 + a12), a22]]


def _cholesky_inv3(A):
    """Inverse of the lower-triangular Cholesky factor of a 3x3 PSD matrix."""
    a00 = A[0][0]
    if not (a00 > 0.0):
        return None
    l00 = math.sqrt(a00)
    l10 = A[1][0] / l00
    l20 = A[2][0] / l00
    d1 = A[1][1] - l10 * l10
    if not (d1 > 0.0):
        return None
    l11 = math.sqrt(d1)
    l21 = (A[2][1] - l20 * l10) / l11
    d2 = A[2][2] - l20 * l20 - l21 * l21
    if not (d2 > 0.0):
        return None
    l22 = math.sqrt(d2)
    return [[1.0 / l00, 0.0, 0.0],
            [-l10 / (l00 * l11), 1.0 / l11, 0.0],
            [(l10 * l21 - l20 * l11) / (l00 * l11 * l22),
             -l21 / (l11 * l22), 1.0 / l22]]


def _congruence(A2, Linv):
    """M = Linv @ A2 @ Linv^T."""
    t = [[sum(A2[i][k] * Linv[j][k] for k in range(3)) for j in range(3)]
         for i in range(3)]
    return [[sum(Linv[i][k] * t[k][j] for k in range(3)) for j in range(3)]
            for i in range(3)]


def _eigvals_sym3(a):
    """Eigenvalues (descending) of a real symmetric 3x3 matrix (Jacobi)."""
    m = [row[:] for row in a]
    for _ in range(50):
        p, q = max(((0, 1), (0, 2), (1, 2)),
                   key=lambda pr: abs(m[pr[0]][pr[1]]))
        off = m[p][q]
        scale = max(max(abs(x) for row in m for x in row), 1.0e-300)
        if abs(off) <= 1.0e-15 * scale:
            break
        tau = (m[q][q] - m[p][p]) / (2.0 * off)
        t = (math.copysign(1.0 / (abs(tau) + math.sqrt(1.0 + tau * tau)), tau)
             if tau != 0.0 else 1.0)
        c = 1.0 / math.sqrt(1.0 + t * t)
        s = t * c
        app, aqq = m[p][p], m[q][q]
        m[p][p] = c * c * app - 2.0 * s * c * off + s * s * aqq
        m[q][q] = s * s * app + 2.0 * s * c * off + c * c * aqq
        m[p][q] = m[q][p] = 0.0
        for k in range(3):
            if k == p or k == q:
                continue
            akp, akq = m[k][p], m[k][q]
            m[k][p] = m[p][k] = c * akp - s * akq
            m[k][q] = m[q][k] = s * akp + c * akq
    return sorted((m[0][0], m[1][1], m[2][2]), reverse=True)


def generalized_check(inner_cov, outer_cov):
    """Return (lambda_min, volume_factor) for the (inner, outer) pair, or None.

    lambda_min = smallest generalized eigenvalue of (Sigma_outer, Sigma_inner).
    volume_factor = sqrt(det Sigma_outer / det Sigma_inner).
    """
    linv = _cholesky_inv3(inner_cov)
    if linv is None:
        return None
    m = _congruence(outer_cov, linv)
    eig = _eigvals_sym3(m)
    lam_min = eig[2]
    if lam_min < -1.0e-10:
        return None  # outer covariance is not positive semidefinite
    prod = eig[0] * eig[1] * max(eig[2], 0.0)
    if prod < 0.0:
        return None
    return lam_min, math.sqrt(prod)


# ---------------------------------------------------------------------------
# Data loading / matching
# ---------------------------------------------------------------------------
def scan_config(root: Path) -> dict:
    """Return {(orbit_id, phase_index): result_json} over <orbit>/p*.json."""
    cases = {}
    if not root.is_dir():
        return cases
    for orbit_dir in sorted(root.iterdir()):
        if not orbit_dir.is_dir():
            continue
        oid = orbit_dir.name
        for f in sorted(orbit_dir.glob("p*.json")):
            try:
                res = json.loads(f.read_text(encoding="utf-8-sig"))
            except Exception:
                continue
            try:
                pi = int(f.stem.split("__")[0][1:])
            except ValueError:
                continue
            cases[(oid, pi)] = res
    return cases


def is_diverged(res: dict) -> bool:
    if res.get("case_diverged") is True:
        return True
    st = res.get("status")
    return st is not None and st < 0


def compare_case(inner_res, outer_res):
    """Aggregate containment over checkpoints. Returns None if either diverged."""
    if is_diverged(inner_res) or is_diverged(outer_res):
        return None
    inner_cp = inner_res.get("checkpoints") or []
    outer_cp = outer_res.get("checkpoints") or []
    n = min(len(inner_cp), len(outer_cp))
    n_valid = 0
    n_contained = 0
    vols = []
    for i in range(n):
        icov = inner_cp[i].get("cov")
        ocov = outer_cp[i].get("cov")
        if not (isinstance(icov, list) and len(icov) == 36
                and isinstance(ocov, list) and len(ocov) == 36):
            continue
        inner3 = pos_cov_3x3(icov)
        outer3 = pos_cov_3x3(ocov)
        if inner3 is None or outer3 is None:
            continue
        r = generalized_check(inner3, outer3)
        if r is None:
            continue
        lam_min, vol = r
        n_valid += 1
        if lam_min >= 1.0 - 1.0e-6:
            n_contained += 1
            vols.append(vol)
    if n_valid == 0:
        return {"n_valid": 0, "n_contained": 0, "rate": None, "vols": []}
    return {"n_valid": n_valid, "n_contained": n_contained,
            "rate": n_contained / n_valid, "vols": vols}


def _agg(vals):
    if not vals:
        return None, None, None, None
    s = sorted(vals)
    n = len(s)
    mean = sum(s) / n
    median = s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])
    return mean, median, s[0], s[-1]


# ---------------------------------------------------------------------------
# CSV
# ---------------------------------------------------------------------------
CSV_HEADER = ("family", "orbit", "amplitude_km", "phase_index", "status",
              "n_checkpoints", "n_contained", "containment_rate",
              "vol_mean", "vol_median", "vol_min", "vol_max")


def _fmt(v):
    if v is None or v == "":
        return ""
    try:
        return f"{float(v):.6g}"
    except (TypeError, ValueError):
        return v


def write_csv(rows, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=CSV_HEADER, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({
                "family": r["family"], "orbit": r["orbit"],
                "amplitude_km": r["amplitude_km"], "phase_index": r["phase_index"],
                "status": r["status"],
                "n_checkpoints": _fmt(r["n_checkpoints"]),
                "n_contained": _fmt(r["n_contained"]),
                "containment_rate": _fmt(r["containment_rate"]),
                "vol_mean": _fmt(r["vol_mean"]), "vol_median": _fmt(r["vol_median"]),
                "vol_min": _fmt(r["vol_min"]), "vol_max": _fmt(r["vol_max"]),
            })


# ---------------------------------------------------------------------------
# Heatmaps
# ---------------------------------------------------------------------------
def _render_heatmap(amps, phases, matrix, title, cmap_name, fmt_spec, path):
    arr = np.array(matrix, dtype=float)
    fig, ax = plt.subplots(figsize=(max(10, 0.45 * len(phases)),
                                    max(6, 0.35 * len(amps))))
    cmap = plt.get_cmap(cmap_name).copy()
    cmap.set_bad(color="0.45")   # missing -> dark gray
    vmax = np.nanmax(arr)
    if not np.isfinite(vmax) or vmax <= 0.0:
        vmax = 1.0
    im = ax.imshow(arr, aspect="auto", origin="lower", cmap=cmap, vmin=0.0, vmax=vmax)
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
                ax.text(j, i, fmt_spec.format(v), ha="center", va="center", fontsize=5.5)
    cbar = fig.colorbar(im, ax=ax, extend="max")
    cbar.set_label(title.split(":")[0])
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def render_family(rows, prefix, outdir, png):
    amps = sorted({r["amplitude_km"] for r in rows})
    phases = sorted({r["phase_index"] for r in rows})
    grid = {(r["amplitude_km"], r["phase_index"]): r for r in rows}

    metrics = [
        ("containment_rate", "containment rate", "YlOrRd", "{:.2f}",
         "heatmap_containment_rate"),
        ("vol_mean", "volume expansion mean", "cividis", "{:.2f}", "heatmap_vol_mean"),
        ("vol_median", "volume expansion median", "cividis", "{:.2f}", "heatmap_vol_median"),
        ("vol_min", "volume expansion min", "cividis", "{:.2f}", "heatmap_vol_min"),
        ("vol_max", "volume expansion max", "cividis", "{:.2f}", "heatmap_vol_max"),
    ]
    for key, label, cmap, spec, fname in metrics:
        matrix = []
        for a in amps:
            row = []
            for p in phases:
                r = grid.get((a, p))
                v = r.get(key) if r else None
                row.append(v if v is not None else float("nan"))
            matrix.append(row)
        title = f"{label}: {prefix}"
        path = outdir / f"{fname}_{prefix}.png"
        if png and HAS_MPL:
            _render_heatmap(amps, phases, matrix, title, cmap, spec, path)
            print(f"  heatmap -> {path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def parser():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--inner-results", type=Path, required=True,
                   help="config dir (Q1) whose covariance should be contained")
    p.add_argument("--outer-results", type=Path, required=True,
                   help="config dir (Q2) that should contain the inner covariance")
    p.add_argument("--orbit-prefix", type=str, default=None,
                   help="restrict to this orbit_id prefix (required for heatmaps)")
    p.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    p.add_argument("--no-png", action="store_true")
    return p


def main(argv=None) -> int:
    args = parser().parse_args(argv)
    inner = scan_config(args.inner_results)
    outer = scan_config(args.outer_results)
    if not inner or not outer:
        print("one or both results dirs are empty / missing", file=sys.stderr)
        return 1

    common = sorted(set(inner) & set(outer))
    if not common:
        print("no orbit+phase cases shared by the two configs", file=sys.stderr)
        return 1

    rows = []
    for oid, pi in common:
        if args.orbit_prefix and not oid.startswith(args.orbit_prefix):
            continue
        r = compare_case(inner[(oid, pi)], outer[(oid, pi)])
        if r is None:
            rows.append({"family": family_of(oid), "orbit": oid,
                         "amplitude_km": parse_amplitude(oid), "phase_index": pi,
                         "status": "diverged", "n_checkpoints": None,
                         "n_contained": None, "containment_rate": None,
                         "vol_mean": None, "vol_median": None,
                         "vol_min": None, "vol_max": None})
            continue
        mean, median, vmin, vmax = _agg(r["vols"])
        rows.append({"family": family_of(oid), "orbit": oid,
                     "amplitude_km": parse_amplitude(oid), "phase_index": pi,
                     "status": "ok", "n_checkpoints": r["n_valid"],
                     "n_contained": r["n_contained"],
                     "containment_rate": r["rate"],
                     "vol_mean": mean, "vol_median": median,
                     "vol_min": vmin, "vol_max": vmax})

    rows.sort(key=lambda r: (FAMILY_ORDER.index(r["family"])
                             if r["family"] in FAMILY_ORDER else 99,
                             r["amplitude_km"], r["phase_index"]))

    args.outdir.mkdir(parents=True, exist_ok=True)
    tag = args.orbit_prefix or "all"
    csv_path = args.outdir / f"cov_containment_{tag}.csv"
    write_csv(rows, csv_path)
    print(f"per-case table -> {csv_path}")

    if args.orbit_prefix:
        render_family(rows, args.orbit_prefix, args.outdir, png=not args.no_png)
    else:
        print("\n(heatmap skipped: pass --orbit-prefix to render one family's "
              "amplitude x phase heatmaps)")

    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(2)
