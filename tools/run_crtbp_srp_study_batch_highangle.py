#!/usr/bin/env python3
"""High-tilt Q2 extension of the SRP sensitivity study.

Reuses the machinery of run_crtbp_srp_study_batch.py but restricts to the
large-tilt (30/45/60/85 deg) extended-beta ([0, 0.2]) flat-plate configs and
raises the default ADS depth to n_split_max=16.

Why 16: the DA map's mu = 1/sqrt(1+u^2+v^2) has unit radius of convergence in
(u,v), while the 85-deg box has corner radius tan(85 deg) ~ 11.4, so ADS must
split the u,v box far more than the 5/10/20-deg runs. Smoke-test one 85-deg
case first (--limit 1 --configs Q2_L2_t85) and inspect n_patches / wall_time_s
before launching the full orbit x phase matrix.

Usage (same flags as the main batch):
    python3 tools/run_crtbp_srp_study_batch_highangle.py --jobs 16 --limit 1 \
        --configs Q2_L2_t85 --outdir results/hi_smoke
"""
from __future__ import annotations

import sys

import run_crtbp_srp_study_batch as base

HIGH_ANGLE_CONFIGS = {
    "Q2_L2_t30": dict(layer="L2", model="flatplate", beta=(0.00, 0.20),
                      tilt_deg=30.0, chi=(0.6, 0.6)),
    "Q2_L2_t45": dict(layer="L2", model="flatplate", beta=(0.00, 0.20),
                      tilt_deg=45.0, chi=(0.6, 0.6)),
    "Q2_L2_t60": dict(layer="L2", model="flatplate", beta=(0.00, 0.20),
                      tilt_deg=60.0, chi=(0.6, 0.6)),
    "Q2_L2_t85": dict(layer="L2", model="flatplate", beta=(0.00, 0.20),
                      tilt_deg=85.0, chi=(0.6, 0.6)),
}


def _has(argv, flag):
    return any(a == flag or a.startswith(flag + "=") for a in argv)


def main() -> int:
    argv = list(sys.argv[1:])
    if not _has(argv, "--outdir"):
        argv += ["--outdir", str(base.ROOT / "results" / "srp_sensitivity_highangle")]
    if not _has(argv, "--n-split-max"):
        argv += ["--n-split-max", "16"]
    base.CONFIGS = HIGH_ANGLE_CONFIGS
    return base.main(argv)


if __name__ == "__main__":
    sys.exit(main())
