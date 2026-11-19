"""Estimate intersection volumes of two endpoint 1-sigma position ellipsoids."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

LU_KM = 384400.0


def endpoint(path):
    record = json.loads(path.read_text(encoding="utf-8-sig"))
    if record.get("status") != 0 or record.get("case_diverged"):
        raise ValueError(f"{path}: endpoint failed or marked diverged")
    return ellipsoid(record["mean_state"], record["cov_state"])


def ellipsoid(mean_state, cov_state):
    mean = np.asarray(mean_state, dtype=float).reshape(6)[:3]
    cov = np.asarray(cov_state, dtype=float).reshape(6, 6)[:3, :3]
    if not (np.isfinite(mean).all() and np.isfinite(cov).all()):
        raise ValueError("nonfinite position mean or covariance")
    try:
        factor = np.linalg.cholesky((cov + cov.T) * 0.5)
    except np.linalg.LinAlgError as exc:
        raise ValueError("3D position covariance is singular or indefinite; "
                         "cannot estimate a full-dimensional volume") from exc
    volume = (4 * math.pi / 3) * np.prod(np.diag(factor) * LU_KM)
    return mean, factor, volume


def ball_points(rng, n):
    unit = rng.normal(size=(n, 3))
    unit *= (np.cbrt(rng.random(n)) / np.linalg.norm(unit, axis=1))[:, None]
    return unit


def overlap_hits(inner, outer, unit):
    small, large = sorted((inner, outer), key=lambda item: item[2])
    transform = np.linalg.solve(large[1], small[1])
    offset = np.linalg.solve(large[1], small[0] - large[0])
    mapped = unit @ transform.T
    return {name: int(np.count_nonzero(np.einsum("ij,ij->i", points, points) <= 1))
            for name, points in (("centered", mapped), ("actual_means", mapped + offset))}


def estimate(inner, outer, samples, seed):
    small = min((inner, outer), key=lambda item: item[2])
    rng = np.random.default_rng(seed)
    hits = dict(centered=0, actual_means=0)
    for start in range(0, samples, 65536):
        n = min(samples - start, 65536)
        # Uniform directions and cube-root radii give uniform volume samples.
        for name, count in overlap_hits(inner, outer, ball_points(rng, n)).items():
            hits[name] += count
    for name, count in hits.items():
        p = count / samples
        volume = small[2] * p
        # Wilson interval for the hit fraction, including zero/all hits.
        z = 1.96
        denom = 1 + z*z/samples
        mid = (p + z*z/(2*samples)) / denom
        half = z * math.sqrt((p*(1-p) + z*z/(4*samples))/samples) / denom
        print(f"{name}: overlap={volume:.8g} km^3, "
              f"MC 95% CI=[{max(0, mid-half)*small[2]:.8g}, "
              f"{min(1, mid+half)*small[2]:.8g}] km^3, "
              f"overlap/Q1={volume/inner[2]:.4%}, overlap/Q2={volume/outer[2]:.4%}, "
              f"IoU={volume/(inner[2]+outer[2]-volume):.4%}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("q1", type=Path)
    parser.add_argument("q2", type=Path)
    parser.add_argument("--samples", type=int, default=1000000)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    if args.samples < 1:
        parser.error("--samples must be positive")
    try:
        q1, q2 = endpoint(args.q1), endpoint(args.q2)
        print(f"V_Q1={q1[2]:.8g} km^3, V_Q2={q2[2]:.8g} km^3")
        estimate(q1, q2, args.samples, args.seed)
    except (OSError, ValueError, KeyError) as exc:
        parser.exit(1, f"error: {exc}\n")
