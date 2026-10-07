#!/usr/bin/env python3
"""Generate the fixed paired SRP Monte Carlo design.

The seven independent standard-normal coordinates are truncated componentwise
at +/- sigma_limit. ADS normalized coordinates are u=z/sigma_limit.
Treatment eta is eta_sigma*z7; the control group always uses eta=0.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import shutil
from pathlib import Path
from typing import Any

import numpy as np

PREFIX_COUNTS = [10_000, 20_000, 50_000]
GENERATOR_NAME = "numpy.random.PCG64DXSM"
FORMAT_VERSION = "srp-mc-samples-v1"


def generate_samples(
    count: int, seed: int, sigma_limit: float
) -> tuple[np.ndarray, int]:
    """Return deterministic (count, 7) componentwise-truncated normal samples."""
    if count <= 0:
        raise ValueError("count must be positive")
    if not math.isfinite(sigma_limit) or sigma_limit <= 0.0:
        raise ValueError("sigma_limit must be finite and positive")

    rng = np.random.Generator(np.random.PCG64DXSM(seed))
    accepted: list[np.ndarray] = []
    remaining = count
    rejected = 0
    chunk_size = 8192

    while remaining:
        candidates = rng.standard_normal((chunk_size, 7))
        mask = np.all(np.abs(candidates) <= sigma_limit, axis=1)
        indices = np.flatnonzero(mask)

        if indices.size >= remaining:
            last = int(indices[remaining - 1])
            accepted.append(candidates[indices[:remaining]])
            rejected += (last + 1) - remaining
            remaining = 0
        else:
            accepted.append(candidates[indices])
            remaining -= int(indices.size)
            rejected += chunk_size - int(indices.size)

    return np.vstack(accepted), rejected


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _atomic_replace_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".tmp")
    temp.write_text(text, encoding="utf-8", newline="\n")
    os.replace(temp, path)


def write_sample_bundle(
    csv_path: str | Path,
    manifest_path: str | Path,
    *,
    count: int = 50_000,
    seed: int = 20_261_005,
    sigma_limit: float = 3.0,
    eta_sigma: float = 0.3,
    npz_path: str | Path | None = None,
    gzip_path: str | Path | None = None,
) -> dict[str, Any]:
    """Write the canonical CSV and its integrity manifest atomically."""
    if not math.isfinite(eta_sigma) or eta_sigma < 0.0:
        raise ValueError("eta_sigma must be finite and nonnegative")

    csv_path = Path(csv_path)
    manifest_path = Path(manifest_path)
    samples, rejected = generate_samples(count, seed, sigma_limit)
    normalized = samples / sigma_limit
    eta = eta_sigma * samples[:, 6]

    csv_path.parent.mkdir(parents=True, exist_ok=True)
    temp_csv = csv_path.with_name(csv_path.name + ".tmp")
    header = (
        ["sample_id"]
        + [f"z{i}" for i in range(1, 8)]
        + [f"u{i}" for i in range(1, 8)]
        + ["eta_treatment"]
    )
    with temp_csv.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(header)
        for index in range(count):
            values = (
                [index + 1]
                + [format(value, ".17g") for value in samples[index]]
                + [format(value, ".17g") for value in normalized[index]]
                + [format(eta[index], ".17g")]
            )
            writer.writerow(values)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temp_csv, csv_path)

    if npz_path is not None:
        npz_path = Path(npz_path)
        npz_path.parent.mkdir(parents=True, exist_ok=True)
        temp_npz = npz_path.with_name(npz_path.name + ".tmp")
        with temp_npz.open("wb") as stream:
            np.savez_compressed(
                stream,
                sample_id=np.arange(1, count + 1, dtype=np.int64),
                z=samples,
                u=normalized,
                eta_treatment=eta,
            )
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temp_npz, npz_path)

    if gzip_path is not None:
        gzip_path = Path(gzip_path)
        gzip_path.parent.mkdir(parents=True, exist_ok=True)
        temp_gzip = gzip_path.with_name(gzip_path.name + ".tmp")
        with csv_path.open("rb") as source, temp_gzip.open("wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
                shutil.copyfileobj(source, compressed, length=1024 * 1024)
            raw.flush()
            os.fsync(raw.fileno())
        os.replace(temp_gzip, gzip_path)

    manifest: dict[str, Any] = {
        "format_version": FORMAT_VERSION,
        "generator": GENERATOR_NAME,
        "seed": int(seed),
        "sample_count": int(count),
        "dimension": 7,
        "distribution": "independent standard normal, componentwise truncated",
        "sigma_limit": float(sigma_limit),
        "coordinate_rule": "u=z/3",
        "eta_treatment_rule": "eta=eta_sigma*z7",
        "eta_sigma": float(eta_sigma),
        "control_eta": 0.0,
        "prefix_counts": PREFIX_COUNTS,
        "rejected_candidate_count": int(rejected),
        "csv_file": csv_path.name,
        "csv_sha256": _sha256(csv_path),
    }
    if npz_path is not None:
        manifest.update({"npz_file": npz_path.name, "npz_sha256": _sha256(npz_path)})
    if gzip_path is not None:
        manifest.update({
            "gzip_csv_file": gzip_path.name,
            "gzip_csv_sha256": _sha256(gzip_path),
        })
    _atomic_replace_text(
        manifest_path,
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
    )
    return manifest


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--count", type=int, default=50_000)
    parser.add_argument("--seed", type=int, default=20_261_005)
    parser.add_argument("--sigma-limit", type=float, default=3.0)
    parser.add_argument("--eta-sigma", type=float, default=0.3)
    parser.add_argument("--npz", type=Path)
    parser.add_argument("--csv-gz", type=Path)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    manifest = write_sample_bundle(
        args.csv,
        args.manifest,
        count=args.count,
        seed=args.seed,
        sigma_limit=args.sigma_limit,
        eta_sigma=args.eta_sigma,
        npz_path=args.npz,
        gzip_path=args.csv_gz,
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
