#!/usr/bin/env python3
"""Diagnose why SRP analysis scripts may not read every result file.

Walks results-root/<config>/<orbit>/p*.json and, per config, counts how many
files exist and how many would be dropped by each filter used in
analyze_srp_cov_axis.scan_results: JSON parse failure, cov_state length != 36,
and status != 0 (broken down by the actual status value).

Usage:
    python3 tools/diagnose_srp_results.py \
        --results-root results/srp_sensitivity/results \
        --orbit-prefix H1N \
        --configs Q1_L1_small Q1_L1_main Q1_L1_ext
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--results-root", type=Path, required=True)
    p.add_argument("--orbit-prefix", type=str, default=None)
    p.add_argument("--configs", nargs="+", default=[])
    args = p.parse_args(argv)

    root = args.results_root
    if not root.is_dir():
        print(f"results-root not found: {root}", file=sys.stderr)
        return 1

    per_cfg = defaultdict(Counter)
    for cfg_dir in sorted(root.iterdir()):
        if not cfg_dir.is_dir():
            continue
        if args.configs and cfg_dir.name not in args.configs:
            continue
        for orbit_dir in sorted(cfg_dir.iterdir()):
            if not orbit_dir.is_dir():
                continue
            if args.orbit_prefix and not orbit_dir.name.startswith(args.orbit_prefix):
                continue
            for f in sorted(orbit_dir.glob("p*.json")):
                c = per_cfg[cfg_dir.name]
                c["files"] += 1
                try:
                    res = json.loads(f.read_text(encoding="utf-8-sig"))
                except Exception:
                    c["json_bad"] += 1
                    continue
                cov = res.get("cov_state")
                if not isinstance(cov, list) or len(cov) != 36:
                    c["cov_bad"] += 1
                st = res.get("status")
                c[f"status_{st}"] += 1

    if not per_cfg:
        print("no matching result files found", file=sys.stderr)
        return 1

    for cfg in sorted(per_cfg):
        c = per_cfg[cfg]
        parts = [f"{cfg}: files={c['files']}",
                 f"json_bad={c['json_bad']}",
                 f"cov_bad={c['cov_bad']}"]
        for k in sorted(c):
            if k.startswith("status_"):
                parts.append(f"{k}={c[k]}")
        print("  ".join(parts))

    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(2)
