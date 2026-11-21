#!/usr/bin/env python3
"""Compare Q1/Q2 synodic position boxes at each saved checkpoint.

Assume a locally untwisted (affine over the full parameter domain) response
to one uniform parameter in Q1 and three independent uniform parameters in
Q2. Along a single coordinate the true affine half-width is the sum of the
absolute scaled response coefficients, whereas sqrt(3)*sigma is their
Euclidean norm. Hence the Q1 half-width is sqrt(3)*sigma, and the Q2
half-width lies between sqrt(3)*sigma and 3*sigma. Sample moments make
these estimates, not rigorous bounds on the propagated trajectories.
"""
from __future__ import annotations

import argparse
from array import array
from collections import Counter, defaultdict
import csv
import json
import math
from pathlib import Path
import sys
from statistics import median
from tempfile import TemporaryDirectory

from analyze_srp_cov_containment import is_diverged, matched_times, scan_config
from srp_covariance_metrics import LU_KM

ND_TO_DAYS = 375237.080249162 / 86400.0
AXES = ("x", "y", "z")
BOXES = (("small", math.sqrt(3)), ("large", 3.0))
FIELDS = ["orbit", "phase_index", "alpha0_rad", "t_nd", "t_days", "status", "judgment"]
for box, _ in BOXES:
    FIELDS.extend([f"{box}_min_margin_km", f"{box}_limiting_face"])
    for axis in AXES:
        FIELDS.extend([f"{box}_{axis}_minus_margin_km", f"{box}_{axis}_plus_margin_km",
                       f"{box}_{axis}_width_gain_km"])


def position_stats(cp):
    mean, cov = cp.get("mean"), cp.get("cov")
    if not isinstance(mean, list) or len(mean) != 6 or not isinstance(cov, list) or len(cov) != 36:
        raise ValueError("expected mean[6] and cov[36]")
    center = [float(x) for x in mean[:3]]
    # The study serializes the 6x6 covariance in row order.
    variance = [float(cov[i]) for i in (0, 7, 14)]
    if not all(map(math.isfinite, center + variance)) or min(variance) < 0:
        raise ValueError("invalid position mean or variance")
    return center, [math.sqrt(v) * LU_KM for v in variance]


def compare_point(inner, outer, tol_km):
    for label, cp in (("q1", inner), ("q2", outer)):
        if cp is None:
            return {"status": f"missing_{label}_checkpoint"}
        if cp.get("n_diverged", 0) or cp.get("n_patches") == 0:
            return {"status": f"unusable_{label}_patches"}
    c1, s1 = position_stats(inner)
    c2, s2 = position_stats(outer)
    result = {"status": "ok"}
    for box, factor in BOXES:
        margins = {}
        for i, axis in enumerate(AXES):
            # Subtract centers before converting to km. A positive margin
            # is clearance inside Q2 on this face, a negative one is excess.
            shift = (c1[i] - c2[i]) * LU_KM
            gain = factor * s2[i] - math.sqrt(3) * s1[i]
            margins[f"{axis}-"] = gain + shift
            margins[f"{axis}+"] = gain - shift
            result[f"{box}_{axis}_minus_margin_km"] = gain + shift
            result[f"{box}_{axis}_plus_margin_km"] = gain - shift
            result[f"{box}_{axis}_width_gain_km"] = 2 * gain
        face = min(margins, key=margins.get)
        result[f"{box}_limiting_face"] = face
        result[f"{box}_min_margin_km"] = margins[face]
    if result["small_min_margin_km"] >= -tol_km:
        result["judgment"] = "contained"
    elif result["large_min_margin_km"] < -tol_km:
        result["judgment"] = "not_contained"
    else:
        result["judgment"] = "indeterminate"
    return result


def case_rows(inner_path, outer_path, tol_km):
    records = []
    for label, path in (("q1", inner_path), ("q2", outer_path)):
        try:
            records.append(json.loads(path.read_text(encoding="utf-8-sig")))
        except (OSError, ValueError) as exc:
            yield {"status": f"invalid_{label}_json: {exc}"}
            return
    inner, outer = records
    for label, record, limit in (("q1", inner, 1), ("q2", outer, 3)):
        indices = record.get("parameter_da_index")
        if indices is not None and sum(i > 0 for i in indices) > limit:
            yield {"status": f"unsupported_{label}_parameter_count"}
            return
    pairs = matched_times(inner, outer)
    if not pairs:
        yield {"status": "no_checkpoints"}
        return
    for t, a, b in pairs:
        row = {"t_nd": t, "t_days": t * ND_TO_DAYS}
        if any(is_diverged(record) and (cp is None or "n_diverged" not in cp)
               for record, cp in zip(records, (a, b))):
            row["status"] = "failed_case_without_checkpoint_quality"
        else:
            try:
                row.update(compare_point(a, b, tol_km))
            except (TypeError, ValueError) as exc:
                row["status"] = f"invalid_statistics: {exc}"
        yield row


JUDGMENTS = ("contained", "indeterminate", "not_contained")
METRICS = [f"{box}_min_margin_km" for box, _ in BOXES] + [
    f"{box}_{axis}_width_gain_km" for box, _ in BOXES for axis in AXES
]


class Summary:
    """Checkpoint-weighted statistics. Keep invalid rows out of percentage denominators."""

    def __init__(self):
        self.total = 0
        self.counts = Counter()
        self.values = {name: array("d") for name in METRICS}
        self.faces = Counter()

    def add(self, judgment, values, face):
        self.total += 1
        if judgment is None:
            return
        self.counts[judgment] += 1
        for name, value in zip(METRICS, values):
            self.values[name].append(value)
        if judgment == "not_contained":
            self.faces[face] += 1

    def result(self):
        valid = sum(self.counts.values())
        out = dict(n_rows=self.total, n_valid=valid, n_invalid=self.total - valid)
        for label in JUDGMENTS:
            out[f"n_{label}"] = self.counts[label]
            out[f"pct_{label}"] = 100 * self.counts[label] / valid if valid else None
        for name, values in self.values.items():
            out[f"median_{name}"] = median(values) if values else None
        for axis in AXES:
            for sign in ("-", "+"):
                out[f"n_not_contained_{axis}{sign}"] = self.faces[axis + sign]
        return out


def summarize(path, outdir):
    """Read an existing detail CSV once. Day d groups checkpoints in [d, d+1)."""
    overall = Summary()
    cases = defaultdict(Summary)
    days = defaultdict(Summary)
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        required = {"orbit", "phase_index", "alpha0_rad", "t_days", "status",
                    "judgment", "large_limiting_face", *METRICS}
        missing = required.difference(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{path}: missing columns {sorted(missing)}")
        for row in reader:
            key = (row["orbit"], int(row["phase_index"]), float(row["alpha0_rad"]))
            groups = [overall, cases[key]]
            if row["t_days"]:
                groups.append(days[math.floor(float(row["t_days"]))])
            judgment, values = None, ()
            if row["status"] == "ok":
                judgment = row["judgment"]
                if judgment not in JUDGMENTS or not row["t_days"]:
                    raise ValueError(f"{path}:{reader.line_num}: invalid judgment or missing time")
                values = tuple(float(row[name]) for name in METRICS)
                if not all(map(math.isfinite, values)):
                    raise ValueError(f"{path}:{reader.line_num}: nonfinite metrics")
            for group in groups:
                group.add(judgment, values, row["large_limiting_face"])

    result = overall.result()
    tables = [
        ("overall", ["n_cases"], [((len(cases),), overall)]),
        ("by_case", ["orbit", "phase_index", "alpha0_rad"], sorted(cases.items())),
        ("by_day", ["day_start"], [((day,), days[day]) for day in sorted(days)]),
    ]
    outdir.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=outdir) as tmp:
        for suffix, labels, groups in tables:
            output = outdir / f"{path.stem}_{suffix}.csv"
            temp = Path(tmp) / output.name
            with temp.open("w", newline="", encoding="utf-8") as stream:
                writer = csv.DictWriter(stream, fieldnames=[*labels, *result])
                writer.writeheader()
                for key, group in groups:
                    writer.writerow({**dict(zip(labels, key)), **group.result()})
            temp.replace(output)
            print(output)
    print(f"{len(cases)} cases, {result['n_valid']}/{overall.total} valid rows")
    for label in JUDGMENTS:
        if result["n_valid"]:
            print(f"  {label}: {result['n_' + label]} ({result['pct_' + label]:.2f}%)")
    return 0 if result["n_valid"] else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inner-results", type=Path, help="Q1 config directory")
    parser.add_argument("--outer-results", type=Path, help="Q2 config directory")
    parser.add_argument("--summary-from-csv", type=Path, help="summarize saved judgments directly")
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--orbit-prefix")
    parser.add_argument("--solar-phase", type=float, help="solar phase in radians")
    parser.add_argument("--tol-km", type=float, help="absolute face tolerance in km (default 0)")
    args = parser.parse_args()
    if args.summary_from_csv is not None:
        if any(x is not None for x in (args.inner_results, args.outer_results, args.orbit_prefix,
                                       args.solar_phase, args.tol_km)):
            parser.error("--summary-from-csv uses the cases and judgments already saved in the CSV")
        return summarize(args.summary_from_csv, args.outdir)
    if args.inner_results is None or args.outer_results is None:
        parser.error("provide both config directories or --summary-from-csv")
    if args.tol_km is None:
        args.tol_km = 0.0
    if not math.isfinite(args.tol_km) or args.tol_km < 0:
        parser.error("--tol-km must be finite and nonnegative")
    inner = scan_config(args.inner_results, args.orbit_prefix, args.solar_phase)
    outer = scan_config(args.outer_results, args.orbit_prefix, args.solar_phase)
    keys = sorted(inner.keys() & outer.keys())
    if not keys:
        parser.error("no common orbit, phase index and solar phase in the two config directories")
    if inner.keys() != outer.keys():
        print(f"Unpaired cases: Q1 only {len(inner.keys() - outer.keys())}, "
              f"Q2 only {len(outer.keys() - inner.keys())}", file=sys.stderr)

    args.outdir.mkdir(parents=True, exist_ok=True)
    output = args.outdir / f"box_containment_checkpoints_{args.orbit_prefix or 'all'}.csv"
    rows = valid = 0
    with TemporaryDirectory(dir=args.outdir) as tmpdir:
        temp = Path(tmpdir) / output.name
        with temp.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=FIELDS)
            writer.writeheader()
            for key in keys:
                orbit, phase, alpha = key
                base = dict(orbit=orbit, phase_index=phase, alpha0_rad=alpha)
                try:
                    for row in case_rows(inner[key], outer[key], args.tol_km):
                        writer.writerow({**base, **row})
                        rows += 1
                        valid += row["status"] == "ok"
                except (AttributeError, KeyError, TypeError, ValueError) as exc:
                    writer.writerow({**base, "status": f"invalid_case: {exc}"})
                    rows += 1
        temp.replace(output)
    print(f"{valid} classified / {rows} rows across {len(keys)} paired cases -> {output}")
    return summarize(output, args.outdir)


if __name__ == "__main__":
    sys.exit(main())
