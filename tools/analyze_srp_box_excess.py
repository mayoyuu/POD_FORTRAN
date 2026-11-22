#!/usr/bin/env python3
"""Summarize Q1 excursions beyond the Q2 large box from saved checkpoint CSVs.

Depths refer to the face with the greatest absolute excess. Relative depth
uses the large Q2 half-width normal to that same face. Percentiles describe
only checkpoints previously classified as not_contained.
"""
import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


class Stats:
    def __init__(self):
        self.total = 0
        self.valid = 0
        self.depth = []
        self.percent = []

    def add(self, valid, depth, percent):
        self.total += 1
        self.valid += int(valid)
        if depth is not None:
            self.depth.append(depth)
        if percent is not None:
            self.percent.append(percent)

    def result(self):
        n = len(self.depth)
        out = dict(n_rows=self.total, n_valid=self.valid, n_outside=n,
                   outside_pct=100*n/self.valid if self.valid else None,
                   n_relative=len(self.percent))
        for unit, sample in (("km", self.depth), ("pct", self.percent)):
            ordered = sorted(sample)
            for name, p in (("median", 0.5), ("p90", 0.9)):
                value = None
                if ordered:
                    j = (len(ordered)-1)*p
                    lo = math.floor(j)
                    value = ordered[lo] + (j-lo)*(ordered[math.ceil(j)]-ordered[lo])
                out[f"{name}_depth_{unit}"] = value
        return out


def depths(row):
    if row["judgment"] != "not_contained":
        return None, None
    face = row["large_limiting_face"]
    if face not in ("x-", "x+", "y-", "y+", "z-", "z+"):
        raise ValueError(f"unknown limiting face {face!r}")
    axis = face[0]
    depth = -float(row["large_min_margin_km"])
    # Width gains use full widths. Subtraction cancels the Q1 contribution.
    wide = float(row[f"large_{axis}_width_gain_km"])
    narrow = float(row[f"small_{axis}_width_gain_km"])
    half = 3*(wide-narrow)/(2*(3-math.sqrt(3)))
    if not all(map(math.isfinite, (depth, half))) or depth <= 0 or half < 0:
        raise ValueError("inconsistent margin or box widths")
    return depth, 100*depth/half if half > 0 else None


def write_rows(path, rows):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)
    print(path)


def plot(path, rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, panels = plt.subplots(3, 1, figsize=(10, 9), sharex=True, layout="constrained")
    fields = ("outside_pct", "median_depth_km", "median_depth_pct")
    labels = ("Outside / valid (%)", "Median excess (km)", "Median excess / Q2 half-width (%)")
    for family in sorted({row["family"] for row in rows}):
        by_day = {row["day_start"]: row for row in rows if row["family"] == family}
        days = range(min(by_day), max(by_day)+1)
        for ax, field in zip(panels, fields):
            ax.plot([d+0.5 for d in days], [by_day.get(d, {}).get(field) for d in days],
                    marker=".", label=family)
    for ax, label in zip(panels, labels):
        ax.set_ylabel(label)
        ax.grid(axis="y", alpha=0.25)
    panels[0].legend(ncols=4, frameon=False)
    panels[-1].set_xlabel("Propagation time (days, daily bins)")
    fig.suptitle("Q1 excursions beyond the estimated Q2 large box")
    fig.savefig(path, dpi=170)
    plt.close(fig)
    print(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path, nargs="+", help="checkpoint detail CSVs from one Q1/Q2 pair of configs")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--plots", action="store_true")
    args = parser.parse_args()
    families, daily = defaultdict(Stats), defaultdict(Stats)
    seen = set()
    for path in args.csv:
        with path.open(encoding="utf-8-sig", newline="") as stream:
            reader = csv.DictReader(stream)
            required = {"orbit", "phase_index", "alpha0_rad", "t_nd", "t_days", "status",
                        "judgment", "large_limiting_face", "large_min_margin_km"}
            required.update(f"{box}_{axis}_width_gain_km"
                            for box in ("small", "large") for axis in "xyz")
            missing = required.difference(reader.fieldnames or ())
            if missing:
                parser.error(f"{path}: missing checkpoint columns {sorted(missing)}")
            for row in reader:
                family = row["orbit"].split("_A")[0]
                valid = row["status"] == "ok"
                try:
                    day = math.floor(float(row["t_days"])) if row["t_days"] else None
                    if valid and (day is None or row["judgment"] not in
                                  ("contained", "indeterminate", "not_contained")):
                        raise ValueError("missing time or invalid judgment")
                    if row["t_nd"]:
                        key = (row["orbit"], int(row["phase_index"]), float(row["alpha0_rad"]), float(row["t_nd"]))
                        if key in seen:
                            raise ValueError(f"duplicate checkpoint {key}; supply nonoverlapping CSVs")
                        seen.add(key)
                    depth, pct = depths(row) if valid else (None, None)
                except (TypeError, ValueError, OverflowError) as exc:
                    parser.error(f"{path}:{reader.line_num}: {exc}")
                families[family].add(valid, depth, pct)
                if day is not None:
                    daily[family, day].add(valid, depth, pct)
    if not daily:
        parser.error("no timed checkpoint rows found")
    args.outdir.mkdir(parents=True, exist_ok=True)
    write_rows(args.outdir / "excess_by_family.csv",
               [dict(family=f, **s.result()) for f, s in sorted(families.items())])
    rows = [dict(family=f, day_start=d, **s.result()) for (f, d), s in sorted(daily.items())]
    write_rows(args.outdir / "excess_by_family_day.csv", rows)
    if args.plots:
        plot(args.outdir / "excess_by_family_day.png", rows)


if __name__ == "__main__":
    main()
