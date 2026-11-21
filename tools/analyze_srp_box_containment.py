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


def plot_results(path, outdir, days, overall, case_key, snapshots):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    colors = ("#138072", "#c6cdd1", "#ce5a35")
    day = sorted(days)
    stats = [days[d].result() for d in day]
    fig = plt.figure(figsize=(13, 7.4), layout="constrained")
    grid = fig.add_gridspec(2, 6)
    freq = fig.add_subplot(grid[0, :4])
    bottom = [0.] * len(day)
    for label, title, color in zip(JUDGMENTS, ("Contained", "Indeterminate", "Outside"), colors):
        pct = [r[f"pct_{label}"] or 0 for r in stats]
        freq.bar(day, pct, bottom=bottom, width=0.94, align="edge", color=color, label=title)
        bottom = [a + b for a, b in zip(bottom, pct)]
    freq.set(title="Daily box judgments", xlabel="Propagation day (bin start)",
             ylabel="Valid checkpoints (%)", ylim=(0, 101))
    freq.legend(ncols=3, loc="upper center", bbox_to_anchor=(0.5, -0.18), frameon=False)

    faces = [axis + sign for axis in AXES for sign in ("-", "+")]
    ax = fig.add_subplot(grid[0, 4:])
    counts = [overall.faces[f] for f in faces]
    bars = ax.barh(faces, counts, color=colors[2])
    ax.bar_label(bars, padding=3)
    ax.set_xlim(0, max(1, max(counts)) * 1.25)
    ax.set(title="Limiting face when outside", xlabel="Checkpoint count")
    ax.invert_yaxis()

    for i, axis in enumerate(AXES):
        ax = fig.add_subplot(grid[1, 2*i:2*i+2])
        # Gaps in the available daily bins should also be gaps in the curves.
        bins = range(min(day), max(day) + 1)
        daily = dict(zip(day, stats))
        for box, style in (("small", "--"), ("large", "-")):
            gains = [daily.get(d, {}).get(f"median_{box}_{axis}_width_gain_km") for d in bins]
            ax.plot([d + 0.5 for d in bins], gains, style, color="#215c84",
                    label=f"Q2 {box} box")
        ax.axhline(0, color="#666666", linewidth=0.7)
        ax.set(title=f"{axis.upper()} full-width change", xlabel="Propagation day",
               ylabel="Median Q2 minus Q1 (km)")
        if i == 0:
            ax.legend(frameon=False)
    fig.suptitle(f"Q2 containment of Q1  |  {overall.total:,} records", fontsize=16)
    output = outdir / f"{path.stem}_overview.png"
    fig.savefig(output, dpi=160)
    plt.close(fig)
    print(output)

    # Width gains in the two conventions differ only by Q2's scale factor.
    # Recover Q1 and Q2 half-widths and the center offset directly from the CSV.
    snapshots.sort(key=lambda r: float(r["t_days"]))
    selected = sorted({0, (len(snapshots)-1)//2, len(snapshots)-1})
    fig, panels = plt.subplots(2, len(selected), squeeze=False, figsize=(4.2*len(selected), 7.1),
                               layout="constrained")
    for col, index in enumerate(selected):
        row = snapshots[index]
        widths, centers = {}, {}
        for axis in AXES:
            small = float(row[f"small_{axis}_width_gain_km"]) / 2
            large = float(row[f"large_{axis}_width_gain_km"]) / 2
            sigma2 = (large-small) / (3-math.sqrt(3))
            widths[axis] = (max(0., math.sqrt(3)*sigma2-small),
                            math.sqrt(3)*sigma2, 3*sigma2)
            centers[axis] = small - float(row[f"small_{axis}_minus_margin_km"])
        for k, (u, v) in enumerate((("x", "y"), ("y", "z"))):
            ax = panels[k, col]
            for b, label, color, style in ((2, "Q2 large", "#215c84", "--"),
                                            (1, "Q2 small", "#215c84", "-"),
                                            (0, "Q1", colors[2], "-")):
                x, y = (centers[u], centers[v]) if b else (0., 0.)
                w, h = widths[u][b], widths[v][b]
                ax.add_patch(Rectangle((x-w, y-h), 2*w, 2*h, fill=False,
                                       edgecolor=color, lw=1.8, linestyle=style, label=label))
            ax.autoscale_view()
            ax.set_aspect("equal", adjustable="datalim")
            ax.set(xlabel=f"{u} offset (km)", ylabel=f"{v} offset (km)")
            if k == 0:
                ax.set_title(f"Day {float(row['t_days']):.2f}  |  {row['judgment'].replace('_',' ')}")
    orbit, phase, alpha = case_key
    fig.suptitle(f"{orbit}  |  phase {phase}  |  solar phase {alpha:g} rad\n"
                 "Estimated boxes relative to the Q1 center")
    panels[-1, 0].legend(loc="upper center", bbox_to_anchor=(0.5, -0.20), ncols=3, frameon=False)
    output = outdir / f"{path.stem}_{orbit}_p{phase:02d}_a{alpha:g}_boxes.png"
    fig.savefig(output, dpi=160)
    plt.close(fig)
    print(output)


def summarize(path, outdir, plots=False):
    """Read an existing detail CSV once. Day d groups checkpoints in [d, d+1)."""
    overall = Summary()
    cases = defaultdict(Summary)
    days = defaultdict(Summary)
    case_key = None
    snapshots = []
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        required = {"orbit", "phase_index", "alpha0_rad", "t_days", "status",
                    "judgment", "large_limiting_face", *METRICS}
        if plots:
            required.update(f"small_{axis}_minus_margin_km" for axis in AXES)
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
                if plots and (case_key is None or key <= case_key):
                    if key != case_key:
                        case_key = key
                        snapshots.clear()
                    snapshots.append(row)
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
    if plots and snapshots:
        plot_results(path, outdir, days, overall, case_key, snapshots)
    return 0 if result["n_valid"] else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inner-results", type=Path, help="Q1 config directory")
    parser.add_argument("--outer-results", type=Path, help="Q2 config directory")
    parser.add_argument("--summary-from-csv", type=Path, help="summarize saved judgments directly")
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--plots", action="store_true",
                        help="write overview and projected boxes for the first valid case")
    parser.add_argument("--orbit-prefix")
    parser.add_argument("--solar-phase", type=float, help="solar phase in radians")
    parser.add_argument("--tol-km", type=float, help="absolute face tolerance in km (default 0)")
    args = parser.parse_args()
    if args.summary_from_csv is not None:
        if any(x is not None for x in (args.inner_results, args.outer_results, args.orbit_prefix,
                                       args.solar_phase, args.tol_km)):
            parser.error("--summary-from-csv uses the cases and judgments already saved in the CSV")
        return summarize(args.summary_from_csv, args.outdir, args.plots)
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
    return summarize(output, args.outdir, args.plots)


if __name__ == "__main__":
    sys.exit(main())
