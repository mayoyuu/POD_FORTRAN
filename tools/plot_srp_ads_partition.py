#!/usr/bin/env python3
"""Plot CRTBP SRP ADS initial-parameter boxes at a requested propagation time.

Replay binary splits in the per-case JSON (call = parent center, wall = parent
FULL width, in normalized root coordinates). This is a parameter-domain plot,
not the spacecraft position distribution. No DA evaluation or propagation occurs.

Angles in --slice and figures are degrees. CSV physical angle bounds are radians.
After the first dropped/diverged patch, only the historical split partition can
be recovered: the JSON does not contain the geometry of deleted patches.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import sys

ND_TO_DAYS = 375237.080249162 / 86400.0
ROOT = Path(__file__).resolve().parents[1]


def finite(value):
    value = float(value)
    if not math.isfinite(value):
        raise ValueError("non-finite number in result or command line")
    return value


def domain_metadata(result):
    names = result.get("parameter_names")
    lower = result.get("parameter_lower")
    upper = result.get("parameter_upper")
    if names is None or lower is None or upper is None:
        if result.get("attitude_coordinates") == "theta_phi":
            raise ValueError("theta_phi result is missing parameter domain metadata")
        # Compatibility with legacy u/v and cannonball outputs only.
        q = math.tan(math.radians(finite(result.get("tilt_max_deg", 0)))) / math.sqrt(2)
        names = ["beta", "chi", "u", "v"]
        lower = [result["beta_lower"], result.get("chi_lower", .6), -q, -q]
        upper = [result["beta_upper"], result.get("chi_upper", .6), q, q]
    if len(names) != 4 or len(lower) != 4 or len(upper) != 4:
        raise ValueError("expected four physical parameter slots")
    if list(names) not in (["beta", "chi", "u", "v"], ["beta", "chi", "theta", "phi"]):
        raise ValueError(f"unsupported parameter order: {names}")
    lower, upper = list(map(finite, lower)), list(map(finite, upper))
    if any(a > b for a, b in zip(lower, upper)):
        raise ValueError("reversed parameter bounds")
    units = result.get("parameter_units", ["m2/kg", "1", "1", "1"])
    if "theta" in names and units != ["m2/kg", "1", "rad", "rad"]:
        raise ValueError("theta/phi JSON bounds must be in radians")
    return list(names), lower, upper


def box_key(center, width, active):
    """Canonical dyadic cells tolerate the Fortran ES20.12 serialization.

    A root cell has level=0, index=0. Index/level matching avoids an O(N^2)
    search over leaves and does not use rounded printed centers as exact keys.
    """
    if len(center) != 4 or len(width) != 4:
        raise ValueError("call/wall must contain four physical slots")
    key = []
    for c, w, enabled in zip(center, width, active):
        c, w = finite(c), finite(w)
        if not enabled:
            if abs(c) > 1e-11 or abs(w) > 1e-11:
                raise ValueError("nonzero geometry for a fixed parameter")
            key.append((-1, 0))
            continue
        if w <= 0:
            raise ValueError("non-positive width for an active parameter")
        level = round(math.log2(2.0 / w))
        if not 0 <= level <= 50:
            raise ValueError("unsupported dyadic box depth")
        exact_width = 2.0 ** (1 - level)
        index = round((c + 1) / exact_width - .5)
        exact_center = -1 + (index + .5) * exact_width
        if (not 0 <= index < 2 ** level
                or not math.isclose(w, exact_width, rel_tol=2e-11, abs_tol=0)
                or abs(c - exact_center) > min(2e-11, exact_width / 8)):
            raise ValueError("invalid or insufficiently precise dyadic call/wall geometry")
        key.append((level, index))
    return tuple(key)


def make_box(key, birth_days):
    return {
        "center": [0.0 if level < 0 else -1 + (index + .5) * 2.0 ** (1 - level)
                   for level, index in key],
        "width": [0.0 if level < 0 else 2.0 ** (1 - level) for level, _ in key],
        "depth": sum(max(0, level) for level, _ in key),
        "birth_days": birth_days,
    }


def reconstruct(result, time_days):
    """Return the split partition at t, including all events recorded at t."""
    names, lower, upper = domain_metadata(result)
    active = [b > a for a, b in zip(lower, upper)]
    end_days = finite(result["t_end_days"])
    time_days = finite(time_days)
    if time_days < 0 or time_days > end_days + 1e-10:
        raise ValueError(f"--time-days must be in [0, {end_days:g}]")
    nd_days = end_days / finite(result["t_end_nd"]) if result.get("t_end_nd", 0) > 0 else ND_TO_DAYS
    # Event td uses only seven significant digits. Include equality to that
    # precision; the file cannot resolve sub-second ordering near a boundary.
    time_tolerance = 5e-7 * max(1.0, time_days)
    events = result.get("split_events")
    expected = int(result["n_splits_total"])
    if not isinstance(events, list) or len(events) != expected:
        raise ValueError("incomplete split_events: count differs from n_splits_total")
    events = sorted(events, key=lambda e: int(e["i"]))
    if [int(e["i"]) for e in events] != list(range(1, expected + 1)):
        raise ValueError("missing or duplicate split event indices")
    n_diverged = int(result.get("n_diverged", 0))
    if n_diverged < 0:
        raise ValueError("negative n_diverged")
    if "n_patches" in result and int(result["n_patches"]) + n_diverged != 1 + expected:
        raise ValueError("n_patches + n_diverged != 1 + n_splits_total; incomplete result")
    root_key = tuple((0, 0) if a else (-1, 0) for a in active)
    leaves = {root_key: make_box(root_key, 0.0)}
    applied = 0
    previous_time = -1.0
    for event in events:
        if "call" not in event or "wall" not in event:
            raise ValueError("split log lacks call/wall; full parameter boxes cannot be recovered")
        td = finite(event["td"]) if "td" in event else finite(event["t"]) * nd_days
        if td < 0 or td < previous_time or td > end_days + 5e-7 * max(1, end_days):
            raise ValueError("event times are out of order or outside the propagation arc")
        previous_time = td
        axis = int(event["ax"]) - 1
        if axis not in range(4) or not active[axis]:
            raise ValueError("event split axis is not an active physical parameter")
        key = box_key(event["call"], event["wall"], active)
        if td > time_days + time_tolerance:
            continue
        if key not in leaves:
            raise ValueError(f"event {event['i']}: parent is not a current leaf; incomplete/inconsistent log")
        parent = leaves.pop(key)
        if "d" in event and int(event["d"]) != parent["depth"]:
            raise ValueError(f"event {event['i']}: parent depth disagrees with call/wall")
        level, index = key[axis]
        for side in (0, 1):
            child = list(key)
            child[axis] = (level + 1, 2 * index + side)
            child = tuple(child)
            leaves[child] = make_box(child, td)
        applied += 1
    first_div = result.get("first_divergence_time")
    first_div_days = finite(first_div) * nd_days if first_div is not None and finite(first_div) >= 0 else None
    historical = n_diverged > 0 and (
        first_div_days is None or first_div_days <= time_days + time_tolerance)
    return {
        "parameter_names": names, "parameter_lower": lower, "parameter_upper": upper,
        "time_days": time_days, "time_tolerance_days": time_tolerance,
        "last_split_days": max((b["birth_days"] for b in leaves.values()), default=0),
        "events_applied": applied, "partition_kind": "historical" if historical else "active",
        "first_divergence_days": first_div_days, "n_diverged_end": n_diverged,
        "boxes": [leaves[key] for key in sorted(leaves)],
    }


def physical_bounds(snapshot, box):
    middle = [(a + b) / 2 for a, b in zip(snapshot["parameter_lower"], snapshot["parameter_upper"])]
    half = [(b - a) / 2 for a, b in zip(snapshot["parameter_lower"], snapshot["parameter_upper"])]
    lower = [m + h * (c - w / 2) for m, h, c, w in zip(middle, half, box["center"], box["width"])]
    upper = [m + h * (c + w / 2) for m, h, c, w in zip(middle, half, box["center"], box["width"])]
    return lower, upper


def select_slice(snapshot, axes, values):
    """Select hidden-coordinate slice; boundary ties belong to the right child.

    All unspecified hidden parameters use the root midpoint. Angle values enter
    in degrees and are returned in radians. A global maximum is included.
    """
    names = snapshot["parameter_names"]
    lower, upper = snapshot["parameter_lower"], snapshot["parameter_upper"]
    if len(axes) != 2 or axes[0] == axes[1] or any(a not in names for a in axes):
        raise ValueError("--axes requires two distinct parameter names in the result")
    if any(lower[names.index(a)] == upper[names.index(a)] for a in axes):
        raise ValueError("plot axes must both be uncertain (nonzero domain width)")
    if any(n not in names or n in axes for n in values):
        raise ValueError("--slice accepts only hidden parameter names")
    fixed = {}
    normalized = {}
    for j, name in enumerate(names):
        if name in axes:
            continue
        value = (lower[j] + upper[j]) / 2
        if name in values:
            value = finite(values[name])
            if name in ("theta", "phi"):
                value = math.radians(value)
        if value < lower[j] - 1e-13 or value > upper[j] + 1e-13:
            raise ValueError(f"slice {name} is outside the root parameter domain")
        value = min(upper[j], max(lower[j], value))
        fixed[name] = value
        if upper[j] > lower[j]:
            normalized[j] = 2 * (value - lower[j]) / (upper[j] - lower[j]) - 1
    selected = []
    for box in snapshot["boxes"]:
        inside = True
        for j, xi in normalized.items():
            lo = box["center"][j] - box["width"][j] / 2
            hi = box["center"][j] + box["width"][j] / 2
            # Snap tiny conversion noise to the nearest boundary, then use
            # half-open intervals (except the global maximum).
            if abs(xi - lo) < 1e-12:
                xi = lo
            if abs(xi - hi) < 1e-12:
                xi = hi
            if not (lo <= xi < hi or xi == hi == 1.0):
                inside = False
                break
        if inside:
            selected.append(box)
    return selected, fixed


def write_boxes(snapshot, selected, path):
    selected_ids = {id(b) for b in selected}
    rows = []
    for i, box in enumerate(snapshot["boxes"], 1):
        lo, hi = physical_bounds(snapshot, box)
        row = {"box_id": i, "depth": box["depth"], "birth_days": box["birth_days"],
               "partition_kind": snapshot["partition_kind"], "in_slice": int(id(box) in selected_ids)}
        for j, name in enumerate(snapshot["parameter_names"]):
            suffix = "_rad" if name in ("theta", "phi") else ""
            row[f"{name}_lower{suffix}"] = lo[j]
            row[f"{name}_upper{suffix}"] = hi[j]
            row[f"{name}_center_xi"] = box["center"][j]
            row[f"{name}_width_xi"] = box["width"][j]
        rows.append(row)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def plot_boxes(snapshot, selected, axes, fixed, path, case_label, projection=False):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.collections import PatchCollection
        from matplotlib.patches import Rectangle
    except ImportError as exc:
        raise ValueError("PNG requires matplotlib: python3 -m pip install matplotlib; "
                         "or pass --no-png for CSV/JSON only") from exc
    indices = [snapshot["parameter_names"].index(a) for a in axes]
    scale = [180 / math.pi if a in ("theta", "phi") else 1 for a in axes]
    rectangles = []
    for box in selected:
        lo, hi = physical_bounds(snapshot, box)
        x, y = [lo[j] * s for j, s in zip(indices, scale)]
        dx, dy = [(hi[j] - lo[j]) * s for j, s in zip(indices, scale)]
        rectangles.append(Rectangle((x, y), dx, dy))
    fig, ax = plt.subplots(figsize=(10, 7), layout="constrained")
    collection = PatchCollection(rectangles, facecolors="none", edgecolors="#174c78", linewidths=.7)
    ax.add_collection(collection)
    for setter, j, s in zip((ax.set_xlim, ax.set_ylim), indices, scale):
        setter(snapshot["parameter_lower"][j] * s, snapshot["parameter_upper"][j] * s)
    for setter, name in zip((ax.set_xlabel, ax.set_ylabel), axes):
        unit = "deg" if name in ("theta", "phi") else ("m2/kg" if name == "beta" else "1")
        setter(f"{name} [{unit}]")
    fixed_text = ", ".join(f"{n}={math.degrees(v) if n in ('theta', 'phi') else v:.6g}"
                           + (" deg" if n in ("theta", "phi") else "") for n, v in fixed.items())
    mode = "2-D projection (boxes may overlap)" if projection else f"Slice: {fixed_text}"
    state = ("HISTORICAL split partition; survival unknown"
             if snapshot["partition_kind"] == "historical" else "Active parameter partition")
    ax.set_title(f"{case_label}\nt = {snapshot['time_days']:g} days | {state}\n"
                 f"{mode} | {len(selected)} displayed / {len(snapshot['boxes'])} boxes", fontsize=10)
    ax.grid(alpha=.15)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--json", required=True, type=Path, help="one per-case CRTBP study result")
    parser.add_argument("--time-days", required=True, type=float)
    parser.add_argument("--axes", nargs=2, default=["theta", "phi"])
    parser.add_argument("--slice", nargs="*", default=[], metavar="NAME=VALUE",
                        help="hidden parameters; angles in degrees; default=root midpoint")
    parser.add_argument("--projection", action="store_true", help="project all boxes; overlaps are expected")
    parser.add_argument("--require-active", action="store_true", help="reject queries after any patch deletion")
    parser.add_argument("--outdir", type=Path, default=ROOT / "results/srp_sensitivity/analysis")
    parser.add_argument("--no-png", action="store_true", help="write CSV/JSON without matplotlib")
    args = parser.parse_args(argv)
    values = {}
    for item in args.slice:
        name, separator, value = item.partition("=")
        if not separator or name in values:
            raise ValueError("use unique --slice NAME=VALUE entries")
        values[name] = finite(value)
    if args.projection and values:
        raise ValueError("--projection and --slice cannot be combined")
    result = json.loads(args.json.read_text(encoding="utf-8-sig"))
    snapshot = reconstruct(result, args.time_days)
    selected, fixed = select_slice(snapshot, args.axes, values)
    if args.projection:
        selected, fixed = snapshot["boxes"], {}
    if snapshot["partition_kind"] == "historical":
        warning = ("Only the historical split partition is recoverable: deleted patch bounds "
                   "were not logged. These boxes must not be interpreted as surviving patches.")
        if args.require_active:
            raise ValueError(warning)
        print("WARNING: " + warning, file=sys.stderr)
    stem = f"partition_{args.json.stem}_t{args.time_days:g}d_{'_'.join(args.axes)}"
    if args.projection:
        stem += "_projection"
    else:
        stem += "_slice_" + "_".join(
            f"{n}{math.degrees(v) if n in ('theta', 'phi') else v:.12g}" for n, v in fixed.items())
    # Keep separate cases in separate directories, even when --outdir is shared.
    out = args.outdir / args.json.parent.parent.name / args.json.parent.name
    out.mkdir(parents=True, exist_ok=True)
    write_boxes(snapshot, selected, out / f"{stem}_boxes.csv")
    report = {k: v for k, v in snapshot.items() if k != "boxes"}
    report.update(source_json=str(args.json.resolve()), axes=args.axes, slice_physical=fixed,
                  mode="projection" if args.projection else "slice",
                  n_boxes=len(snapshot["boxes"]), n_displayed=len(selected),
                  physical_angle_units="rad", plot_angle_units="deg")
    (out / f"{stem}_summary.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    if not args.no_png:
        label = f"{args.json.parent.parent.name}\n{args.json.parent.name} / {args.json.stem}"
        plot_boxes(snapshot, selected, args.axes, fixed, out / f"{stem}.png", label, args.projection)
    print(f"{snapshot['partition_kind']}: {len(snapshot['boxes'])} boxes, "
          f"{len(selected)} displayed, {snapshot['events_applied']} split events -> {out}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, KeyError, TypeError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(2)
