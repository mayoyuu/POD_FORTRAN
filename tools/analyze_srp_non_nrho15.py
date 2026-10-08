#!/usr/bin/env python3
"""Error-cloud statistics and first-split comparison for the 15 non-L2NRHO orbits.

Implements sections 1-5 of 2026-09-28-srp-cloud-statistics-design.md:

  * Per-group, per-hour error-cloud moments over the 512 fixed design points:
    mean offset from the origin orbit, longest covariance RMS axis, maximum
    absolute skewness, maximum absolute excess kurtosis, and the maximum paired
    treatment-control point difference relative to each cloud's longest RMS axis.
  * First-split comparison derived from the cumulative ``split_dim*`` counters.
  * Point-to-point SRP increments re-reported in km and km/s.

Only the continuous prefix where BOTH the treatment and control groups pass ADS
validation (no ``validation_failed`` / ``precision_unconfirmed``) enters the
figures, peaks, population quantiles and integrals; the first invalid checkpoint
and every later one are excluded. Raw values are retained for auditability.

Position quantities are in km and velocity quantities in km/s throughout; no
m / mm/s fields are emitted.
"""
from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ORBIT_NAMES = [
    "DRO-1", "DRO-2", "DRO-3",
    "L1Halo-1", "L1Halo-2", "L1Halo-3",
    "L1PL-1", "L1PL-2", "L1PL-3",
    "L2Halo-1", "L2Halo-2", "L2Halo-3",
    "L2PL-1", "L2PL-2", "L2PL-3",
]
COMPONENTS = ("dx", "dy", "dz", "dvx", "dvy", "dvz")
POSITION = ("dx", "dy", "dz")
VELOCITY = ("dvx", "dvy", "dvz")
ALIASES = {
    "hour": ("hour", "time_hour", "time_hours", "elapsed_hour", "elapsed_hours", "propagation_hour"),
    "point_id": ("point_id", "shape_point_id", "sample_id", "design_point_id", "id"),
    "dx": ("dx_km", "err_x_km", "error_x_km", "x_error_km", "delta_x_km", "x_err_km", "dx", "err_x", "x_error", "delta_x", "x_km"),
    "dy": ("dy_km", "err_y_km", "error_y_km", "y_error_km", "delta_y_km", "y_err_km", "dy", "err_y", "y_error", "delta_y", "y_km"),
    "dz": ("dz_km", "err_z_km", "error_z_km", "z_error_km", "delta_z_km", "z_err_km", "dz", "err_z", "z_error", "delta_z", "z_km"),
    "dvx": ("dvx_kms", "dvx_km_s", "err_vx_km_s", "error_vx_km_s", "vx_error_km_s", "delta_vx_km_s", "vx_err_km_s", "dvx", "err_vx", "vx_error", "delta_vx", "vx_km_s"),
    "dvy": ("dvy_kms", "dvy_km_s", "err_vy_km_s", "error_vy_km_s", "vy_error_km_s", "delta_vy_km_s", "vy_err_km_s", "dvy", "err_vy", "vy_error", "delta_vy", "vy_km_s"),
    "dvz": ("dvz_kms", "dvz_km_s", "err_vz_km_s", "error_vz_km_s", "vz_error_km_s", "delta_vz_km_s", "vz_err_km_s", "dvz", "err_vz", "vz_error", "delta_vz", "vz_km_s"),
}
HISTORY_COLUMNS = (
    "split_dim1", "split_dim2", "split_dim3", "split_dim4", "split_dim5",
    "split_dim6", "split_dim7", "validation_failed", "precision_unconfirmed",
)
CLOUD_METRIC_COLUMNS = (
    "pos_mean_offset_km", "vel_mean_offset_km_s",
    "pos_longest_rms_axis_km", "vel_longest_rms_axis_km_s",
    "skew_dx", "skew_dy", "skew_dz", "skew_dvx", "skew_dvy", "skew_dvz",
    "skew_max_abs", "skew_max_component", "skew_max_signed",
    "kurt_dx", "kurt_dy", "kurt_dz", "kurt_dvx", "kurt_dvy", "kurt_dvz",
    "kurt_max_abs", "kurt_max_component", "kurt_max_signed",
    "rel_max_pos", "rel_max_pos_num_km", "rel_max_pos_den_km",
    "rel_max_vel", "rel_max_vel_num_km_s", "rel_max_vel_den_km_s",
)
SUMMARY_METRICS = (
    "pos_mean_offset_km", "vel_mean_offset_km_s",
    "pos_longest_rms_axis_km", "vel_longest_rms_axis_km_s",
    "skew_max_abs", "kurt_max_abs",
    "rel_max_pos", "rel_max_vel",
)
DELTA_METRICS = (
    "delta_r_max_km", "delta_r_rms_km", "delta_r_p95_km",
    "delta_v_max_km_s", "delta_v_rms_km_s", "delta_v_p95_km_s",
)

# A variance below this threshold (std < 1e-12) is treated as "near zero", so
# the corresponding skewness / kurtosis / relative-error ratio is not computable.
NEAR_ZERO_VARIANCE = 1e-24


class AnalysisError(RuntimeError):
    pass


def norm_name(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", name.strip().lower()).strip("_")


def find_column(columns: Iterable[str], logical: str) -> str:
    lookup = {norm_name(c): c for c in columns}
    for alias in ALIASES[logical]:
        if norm_name(alias) in lookup:
            return lookup[norm_name(alias)]
    raise AnalysisError(
        f"Cannot identify {logical!r}. Available columns: " + ", ".join(map(str, columns))
    )


def find_one(directory: Path, suffix: str) -> Path:
    matches = sorted(directory.glob(f"*{suffix}"))
    if len(matches) != 1:
        raise AnalysisError(
            f"Expected one *{suffix} in {directory}, found {len(matches)}: "
            + ", ".join(p.name for p in matches)
        )
    return matches[0]


def read_shape(path: Path) -> pd.DataFrame:
    header = pd.read_csv(path, nrows=0).columns.tolist()
    schema = {name: find_column(header, name) for name in ("hour", "point_id", *COMPONENTS)}
    frame = pd.read_csv(path, usecols=list(schema.values()))
    frame = frame.rename(columns={actual: logical for logical, actual in schema.items()})
    frame["hour"] = pd.to_numeric(frame["hour"], errors="raise")
    frame["point_id"] = pd.to_numeric(frame["point_id"], errors="raise").astype(np.int64)
    for name in COMPONENTS:
        frame[name] = pd.to_numeric(frame[name], errors="raise")
    return frame


def read_history(path: Path, prefix: str) -> pd.DataFrame:
    header = pd.read_csv(path, nrows=0).columns.tolist()
    hour = find_column(header, "hour")
    normalized = {norm_name(c): c for c in header}
    selected = [hour]
    for candidate in HISTORY_COLUMNS:
        actual = normalized.get(norm_name(candidate))
        if actual is not None:
            selected.append(actual)
    frame = pd.read_csv(path, usecols=list(dict.fromkeys(selected))).rename(columns={hour: "hour"})
    frame = frame.rename(columns={c: prefix + norm_name(c) for c in frame.columns if c != "hour"})
    if frame["hour"].duplicated().any():
        raise AnalysisError(f"Duplicate hours in {path}")
    return frame


def read_initial_design(path: Path) -> pd.DataFrame:
    header = pd.read_csv(path, nrows=0).columns.tolist()
    point = find_column(header, "point_id")
    normalized = {norm_name(c): c for c in header}
    ucols = [normalized[f"u{i}"] for i in range(1, 7) if f"u{i}" in normalized]
    if len(ucols) != 6:
        raise AnalysisError(f"{path} does not contain all u1...u6 columns")
    frame = pd.read_csv(path, usecols=[point, *ucols])
    frame = frame.rename(columns={point: "point_id", **{c: f"u{i+1}" for i, c in enumerate(ucols)}})
    frame["point_id"] = pd.to_numeric(frame["point_id"], errors="raise").astype(np.int64)
    return frame.sort_values("point_id").reset_index(drop=True)


def verify_initial_design(treatment: Path, control: Path, expected_points: int) -> float:
    left, right = read_initial_design(treatment), read_initial_design(control)
    if len(left) != expected_points or len(right) != expected_points:
        raise AnalysisError(
            f"Initial design must have {expected_points} points; treatment={len(left)}, control={len(right)}"
        )
    merged = left.merge(right, on="point_id", suffixes=("_t", "_c"), validate="one_to_one")
    if len(merged) != expected_points:
        raise AnalysisError("Treatment/control initial point_id sets differ")
    maximum = max(float(np.max(np.abs(merged[f"u{i}_t"] - merged[f"u{i}_c"]))) for i in range(1, 7))
    if maximum > 1.0e-12:
        raise AnalysisError(f"Treatment/control u1...u6 mismatch: {maximum:.6e}")
    return maximum


def rms(values: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.square(values))))


def safe_ratio(numerator: float, denominator: float) -> float:
    if not np.isfinite(denominator) or abs(denominator) <= 1e-12:
        return math.nan
    return float(numerator / denominator)


def longest_rms_axis(values: np.ndarray) -> float:
    """Longest RMS semiaxis of a fixed-design point cloud (N denominator)."""
    centered = values - np.mean(values, axis=0)
    covariance = centered.T @ centered / max(len(values), 1)
    eigenvalues = np.linalg.eigvalsh(covariance)
    return float(np.sqrt(max(float(np.max(eigenvalues)), 0.0)))


def skewness(x: np.ndarray) -> float:
    mu = float(np.mean(x))
    centered = x - mu
    m2 = float(np.mean(np.square(centered)))
    if m2 <= NEAR_ZERO_VARIANCE:
        return math.nan
    m3 = float(np.mean(np.power(centered, 3)))
    return m3 / (m2 ** 1.5)


def excess_kurtosis(x: np.ndarray) -> float:
    mu = float(np.mean(x))
    centered = x - mu
    m2 = float(np.mean(np.square(centered)))
    if m2 <= NEAR_ZERO_VARIANCE:
        return math.nan
    m4 = float(np.mean(np.power(centered, 4)))
    return m4 / (m2 ** 2) - 3.0


def max_abs_component(moments: dict[str, float]) -> tuple[float, str, float]:
    finite = {k: v for k, v in moments.items() if np.isfinite(v)}
    if not finite:
        return math.nan, "", math.nan
    component = max(finite, key=lambda k: abs(finite[k]))
    return abs(finite[component]), component, finite[component]


def attach_paired_relative_max(
    treatment_cloud: pd.DataFrame,
    control_cloud: pd.DataFrame,
    delta: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Use paired treatment-control maxima as the numerator for both clouds.

    At each hour the shared position numerator is
    max_i ||r_treatment_i - r_control_i|| over matching point_id values.
    It is divided separately by the treatment and control position longest RMS
    axes. Velocity uses the identical construction in km/s.
    """
    numerator = delta[["hour", "delta_r_max_km", "delta_v_max_km_s"]]

    def attach(cloud: pd.DataFrame) -> pd.DataFrame:
        result = cloud.merge(numerator, on="hour", how="left", validate="one_to_one")
        result["rel_max_pos_num_km"] = result["delta_r_max_km"]
        result["rel_max_pos_den_km"] = result["pos_longest_rms_axis_km"]
        result["rel_max_pos"] = [
            safe_ratio(num, den)
            for num, den in zip(result["rel_max_pos_num_km"], result["rel_max_pos_den_km"])
        ]
        result["rel_max_vel_num_km_s"] = result["delta_v_max_km_s"]
        result["rel_max_vel_den_km_s"] = result["vel_longest_rms_axis_km_s"]
        result["rel_max_vel"] = [
            safe_ratio(num, den)
            for num, den in zip(result["rel_max_vel_num_km_s"], result["rel_max_vel_den_km_s"])
        ]
        return result.drop(columns=["delta_r_max_km", "delta_v_max_km_s"])

    return attach(treatment_cloud), attach(control_cloud)


def cloud_hour(hour: float, group: pd.DataFrame) -> dict:
    pos = group[list(POSITION)].to_numpy(float)
    vel = group[list(VELOCITY)].to_numpy(float)
    row: dict = {"hour": float(hour)}
    row["pos_mean_offset_km"] = float(np.linalg.norm(np.mean(pos, axis=0)))
    row["vel_mean_offset_km_s"] = float(np.linalg.norm(np.mean(vel, axis=0)))
    row["pos_longest_rms_axis_km"] = longest_rms_axis(pos)
    row["vel_longest_rms_axis_km_s"] = longest_rms_axis(vel)
    skew: dict[str, float] = {}
    kurt: dict[str, float] = {}
    for name, column in zip(COMPONENTS, group[list(COMPONENTS)].to_numpy(float).T):
        skew[name] = skewness(column)
        kurt[name] = excess_kurtosis(column)
        row[f"skew_{name}"] = skew[name]
        row[f"kurt_{name}"] = kurt[name]
    row["skew_max_abs"], row["skew_max_component"], row["skew_max_signed"] = max_abs_component(skew)
    row["kurt_max_abs"], row["kurt_max_component"], row["kurt_max_signed"] = max_abs_component(kurt)
    # Paired treatment-control numerators are attached after both clouds have
    # been aligned by (hour, point_id); they cannot be computed from one group.
    return row


def cloud_statistics(shape: pd.DataFrame) -> pd.DataFrame:
    records = [cloud_hour(hour, group) for hour, group in shape.groupby("hour", sort=True)]
    return pd.DataFrame.from_records(records)


def compute_validity(t_history: pd.DataFrame, c_history: pd.DataFrame) -> pd.DataFrame:
    merged = t_history.merge(c_history, on="hour", how="left", validate="one_to_one")
    invalid = pd.Series(False, index=merged.index)
    for column in ("t_validation_failed", "t_precision_unconfirmed",
                   "c_validation_failed", "c_precision_unconfirmed"):
        if column in merged:
            invalid = invalid | pd.to_numeric(merged[column], errors="coerce").fillna(0).ne(0)
    valid = (~invalid.cummax()).astype(np.int8)
    return pd.DataFrame({"hour": merged["hour"], "valid": valid})


def pair_shapes(treatment: pd.DataFrame, control: pd.DataFrame) -> pd.DataFrame:
    keys = ["hour", "point_id"]
    if treatment.duplicated(keys).any() or control.duplicated(keys).any():
        raise AnalysisError("Duplicate (hour, point_id) rows in a shape file")
    left = treatment[keys + list(COMPONENTS)]
    right = control[keys + list(COMPONENTS)]
    paired = left.merge(right, on=keys, suffixes=("_t", "_c"), validate="one_to_one")
    if len(paired) != len(left) or len(paired) != len(right):
        raise AnalysisError(
            f"Shape key sets differ: treatment={len(left)}, control={len(right)}, paired={len(paired)}"
        )
    for component in COMPONENTS:
        paired[f"delta_{component}"] = paired[f"{component}_t"] - paired[f"{component}_c"]
    paired["delta_r_km"] = np.linalg.norm(paired[["delta_dx", "delta_dy", "delta_dz"]], axis=1)
    paired["delta_v_km_s"] = np.linalg.norm(paired[["delta_dvx", "delta_dvy", "delta_dvz"]], axis=1)
    return paired


def delta_by_hour(paired: pd.DataFrame) -> pd.DataFrame:
    records = []
    for hour, group in paired.groupby("hour", sort=True):
        dr = group["delta_r_km"].to_numpy(float)
        dv = group["delta_v_km_s"].to_numpy(float)
        records.append({
            "hour": float(hour),
            "delta_r_max_km": float(np.max(dr)),
            "delta_r_rms_km": rms(dr),
            "delta_r_p95_km": float(np.quantile(dr, 0.95)),
            "delta_v_max_km_s": float(np.max(dv)),
            "delta_v_rms_km_s": rms(dv),
            "delta_v_p95_km_s": float(np.quantile(dv, 0.95)),
        })
    return pd.DataFrame.from_records(records)


def dim_number(column: str) -> int:
    match = re.search(r"split_dim(\d+)$", column)
    return int(match.group(1)) if match else -1


def first_split(history: pd.DataFrame, prefix: str, dimensions: Iterable[int]) -> dict:
    columns = [f"{prefix}split_dim{i}" for i in dimensions if f"{prefix}split_dim{i}" in history]
    if not columns:
        return {"hour": math.nan, "interval_start": math.nan,
                "interval_end": math.nan, "directions": [], "ever_nonzero": False}
    values = history[columns].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    if not values.gt(0.0).any(axis=1).any():
        return {"hour": math.nan, "interval_start": math.nan,
                "interval_end": math.nan, "directions": [], "ever_nonzero": False}
    first_index = int(np.flatnonzero(values.gt(0.0).any(axis=1).to_numpy())[0])
    hour = float(history["hour"].to_numpy()[first_index])
    first_row = values.iloc[first_index]
    directions = sorted(dim_number(c) for c in columns if float(first_row[c]) > 0.0)
    interval_start = 0.0 if hour == 0.0 else hour - 1.0
    return {"hour": hour, "interval_start": interval_start,
            "interval_end": hour, "directions": directions, "ever_nonzero": True}


def dim7_ever_nonzero(history: pd.DataFrame) -> bool:
    column = "c_split_dim7"
    if column not in history:
        return False
    return bool(pd.to_numeric(history[column], errors="coerce").fillna(0).gt(0).any())


def summarize_orbit(orbit: str, cloud: pd.DataFrame, validity: pd.DataFrame) -> dict:
    valid_cloud = cloud.loc[cloud["valid"].eq(1)].copy()
    all_hours = sorted(cloud["hour"].unique())
    invalid_hours = sorted(validity.loc[validity["valid"].eq(0), "hour"].unique())
    first_invalid = float(invalid_hours[0]) if invalid_hours else None
    last_valid = float(all_hours[-1]) if not invalid_hours else float(invalid_hours[0] - 1.0)
    row: dict = {
        "orbit": orbit,
        "family": orbit.rsplit("-", 1)[0],
        "total_checkpoint_count": int(len(all_hours)),
        "valid_checkpoint_count": int(valid_cloud["hour"].nunique()),
        "first_invalid_hour": first_invalid,
        "last_valid_hour": last_valid,
    }
    for group_name in ("treatment", "control"):
        group = valid_cloud.loc[valid_cloud["group"].eq(group_name)]
        hours = group["hour"].to_numpy(float)
        for metric in SUMMARY_METRICS:
            values = pd.to_numeric(group[metric], errors="coerce").to_numpy(float)
            finite_mask = np.isfinite(values)
            finite_values = values[finite_mask]
            finite_hours = hours[finite_mask]
            prefix = f"{group_name}_{metric}"
            if len(finite_values):
                peak_index = int(np.argmax(finite_values))
                row[f"{prefix}_endpoint"] = float(finite_values[-1])
                row[f"{prefix}_peak"] = float(finite_values[peak_index])
                row[f"{prefix}_peak_hour"] = float(finite_hours[peak_index])
                row[f"{prefix}_mean"] = float(np.mean(finite_values))
            else:
                row[f"{prefix}_endpoint"] = None
                row[f"{prefix}_peak"] = None
                row[f"{prefix}_peak_hour"] = None
                row[f"{prefix}_mean"] = None
    return row


def analyse_orbit(input_root: Path, orbit: str, expected_hours: int,
                  expected_points: int) -> dict:
    treatment_dir, control_dir = input_root / "treatment" / orbit, input_root / "control" / orbit
    if not treatment_dir.is_dir() or not control_dir.is_dir():
        raise AnalysisError(f"Missing treatment/control directory for {orbit}")
    t_shape = read_shape(find_one(treatment_dir, "_shape.csv"))
    c_shape = read_shape(find_one(control_dir, "_shape.csv"))
    t_history = read_history(find_one(treatment_dir, "_history.csv"), "t_")
    c_history = read_history(find_one(control_dir, "_history.csv"), "c_")
    max_u_difference = verify_initial_design(
        find_one(treatment_dir, "_initial_points.csv"),
        find_one(control_dir, "_initial_points.csv"),
        expected_points,
    )

    hours = np.sort(t_shape["hour"].unique())
    if len(hours) != expected_hours:
        raise AnalysisError(f"{orbit}: expected {expected_hours} hours, found {len(hours)}")
    for frame, label in ((t_shape, "treatment"), (c_shape, "control")):
        counts = frame.groupby("hour")["point_id"].count()
        if not counts.eq(expected_points).all():
            raise AnalysisError(
                f"{orbit} {label}: expected {expected_points} points/hour; "
                f"bad={counts[counts.ne(expected_points)].to_dict()}"
            )

    validity = compute_validity(t_history, c_history)

    # Pair first so the relative-maximum numerator is the maximum difference
    # between corresponding treatment/control points, rather than either
    # cloud's maximum absolute distance from the origin.
    delta = delta_by_hour(pair_shapes(t_shape, c_shape))

    t_cloud = cloud_statistics(t_shape)
    c_cloud = cloud_statistics(c_shape)
    t_cloud, c_cloud = attach_paired_relative_max(t_cloud, c_cloud, delta)
    t_cloud.insert(0, "group", "treatment")
    c_cloud.insert(0, "group", "control")
    cloud = pd.concat([t_cloud, c_cloud], ignore_index=True)
    cloud = cloud.merge(validity, on="hour", how="left", validate="many_to_one")

    delta = delta.merge(validity, on="hour", how="left", validate="one_to_one")

    t_split = first_split(t_history, "t_", range(1, 8))
    c_split = first_split(c_history, "c_", range(1, 7))
    control_dim7_nonzero = dim7_ever_nonzero(c_history)

    return {
        "orbit": orbit,
        "cloud": cloud,
        "delta": delta,
        "summary": summarize_orbit(orbit, cloud, validity),
        "first_split": make_first_split_row(orbit, t_split, c_split, control_dim7_nonzero),
        "max_u_difference": max_u_difference,
    }


def make_first_split_row(orbit: str, t_split: dict, c_split: dict,
                         control_dim7_nonzero: bool) -> dict:
    t_hour, c_hour = t_split["hour"], c_split["hour"]
    both_finite = np.isfinite(t_hour) and np.isfinite(c_hour)
    return {
        "orbit": orbit,
        "treatment_first_split_hour": t_hour,
        "treatment_first_split_interval_start": t_split["interval_start"],
        "treatment_first_split_interval_end": t_split["interval_end"],
        "treatment_first_nonzero_directions": ",".join(map(str, t_split["directions"])),
        "control_first_split_hour": c_hour,
        "control_first_split_interval_start": c_split["interval_start"],
        "control_first_split_interval_end": c_split["interval_end"],
        "control_first_nonzero_directions": ",".join(map(str, c_split["directions"])),
        "control_dim7_ever_nonzero": control_dim7_nonzero,
        "split_hour_difference": float(t_hour - c_hour) if both_finite else math.nan,
    }


def population_delta(deltas: pd.DataFrame) -> pd.DataFrame:
    valid = deltas.loc[deltas["valid"].eq(1)]
    orbit_count = int(deltas["orbit"].nunique())
    records = []
    for hour in sorted(deltas["hour"].unique()):
        group = valid.loc[valid["hour"].eq(hour)]
        row: dict = {
            "hour": float(hour),
            "valid_orbit_count": int(group["orbit"].nunique()),
            "excluded_orbit_count": int(orbit_count - group["orbit"].nunique()),
        }
        for metric in DELTA_METRICS:
            values = pd.to_numeric(group[metric], errors="coerce").dropna().to_numpy(float)
            if len(values):
                row[f"{metric}_min"] = float(np.min(values))
                row[f"{metric}_q25"] = float(np.quantile(values, 0.25))
                row[f"{metric}_median"] = float(np.median(values))
                row[f"{metric}_q75"] = float(np.quantile(values, 0.75))
                row[f"{metric}_max"] = float(np.max(values))
        records.append(row)
    return pd.DataFrame.from_records(records)


# --------------------------------------------------------------------------- #
# Figures
# --------------------------------------------------------------------------- #
def _plot_overlaid(ax, cloud: pd.DataFrame, metric: str, ylabel: str) -> None:
    for group, color in (("treatment", "#d62728"), ("control", "#1f77b4")):
        frame = cloud.loc[(cloud["group"].eq(group)) & (cloud["valid"].eq(1))]
        ax.plot(frame["hour"], frame[metric], lw=1.0, color=color, label=group)
    ax.set(xlabel="Propagation time (h)", ylabel=ylabel)
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)


def make_cloud_moments(cloud: pd.DataFrame, orbit: str, output_dir: Path) -> None:
    panels = (
        ("pos_mean_offset_km", "Position mean offset (km)"),
        ("vel_mean_offset_km_s", "Velocity mean offset (km/s)"),
        ("pos_longest_rms_axis_km", "Position longest RMS axis (km)"),
        ("vel_longest_rms_axis_km_s", "Velocity longest RMS axis (km/s)"),
        ("skew_max_abs", "Max |skewness| (dimensionless)"),
        ("kurt_max_abs", "Max |excess kurtosis| (dimensionless)"),
    )
    fig, axes = plt.subplots(3, 2, figsize=(13, 13), constrained_layout=True)
    for ax, (metric, ylabel) in zip(axes.flat, panels):
        _plot_overlaid(ax, cloud, metric, ylabel)
    fig.suptitle(f"{orbit}: error-cloud moments (treatment vs control)")
    fig.savefig(output_dir / f"{orbit}_cloud_moments.png", dpi=180)
    plt.close(fig)


def make_relative_max_error(cloud: pd.DataFrame, orbit: str, output_dir: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)
    _plot_overlaid(
        axes[0], cloud, "rel_max_pos",
        "Max paired treatment-control position difference / group RMS axis",
    )
    _plot_overlaid(
        axes[1], cloud, "rel_max_vel",
        "Max paired treatment-control velocity difference / group RMS axis",
    )
    fig.suptitle(f"{orbit}: paired SRP effect relative to each cloud's longest RMS axis")
    fig.savefig(output_dir / f"{orbit}_relative_max_error.png", dpi=180)
    plt.close(fig)


def make_first_split_figure(rows: list[dict], output_dir: Path) -> None:
    orbits = [row["orbit"] for row in rows]
    t_hours = [row["treatment_first_split_hour"] for row in rows]
    c_hours = [row["control_first_split_hour"] for row in rows]
    x = np.arange(len(orbits))
    finite = [h for h in t_hours + c_hours if np.isfinite(h)]
    top = max(finite) if finite else 1.0

    fig, ax = plt.subplots(figsize=(12, 6), constrained_layout=True)
    labeled_treatment = labeled_control = False
    for i, (th, ch) in enumerate(zip(t_hours, c_hours)):
        if np.isfinite(th):
            ax.plot(i - 0.18, th, "o", color="#d62728", ms=7,
                    label="treatment" if not labeled_treatment else None)
            labeled_treatment = True
        if np.isfinite(ch):
            ax.plot(i + 0.18, ch, "s", color="#1f77b4", ms=7,
                    label="control" if not labeled_control else None)
            labeled_control = True
        if not np.isfinite(th) and not np.isfinite(ch):
            ax.text(i, top * 1.02, "no split", ha="center", va="bottom",
                    fontsize=8, color="gray")
    ax.set_ylim(-top * 0.06, top * 1.18)
    ax.set_xticks(x, orbits, rotation=45, ha="right")
    ax.set(xlabel="Orbit", ylabel="First split saved hour (h)",
           title="First cumulative split hour: treatment vs control (no-split orbits marked)")
    ax.grid(True, axis="y", alpha=0.25)
    if labeled_treatment or labeled_control:
        ax.legend(fontsize=8)
    fig.savefig(output_dir / "first_split_comparison.png", dpi=180)
    plt.close(fig)


def make_delta_figures(deltas: pd.DataFrame, population: pd.DataFrame, output_dir: Path) -> None:
    valid = deltas.loc[deltas["valid"].eq(1)]
    for metric, ylabel, base in (
        ("delta_r_max_km", "Max paired SRP position increment (km)", "delta_r_max"),
        ("delta_v_max_km_s", "Max paired SRP velocity increment (km/s)", "delta_v_max"),
    ):
        fig, ax = plt.subplots(figsize=(11, 6), constrained_layout=True)
        for orbit, group in valid.groupby("orbit", sort=False):
            ax.plot(group["hour"], group[metric], lw=1.1, label=orbit)
        ax.set(xlabel="Propagation time (h)", ylabel=ylabel,
               title="Validated checkpoints only; each line stops at first invalid checkpoint")
        ax.grid(True, alpha=0.25)
        ax.legend(ncol=3, fontsize=8)
        fig.savefig(output_dir / f"{base}_by_orbit.png", dpi=180)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(10, 5.5), constrained_layout=True)
        x = population["hour"].to_numpy(float)
        ax.plot(x, population[f"{metric}_median"], color="#1f77b4", label="Median of valid orbits")
        ax.fill_between(x, population[f"{metric}_q25"], population[f"{metric}_q75"],
                        color="#1f77b4", alpha=0.22, label="25%-75% range")
        ax.set(xlabel="Propagation time (h)", ylabel=ylabel,
               title="Population statistics use only ADS-validated orbit results")
        ax.grid(True, alpha=0.25)
        ax.legend()
        fig.savefig(output_dir / f"{base}_population.png", dpi=180)
        plt.close(fig)


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    script_dir = Path(__file__).resolve().parent
    parser.add_argument("--input-root", type=Path, default=script_dir / "non_nrho15")
    parser.add_argument("--output-root", type=Path, default=None)
    parser.add_argument("--expected-hours", type=int, default=361)
    parser.add_argument("--expected-points", type=int, default=512)
    parser.add_argument("--no-plots", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    input_root = args.input_root.expanduser().resolve()
    output_root = (args.output_root or (input_root.parent / "analysis_non_nrho15")).expanduser().resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    cloud_dir = output_root / "cloud_statistics"
    cloud_dir.mkdir(parents=True, exist_ok=True)
    figures_dir = output_root / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    clouds, deltas, summaries, split_rows, failures = [], [], [], [], []
    for index, orbit in enumerate(ORBIT_NAMES, 1):
        print(f"[{index:02d}/15] {orbit}", flush=True)
        try:
            result = analyse_orbit(input_root, orbit, args.expected_hours, args.expected_points)
            cloud = result["cloud"]
            cloud.insert(0, "family", orbit.rsplit("-", 1)[0])
            cloud.insert(0, "orbit", orbit)
            clouds.append(cloud)
            delta = result["delta"].copy()
            delta.insert(0, "orbit", orbit)
            deltas.append(delta)
            summaries.append(result["summary"])
            split_rows.append(result["first_split"])
            cloud[["orbit", "family", "group", "hour", "valid", *CLOUD_METRIC_COLUMNS]].to_csv(
                cloud_dir / f"{orbit}_cloud_statistics.csv", index=False, float_format="%.15g"
            )
            if not args.no_plots:
                make_cloud_moments(cloud, orbit, figures_dir)
                make_relative_max_error(cloud, orbit, figures_dir)
        except Exception as exc:
            failures.append({"orbit": orbit, "error": f"{type(exc).__name__}: {exc}"})
            print(f"  ERROR: {exc}", file=sys.stderr, flush=True)

    figure_warnings: list[str] = []
    if summaries:
        pd.DataFrame.from_records(summaries).to_csv(
            output_root / "cloud_orbit_summary.csv", index=False, float_format="%.15g"
        )
        split_frame = pd.DataFrame.from_records(split_rows)
        split_frame.to_csv(output_root / "first_split_comparison.csv", index=False, float_format="%.15g")
        all_deltas = pd.concat(deltas, ignore_index=True)
        population = population_delta(all_deltas)
        population.to_csv(output_root / "paired_population_history.csv", index=False, float_format="%.15g")
        if not args.no_plots:
            make_first_split_figure(split_rows, figures_dir)
            make_delta_figures(all_deltas, population, figures_dir)

    control_dim7_violations = [r["orbit"] for r in split_rows if r["control_dim7_ever_nonzero"]]
    manifest = {
        "analysis_format_version": 1,
        "status": "complete" if not failures and len(summaries) == len(ORBIT_NAMES) else "incomplete",
        "input_root": str(input_root),
        "output_root": str(output_root),
        "included_orbits": [r["orbit"] for r in summaries],
        "explicitly_excluded_orbits": ["L2NRHO-1", "L2NRHO-2"],
        "failed_orbits": failures,
        "expected_checkpoint_count": args.expected_hours,
        "expected_shape_point_count": args.expected_points,
        "units": {"position": "km", "velocity": "km/s",
                  "skewness_kurtosis_relative_max_error": "dimensionless", "time": "h"},
        "exclusion_policy": {
            "rule": "exclude first invalid checkpoint and all later checkpoints when either group has validation_failed or precision_unconfirmed",
            "figures_and_population_statistics": "validated checkpoints only",
        },
        "relative_max_definition": {
            "shared_numerator_position": "max over matching point_id of norm(treatment position error - control position error), km",
            "shared_numerator_velocity": "max over matching point_id of norm(treatment velocity error - control velocity error), km/s",
            "treatment_denominator": "treatment cloud longest covariance RMS axis",
            "control_denominator": "control cloud longest covariance RMS axis",
        },
        "near_zero_variance_threshold": NEAR_ZERO_VARIANCE,
        "control_dim7_nonzero_orbits": control_dim7_violations,
        "first_split_limitations": [
            "Treatment uses a legacy schema without per-event split criteria.",
            "Split hours are derived from cumulative counters at hourly save points.",
        ],
        "warnings": figure_warnings,
    }
    with (output_root / "cloud_analysis_manifest.json").open("w", encoding="utf-8", newline="\n") as stream:
        json.dump(manifest, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write("\n")
    print(f"Completed {len(summaries)}/15 orbits; results: {output_root}")
    return 0 if not failures and len(summaries) == len(ORBIT_NAMES) else 1


if __name__ == "__main__":
    raise SystemExit(main())
