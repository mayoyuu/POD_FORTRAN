#!/usr/bin/env python3
"""Compare paired SRP-uncertain and eta=0 ADS campaigns without a delta point cloud."""
from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
import statistics
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Iterator

PAIR_SCHEMA_VERSION = 1
COMPONENTS = (
    ("dx_km", "x", 0),
    ("dy_km", "y", 1),
    ("dz_km", "z", 2),
    ("dvx_kms", "vx", 3),
    ("dvy_kms", "vy", 4),
    ("dvz_kms", "vz", 5),
)
NOMINAL_COLUMNS = (
    "nominal_x_km",
    "nominal_y_km",
    "nominal_z_km",
    "nominal_vx_kms",
    "nominal_vy_kms",
    "nominal_vz_kms",
)
INITIAL_STATE_COLUMNS = (
    "x0_km",
    "y0_km",
    "z0_km",
    "vx0_kms",
    "vy0_kms",
    "vz0_kms",
)
COMMON_FLOAT_PARAMETERS = (
    "duration_hours",
    "save_hours",
    "position_tolerance_km",
    "velocity_tolerance_kms",
)
COMMON_INTEGER_PARAMETERS = (
    "da_order",
    "max_split_depth",
    "shape_point_count",
    "validation_point_count",
)
EVENT_FIELDS = (
    "first_split_hour",
    "first_split_interval_start_hour",
    "first_split_component",
    "first_split_direction",
    "first_split_parent_depth",
    "first_split_bfs_iteration",
    "first_split_truncation_estimates",
    "first_split_component_tolerances",
    "first_split_positive_excess",
    "first_split_direction_scores",
    "first_split_parent_unit_center",
    "first_split_parent_unit_width",
    "first_unconfirmed_hour",
    "first_failed_hour",
    "previous_passing_hour",
    "worst_point_id",
    "worst_point_position_error_km",
    "worst_point_velocity_error_kms",
)
ETA_EVENT_FIELDS = (
    "first_eta_split_hour",
    "first_eta_split_interval_start_hour",
    "first_eta_split_component",
    "first_eta_split_parent_depth",
    "first_eta_split_bfs_iteration",
    "first_eta_split_truncation_estimates",
    "first_eta_split_component_tolerances",
    "first_eta_split_positive_excess",
    "first_eta_split_direction_scores",
    "first_eta_split_parent_unit_center",
    "first_eta_split_parent_unit_width",
)


class PairingError(RuntimeError):
    pass


def _float(row: dict, key: str) -> float:
    try:
        value = float(row[key])
    except (KeyError, TypeError, ValueError) as exc:
        raise PairingError(f"missing or invalid numeric field {key}") from exc
    if not math.isfinite(value):
        raise PairingError(f"non-finite numeric field {key}")
    return value


def _int(row: dict, key: str) -> int:
    try:
        return int(row[key])
    except (KeyError, TypeError, ValueError) as exc:
        raise PairingError(f"missing or invalid integer field {key}") from exc


def _close(a: float, b: float, tolerance: float = 1e-11) -> bool:
    return math.isclose(a, b, rel_tol=tolerance, abs_tol=tolerance)


def _load_json(path: Path) -> dict:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise PairingError(f"cannot read JSON {path}: {exc}") from exc


def _load_csv(path: Path) -> list[dict]:
    try:
        with path.open(newline="", encoding="utf-8") as handle:
            return list(csv.DictReader(handle))
    except OSError as exc:
        raise PairingError(f"cannot read CSV {path}: {exc}") from exc


def _write_json(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    temp.replace(path)


def _write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise PairingError(f"cannot write empty CSV {path}")
    temp = path.with_suffix(path.suffix + ".tmp")
    with temp.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        for row in rows:
            writer.writerow({key: "" if value is None else value for key, value in row.items()})
    temp.replace(path)


def design_quantile(values: Iterable[float], probability: float) -> float:
    ordered = sorted(float(value) for value in values)
    if not ordered:
        raise PairingError("quantile requires at least one value")
    if probability < 0.0 or probability > 1.0:
        raise PairingError("quantile probability must lie in [0,1]")
    position = (len(ordered) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    fraction = position - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def _validate_reports(treatment: dict, control: dict) -> None:
    if int(treatment.get("output_schema_version", 0)) < 2:
        raise PairingError("treatment output_schema_version is incompatible")
    if int(control.get("output_schema_version", 0)) < 2:
        raise PairingError("control output_schema_version is incompatible")
    if treatment.get("n_variables") != 7:
        raise PairingError("treatment n_variables must be 7")
    if control.get("n_variables") != 6:
        raise PairingError("control n_variables must be 6")
    if treatment.get("nominal_srp_enabled") is not True:
        raise PairingError("treatment nominal SRP is not enabled")
    if control.get("nominal_srp_enabled") is not True:
        raise PairingError("control nominal SRP is not enabled")
    if treatment.get("srp_uncertainty_enabled") is not True:
        raise PairingError("treatment SRP uncertainty is not enabled")
    if control.get("srp_uncertainty_enabled") is not False:
        raise PairingError("control SRP uncertainty must be disabled")
    if float(treatment.get("eta_sigma", 0.0)) <= 0.0:
        raise PairingError("treatment eta_sigma must be positive")
    if not _close(float(control.get("eta_sigma", math.nan)), 0.0):
        raise PairingError("control eta_sigma must be zero")

    for key in COMMON_FLOAT_PARAMETERS:
        if key not in treatment or key not in control:
            raise PairingError(f"missing paired parameter {key}")
        if not _close(float(treatment[key]), float(control[key])):
            raise PairingError(f"paired parameter mismatch: {key}")
    for key in COMMON_INTEGER_PARAMETERS:
        if treatment.get(key) != control.get(key):
            raise PairingError(f"paired parameter mismatch: {key}")
    for key in ("opm_sha256", "config_sha256"):
        if not treatment.get(key) or treatment.get(key) != control.get(key):
            raise PairingError(f"paired hash mismatch: {key}")
    for key in ("executable_sha256", "source_revision"):
        if treatment.get(key) is not None and control.get(key) is not None:
            if treatment.get(key) != control.get(key):
                raise PairingError(f"paired provenance mismatch: {key}")


def _validate_initial_points(
    treatment_prefix: Path, control_prefix: Path, point_count: int
) -> dict:
    treatment = _load_csv(Path(str(treatment_prefix) + "_initial_points.csv"))
    control = _load_csv(Path(str(control_prefix) + "_initial_points.csv"))
    if len(treatment) != point_count or len(control) != point_count:
        raise PairingError("initial point count mismatch")
    max_initial_state_difference = 0.0
    for index, (left, right) in enumerate(zip(treatment, control), 1):
        if _int(left, "point_id") != index or _int(right, "point_id") != index:
            raise PairingError(f"initial point_id mismatch at row {index}")
        for key in tuple(f"u{k}" for k in range(1, 8)):
            if not _close(_float(left, key), _float(right, key)):
                raise PairingError(f"paired point-set mismatch: point_id={index} {key}")
        difference = math.sqrt(
            sum(
                (_float(left, key) - _float(right, key)) ** 2
                for key in INITIAL_STATE_COLUMNS[:3]
            )
        )
        velocity_difference = math.sqrt(
            sum(
                (_float(left, key) - _float(right, key)) ** 2
                for key in INITIAL_STATE_COLUMNS[3:]
            )
        )
        max_initial_state_difference = max(
            max_initial_state_difference, difference, velocity_difference
        )
        if difference > 1e-10 or velocity_difference > 1e-12:
            raise PairingError(f"initial physical state mismatch at point_id={index}")
        if abs(_float(right, "eta")) > 1e-15:
            raise PairingError(f"control eta is nonzero at point_id={index}")
        if _int(left, "validation") != _int(right, "validation"):
            raise PairingError(f"validation point-set mismatch at point_id={index}")

    if point_count >= 128:
        for point_id in range(1, 65):
            minus = control[point_id - 1]
            plus = control[point_id + 63]
            for key in tuple(f"u{k}" for k in range(1, 7)) + INITIAL_STATE_COLUMNS:
                if not _close(_float(minus, key), _float(plus, key)):
                    raise PairingError(
                        f"control u7 projection mismatch for point pair "
                        f"{point_id},{point_id + 64}"
                    )
    return {
        "point_count": point_count,
        "max_initial_state_difference": max_initial_state_difference,
    }


def _load_histories(
    treatment_prefix: Path,
    control_prefix: Path,
    report: dict,
) -> tuple[list[dict], list[dict]]:
    treatment = _load_csv(Path(str(treatment_prefix) + "_history.csv"))
    control = _load_csv(Path(str(control_prefix) + "_history.csv"))
    expected = math.ceil(
        float(report["duration_hours"]) / float(report["save_hours"]) - 1e-14
    ) + 1
    if len(treatment) != expected or len(control) != expected:
        raise PairingError(
            f"checkpoint count mismatch: expected {expected}, "
            f"treatment {len(treatment)}, control {len(control)}"
        )
    for index, (left, right) in enumerate(zip(treatment, control)):
        hour_left, hour_right = _float(left, "hour"), _float(right, "hour")
        epoch_left, epoch_right = _float(left, "epoch_et"), _float(right, "epoch_et")
        if not _close(hour_left, hour_right):
            raise PairingError(f"hour mismatch at checkpoint {index}")
        if not _close(epoch_left, epoch_right):
            raise PairingError(f"epoch mismatch at checkpoint {index}")
    return treatment, control


def _shape_groups(path: Path, point_count: int) -> Iterator[tuple[float, list[list[float]]]]:
    try:
        handle = path.open(newline="", encoding="utf-8")
    except OSError as exc:
        raise PairingError(f"cannot read shape CSV {path}: {exc}") from exc
    with handle:
        reader = csv.DictReader(handle)
        while True:
            first = next(reader, None)
            if first is None:
                return
            hour = _float(first, "hour")
            rows = [first]
            for _ in range(point_count - 1):
                row = next(reader, None)
                if row is None:
                    raise PairingError(f"incomplete shape checkpoint at hour {hour}")
                rows.append(row)
            values: list[list[float]] = []
            for expected_id, row in enumerate(rows, 1):
                if _int(row, "point_id") != expected_id:
                    raise PairingError(
                        f"shape point_id mismatch at hour {hour}: "
                        f"expected {expected_id}, got {row.get('point_id')}"
                    )
                if not _close(_float(row, "hour"), hour):
                    raise PairingError(f"shape hour changed within checkpoint {hour}")
                values.append([_float(row, name) for name, _, _ in COMPONENTS])
            yield hour, values


def _norm(values: list[float]) -> float:
    return math.sqrt(sum(value * value for value in values))


def _ratio(numerator: float, denominator: float) -> float | None:
    return numerator / denominator if abs(denominator) > 1e-30 else None


def _paired_hour(
    treatment_history: dict,
    control_history: dict,
    treatment_shape: list[list[float]],
    control_shape: list[list[float]],
) -> tuple[dict, float, float, float]:
    point_count = len(treatment_shape)
    if len(control_shape) != point_count:
        raise PairingError("shape point count differs between groups")
    hour = _float(treatment_history, "hour")
    row: dict = {"hour": hour, "epoch_et": _float(treatment_history, "epoch_et")}

    for component, _, _ in COMPONENTS:
        t_lower = _float(treatment_history, f"lower_{component}")
        t_upper = _float(treatment_history, f"upper_{component}")
        c_lower = _float(control_history, f"lower_{component}")
        c_upper = _float(control_history, f"upper_{component}")
        t_width, c_width = t_upper - t_lower, c_upper - c_lower
        row.update(
            {
                f"treatment_lower_{component}": t_lower,
                f"treatment_upper_{component}": t_upper,
                f"treatment_width_{component}": t_width,
                f"control_lower_{component}": c_lower,
                f"control_upper_{component}": c_upper,
                f"control_width_{component}": c_width,
                f"width_abs_increment_{component}": abs(t_width - c_width),
                f"width_amplification_{component}": _ratio(t_width, c_width),
            }
        )

    for group_name, history in (
        ("treatment", treatment_history),
        ("control", control_history),
    ):
        for key in ("pos_max_km", "vel_max_kms"):
            row[f"{group_name}_{key}"] = _float(history, key)
        for kind, unit in (("pos", "km"), ("vel", "kms")):
            axes = [_float(history, f"{kind}_axis{k}_{unit}") for k in range(1, 4)]
            for index, value in enumerate(axes, 1):
                row[f"{group_name}_{kind}_axis{index}_{unit}"] = value
            volume_unit = "km3" if kind == "pos" else "kms3"
            row[f"{group_name}_{kind}_volume_proxy_{volume_unit}"] = math.prod(axes)

    for key in ("pos_max_km", "vel_max_kms"):
        treatment_value = row[f"treatment_{key}"]
        control_value = row[f"control_{key}"]
        row[f"{key}_abs_increment"] = abs(treatment_value - control_value)
        row[f"{key}_amplification"] = _ratio(treatment_value, control_value)
    for kind, unit, volume_unit in (
        ("pos", "km", "km3"),
        ("vel", "kms", "kms3"),
    ):
        for axis in range(1, 4):
            t_key = f"treatment_{kind}_axis{axis}_{unit}"
            c_key = f"control_{kind}_axis{axis}_{unit}"
            row[f"{kind}_axis{axis}_{unit}_abs_increment"] = abs(row[t_key] - row[c_key])
            row[f"{kind}_axis{axis}_{unit}_amplification"] = _ratio(row[t_key], row[c_key])
        t_key = f"treatment_{kind}_volume_proxy_{volume_unit}"
        c_key = f"control_{kind}_volume_proxy_{volume_unit}"
        row[f"{kind}_volume_proxy_{volume_unit}_abs_increment"] = abs(row[t_key] - row[c_key])
        row[f"{kind}_volume_proxy_{volume_unit}_amplification"] = _ratio(row[t_key], row[c_key])

    nominal_t = [_float(treatment_history, key) for key in NOMINAL_COLUMNS]
    nominal_c = [_float(control_history, key) for key in NOMINAL_COLUMNS]
    nominal_delta = [left - right for left, right in zip(nominal_t, nominal_c)]
    nominal_position_difference = _norm(nominal_delta[:3])
    nominal_velocity_difference = _norm(nominal_delta[3:])

    paired_delta: list[list[float]] = []
    for left, right in zip(treatment_shape, control_shape):
        paired_delta.append(
            [
                nominal_delta[index] + left[index] - right[index]
                for index in range(6)
            ]
        )

    for component, _, index in COMPONENTS:
        treatment_centroid = statistics.fmean(point[index] for point in treatment_shape)
        control_centroid = statistics.fmean(point[index] for point in control_shape)
        values = [point[index] for point in paired_delta]
        row[f"treatment_centroid_{component}"] = treatment_centroid
        row[f"control_centroid_{component}"] = control_centroid
        row[f"centroid_delta_{component}"] = statistics.fmean(values)
        row[f"srp_delta_{component}_min"] = min(values)
        row[f"srp_delta_{component}_max"] = max(values)

    delta_r = [_norm(point[:3]) for point in paired_delta]
    delta_v = [_norm(point[3:]) for point in paired_delta]
    worst_r = max(range(point_count), key=delta_r.__getitem__)
    worst_v = max(range(point_count), key=delta_v.__getitem__)
    row.update(
        {
            "delta_r_max_km": delta_r[worst_r],
            "delta_r_rms_km": math.sqrt(statistics.fmean(value * value for value in delta_r)),
            "delta_r_design_p95_km": design_quantile(delta_r, 0.95),
            "delta_r_worst_point_id": worst_r + 1,
            "delta_v_max_kms": delta_v[worst_v],
            "delta_v_rms_kms": math.sqrt(statistics.fmean(value * value for value in delta_v)),
            "delta_v_design_p95_kms": design_quantile(delta_v, 0.95),
            "delta_v_worst_point_id": worst_v + 1,
        }
    )

    if point_count < 128:
        raise PairingError("at least 128 corner points are required for odd/even metrics")
    odd_values, even_values = [], []
    for index in range(64):
        minus = paired_delta[index]
        plus = paired_delta[index + 64]
        odd_values.append([(p - m) / 2.0 for p, m in zip(plus, minus)])
        even_values.append([(p + m) / 2.0 for p, m in zip(plus, minus)])
    for component, _, index in COMPONENTS:
        odd_component = [point[index] for point in odd_values]
        even_component = [point[index] for point in even_values]
        row[f"odd_{component}_min"] = min(odd_component)
        row[f"odd_{component}_max"] = max(abs(value) for value in odd_component)
        row[f"even_{component}_min"] = min(even_component)
        row[f"even_{component}_max"] = max(abs(value) for value in even_component)
    row["odd_delta_r_max_km"] = max(_norm(point[:3]) for point in odd_values)
    row["even_delta_r_max_km"] = max(_norm(point[:3]) for point in even_values)
    row["odd_delta_v_max_kms"] = max(_norm(point[3:]) for point in odd_values)
    row["even_delta_v_max_kms"] = max(_norm(point[3:]) for point in even_values)

    duplicate_difference = 0.0
    for index in range(64):
        duplicate_difference = max(
            duplicate_difference,
            max(
                abs(a - b)
                for a, b in zip(control_shape[index], control_shape[index + 64])
            ),
        )

    for group_name, history in (
        ("treatment", treatment_history),
        ("control", control_history),
    ):
        for key in (
            "patches",
            "depth_limited",
            "bfs_iterations",
            "max_queue_size",
            "step_elapsed_seconds",
            "max_validation_pos_km",
            "max_validation_vel_kms",
        ):
            row[f"{group_name}_{key}"] = _float(history, key)
        for direction in range(1, 8):
            row[f"{group_name}_split_dim{direction}"] = _float(
                history, f"split_dim{direction}"
            )
    if row["control_split_dim7"] != 0.0:
        raise PairingError(f"control split_dim7 is nonzero at hour {hour}")
    return (
        row,
        nominal_position_difference,
        nominal_velocity_difference,
        duplicate_difference,
    )


def _aggregate_metric(rows: list[dict], key: str) -> dict | None:
    samples = [
        (float(row["hour"]), row[key])
        for row in rows
        if row.get(key) is not None and isinstance(row.get(key), (int, float))
    ]
    if not samples:
        return None
    peak_hour, peak_value = max(samples, key=lambda item: abs(item[1]))
    integral = 0.0
    for (hour0, value0), (hour1, value1) in zip(samples, samples[1:]):
        integral += 0.5 * (value0 + value1) * (hour1 - hour0)
    duration = samples[-1][0] - samples[0][0]
    return {
        "endpoint": samples[-1][1],
        "peak": peak_value,
        "peak_hour": peak_hour,
        "time_average": integral / duration if duration > 0.0 else samples[-1][1],
        "trapezoid_integral": integral,
    }


def _is_effect_metric(key: str) -> bool:
    return (
        "abs_increment" in key
        or "amplification" in key
        or key.startswith("centroid_delta_")
        or key.startswith("srp_delta_")
        or key.startswith("delta_")
        or key.startswith("odd_")
        or key.startswith("even_")
    ) and not key.endswith("point_id")


def _event_data(report: dict) -> dict:
    return {key: report.get(key) for key in EVENT_FIELDS}


def _eta_event_data(report: dict) -> dict:
    return {key: report.get(key) for key in ETA_EVENT_FIELDS}


def _event_difference(treatment: dict, control: dict, key: str) -> float | None:
    left, right = treatment.get(key), control.get(key)
    if left is None or right is None:
        return None
    return float(left) - float(right)


def compare_orbit(
    treatment_prefix: Path | str,
    control_prefix: Path | str,
    output_prefix: Path | str,
) -> dict:
    treatment_prefix = Path(treatment_prefix)
    control_prefix = Path(control_prefix)
    output_prefix = Path(output_prefix)
    treatment_report_path = Path(str(treatment_prefix) + "_report.json")
    control_report_path = Path(str(control_prefix) + "_report.json")
    treatment_report = _load_json(treatment_report_path)
    control_report = _load_json(control_report_path)
    _validate_reports(treatment_report, control_report)

    point_count = int(treatment_report["shape_point_count"])
    initial_integrity = _validate_initial_points(
        treatment_prefix, control_prefix, point_count
    )
    treatment_history, control_history = _load_histories(
        treatment_prefix, control_prefix, treatment_report
    )

    shape_t = _shape_groups(Path(str(treatment_prefix) + "_shape.csv"), point_count)
    shape_c = _shape_groups(Path(str(control_prefix) + "_shape.csv"), point_count)
    paired_rows: list[dict] = []
    nominal_position_max = 0.0
    nominal_velocity_max = 0.0
    duplicate_max = 0.0
    sentinel = object()
    for checkpoint, groups in enumerate(
        itertools.zip_longest(shape_t, shape_c, fillvalue=sentinel)
    ):
        group_t, group_c = groups
        if group_t is sentinel or group_c is sentinel:
            raise PairingError("shape checkpoint count mismatch")
        hour_t, values_t = group_t
        hour_c, values_c = group_c
        if checkpoint >= len(treatment_history):
            raise PairingError("shape file has more checkpoints than history")
        expected_hour = _float(treatment_history[checkpoint], "hour")
        if not _close(hour_t, hour_c) or not _close(hour_t, expected_hour):
            raise PairingError(f"shape/history hour mismatch at checkpoint {checkpoint}")
        row, nominal_position, nominal_velocity, duplicate_difference = _paired_hour(
            treatment_history[checkpoint],
            control_history[checkpoint],
            values_t,
            values_c,
        )
        paired_rows.append(row)
        nominal_position_max = max(nominal_position_max, nominal_position)
        nominal_velocity_max = max(nominal_velocity_max, nominal_velocity)
        duplicate_max = max(duplicate_max, duplicate_difference)

    if len(paired_rows) != len(treatment_history):
        raise PairingError(
            f"shape checkpoint count mismatch: expected {len(treatment_history)}, "
            f"got {len(paired_rows)}"
        )
    _write_csv(Path(str(output_prefix) + "_paired_history.csv"), paired_rows)

    metrics = {}
    for key in paired_rows[0]:
        if _is_effect_metric(key):
            aggregate = _aggregate_metric(paired_rows, key)
            if aggregate is not None:
                metrics[key] = aggregate

    event_differences = {
        key + "_difference_hours": _event_difference(
            treatment_report, control_report, key
        )
        for key in (
            "first_split_hour",
            "first_split_interval_start_hour",
            "first_unconfirmed_hour",
            "first_failed_hour",
            "previous_passing_hour",
        )
    }
    report = {
        "output_schema_version": PAIR_SCHEMA_VERSION,
        "method": "streamed fixed-design paired SRP ADS comparison",
        "sampling_interpretation": (
            "principal axes, RMS, and 95th percentile are fixed-design sample metrics, "
            "not probability confidence intervals"
        ),
        "treatment_prefix": str(treatment_prefix.resolve()),
        "control_prefix": str(control_prefix.resolve()),
        "paired_history": str(
            Path(str(output_prefix) + "_paired_history.csv").resolve()
        ),
        "checkpoint_count": len(paired_rows),
        "point_count_per_checkpoint": point_count,
        "parameters": {
            key: treatment_report.get(key)
            for key in COMMON_FLOAT_PARAMETERS + COMMON_INTEGER_PARAMETERS
        }
        | {
            "treatment_eta_sigma": treatment_report["eta_sigma"],
            "control_eta_sigma": control_report["eta_sigma"],
        },
        "sources": {
            "treatment": {
                "report": str(treatment_report_path.resolve()),
                "source_opm": treatment_report.get("source_opm"),
                "source_config": treatment_report.get("source_config"),
                "opm_sha256": treatment_report.get("opm_sha256"),
                "config_sha256": treatment_report.get("config_sha256"),
                "executable_sha256": treatment_report.get("executable_sha256"),
                "source_revision": treatment_report.get("source_revision"),
            },
            "control": {
                "report": str(control_report_path.resolve()),
                "source_opm": control_report.get("source_opm"),
                "source_config": control_report.get("source_config"),
                "opm_sha256": control_report.get("opm_sha256"),
                "config_sha256": control_report.get("config_sha256"),
                "executable_sha256": control_report.get("executable_sha256"),
                "source_revision": control_report.get("source_revision"),
            },
        },
        "pairing_integrity": {
            "complete": True,
            **initial_integrity,
            "stable_point_ids": True,
            "matching_hours_and_epochs": True,
            "matching_parameters_and_hashes": True,
            "control_u7_projection_duplicate_max_abs_state_difference": duplicate_max,
        },
        "nominal_center_max_position_difference_km": nominal_position_max,
        "nominal_center_max_velocity_difference_kms": nominal_velocity_max,
        "events": {
            "treatment": _event_data(treatment_report),
            "control": _event_data(control_report),
            "treatment_first_eta_split": _eta_event_data(treatment_report),
            "time_differences": event_differences,
        },
        "metrics": metrics,
    }
    _write_json(Path(str(output_prefix) + "_paired_report.json"), report)
    return report


def _population_rows(history_paths: list[Path]) -> list[dict]:
    by_hour: dict[float, dict[str, list[float]]] = defaultdict(
        lambda: defaultdict(list)
    )
    metric_names = (
        "delta_r_max_km",
        "delta_r_rms_km",
        "delta_r_design_p95_km",
        "delta_v_max_kms",
        "delta_v_rms_kms",
        "delta_v_design_p95_kms",
        "pos_max_km_abs_increment",
        "vel_max_kms_abs_increment",
        "odd_delta_r_max_km",
        "even_delta_r_max_km",
        "odd_delta_v_max_kms",
        "even_delta_v_max_kms",
    )
    for path in history_paths:
        for row in _load_csv(path):
            hour = _float(row, "hour")
            for key in metric_names:
                by_hour[hour][key].append(_float(row, key))
    rows = []
    for hour in sorted(by_hour):
        row: dict = {"hour": hour, "orbit_count": len(next(iter(by_hour[hour].values())))}
        for key, values in by_hour[hour].items():
            row[f"{key}_min"] = min(values)
            row[f"{key}_q25"] = design_quantile(values, 0.25)
            row[f"{key}_median"] = design_quantile(values, 0.5)
            row[f"{key}_q75"] = design_quantile(values, 0.75)
            row[f"{key}_max"] = max(values)
        rows.append(row)
    return rows


def compare_campaign(
    treatment_root: Path,
    control_root: Path,
    output_root: Path,
) -> dict:
    treatment_root = treatment_root.resolve()
    control_root = control_root.resolve()
    output_root = output_root.resolve()
    reports = sorted(treatment_root.rglob("*_report.json"))
    orbit_rows = []
    histories = []
    errors = []
    for report_path in reports:
        treatment_prefix = Path(str(report_path)[: -len("_report.json")])
        relative = treatment_prefix.relative_to(treatment_root)
        control_prefix = control_root / relative
        output_prefix = output_root / relative
        try:
            report = compare_orbit(treatment_prefix, control_prefix, output_prefix)
            key_metrics = report["metrics"]
            orbit_rows.append(
                {
                    "relative_prefix": str(relative),
                    "status": "complete",
                    "checkpoint_count": report["checkpoint_count"],
                    "first_treatment_split_hour": report["events"]["treatment"].get(
                        "first_split_hour"
                    ),
                    "first_control_split_hour": report["events"]["control"].get(
                        "first_split_hour"
                    ),
                    "first_eta_split_hour": report["events"][
                        "treatment_first_eta_split"
                    ].get("first_eta_split_hour"),
                    "first_treatment_failure_hour": report["events"]["treatment"].get(
                        "first_failed_hour"
                    ),
                    "first_control_failure_hour": report["events"]["control"].get(
                        "first_failed_hour"
                    ),
                    "nominal_center_max_position_difference_km": report[
                        "nominal_center_max_position_difference_km"
                    ],
                    "nominal_center_max_velocity_difference_kms": report[
                        "nominal_center_max_velocity_difference_kms"
                    ],
                    "delta_r_max_endpoint_km": key_metrics["delta_r_max_km"]["endpoint"],
                    "delta_r_max_peak_km": key_metrics["delta_r_max_km"]["peak"],
                    "delta_r_max_peak_hour": key_metrics["delta_r_max_km"]["peak_hour"],
                    "delta_v_max_endpoint_kms": key_metrics["delta_v_max_kms"]["endpoint"],
                    "delta_v_max_peak_kms": key_metrics["delta_v_max_kms"]["peak"],
                    "delta_v_max_peak_hour": key_metrics["delta_v_max_kms"]["peak_hour"],
                    "error": "",
                }
            )
            histories.append(Path(str(output_prefix) + "_paired_history.csv"))
        except (PairingError, OSError, ValueError, KeyError) as exc:
            errors.append({"relative_prefix": str(relative), "error": str(exc)})
            orbit_rows.append(
                {
                    "relative_prefix": str(relative),
                    "status": "rejected",
                    "checkpoint_count": "",
                    "first_treatment_split_hour": "",
                    "first_control_split_hour": "",
                    "first_eta_split_hour": "",
                    "first_treatment_failure_hour": "",
                    "first_control_failure_hour": "",
                    "nominal_center_max_position_difference_km": "",
                    "nominal_center_max_velocity_difference_kms": "",
                    "delta_r_max_endpoint_km": "",
                    "delta_r_max_peak_km": "",
                    "delta_r_max_peak_hour": "",
                    "delta_v_max_endpoint_kms": "",
                    "delta_v_max_peak_kms": "",
                    "delta_v_max_peak_hour": "",
                    "error": str(exc),
                }
            )

    if not reports:
        raise PairingError(f"no treatment reports found below {treatment_root}")
    _write_csv(output_root / "paired_orbit_summary.csv", orbit_rows)
    if histories:
        _write_csv(
            output_root / "paired_population_history.csv",
            _population_rows(histories),
        )
    summary = {
        "output_schema_version": PAIR_SCHEMA_VERSION,
        "status": "complete" if not errors else "completed_with_rejections",
        "treatment_root": str(treatment_root),
        "control_root": str(control_root),
        "output_root": str(output_root),
        "discovered": len(reports),
        "completed": len(histories),
        "rejected": len(errors),
        "errors": errors,
        "orbit_summary": str(output_root / "paired_orbit_summary.csv"),
        "population_history": (
            str(output_root / "paired_population_history.csv") if histories else None
        ),
    }
    _write_json(output_root / "paired_summary.json", summary)
    return summary


def make_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--treatment-root",
        type=Path,
        default=root / "SRP" / "SRP_UCERTAINTY_PROP_WITHINIT",
    )
    parser.add_argument(
        "--control-root",
        type=Path,
        default=root / "SRP" / "SRP_UCERTAINTY_PROP_WITHINIT_ETA0",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=root / "SRP" / "SRP_UCERTAINTY_EFFECT_WITHINIT",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = make_parser().parse_args(argv)
    try:
        summary = compare_campaign(
            args.treatment_root, args.control_root, args.output_root
        )
    except PairingError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    print(f"Paired summary: {args.output_root.resolve() / 'paired_summary.json'}")
    return 0 if summary["rejected"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
