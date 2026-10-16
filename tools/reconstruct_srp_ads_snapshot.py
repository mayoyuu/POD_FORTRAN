#!/usr/bin/env python3
"""Reconstruct dense point clouds from sparse SRP ADS patch snapshots."""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import math
import os
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import NamedTuple

FORMAT_VERSION = 1
POSITION_VERIFY_TOL_KM = 1.0e-8
VELOCITY_VERIFY_TOL_KMS = 1.0e-11
STATE_NAMES = ("x_km", "y_km", "z_km", "vx_kms", "vy_kms", "vz_kms")
ERROR_NAMES = ("dx_km", "dy_km", "dz_km", "dvx_kms", "dvy_kms", "dvz_kms")
PRIMES = (2, 3, 5, 7, 11, 13, 17, 19)


class SnapshotError(RuntimeError):
    """A snapshot is incomplete or internally inconsistent."""


class Sample(NamedTuple):
    kind: str
    coordinates: tuple[float, ...]


class Patch:
    def __init__(self, patch_id, time_seconds, center, width, history, terms):
        self.patch_id = patch_id
        self.time_seconds = time_seconds
        self.center = tuple(center)
        self.width = tuple(width)
        self.split_history = tuple(history)
        self.terms = terms

    def contains(self, point, n_variables):
        return all(
            abs(point[i] - self.center[i]) <= 0.5 * self.width[i] + 1.0e-13
            for i in range(n_variables)
        )

    def global_to_local(self, point, n_variables):
        local = [0.0] * len(self.center)
        for i in range(n_variables):
            if self.width[i] <= 0.0:
                raise SnapshotError(
                    f"patch {self.patch_id} has non-positive width{i + 1}"
                )
            local[i] = 2.0 * (point[i] - self.center[i]) / self.width[i]
        return tuple(local)

    def local_to_global(self, local, n_variables):
        point = [0.0] * len(self.center)
        for i in range(n_variables):
            point[i] = self.center[i] + 0.5 * self.width[i] * local[i]
        return tuple(point)

    def evaluate_local(self, local):
        state = []
        for component in range(1, 7):
            value = 0.0
            for powers, coefficient in self.terms[component]:
                term = coefficient
                for i, power in enumerate(powers):
                    if power:
                        term *= local[i] ** power
                value += term
            state.append(value)
        return tuple(state)


class Snapshot:
    def __init__(self, prefix, n_variables, da_order, domain_sigma,
                 eta_sigma, time_seconds, patches, source_files, schema_version=1,
                 parameterization=None, mode_hashes=None):
        self.prefix = prefix
        self.n_variables = n_variables
        self.da_order = da_order
        self.domain_sigma = domain_sigma
        self.eta_sigma = eta_sigma
        self.time_seconds = time_seconds
        self.patches = tuple(sorted(patches, key=lambda patch: patch.patch_id))
        self.source_files = source_files
        self.schema_version=schema_version
        self.storage_variables=len(self.patches[0].center)
        self.parameterization=parameterization or ('SRP_ETA' if n_variables==7 else 'NONE')
        self.mode_hashes=mode_hashes or {}

    def evaluate_global(self, point):
        if len(point) not in (self.n_variables, self.storage_variables):
            raise SnapshotError(
                f"point has {len(point)} coordinates; expected "
                f"{self.n_variables} or {self.storage_variables}"
            )
        for patch in self.patches:
            if patch.contains(point, self.n_variables):
                local = patch.global_to_local(point, self.n_variables)
                return patch.patch_id, patch.evaluate_local(local)
        raise SnapshotError("point is outside every saved ADS patch")

    def physical_eta(self, normalized_eta):
        if self.parameterization=='SRP_MODES':
            raise SnapshotError('mode coordinates are not scalar eta')
        if self.n_variables == 6:
            return 0.0
        return self.domain_sigma * self.eta_sigma * normalized_eta

    def physical_modes(self,q7,q8):
        if self.parameterization!='SRP_MODES' or not all(math.isfinite(q) and abs(q)<=1+1e-13 for q in (q7,q8)):
            raise SnapshotError('invalid mode unit coordinates')
        return tuple(0.5*(1+max(-1,min(1,q))) for q in (q7,q8))


def _finite_float(value, label):
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise SnapshotError(f"invalid {label}: {value!r}") from exc
    if not math.isfinite(parsed):
        raise SnapshotError(f"non-finite {label}: {value!r}")
    return parsed


def _integer(value, label):
    try:
        return int(value)
    except (TypeError, ValueError) as exc:
        raise SnapshotError(f"invalid {label}: {value!r}") from exc


def _read_csv(path):
    if not path.is_file():
        raise SnapshotError(f"missing file: {path}")
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None:
            raise SnapshotError(f"CSV has no header: {path}")
        return list(reader)


def replay_split_history(entries, n_variables, storage_variables=None):
    """Return global unit-box center and width after signed splits."""
    if n_variables not in (6, 7, 8):
        raise SnapshotError(f"unsupported ADS dimension {n_variables}")
    columns=storage_variables or max(7,n_variables)
    center = [0.0] * columns
    width = [2.0] * columns
    if n_variables == 6 and storage_variables is None:
        width[6] = 0.0
    for entry in entries:
        direction = abs(int(entry))
        if entry == 0 or direction > n_variables:
            raise SnapshotError(
                f"invalid split direction {entry} for {n_variables}D snapshot"
            )
        index = direction - 1
        width[index] *= 0.5
        center[index] += (1.0 if entry > 0 else -1.0) * 0.5 * width[index]
    return tuple(center), tuple(width)


def _map_values(rows):
    values = {}
    for row in rows:
        key = (row.get("key") or "").strip()
        if not key:
            raise SnapshotError("empty key in snapshot map")
        if key in values:
            if key in {"da_order", "domain_sigma", "eta_sigma", "schema_version", "n_variables", "da_variables",
                       "parameterization", "mode_csv_sha256", "mode_metadata_sha256"}:
                raise SnapshotError(f"duplicate scalar snapshot map key {key}")
            continue
        values[key] = (row.get("value") or "").strip()
    return values


def load_snapshot(prefix):
    prefix = Path(prefix)
    files = {
        "meta": Path(f"{prefix}_meta.csv"),
        "coeff": Path(f"{prefix}_coeff.csv"),
        "splits": Path(f"{prefix}_splits.csv"),
        "map": Path(f"{prefix}_map.csv"),
    }
    meta_rows = _read_csv(files["meta"])
    coeff_rows = _read_csv(files["coeff"])
    split_rows = _read_csv(files["splits"])
    map_rows=_read_csv(files['map'])
    mapping = _map_values(map_rows)
    if not meta_rows:
        raise SnapshotError("snapshot contains no patches")
    for key in ("da_order", "domain_sigma", "eta_sigma"):
        if key not in mapping:
            raise SnapshotError(f"snapshot map is missing {key}")

    da_order_value = _finite_float(mapping["da_order"], "map value da_order")
    da_order = int(round(da_order_value))
    if abs(da_order_value - da_order) > 1.0e-12 or da_order < 0:
        raise SnapshotError(f"invalid da_order {da_order_value}")
    domain_sigma = _finite_float(
        mapping["domain_sigma"], "map value domain_sigma"
    )
    eta_sigma = _finite_float(mapping["eta_sigma"], "map value eta_sigma")
    if domain_sigma <= 0.0 or eta_sigma < 0.0:
        raise SnapshotError(
            "domain_sigma must be positive and eta_sigma must be non-negative"
        )
    n_variables = 7 if eta_sigma > 0.0 else 6
    version=_integer(mapping.get('schema_version','1'),'schema_version')
    columns=7
    parameterization='SRP_ETA' if n_variables==7 else 'NONE'
    mode_hashes={}
    if version not in (1,2):
        raise SnapshotError('unsupported snapshot schema')
    if version==2:
        n_variables=_integer(mapping.get('n_variables',''),'n_variables')
        columns=_integer(mapping.get('da_variables',''),'da_variables')
        parameterization=mapping.get('parameterization')
        expected_dim={'NONE':6,'SRP_ETA':7,'SRP_MODES':8}.get(parameterization)
        if n_variables!=expected_dim or not n_variables<=columns<=8:
            raise SnapshotError('v2 parameterization/dimensions disagree')
        if (parameterization=='SRP_ETA') != (eta_sigma>0):
            raise SnapshotError('v2 eta scale and parameterization disagree')
        roles=[row for row in map_rows if row['key']=='variable_role']
        expected_roles=['x','y','z','vx','vy','vz']+({'NONE':[], 'SRP_ETA':['eta'], 'SRP_MODES':['s','u']}[parameterization])
        expected_roles+=['inactive']*(columns-n_variables)
        if len(roles)!=columns or sorted((_integer(row['i'],'role index'),row['value']) for row in roles)!=list(enumerate(expected_roles,1)):
            raise SnapshotError('v2 variable roles do not match parameterization')
        if parameterization=='SRP_MODES':
            for key in ('mode_csv_sha256','mode_metadata_sha256'):
                value=mapping.get(key,'')
                if re.fullmatch('[0-9a-fA-F]{64}',value) is None:
                    raise SnapshotError(f'invalid v2 {key}')
                mode_hashes[key]=value.lower()

    meta = {}
    common_time = None
    for row in meta_rows:
        patch_id = _integer(row.get("patch_id", ""), "patch_id")
        if patch_id <= 0 or patch_id in meta:
            raise SnapshotError(f"invalid or duplicate patch_id {patch_id}")
        time_seconds = _finite_float(
            row.get("time_seconds", ""), "time_seconds"
        )
        center = tuple(
            _finite_float(row.get(f"center{i}", ""), f"center{i}")
            for i in range(1, columns+1)
        )
        width = tuple(
            _finite_float(row.get(f"width{i}", ""), f"width{i}")
            for i in range(1, columns+1)
        )
        split_count = _integer(row.get("split_count", ""), "split_count")
        if split_count < 0:
            raise SnapshotError(f"negative split_count for patch {patch_id}")
        if any(width[i] <= 0.0 for i in range(n_variables)):
            raise SnapshotError(f"patch {patch_id} has non-positive active width")
        if version==1 and n_variables == 6 and (
            abs(center[6]) > 1.0e-14 or abs(width[6]) > 1.0e-14
        ):
            raise SnapshotError(
                "6D control snapshot must have center7=width7=0"
            )
        if common_time is None:
            common_time = time_seconds
        elif abs(time_seconds - common_time) > 1.0e-8:
            raise SnapshotError("patch metadata contains inconsistent times")
        meta[patch_id] = (time_seconds, center, width, split_count)

    raw_histories = defaultdict(list)
    for row in split_rows:
        patch_id = _integer(row.get("patch_id", ""), "split patch_id")
        if patch_id not in meta:
            raise SnapshotError(
                f"split history references unknown patch {patch_id}"
            )
        raw_histories[patch_id].append(
            (
                _integer(row.get("sequence", ""), "split sequence"),
                _integer(row.get("split_direction", ""), "split direction"),
            )
        )

    histories = {}
    for patch_id, (_, center, width, split_count) in meta.items():
        ordered = sorted(raw_histories.get(patch_id, []))
        if (
            [item[0] for item in ordered] != list(range(1, split_count + 1))
            or len(ordered) != split_count
        ):
            raise SnapshotError(
                f"split history count/sequence mismatch for patch {patch_id}"
            )
        entries = tuple(item[1] for item in ordered)
        replay_center, replay_width = replay_split_history(
            entries, n_variables,columns if version==2 else None
        )
        if any(
            abs(center[i] - replay_center[i]) > 1.0e-12
            or abs(width[i] - replay_width[i]) > 1.0e-12
            for i in range(columns)
        ):
            raise SnapshotError(
                f"split history does not reproduce metadata for patch {patch_id}"
            )
        histories[patch_id] = entries

    terms = {patch_id: defaultdict(list) for patch_id in meta}
    for row in coeff_rows:
        patch_id = _integer(
            row.get("patch_id", ""), "coefficient patch_id"
        )
        if patch_id not in meta:
            raise SnapshotError(
                f"coefficient references unknown patch {patch_id}"
            )
        component = _integer(row.get("component", ""), "component")
        if not 1 <= component <= 6:
            raise SnapshotError(f"invalid state component {component}")
        powers = tuple(
            _integer(row.get(f"power{i}", ""), f"power{i}")
            for i in range(1, columns+1)
        )
        if any(power < 0 for power in powers) or sum(powers) > da_order:
            raise SnapshotError(
                f"polynomial power exceeds DA order for patch {patch_id}"
            )
        if any(powers[n_variables:]):
            raise SnapshotError(
                "6D control snapshot contains a non-zero power7 term"
            )
        terms[patch_id][component].append(
            (powers, _finite_float(row.get("coefficient", ""), "coefficient"))
        )

    patches = []
    for patch_id in sorted(meta):
        time_seconds, center, width, _ = meta[patch_id]
        missing = [
            component for component in range(1, 7)
            if component not in terms[patch_id]
        ]
        if missing:
            raise SnapshotError(
                f"patch {patch_id} has no coefficients for components {missing}"
            )
        patches.append(
            Patch(patch_id, time_seconds, center, width,
                  histories[patch_id], dict(terms[patch_id]))
        )
    return Snapshot(
        prefix, n_variables, da_order, domain_sigma, eta_sigma,
        common_time, patches, files,version,parameterization,mode_hashes
    )


def _halton(index, base):
    value = 0.0
    factor = 1.0 / base
    while index:
        value += factor * (index % base)
        index //= base
        factor /= base
    return value


def _padded(values, n_variables):
    return tuple(values[:n_variables]) + (0.0,) * (max(7,n_variables) - n_variables)


def generate_local_samples(n_variables, count, boundary_fraction):
    """Make deterministic local points including every corner and axis."""
    if n_variables not in (6, 7, 8):
        raise SnapshotError(f"unsupported ADS dimension {n_variables}")
    minimum = 2**n_variables + 2 * n_variables + 1
    if count < minimum:
        raise SnapshotError(
            f"points-per-patch must be at least {minimum} for {n_variables}D"
        )
    if (
        not math.isfinite(boundary_fraction)
        or not 0.0 <= boundary_fraction <= 1.0
    ):
        raise SnapshotError("boundary-fraction must lie in [0,1]")

    samples = [
        Sample("corner", _padded(bits, n_variables))
        for bits in itertools.product((-1.0, 1.0), repeat=n_variables)
    ]
    for dimension in range(n_variables):
        for sign in (-1.0, 1.0):
            values = [0.0] * n_variables
            values[dimension] = sign
            samples.append(Sample("axis", _padded(values, n_variables)))
    samples.append(Sample("center", (0.0,) * max(7,n_variables)))

    extra = count - len(samples)
    face_count = int(round(extra * boundary_fraction))
    if 0.0 < boundary_fraction < 1.0 and extra >= 2:
        face_count = max(1, min(extra - 1, face_count))
    for sample_index in range(face_count):
        index = sample_index + 1
        values = [
            2.0 * _halton(index, PRIMES[i]) - 1.0
            for i in range(n_variables)
        ]
        dimension = sample_index % n_variables
        values[dimension] = (
            -1.0 if (sample_index // n_variables) % 2 == 0 else 1.0
        )
        samples.append(Sample("face", _padded(values, n_variables)))
    for sample_index in range(extra - face_count):
        index = face_count + sample_index + 1
        values = [
            2.0 * _halton(index, PRIMES[i]) - 1.0
            for i in range(n_variables)
        ]
        samples.append(Sample("interior", _padded(values, n_variables)))
    return samples


def _sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _infer_run_prefix(snapshot_prefix):
    match = re.search(r"_patch_s\d{7}$", snapshot_prefix.name)
    if match is None:
        return None
    return snapshot_prefix.with_name(snapshot_prefix.name[: match.start()])


def _nominal_state(history_path, hour):
    rows = _read_csv(history_path)
    matches = [
        row for row in rows
        if abs(_finite_float(row.get("hour", ""), "history hour") - hour)
        <= 1.0e-8
    ]
    if len(matches) != 1:
        raise SnapshotError(
            f"history must contain exactly one row for hour {hour:.12g}"
        )
    row = matches[0]
    names = (
        "nominal_x_km", "nominal_y_km", "nominal_z_km",
        "nominal_vx_kms", "nominal_vy_kms", "nominal_vz_kms",
    )
    return tuple(_finite_float(row.get(name, ""), name) for name in names)


def _write_csv_atomic(path, header, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    try:
        with temporary.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(header)
            writer.writerows(rows)
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def _write_json_atomic(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    try:
        with temporary.open("w", encoding="utf-8") as stream:
            json.dump(
                payload, stream, indent=2, sort_keys=True, allow_nan=False
            )
            stream.write("\n")
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def _verify_saved_cloud(snapshot, run_prefix, nominal):
    initial_path = Path(f"{run_prefix}_initial_points.csv")
    shape_path = Path(f"{run_prefix}_shape.csv")
    points = {}
    for row in _read_csv(initial_path):
        point_id = _integer(row.get("point_id", ""), "initial point_id")
        if point_id in points:
            raise SnapshotError(f"duplicate initial point_id {point_id}")
        points[point_id] = tuple(
            _finite_float(row.get(f"u{i}", ""), f"u{i}")
            for i in range(1, snapshot.storage_variables+1)
        )

    hour = snapshot.time_seconds / 3600.0
    expected = {}
    for row in _read_csv(shape_path):
        if (
            abs(_finite_float(row.get("hour", ""), "shape hour") - hour)
            > 1.0e-8
        ):
            continue
        point_id = _integer(row.get("point_id", ""), "shape point_id")
        if point_id in expected:
            raise SnapshotError(
                f"duplicate shape point_id {point_id} at hour {hour:.12g}"
            )
        deviation = tuple(
            _finite_float(row.get(name, ""), name) for name in ERROR_NAMES
        )
        expected[point_id] = tuple(
            nominal[i] + deviation[i] for i in range(6)
        )
    if not points or set(points) != set(expected):
        raise SnapshotError(
            "initial point IDs and saved shape point IDs do not match"
        )

    max_component = [0.0] * 6
    max_position = 0.0
    max_velocity = 0.0
    worst_point = None
    worst_score = -1.0
    for point_id in sorted(points):
        _, actual = snapshot.evaluate_global(points[point_id])
        errors = [
            abs(actual[i] - expected[point_id][i]) for i in range(6)
        ]
        max_component = [
            max(max_component[i], errors[i]) for i in range(6)
        ]
        position = math.sqrt(sum(errors[i] ** 2 for i in range(3)))
        velocity = math.sqrt(sum(errors[i] ** 2 for i in range(3, 6)))
        score = max(
            position / POSITION_VERIFY_TOL_KM,
            velocity / VELOCITY_VERIFY_TOL_KMS,
        )
        if score > worst_score:
            worst_score = score
            worst_point = point_id
        max_position = max(max_position, position)
        max_velocity = max(max_velocity, velocity)

    passed = (
        max_position <= POSITION_VERIFY_TOL_KM
        and max_velocity <= VELOCITY_VERIFY_TOL_KMS
    )
    return {
        "status": "passed" if passed else "failed",
        "point_count": len(points),
        "worst_point_id": worst_point,
        "max_component_absolute_error": dict(
            zip(STATE_NAMES, max_component)
        ),
        "max_position_error_km": max_position,
        "max_velocity_error_kms": max_velocity,
        "position_tolerance_km": POSITION_VERIFY_TOL_KM,
        "velocity_tolerance_kms": VELOCITY_VERIFY_TOL_KMS,
        "initial_points_file": str(initial_path),
        "shape_file": str(shape_path),
    }


def reconstruct_snapshot(
    prefix,
    output=None,
    report_path=None,
    points_per_patch=4096,
    boundary_fraction=0.5,
    history_path=None,
    absolute_only=False,
    verify=True,
):
    snapshot = load_snapshot(prefix)
    prefix = Path(prefix)
    output_path = (
        Path(output)
        if output is not None
        else Path(f"{prefix}_reconstructed.csv")
    )
    json_path = (
        Path(report_path)
        if report_path is not None
        else Path(f"{prefix}_reconstruction.json")
    )
    run_prefix = _infer_run_prefix(prefix)
    history = (
        Path(history_path)
        if history_path is not None
        else (
            Path(f"{run_prefix}_history.csv")
            if run_prefix is not None
            else None
        )
    )
    nominal = None
    if not absolute_only:
        if history is None or not history.is_file():
            raise SnapshotError(
                "history file is required unless --absolute-only is used"
            )
        nominal = _nominal_state(
            history, snapshot.time_seconds / 3600.0
        )

    if verify:
        if run_prefix is None:
            raise SnapshotError(
                "cannot infer parent run prefix for saved-cloud verification"
            )
        verification_nominal = nominal
        if verification_nominal is None:
            inferred_history = Path(f"{run_prefix}_history.csv")
            if not inferred_history.is_file():
                raise SnapshotError(
                    "saved-cloud verification requires parent history"
                )
            verification_nominal = _nominal_state(
                inferred_history, snapshot.time_seconds / 3600.0
            )
        verification = _verify_saved_cloud(
            snapshot, run_prefix, verification_nominal
        )
        if verification["status"] != "passed":
            raise SnapshotError(
                "saved-cloud verification failed: "
                f"position={verification['max_position_error_km']:.6g} km, "
                f"velocity={verification['max_velocity_error_kms']:.6g} km/s"
            )
    else:
        verification = {
            "status": "skipped",
            "reason": "disabled by caller",
        }


    samples = generate_local_samples(
        snapshot.n_variables, points_per_patch, boundary_fraction
    )
    samples = [Sample(sample.kind, sample.coordinates[:snapshot.storage_variables] +
                      (0.0,) * max(0, snapshot.storage_variables-len(sample.coordinates)))
               for sample in samples]
    state_min = [math.inf] * 6
    state_max = [-math.inf] * 6
    error_min = [math.inf] * 6
    error_max = [-math.inf] * 6
    header = ["patch_id", "sample_id", "sample_kind"]
    header += [f"local{i}" for i in range(1, snapshot.storage_variables+1)]
    header += [f"u{i}" for i in range(1, snapshot.storage_variables+1)]
    header += (['s','u'] if snapshot.parameterization=='SRP_MODES' else ['eta']) + list(STATE_NAMES)
    if nominal is not None:
        header += list(ERROR_NAMES)

    def rows():
        for patch in snapshot.patches:
            for sample_id, sample in enumerate(samples, start=1):
                point = patch.local_to_global(
                    sample.coordinates, snapshot.n_variables
                )
                state = patch.evaluate_local(sample.coordinates)
                physical_parameters=snapshot.physical_modes(point[6],point[7]) if snapshot.parameterization=='SRP_MODES' else \
                    (snapshot.physical_eta(point[6] if snapshot.n_variables>6 else 0),)
                for i in range(6):
                    state_min[i] = min(state_min[i], state[i])
                    state_max[i] = max(state_max[i], state[i])
                row = [
                    patch.patch_id,
                    sample_id,
                    sample.kind,
                    *sample.coordinates[:snapshot.storage_variables],
                    *point,
                    *physical_parameters,
                    *state,
                ]
                if nominal is not None:
                    error = tuple(
                        state[i] - nominal[i] for i in range(6)
                    )
                    for i in range(6):
                        error_min[i] = min(error_min[i], error[i])
                        error_max[i] = max(error_max[i], error[i])
                    row.extend(error)
                yield row

    _write_csv_atomic(output_path, header, rows())

    report = {
        "format_version": snapshot.schema_version,
        "output_schema_version": snapshot.schema_version,
        "n_variables": snapshot.n_variables,
        "snapshot_prefix": str(prefix),
        "output_file": str(output_path),
        "report_file": str(json_path),
        "time_seconds": snapshot.time_seconds,
        "hour": snapshot.time_seconds / 3600.0,
        "ads_dimension": snapshot.n_variables,
        "da_order": snapshot.da_order,
        "domain_sigma": snapshot.domain_sigma,
        "eta_sigma": snapshot.eta_sigma,
        "patch_count": len(snapshot.patches),
        "points_per_patch": points_per_patch,
        "total_samples": points_per_patch * len(snapshot.patches),
        "boundary_fraction": boundary_fraction,
        "sampling": {
            "deterministic": True,
            "mandatory_points": "all local corners, signed axes, and center",
            "additional_points": "Halton face and interior points",
        },
        "state_min": dict(zip(STATE_NAMES, state_min)),
        "state_max": dict(zip(STATE_NAMES, state_max)),
        "source_files": {
            key: {"path": str(path), "sha256": _sha256(path)}
            for key, path in snapshot.source_files.items()
        },
        "structural_validation": {
            "status": "passed",
            "split_histories_replayed": True,
            "inactive_eta_dimension_zero": snapshot.n_variables == 6,
        },
        "saved_cloud_verification": verification,
        "verification": verification,
        "interpretation": (
            "Samples reconstruct the saved ADS polynomial approximation; "
            "they are not an independent high-fidelity truth guarantee "
            "between sampled points."
        ),
    }
    if nominal is not None:
        report["history_file"] = str(history)
        report["nominal_state"] = dict(zip(STATE_NAMES, nominal))
        report["deviation_min"] = dict(zip(ERROR_NAMES, error_min))
        report["deviation_max"] = dict(zip(ERROR_NAMES, error_max))
    if snapshot.schema_version==2:
        report.update(parameterization=snapshot.parameterization,mode_hashes=snapshot.mode_hashes)
    _write_json_atomic(json_path, report)
    return report


def _parser():
    parser = argparse.ArgumentParser(
        description=(
            "Reconstruct a dense point cloud from a saved SRP ADS "
            "patch snapshot."
        )
    )
    parser.add_argument(
        "--prefix",
        required=True,
        help=(
            "Snapshot prefix before "
            "_meta.csv/_coeff.csv/_splits.csv/_map.csv"
        ),
    )
    parser.add_argument("--points-per-patch", type=int, default=4096)
    parser.add_argument("--boundary-fraction", type=float, default=0.5)
    parser.add_argument("--output", help="Output reconstructed CSV path")
    parser.add_argument("--report", help="Output JSON report path")
    parser.add_argument(
        "--history", help="Override parent *_history.csv path"
    )
    parser.add_argument(
        "--absolute-only",
        action="store_true",
        help="Write absolute states without requiring nominal history",
    )
    parser.add_argument(
        "--no-verify",
        action="store_true",
        help="Skip verification against the saved 512-point cloud",
    )
    return parser


def main(argv=None):
    args = _parser().parse_args(argv)
    try:
        report = reconstruct_snapshot(
            prefix=args.prefix,
            output=args.output,
            report_path=args.report,
            points_per_patch=args.points_per_patch,
            boundary_fraction=args.boundary_fraction,
            history_path=args.history,
            absolute_only=args.absolute_only,
            verify=not args.no_verify,
        )
    except (SnapshotError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    print(
        f"Wrote {report['total_samples']} reconstructed samples "
        f"to {report['output_file']}"
    )
    print(f"Report: {report['report_file']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
