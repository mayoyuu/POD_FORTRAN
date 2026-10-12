#!/usr/bin/env python3
"""Pack single-group SRP point CSV files into validated compressed HDF5."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import tempfile

import h5py
import numpy as np


SOURCE_SCHEMA = "srp_mc_initial_points_v1"
OUTPUT_SCHEMA = "srp_mc_initial_points_hdf5_v1"
STATE_COLUMNS = ("x_km", "y_km", "z_km", "vx_kms", "vy_kms", "vz_kms")
ERROR_COLUMNS = ("dx_km", "dy_km", "dz_km", "dvx_kms", "dvy_kms", "dvz_kms")
VALUE_COLUMNS = STATE_COLUMNS + ERROR_COLUMNS


class PointFormatError(ValueError):
    """Raised when a point CSV violates the fixed-design output contract."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _rows(path: Path):
    stream = path.open(newline="")
    try:
        schema = stream.readline().rstrip("\r\n")
        if schema != f"schema={SOURCE_SCHEMA}":
            raise PointFormatError("point CSV schema line is missing or unsupported")
        reader = csv.DictReader(stream)
        required = {"hour", "keyframe_kind", "sample_id", *VALUE_COLUMNS}
        if reader.fieldnames is None or set(reader.fieldnames) != required:
            raise PointFormatError(
                "point CSV header does not match srp_mc_initial_points_v1"
            )
        yield from reader
    finally:
        stream.close()


def _parse_row(row: dict[str, str], row_number: int):
    try:
        hour = float(row["hour"])
        sample_id = int(row["sample_id"])
        kind = row["keyframe_kind"]
        state = np.asarray(
            [float(row[name]) for name in STATE_COLUMNS], dtype=np.float64
        )
        error = np.asarray(
            [float(row[name]) for name in ERROR_COLUMNS], dtype=np.float64
        )
    except (KeyError, ValueError, TypeError) as exc:
        raise PointFormatError(f"malformed point row {row_number}") from exc
    if not math.isfinite(hour) or sample_id < 1 or not kind or "," in kind:
        raise PointFormatError(f"invalid point identity at row {row_number}")
    if not np.isfinite(state).all() or not np.isfinite(error).all():
        raise PointFormatError(f"nonfinite point value at row {row_number}")
    return hour, kind, sample_id, state, error


def _scan(path: Path, expected_points: int):
    if expected_points < 1:
        raise PointFormatError("expected point count must be positive")
    frames = []
    finished_keys = set()
    current_key = None
    current_ids = []
    baseline_ids = None
    row_count = 0

    def finish_frame():
        nonlocal baseline_ids
        if current_key is None:
            return
        if len(current_ids) != expected_points:
            raise PointFormatError(
                f"frame {current_key!r} has {len(current_ids)} points; "
                f"expected {expected_points}"
            )
        if len(set(current_ids)) != len(current_ids):
            raise PointFormatError(f"duplicate sample_id in frame {current_key!r}")
        if baseline_ids is None:
            baseline_ids = tuple(current_ids)
        elif tuple(current_ids) != baseline_ids:
            raise PointFormatError(f"sample_id order changed in frame {current_key!r}")
        frames.append(current_key)

    for row_count, row in enumerate(_rows(path), start=1):
        hour, kind, sample_id, _, _ = _parse_row(row, row_count)
        key = (hour, kind)
        if current_key is None:
            current_key = key
        elif key != current_key:
            finish_frame()
            finished_keys.add(current_key)
            if key in finished_keys:
                raise PointFormatError(f"noncontiguous duplicate frame {key!r}")
            current_key = key
            current_ids = []
        current_ids.append(sample_id)
    finish_frame()
    if row_count == 0:
        raise PointFormatError("point CSV contains no data rows")
    if len(set(frames)) != len(frames):
        raise PointFormatError("duplicate frame key")
    return frames, np.asarray(baseline_ids, dtype=np.int64), row_count


def pack_initial_points(
    csv_path: str | Path,
    h5_path: str | Path,
    expected_points: int,
    kind: str,
    report_path: str | Path | None = None,
) -> dict:
    """Validate and atomically pack one visualization or keyframe CSV."""
    source = Path(csv_path)
    output = Path(h5_path)
    if kind not in ("visualization", "keyframe"):
        raise PointFormatError("point role must be visualization or keyframe")
    frames, sample_ids, row_count = _scan(source, expected_points)
    report = {}
    report_sha256 = ""
    if report_path is not None:
        report_file = Path(report_path)
        try:
            report = json.loads(report_file.read_text())
        except (OSError, json.JSONDecodeError) as exc:
            raise PointFormatError("report JSON is missing or malformed") from exc
        report_sha256 = _sha256(report_file)

    output.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temp_name = tempfile.mkstemp(
        prefix=output.name + ".", suffix=".tmp", dir=output.parent
    )
    os.close(descriptor)
    temp_path = Path(temp_name)
    frame_count = len(frames)
    vector_chunks = (1, min(expected_points, 2048), 6)
    id_chunks = (1, min(expected_points, 8192))
    try:
        with h5py.File(temp_path, "w") as h5:
            h5.attrs["output_schema_version"] = OUTPUT_SCHEMA
            h5.attrs["source_schema_version"] = SOURCE_SCHEMA
            h5.attrs["point_role"] = kind
            h5.attrs["source_csv_sha256"] = _sha256(source)
            h5.attrs["source_csv_name"] = source.name
            h5.attrs["report_sha256"] = report_sha256
            h5.attrs["report_schema_version"] = str(
                report.get("output_schema_version", "")
            )
            h5.attrs["srp_model"] = str(report.get("srp_model", ""))
            h5.attrs["srp_attitude_mode"] = str(
                report.get("srp_attitude_mode", "")
            )
            h5.attrs["frame_count"] = frame_count
            h5.attrs["point_count_per_frame"] = expected_points
            h5.attrs["state_units"] = "km,km,km,km/s,km/s,km/s"
            h5.attrs["error_reference"] = "model nominal center at the same hour"

            h5.create_dataset(
                "hour",
                data=np.asarray([frame[0] for frame in frames], dtype=np.float64),
                dtype="f8",
            )
            string_dtype = h5py.string_dtype(encoding="utf-8")
            h5.create_dataset(
                "keyframe_kind",
                data=np.asarray([frame[1] for frame in frames], dtype=object),
                dtype=string_dtype,
            )
            id_dataset = h5.create_dataset(
                "sample_id",
                shape=(frame_count, expected_points),
                dtype="i8",
                chunks=id_chunks,
                compression="gzip",
                compression_opts=4,
                shuffle=True,
            )
            state_dataset = h5.create_dataset(
                "state",
                shape=(frame_count, expected_points, 6),
                dtype="f8",
                chunks=vector_chunks,
                compression="gzip",
                compression_opts=4,
                shuffle=True,
            )
            error_dataset = h5.create_dataset(
                "error",
                shape=(frame_count, expected_points, 6),
                dtype="f8",
                chunks=vector_chunks,
                compression="gzip",
                compression_opts=4,
                shuffle=True,
            )

            state = np.empty((expected_points, 6), dtype=np.float64)
            error = np.empty_like(state)
            frame_index = 0
            point_index = 0
            active_key = frames[0]
            for row_number, row in enumerate(_rows(source), start=1):
                hour, frame_kind, sample_id, row_state, row_error = _parse_row(
                    row, row_number
                )
                key = (hour, frame_kind)
                if key != active_key:
                    id_dataset[frame_index, :] = sample_ids
                    state_dataset[frame_index, :, :] = state
                    error_dataset[frame_index, :, :] = error
                    frame_index += 1
                    active_key = key
                    point_index = 0
                if sample_id != sample_ids[point_index]:
                    raise PointFormatError(
                        "sample order changed between validation passes"
                    )
                state[point_index, :] = row_state
                error[point_index, :] = row_error
                point_index += 1
            id_dataset[frame_index, :] = sample_ids
            state_dataset[frame_index, :, :] = state
            error_dataset[frame_index, :, :] = error
            h5.flush()
        os.replace(temp_path, output)
    except BaseException:
        try:
            temp_path.unlink()
        except FileNotFoundError:
            pass
        raise

    return {
        "output_schema_version": OUTPUT_SCHEMA,
        "source": str(source),
        "output": str(output),
        "source_csv_sha256": _sha256(source),
        "frame_count": frame_count,
        "point_count": expected_points,
        "row_count": row_count,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--expected-per-frame", required=True, type=int)
    parser.add_argument(
        "--point-role", required=True, choices=("visualization", "keyframe")
    )
    parser.add_argument("--report-json", type=Path)
    parser.add_argument("--manifest", type=Path)
    args = parser.parse_args()
    result = pack_initial_points(
        args.input,
        args.output,
        args.expected_per_frame,
        args.point_role,
        args.report_json,
    )
    text = json.dumps(result, indent=2) + "\n"
    if args.manifest:
        args.manifest.write_text(text)
    print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
