"""Shared selection policy for first-split times in SRP analysis tools."""
from __future__ import annotations

import math


def _finite_float(value):
    if value is None or isinstance(value, bool):
        return None
    try:
        converted = float(value)
    except (TypeError, ValueError):
        return None
    return converted if math.isfinite(converted) else None


def select_first_split_time(record: dict):
    """Return ``(time_nd, source)`` using event time before checkpoint time.

    A nonnegative ``first_split_time_accurate`` is the event-detected value.
    When it is unavailable or negative, a nonnegative checkpoint value is used.
    Negative sentinels mean that no split was detected; absent or malformed
    values are reported separately as missing.
    """
    event_time = _finite_float(record.get("first_split_time_accurate"))
    checkpoint_time = _finite_float(record.get("first_split_time"))

    if event_time is not None and event_time >= 0.0:
        return event_time, "event"
    if checkpoint_time is not None and checkpoint_time >= 0.0:
        return checkpoint_time, "checkpoint"
    if ((event_time is not None and event_time < 0.0)
            or (checkpoint_time is not None and checkpoint_time < 0.0)):
        return None, "never"
    return None, "missing"
