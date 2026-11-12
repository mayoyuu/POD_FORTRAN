"""Shared position-covariance metrics for SRP study analysis tools."""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Sequence


LU_KM = 384400.0


@dataclass(frozen=True)
class PositionCovarianceMetrics:
    """One-sigma principal axes and directions of a position covariance."""

    axes_km: tuple[float, float, float]
    directions: tuple[
        tuple[float, float, float],
        tuple[float, float, float],
        tuple[float, float, float],
    ]
    total_rms_km: float

    @property
    def longest_axis_km(self) -> float:
        return self.axes_km[0]

    @property
    def longest_direction(self) -> tuple[float, float, float]:
        return self.directions[0]


def position_covariance_metrics(
        cov_state: Sequence[float],
        *,
        lu_km: float = LU_KM,
        psd_rtol: float = 1.0e-10) -> PositionCovarianceMetrics:
    """Return ordered position-covariance RMS axes and total RMS in kilometres.

    cov_state is the row-major flattened 6x6 state covariance written by the
    study executable. Small serialization asymmetries are removed before the
    symmetric eigenvalue decomposition. Eigenvector signs are canonicalized
    because an axis is physically unoriented.
    """
    if len(cov_state) != 36:
        raise ValueError("cov_state must contain exactly 36 values")
    values = [float(value) for value in cov_state]
    if not all(math.isfinite(value) for value in values):
        raise ValueError("cov_state must contain only finite values")
    if not math.isfinite(lu_km) or lu_km <= 0.0:
        raise ValueError("lu_km must be finite and positive")
    if not math.isfinite(psd_rtol) or psd_rtol < 0.0:
        raise ValueError("psd_rtol must be finite and nonnegative")

    covariance = [
        [0.5 * (values[6 * i + j] + values[6 * j + i])
         for j in range(3)]
        for i in range(3)
    ]
    eigenvalues, eigenvectors = _symmetric_eigenpairs_3x3(covariance)
    order = sorted(range(3), key=lambda index: eigenvalues[index], reverse=True)
    eigenvalues = [eigenvalues[index] for index in order]
    eigenvectors = [
        tuple(eigenvectors[row][index] for row in range(3))
        for index in order
    ]

    scale = max(max(abs(value) for value in eigenvalues), float.fromhex("0x1p-1022"))
    if eigenvalues[-1] < -psd_rtol * scale:
        raise ValueError(
            "position covariance is not positive semidefinite: "
            f"minimum eigenvalue={eigenvalues[-1]:.6e}")
    eigenvalues = [max(value, 0.0) for value in eigenvalues]

    directions = []
    for vector in eigenvectors:
        direction = list(vector)
        pivot = max(range(3), key=lambda index: abs(direction[index]))
        if direction[pivot] < 0.0:
            direction = [-value for value in direction]
        directions.append(tuple(float(value) for value in direction))

    axes_km = tuple(math.sqrt(value) * lu_km for value in eigenvalues)
    total_rms_km = math.sqrt(sum(eigenvalues)) * lu_km
    return PositionCovarianceMetrics(
        axes_km=axes_km,
        directions=tuple(directions),
        total_rms_km=total_rms_km,
    )


def _symmetric_eigenpairs_3x3(matrix):
    """Jacobi decomposition of a real symmetric 3x3 matrix."""
    a = [row[:] for row in matrix]
    vectors = [[1.0 if i == j else 0.0 for j in range(3)]
               for i in range(3)]

    for _ in range(50):
        p, q = max(((0, 1), (0, 2), (1, 2)),
                   key=lambda pair: abs(a[pair[0]][pair[1]]))
        off_diagonal = a[p][q]
        scale = max(max(abs(value) for row in a for value in row),
                    float.fromhex("0x1p-1022"))
        if abs(off_diagonal) <= 1.0e-15 * scale:
            break

        tau = (a[q][q] - a[p][p]) / (2.0 * off_diagonal)
        if tau == 0.0:
            tangent = 1.0
        else:
            tangent = math.copysign(
                1.0 / (abs(tau) + math.sqrt(1.0 + tau * tau)), tau)
        cosine = 1.0 / math.sqrt(1.0 + tangent * tangent)
        sine = tangent * cosine

        app, aqq = a[p][p], a[q][q]
        for k in range(3):
            if k == p or k == q:
                continue
            akp, akq = a[k][p], a[k][q]
            a[k][p] = a[p][k] = cosine * akp - sine * akq
            a[k][q] = a[q][k] = sine * akp + cosine * akq
        a[p][p] = (cosine * cosine * app
                   - 2.0 * sine * cosine * off_diagonal
                   + sine * sine * aqq)
        a[q][q] = (sine * sine * app
                   + 2.0 * sine * cosine * off_diagonal
                   + cosine * cosine * aqq)
        a[p][q] = a[q][p] = 0.0

        for k in range(3):
            vkp, vkq = vectors[k][p], vectors[k][q]
            vectors[k][p] = cosine * vkp - sine * vkq
            vectors[k][q] = sine * vkp + cosine * vkq

    return [a[i][i] for i in range(3)], vectors
