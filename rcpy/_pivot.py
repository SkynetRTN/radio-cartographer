"""Verbatim ports of C++ `Tools::performPivot` / `pivotSystem` /
`lowerTriangleSolver` from src/Tools.cpp.

These are the linear-system solvers used throughout the C++ codebase for
weighted least-squares regressions (BG sub, surface fit, edge fit, etc).
The elimination order is unusual — from the LAST column down to column 1
— and produces a lower-triangular system that is then forward-solved.
On near-singular matrices this gives different numeric "garbage" than
`np.linalg.solve` (LU partial pivoting) or `np.linalg.lstsq` (SVD).
Matching it bit-for-bit is required for parity with the C++ reference
on degenerate inputs (e.g. 2 distinct samples passed to a 3-term
regression, or N-clustered samples passed to M10).
"""
from __future__ import annotations

import numpy as np


def _perform_pivot(cc: int, A_in: np.ndarray, b_in: np.ndarray) -> np.ndarray:
    """Solve A · x = b using the C++ Tools::performPivot algorithm.

    Eliminates from column cc-1 downwards (NOT standard forward LU),
    reducing A to lower-triangular form, then forward-substitutes via
    `lowerTriangleSolver`. Verbatim port of Tools.cpp:171-212 +
    Tools.cpp:157-169.
    """
    A = A_in.flatten().astype(np.float64).copy()
    b = b_in.astype(np.float64).copy()

    with np.errstate(invalid="ignore", divide="ignore"):
        for i in range(cc - 1, 0, -1):
            pivot = abs(A[cc * i + i])
            pivot_index = i
            for j in range(i, -1, -1):
                if abs(A[cc * j + i]) > pivot:
                    pivot = abs(A[cc * j + i])
                    pivot_index = j

            if pivot_index == i:
                for j in range(cc):
                    A[cc * i + j] = A[cc * i + j] / pivot
                b[i] = b[i] / pivot
            else:
                for j in range(cc):
                    tmp = A[cc * i + j]
                    A[cc * i + j] = A[cc * pivot_index + j] / pivot
                    A[cc * pivot_index + j] = tmp
                tmp = b[i]
                b[i] = b[pivot_index] / pivot
                b[pivot_index] = tmp

            for j in range(i - 1, -1, -1):
                multiplier = A[j * cc + i] / A[i * cc + i]
                for k in range(cc):
                    A[j * cc + k] -= multiplier * A[i * cc + k]
                b[j] -= multiplier * b[i]

        x = np.zeros(cc, dtype=np.float64)
        for i in range(cc):
            sub = 0.0
            for j in range(i + 1):
                sub += A[cc * i + j] * x[j]
            x[i] = (b[i] - sub) / A[cc * i + i]
        return x


def _pivot_system_const_term(cc: int, A_in: np.ndarray, b_in: np.ndarray) -> float:
    """Return the constant term (x[0]) of the system A · x = b.

    Matches the C++ pattern in Cartographer::surfaceFit10/6/3:
        finalHold = Tools::pivotSystem(columnCount, A, b);
        finalAnswer = finalHold[1][0] / finalHold[0][0];

    After pivotSystem reduces A to lower-triangular form, the first
    row's solution is just b[0] / A[0][0] (the i=0 case of
    lowerTriangleSolver, where the sub-sum is zero). This is exactly
    equivalent to taking x[0] from _perform_pivot, so we just call
    that and return its first element.
    """
    return float(_perform_pivot(cc, A_in, b_in)[0])
