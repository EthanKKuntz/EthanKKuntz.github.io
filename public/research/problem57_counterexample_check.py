#!/usr/bin/env python3
"""Exact certificate checks for a counterexample to Green's Problem 57.

The unrestricted witness is encoded by exponent matrices for a primitive 48th
root of unity.  The restricted upper bound is certified by three Hermitian
positive-definite matrices over Q(omega), omega^2 + omega + 1 = 0.  The three
cases represent the possible equality patterns of (j_0,j_1,j_2) in Z/3Z.
"""

from __future__ import annotations

from collections import defaultdict
from fractions import Fraction as F


H = (20, -14, 5)
N = 48

# f1 is indexed by (x2,x3), f2 by (x1,x3), f3 by (x1,x2).
E1 = ((43, 35, 38), (1, 12, 9), (47, 20, 4))
E2 = ((13, 40, 30), (6, 5, 36), (7, 16, 35))
E3 = ((41, 8, 37), (21, 31, 46), (46, 42, 12))


def witness_coefficients() -> dict[int, int]:
    """Return c_m for T = sum_m c_m zeta_48^m."""
    coeffs: dict[int, int] = defaultdict(int)
    for x1 in range(3):
        for x2 in range(3):
            for x3 in range(3):
                exponent = (E1[x2][x3] + E2[x1][x3] + E3[x1][x2]) % N
                coeffs[exponent] += H[(x1 + x2 + x3) % 3]
    return dict(sorted(coeffs.items()))


# Elements a+b*w of Q(w), where w^2+w+1=0.
Qw = tuple[F, F]
ZERO: Qw = (F(0), F(0))
ONE: Qw = (F(1), F(0))


def qw_add(x: Qw, y: Qw) -> Qw:
    return x[0] + y[0], x[1] + y[1]


def qw_neg(x: Qw) -> Qw:
    return -x[0], -x[1]


def qw_sub(x: Qw, y: Qw) -> Qw:
    return qw_add(x, qw_neg(y))


def qw_mul(x: Qw, y: Qw) -> Qw:
    a, b = x
    c, d = y
    return a * c - b * d, a * d + b * c - b * d


def qw_inv(x: Qw) -> Qw:
    a, b = x
    norm = a * a - a * b + b * b
    return (a - b) / norm, -b / norm


def qw_div(x: Qw, y: Qw) -> Qw:
    return qw_mul(x, qw_inv(y))


def qw_conj(x: Qw) -> Qw:
    return x[0] - x[1], -x[1]


def qw_rat(x: int | F) -> Qw:
    return F(x), F(0)


def omega_pow(k: int) -> Qw:
    return (ONE, (F(0), F(1)), (F(-1), F(-1)))[k % 3]


def determinant(matrix: list[list[Qw]]) -> Qw:
    """Exact Gaussian-elimination determinant over Q(w)."""
    a = [row[:] for row in matrix]
    out = ONE
    n = len(a)
    for col in range(n):
        pivot = next(row for row in range(col, n) if a[row][col] != ZERO)
        if pivot != col:
            a[col], a[pivot] = a[pivot], a[col]
            out = qw_neg(out)
        q = a[col][col]
        out = qw_mul(out, q)
        for row in range(col + 1, n):
            ratio = qw_div(a[row][col], q)
            for j in range(col + 1, n):
                a[row][j] = qw_sub(a[row][j], qw_mul(ratio, a[col][j]))
    return out


def dual_matrix(js: tuple[int, int, int], y100: tuple[int, ...]) -> list[list[Qw]]:
    """Build D-C, where C=[[0,A*/2],[A/2,0]]."""
    A = [
        [qw_mul(qw_rat(H[(s + z) % 3]), omega_pow(js[z] * s)) for s in range(3)]
        for z in range(3)
    ]
    out = [[ZERO for _ in range(6)] for _ in range(6)]
    for i, value in enumerate(y100):
        out[i][i] = qw_rat(F(value, 100))
    for s in range(3):
        for z in range(3):
            out[s][3 + z] = qw_neg(qw_mul(qw_conj(A[z][s]), qw_rat(F(1, 2))))
            out[3 + z][s] = qw_neg(qw_mul(A[z][s], qw_rat(F(1, 2))))
    return out


CERTIFICATES = (
    ((0, 0, 0), (1478, 1478, 1478, 1478, 1478, 1478)),
    ((0, 0, 1), (1772, 1926, 1448, 1542, 1786, 1817)),
    ((0, 1, 2), (1713, 1713, 1713, 1713, 1713, 1713)),
)


def main() -> None:
    expected = {
        0: 30,
        1: 65,
        2: 20,
        3: 40,
        12: 5,
        13: 10,
        18: -14,
        20: -14,
        22: -42,
        23: -42,
        28: -14,
        38: 15,
        47: 40,
    }
    coeffs = witness_coefficients()
    assert coeffs == expected
    print("witness coefficients:", coeffs)

    for js, y100 in CERTIFICATES:
        assert sum(y100) < 10300
        matrix = dual_matrix(js, y100)
        minors = []
        for k in range(1, 7):
            minor = determinant([row[:k] for row in matrix[:k]])
            assert minor[1] == 0 and minor[0] > 0
            minors.append(minor[0])
        print(f"pattern {js}: sum(y)={F(sum(y100), 100)}")
        print("  leading minors:", minors)

    print("all exact certificate checks passed")


if __name__ == "__main__":
    main()
