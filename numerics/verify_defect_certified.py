# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2"]
# ///
"""RIGOROUS certified lower bounds on the composite defect U_N - A_N (the gap A_N < U_N).

A_N = max over the orbit-closure subtorus of g(psi) = sum_r b_r sin(phi_r), where the free phases are the
independent prefix psi_1..psi_M (M = phi(2N)/2) and each dependent phi_r is an INTEGER combination
phi_r = sum_s C[r,s] psi_s of them (C computed exactly from the cyclotomic coordinates, then verified to be
integer and to satisfy sum_s C[r,s] omega_s = omega_r).

The Lipschitz-grid certificate -- the elementary member of the Lasserre moment-SOS hierarchy -- gives a
RIGOROUS upper bound:  max_psi g <= (max of g over a grid of cell centers) + (1/2) sum_s L_s h,  where
L_s = sum_r b_r |C[r,s]| bounds |d g / d psi_s| and h = 2 pi / m is the grid step.  Hence

    U_N - A_N  >=  U_N - [gridmax + (h/2) sum_s L_s]   (a proof, no SDP solver).

This certifies A_N < U_N with an explicit gap for N = 6, 9, 12, 15 (M up to 4), upgrading the defect lower
bound from numerical to rigorous.  (The grid is exponential in M, so primorials -- where the defect grows --
stay out of reach; their growth is the numerical trend of verify_defect_growth.py.)
"""

from __future__ import annotations

import sys
from fractions import Fraction
from math import gcd
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from large_wave_effect.cyclotomic import sin_coords


def totient(n: int) -> int:
    return sum(1 for k in range(1, n + 1) if gcd(k, n) == 1)


def weights(n: int) -> tuple[float, np.ndarray, np.ndarray]:
    m = n // 2
    r = np.arange(1, m + 1)
    omega = 2.0 * np.sin(np.pi * r / n)
    b = 2.0 / (n * omega)
    if n % 2 == 0:
        b[-1] = 1.0 / (n * omega[-1])
    return float(b.sum()), b, omega


def integer_coord_matrix(n: int) -> np.ndarray:
    """C (m x M), phi_r = sum_s C[r,s] psi_s; exact via Fraction solve, verified integer + consistent."""
    m, big_m = n // 2, totient(2 * n) // 2
    co = [[Fraction(x) for x in row] for row in sin_coords(2 * n, list(range(1, m + 1)))]
    dim = len(co[0])
    # solve, for each r, co[r] = sum_s x_s co[s] (s < M) by Gaussian elimination on the prefix columns
    aug = [[co[s][j] for s in range(big_m)] + [co[r][j] for r in range(m)] for j in range(dim)]
    pr = 0
    for col in range(big_m):
        sel = next((i for i in range(pr, dim) if aug[i][col] != 0), None)
        if sel is None:
            continue
        aug[pr], aug[sel] = aug[sel], aug[pr]
        p = aug[pr][col]
        aug[pr] = [v / p for v in aug[pr]]
        for i in range(dim):
            if i != pr and aug[i][col] != 0:
                f = aug[i][col]
                aug[i] = [a - f * b for a, b in zip(aug[i], aug[pr], strict=True)]
        pr += 1
        if pr == big_m:
            break
    c = np.zeros((m, big_m))
    for s in range(big_m):
        for r in range(m):
            val = aug[s][big_m + r]
            assert val.denominator == 1, f"non-integer coordinate at N={n}"
            c[r, s] = float(val)
    return c.astype(int)


def grid_max(c: np.ndarray, b: np.ndarray, m: int) -> float:
    big_m = c.shape[1]
    h = 2 * np.pi / m
    ctr = (np.arange(m) + 0.5) * h
    cf = c.astype(float)
    best = -1e9
    if big_m == 1:
        ph = cf @ ctr[None, :]
        return float((b @ np.sin(ph)).max())
    for x0 in ctr:  # block over the first axis to bound memory
        mesh = np.meshgrid(*([ctr] * (big_m - 1)), indexing="ij")
        psi = np.stack([np.full(mesh[0].size, x0)] + [mm.ravel() for mm in mesh])
        best = max(best, float((b @ np.sin(cf @ psi)).max()))
    return best


def main() -> int:
    print("Certified gap  A_N < U_N  (rigorous Lipschitz-grid upper bound on A_N):")
    print(f"  {'N':>3} {'M':>2} {'U_N':>7} {'cert. A_N <=':>12} {'defect >=':>10}")
    grids = {2: 2000, 3: 260, 4: 64}
    ok = True
    for n in [6, 9, 12, 15]:
        u, b, omega = weights(n)
        c = integer_coord_matrix(n)
        assert np.allclose(c @ omega[: c.shape[1]], omega, atol=1e-9), f"coord inconsistency N={n}"
        big_m = c.shape[1]
        gm = grid_max(c, b, grids[big_m])
        slack = (np.abs(c) * b[:, None]).sum(axis=0).sum() * (2 * np.pi / grids[big_m]) / 2
        upper = gm + slack
        defect_lb = u - upper
        if defect_lb <= 0:
            ok = False
        print(f"  {n:>3} {big_m:>2} {u:>7.4f} {upper:>12.4f} {defect_lb:>10.4f}")
    print("=" * 52)
    print("RESULT:", "PASS -- A_N < U_N certified with explicit gap" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
