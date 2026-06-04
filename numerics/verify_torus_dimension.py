# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2"]
# ///
"""The large wave is one-dimensional: the ceiling order is ln N (d=1) but O(1) for d>=2.

On the d-dimensional discrete torus (Z_N)^d the Laplacian is a block-circulant-with-circulant-blocks
(BCCB) matrix, diagonalized by the d-dimensional DFT with eigenvalues
    lambda_k = sum_{i=1}^d 4 sin^2(pi k_i / N),   k in (Z_N)^d,   omega_k = sqrt(lambda_k).
The velocity Green's function obeys the UNCONDITIONAL Bohr/triangle bound A_N^{(d)} <= U_N^{(d)} with
    U_N^{(d)} = (1/N^d) sum_{k != 0} 1/omega_k.
Near the origin omega_k ~ (2 pi/N)|k|, so the small-k part of N^d U_N^{(d)} behaves like the lattice
sum sum 1/|k|, i.e. the infrared integral int d^d k / |k|: log-divergent in d=1 (=> U_N ~ (1/pi) ln N)
but CONVERGENT in d>=2 (=> U_N^{(d)} -> const). Hence the amplification is unbounded only in 1D; on
the 2- and 3-torus the velocity Green's function is uniformly bounded -- no large wave.

We verify: (i) d=1 matches (1/pi) ln N + C_0; (ii) d=2,3 saturate (U_{2N}-U_N -> 0), while the d=1
increment -> (1/pi) ln 2; the bound A_N <= U_N is the elementary triangle inequality (no independence).
"""

from __future__ import annotations

import numpy as np


def ceiling(n: int, d: int) -> float:
    """U_N^{(d)} = (1/N^d) sum_{k!=0} 1/omega_k on the d-torus (Z_N)^d."""
    s = 2.0 * np.sin(np.pi * np.arange(n) / n)  # |2 sin(pi k_i/N)| per axis; lambda = sum s^2
    if d == 1:
        w2 = s**2
    elif d == 2:
        w2 = s[:, None] ** 2 + s[None, :] ** 2
    else:
        w2 = s[:, None, None] ** 2 + s[None, :, None] ** 2 + s[None, None, :] ** 2
    w = np.sqrt(w2)
    inv = np.zeros_like(w)
    nz = w > 1e-12
    inv[nz] = 1.0 / w[nz]  # drop the single zero mode k=0
    return float(inv.sum() / n**d)


def main() -> int:
    print("=" * 66)
    print("Higher-dimensional tori: the large wave is one-dimensional")
    print("=" * 66)
    ok = True
    c0 = (np.euler_gamma + np.log(2 / np.pi)) / np.pi

    print("\n(1) d=1: U_N matches (1/pi) ln N + C_0")
    d1_ok = True
    for n in (64, 256, 1024, 4096):
        u = ceiling(n, 1)
        approx = np.log(n) / np.pi + c0
        d1_ok &= abs(u - approx) < 0.01
        print(f"    N={n:>5}: U_N={u:.5f}  (1/pi)lnN+C_0={approx:.5f}  diff={u - approx:+.2e}")
    ok &= d1_ok

    print("\n(2) Saturation: increment U_{2N}-U_N  (d=1 -> (ln2)/pi=0.2206; d>=2 -> 0)")
    print(f"    {'N':>5} {'d=1 incr':>10} {'d=2 incr':>10} {'d=3 incr':>10}")
    incr2, incr3 = [], []
    grid = [8, 16, 32, 64]
    for n in grid:
        i1 = ceiling(2 * n, 1) - ceiling(n, 1)
        i2 = ceiling(2 * n, 2) - ceiling(n, 2)
        i3 = ceiling(2 * n, 3) - ceiling(n, 3)
        incr2.append(i2)
        incr3.append(i3)
        print(f"    {n:>5} {i1:>10.4f} {i2:>10.4f} {i3:>10.4f}")
    # d=1 increment -> ln2/pi; d>=2 increments shrink toward 0
    sat_ok = abs((ceiling(8192, 1) - ceiling(4096, 1)) - np.log(2) / np.pi) < 0.01
    sat_ok &= abs(incr2[-1]) < abs(incr2[0]) and abs(incr2[-1]) < 0.05
    sat_ok &= abs(incr3[-1]) < abs(incr3[0]) and abs(incr3[-1]) < 0.03
    ok &= sat_ok
    print(f"    -> d=1 grows by (ln2)/pi each doubling; d=2,3 increments -> 0: "
          f"{'OK' if sat_ok else 'FAIL'}")

    print("\n(3) Bounded ceilings in d>=2 (A_N^{(d)} <= U_N^{(d)} by triangle ineq, unconditional)")
    u2 = ceiling(128, 2)
    u3 = ceiling(64, 3)
    print(f"    d=2: U_128 = {u2:.4f}   d=3: U_64 = {u3:.4f}   (both O(1), vs d=1 U_128={ceiling(128, 1):.4f})")

    print("\n" + "=" * 66)
    print("RESULT:", "TORUS DIMENSION VERIFIED" if ok else "CHECK FAILED")
    print("Ceiling order: ln N (d=1), O(1) (d>=2). The large-wave amplification is intrinsically 1D;")
    print("on the 2-/3-torus the velocity Green's function is uniformly bounded (no large wave).")
    print("=" * 66)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
