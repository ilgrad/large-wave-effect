# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2"]
# ///
"""Spectral determinant of L_N and the spanning trees of the cycle (Kirchhoff).

The spectral zeta zeta_{L_N}(s) = sum_{r=1}^{N-1} lambda_r^{-s} continues analytically to s=0 with
    zeta'_{L_N}(0) = - sum_{r != 0} log lambda_r = - log (det' L_N),
the (zeta-regularized) spectral determinant. Since lambda_r = 4 sin^2(pi r/N) and the classical product
    prod_{r=1}^{N-1} 2 sin(pi r/N) = N,
we get det' L_N = prod_{r != 0} lambda_r = N^2, hence by the matrix-tree theorem the number of
spanning trees of the cycle C_N is
    tau(C_N) = det'(L_N) / N = N.
So the large-wave spectral object (the zeta of L_N) ties at s=0 to a combinatorial invariant of the
graph. All verified below.
"""

from __future__ import annotations

import numpy as np


def lambdas(n: int) -> np.ndarray:
    return 4.0 * np.sin(np.pi * np.arange(1, n) / n) ** 2


def main() -> int:
    print("=" * 66)
    print("Spectral determinant of L_N and spanning trees of C_N")
    print("=" * 66)
    ok = True

    print(f"\n  {'N':>5} {'prod 2sin':>12} {'det L_N=prod lam':>18} {'N^2':>10} {'tau=det/N':>10}")
    for n in (4, 8, 16, 31, 64, 257, 1024):
        lam = lambdas(n)
        prod2sin = float(np.prod(2.0 * np.sin(np.pi * np.arange(1, n) / n)))
        det_prime = float(np.prod(lam))  # = product of nonzero eigenvalues
        tau = det_prime / n
        good = abs(prod2sin - n) / n < 1e-6 and abs(det_prime - n * n) / (n * n) < 1e-6
        good = good and abs(tau - n) / n < 1e-6
        ok &= good
        print(f"  {n:>5} {prod2sin:>12.4f} {det_prime:>18.6g} {n * n:>10} {tau:>10.4f}  {'OK' if good else 'FAIL'}")

    # zeta'(0) = -log det' = -2 ln N  (spectral determinant via -sum log lambda)
    print("\n  zeta'_{L_N}(0) = -sum log lambda_r = -2 ln N  (spectral determinant = N^2)")
    z_ok = True
    for n in (16, 64, 256):
        zp0 = -float(np.sum(np.log(lambdas(n))))
        target = -2.0 * np.log(n)
        z_ok &= abs(zp0 - target) < 1e-9
        print(f"    N={n:>4}: zeta'(0)={zp0:.5f}  -2 ln N={target:.5f}  {'OK' if abs(zp0 - target) < 1e-9 else 'FAIL'}")
    ok &= z_ok

    print("\n" + "=" * 66)
    print("RESULT:", "SPANNING TREES / SPECTRAL DETERMINANT VERIFIED" if ok else "CHECK FAILED")
    print("det'(L_N)=N^2, tau(C_N)=N: the s=0 value of the large-wave spectral zeta is a graph invariant.")
    print("=" * 66)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
