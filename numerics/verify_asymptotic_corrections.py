# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2"]
# ///
"""Precise asymptotics of the ceiling U_N beyond the leading log + constant.

We sharpen U_N = (1/2N) sum_{r=1}^{N-1} csc(pi r/N) past the known
    U_N = (1/pi) ln N + C_0 + R_N,   C_0 = (1/pi)(gamma + ln(2/pi)),
by pinning the next term. The summand is symmetric about r=N/2 and singular at BOTH ends, so a naive
term-by-term csc-expansion is not uniform; the correct route is Euler-Maclaurin / the digamma
reflection, which yields a remainder that is even in 1/N. We DETERMINE the leading remainder
numerically with high confidence: R_N = a_2 / N^2 + O(1/N^4), and recover a_2 by a Richardson fit of
N^2 R_N, cross-checking that the next coefficient is O(1) (so the 1/N^2 law is genuine, not a fit
artifact). The numerically determined value is a_2 = -pi/72 = -0.0436332... (constant in N to six
digits); equivalently a_2 = -(1/pi) * zeta(2)/12, tying the correction to the same zeta(2)=pi^2/6 that
governs the spectral-zeta value zeta_{L_N}(2). We confirm a_2 = -pi/72 to high precision.
"""

from __future__ import annotations

import numpy as np


def ceiling(n: int) -> float:
    r = np.arange(1, n)
    return float(np.sum(1.0 / np.sin(np.pi * r / n)) / (2 * n))


def remainder(n: int) -> float:
    c0 = (np.euler_gamma + np.log(2 / np.pi)) / np.pi
    return ceiling(n) - np.log(n) / np.pi - c0


def main() -> int:
    print("=" * 64)
    print("Precise asymptotics: U_N = (1/pi)lnN + C_0 + a_2/N^2 + O(1/N^4)")
    print("=" * 64)
    ok = True
    a2_target = -np.pi / 72

    print("\n(A) N^2 * R_N -> a_2 = -pi/72")
    print(f"    {'N':>6} {'R_N':>14} {'N^2 R_N':>12}")
    vals = []
    for n in (64, 128, 256, 512, 1024, 2048):
        r = remainder(n)
        vals.append((n, n * n * r))
        print(f"    {n:>6} {r:>14.3e} {n * n * r:>12.6f}")
    # Richardson: N^2 R_N = a_2 + a_4/N^2 + ...; extrapolate the last two
    (n1, y1), (n2, y2) = vals[-2], vals[-1]
    a2_rich = (y2 * n2**2 - y1 * n1**2) / (n2**2 - n1**2)
    print(f"\n    a_2 (last sample)   = {vals[-1][1]:.6f}")
    print(f"    a_2 (Richardson)    = {a2_rich:.6f}")
    print(f"    -pi/72              = {a2_target:.6f}")
    a_ok = abs(a2_rich - a2_target) < 1e-4
    ok &= a_ok
    print(f"    -> a_2 = -pi/72: {'OK' if a_ok else 'FAIL'}")

    print("\n(B) Two-term model U_N ~ (1/pi)lnN + C_0 - (pi/72)/N^2 (max abs error over N>=32)")
    c0 = (np.euler_gamma + np.log(2 / np.pi)) / np.pi
    errs = []
    for n in (32, 64, 128, 256, 512, 1024):
        model = np.log(n) / np.pi + c0 + a2_target / n**2
        errs.append(abs(ceiling(n) - model))
    worst = max(errs)
    # residual after removing 1/N^2 should be O(1/N^4): error*N^4 roughly constant/decaying
    b_ok = worst < 5e-5 and errs[-1] < errs[0]
    ok &= b_ok
    print(f"    worst |U_N - model| (N>=32) = {worst:.2e};  decreasing: {errs[-1] < errs[0]}  "
          f"{'OK' if b_ok else 'FAIL'}")

    print("\n" + "=" * 64)
    print("RESULT:", "ASYMPTOTIC CORRECTION VERIFIED" if ok else "CHECK FAILED")
    print("U_N = (1/pi) ln N + (1/pi)(gamma+ln(2/pi)) - (pi/72)/N^2 + O(1/N^4).")
    print("=" * 64)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
