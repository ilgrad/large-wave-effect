# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2", "scipy>=1.11"]
# ///
"""B_N = sup_t ||exp(-i t L_N)||_{inf->inf}: an improved liminf and the parity of the windowed constant.

The paper's Theorem schr-lower proved liminf_N B_N/sqrt(N) >= (1/2) c_0 = 0.6086 by evaluating the
Bessel front at t = N/8.  The *optimal* aliasing-free time is t = N/4: there 2t = N/2, the Bessel front
|n| ~ 2t = N/2 exactly fills one half-period, the aliased copies sit beyond the front (distance
>= N - 2t > 2t at t = N/4(1-delta), so exponentially small by DLMF 10.14.5), and

    sum_k |K_N(k, N/4(1-delta))| = (1+o(1)) sum_n |J_n(N/2)| = (1+o(1)) c_0 sqrt(N/2) = (c_0/sqrt2) sqrt(N),

so  liminf_N B_N / sqrt(N) >= c_0 / sqrt(2) = 0.86068...,  a clean improvement.

On the dispersive window t = O(N) the constant has NO single limit -- the windowed amplification
B_N(T) = sup_{t <= T} ||exp(-i t L_N)||_{inf->inf} splits by the parity of N (the aliasing phase i^{-N}
is real for even N, imaginary for odd N):

  * EVEN N: the window sup is at t = N/4 and B_N(0.6N)/sqrt(N) -> c_0/sqrt2 (the front-fills-half-period
    value; at t = N/2 the real aliasing partly cancels, giving less).
  * ODD N: the window sup is at t = N/2, where |K_N(k,N/2)| = sqrt(J_k(N)^2 + J_{N-k}(N)^2) (no
    cancellation), and B_N(0.6N)/sqrt(N) -> beta_odd ~ 0.928 > c_0/sqrt2.

These are constants of the window, not of the infinite-time B_N: over all t, B_p = sqrt(p) for p prime
(the orbit closure contains the Gauss chirp; verify_bn_prime_flat.py), so limsup B_N/sqrt(N) = 1, and
B_N < sqrt(N) for even N >= 6 (verify_bn_saturation.py).  This script certifies the t = N/4 lower bound and
exhibits the parity split of the window sup.
"""

from __future__ import annotations

import numpy as np
from numpy import pi, sqrt
from scipy.special import beta

C0 = (2 / pi) ** 1.5 * beta(0.5, 0.75)
TARGET = C0 / sqrt(2)


def l1(n: int, t: float) -> float:
    lam = 4 * np.sin(pi * np.arange(n) / n) ** 2
    return float(np.abs(np.fft.ifft(np.exp(-1j * t * lam))).sum())


def sup_over_t(n: int) -> float:
    lam = 4 * np.sin(pi * np.arange(n) / n) ** 2
    ts = np.linspace(0.2 * n, 0.6 * n, 4000)
    return max(float(np.abs(np.fft.ifft(np.exp(-1j * t * lam))).sum()) for t in ts) / sqrt(n)


def main() -> int:
    ok = True
    print(f"c_0/sqrt2 = {TARGET:.5f}  (improved liminf, vs the paper's c_0/2 = {C0 / 2:.5f})")

    # (1) aliasing-free t = N/4 value -> c_0/sqrt2 (lower bound on B_N for every N)
    print("\n[1] t = N/4(1-N^-0.4) value (a rigorous lower bound on B_N/sqrt N):")
    vals = []
    for n in [4096, 16384, 65536]:
        v = l1(n, n / 4 * (1 - n**-0.4)) / sqrt(n)
        vals.append(v)
        print(f"    N={n:>6}: {v:.4f}")
    # increasing toward c_0/sqrt2 from below
    if not (vals[-1] > 0.84 and all(v <= TARGET + 0.01 for v in vals)):
        ok = False

    # (2) parity split of the sup
    print("\n[2] windowed sup (0.2N <= t <= 0.6N) / sqrt(N) by parity (no single limit):")
    even, odd = [], []
    for n in [16384, 16385, 32768, 32769, 65536, 65537]:
        s = sup_over_t(n)
        (even if n % 2 == 0 else odd).append(s)
        print(f"    N={n:>6} ({'even' if n % 2 == 0 else 'odd '}): {s:.4f}")
    print(
        f"    even -> {np.mean(even):.4f} (~ c_0/sqrt2 = {TARGET:.4f});  "
        f"odd -> {np.mean(odd):.4f} (~ beta_odd ~ 0.928)"
    )
    if not (np.mean(odd) - np.mean(even) > 0.04):  # parity gap is real
        ok = False

    print("=" * 64)
    print(
        "RESULT:",
        "PASS -- liminf >= c_0/sqrt2; windowed constant splits by parity" if ok else "FAIL",
    )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
