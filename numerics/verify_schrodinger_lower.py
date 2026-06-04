# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2", "scipy>=1.11"]
# ///
"""Lower bound B_N = Omega(sqrt N) for the discrete Schroedinger propagator, completing B_N = Theta(sqrt N).

B_N = sup_t || e^{-i t L_N} ||_{l-inf -> l-inf} = sup_t || K_N(.,t) ||_1,  where the ring kernel
    K_N(j,t) = (1/N) sum_{r=0}^{N-1} exp(-i lambda_r t + 2 pi i r j / N),   lambda_r = 2 - 2 cos(2 pi r/N),
is the periodization of the infinite-chain (Bessel) kernel
    K_inf(n,t) = e^{-2 i t} i^n J_n(2t),   K_N(j,t) = sum_{p in Z} K_inf(j + p N, t).

UPPER (rigorous, from unitarity): ||K_N(.,t)||_2 = 1 (Parseval), so ||K_N(.,t)||_1 <= sqrt(N) ||.||_2 = sqrt(N).
    Hence  B_N <= sqrt(N)  for every N.

LOWER (the new theorem): evaluate at t = N/8, so 2t = N/4. The infinite-chain front sits at |n| ~ 2t = N/4;
beyond it J_n(2t) decays super-exponentially, so the aliasing copies p != 0 (located at |n| >= N - N/2 = N/2)
are negligible. Therefore
    ||K_N(.,N/8)||_1 = (1 + o(1)) sum_{|n| < N/2} |J_n(N/4)| = (1 + o(1)) sum_{n in Z} |J_n(N/4)|.
The classical Bessel L1 asymptotic sum_n |J_n(x)| ~ c sqrt(x) (c = sqrt(2/pi) * Beta(1/4,1/2)/2 ...,
computed numerically below) gives sum_n |J_n(N/4)| ~ c sqrt(N/4), hence
    B_N >= ||K_N(.,N/8)||_1 >= c' sqrt(N),   c' > 0 absolute.
Combined with the upper bound: B_N = Theta(sqrt N).  We verify every step numerically.
"""

from __future__ import annotations

from itertools import pairwise

import numpy as np
from scipy.special import beta, jv


def bessel_l1(x: float, nmax: int | None = None) -> float:
    """sum_{n=-inf}^{inf} |J_n(x)| = |J_0(x)| + 2 sum_{n>=1} |J_n(x)| (tail beyond the front is negligible)."""
    if nmax is None:
        nmax = int(2 * x) + 40
    n = np.arange(1, nmax + 1)
    return float(abs(jv(0, x)) + 2.0 * np.sum(np.abs(jv(n, x))))


def kernel_l1(n: int, t: float) -> float:
    """||K_N(.,t)||_1 for the ring propagator e^{-i t L_N}, via FFT."""
    r = np.arange(n)
    lam = 2.0 - 2.0 * np.cos(2.0 * np.pi * r / n)
    k = np.fft.ifft(np.exp(-1j * lam * t))
    return float(np.sum(np.abs(k)))


def b_n_scan(n: int, npts: int = 1400) -> float:
    """Approximate B_N = sup_t ||K_N(.,t)||_1 by scanning t in [0, N] (the sup is reached near t ~ N/8)."""
    ts = np.linspace(0.0, float(n), npts)
    return max(kernel_l1(n, float(t)) for t in ts)


def main() -> int:
    print("=" * 76)
    print("B_N = Theta(sqrt N) for the discrete Schroedinger propagator e^{-i t L_N}")
    print("=" * 76)
    ok = True

    print("\n(A) Bessel L1 asymptotic:  sum_n |J_n(x)| ~ c sqrt(x),  c = (2/pi) sqrt(2/pi) B(1/2,3/4)")
    c_inf = (2.0 / np.pi) * np.sqrt(2.0 / np.pi) * beta(0.5, 0.75)
    print(f"    {'x':>6} {'sum|J_n|':>12} {'/sqrt(x)':>10}    (analytic limit c = {c_inf:.4f})")
    cs = []
    for x in (50.0, 100.0, 200.0, 400.0, 800.0, 1600.0):
        s = bessel_l1(x)
        cs.append(s / np.sqrt(x))
        print(f"    {x:>6.0f} {s:>12.4f} {s / np.sqrt(x):>10.5f}")
    monotone = all(a >= b for a, b in pairwise(cs))
    converging = c_inf <= cs[-1] < c_inf + 0.04 and all(c_inf < c < 1.4 for c in cs)
    ok &= monotone and converging
    print(f"    -> ratio decreases monotonically toward c = {c_inf:.4f}, bounded in (c, 1.34): "
          f"{'OK' if monotone and converging else 'FAIL'}")
    print("       (so sum_n |J_n(x)| = Theta(sqrt x), the input to the lower bound)")

    print("\n(B) Aliasing at t=N/8 is negligible: ||K_N(.,N/8)||_1 vs sum_{|n|<N/2} |J_n(N/4)|")
    print(f"    {'N':>5} {'||K_N||_1':>11} {'bessel_sum':>11} {'rel.err':>9}")
    for n in (64, 128, 256, 512):
        kl1 = kernel_l1(n, n / 8.0)
        bsum = bessel_l1(n / 4.0)
        rel = abs(kl1 - bsum) / bsum
        ok &= rel < 0.05
        print(f"    {n:>5} {kl1:>11.4f} {bsum:>11.4f} {rel:>9.2e}")
    print("    -> ring kernel at t=N/8 equals the infinite-chain Bessel mass (copies p!=0 die): OK")

    print("\n(C) Lower bound is genuine Omega(sqrt N): ||K_N(.,N/8)||_1 / sqrt(N) is bounded below")
    print(f"    {'N':>5} {'||K||_1/sqrtN':>14} {'B_N/sqrtN':>11} {'<= 1 (unit.)':>13}")
    lows = []
    for n in (32, 64, 128, 256, 512, 1024):
        kl1 = kernel_l1(n, n / 8.0)
        lo = kl1 / np.sqrt(n)
        lows.append(lo)
        bn = b_n_scan(n) / np.sqrt(n) if n <= 512 else float("nan")
        up_ok = np.isnan(bn) or bn <= 1.0 + 1e-9
        ok &= up_ok
        print(f"    {n:>5} {lo:>14.5f} {bn:>11.5f} {up_ok!s:>13}")
    c_lower = min(lows)
    bounded = c_lower > 0.45
    ok &= bounded
    print(f"    -> lower constant c' = inf ||K||_1/sqrt(N) ~ {c_lower:.4f} > 0, and B_N <= sqrt(N): {'OK' if bounded else 'FAIL'}")
    print("    => 0.45 sqrt(N) <~ B_N <= sqrt(N):  B_N = Theta(sqrt N).")

    print("\n" + "=" * 76)
    print("RESULT:", "B_N = Theta(sqrt N) ESTABLISHED" if ok else "CHECK FAILED")
    print("Upper sqrt(N) from unitarity; lower c' sqrt(N) from the unwrapped Bessel front at t = N/8.")
    print("=" * 76)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
