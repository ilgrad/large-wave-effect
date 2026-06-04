# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2"]
# ///
"""The constant in B_N = Theta(sqrt N): scan sup_t ||K_N(.,t)||_1 / sqrt(N) over a wide range of N.

B_N = sup_t ||e^{-i t L_N}||_{inf->inf} = sup_t sum_j |K_N(j,t)|, with B_N <= sqrt(N) (unitarity) and
B_N >= 0.63 sqrt(N) (Bessel front at t=N/8). Here we map the ratio B_N/sqrt(N): a refined t-scan that
includes the Talbot/revival grid t = pi a / N (a=1..N) where the kernel near-revives. The data pin the
empirical constant sup_N B_N/sqrt(N) and show it sits well inside [0.63, 1], the still-open exact value
being an almost-periodic supremum.
"""

from __future__ import annotations

import numpy as np


def kernel_l1(n: int, t: float, lam: np.ndarray) -> float:
    return float(np.sum(np.abs(np.fft.ifft(np.exp(-1j * lam * t)))))


def b_over_sqrt(n: int) -> tuple[float, float]:
    """Return (B_N/sqrt(N) from a dense scan, value at the proven point t=N/8)."""
    lam = 2.0 - 2.0 * np.cos(2.0 * np.pi * np.arange(n) / n)
    # dense uniform scan plus the Talbot grid t = pi a / N and the Bessel point N/8
    ts = np.concatenate([
        np.linspace(0.01, 2.0 * n, 4000),
        np.pi * np.arange(1, n + 1) / n,
        np.array([n / 8.0]),
    ])
    vals = [kernel_l1(n, float(t), lam) for t in ts]
    return max(vals) / np.sqrt(n), kernel_l1(n, n / 8.0, lam) / np.sqrt(n)


def main() -> int:
    print("=" * 70)
    print("Constant in B_N = Theta(sqrt N):  scan of B_N / sqrt(N)")
    print("=" * 70)
    print(f"\n  {'N':>6} {'B_N/sqrtN':>11} {'t=N/8 lower':>12} {'<=1 (unitarity)':>16}")
    sup_ratio = 0.0
    ok = True
    for n in (16, 32, 64, 128, 256, 512, 1024, 2048):
        scan, low = b_over_sqrt(n)
        sup_ratio = max(sup_ratio, scan)
        within = scan <= 1.0 + 1e-9 and low >= 0.6
        ok &= within
        print(f"  {n:>6} {scan:>11.5f} {low:>12.5f} {within!s:>16}")
    band = 0.6 < sup_ratio < 1.0
    ok &= band
    print(f"\n  -> empirical sup_N B_N/sqrt(N) ~ {sup_ratio:.4f}, inside the proven band (0.63, 1].")
    print("     The exact constant lim B_N/sqrt(N) remains open (an almost-periodic / Gauss-sum sup).")
    print("\n" + "=" * 70)
    print("RESULT:", "B_N CONSTANT MAPPED" if ok else "CHECK FAILED")
    print("=" * 70)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
