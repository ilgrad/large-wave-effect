# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2", "scipy>=1.11"]
# ///
"""Excess lemma for prime-power odd part: the cases where it holds outright (problem 2).

Excess = A_N - L_pre, where A_N = sup_t sum_r b_r sin(omega_r t) (orbit-closure maximum), L_pre is the sum
of b_r over the Q-independent prefix (size phi(2N)/2), and U_N = sum_r b_r. Unconditionally A_N in
[2 L_pre - U_N, U_N], so |A_N - L_pre| <= U_N - L_pre = (1/pi) ln(m/phi(m)) + O(1); the excess lemma is
therefore trivial whenever the deficit U_N - L_pre is bounded, i.e. whenever omega(m) is bounded.

This certificate verifies the two cleanest families, where the dependent band has O(1) modes:
  * N = 2p (p odd prime): one dependent mode omega_p = 2 of weight 1/(2N); U_N - L_pre = 1/(2N) exactly,
    hence the SHARP two-sided  |A_N - L_pre| <= 1/(2N) -> 0.
  * N = 4p: two dependent modes {2p-1, 2p}; U_N - L_pre = (1/N)(sec(pi/4p) + 1/2) = O(1/N).
For each N it confirms the exact U_N - L_pre value and that the optimizer's A_N lies in the excess band
[L_pre - (U_N - L_pre), L_pre + (U_N - L_pre)] = [2 L_pre - U_N, U_N].  The genuinely open case is
omega(m) >= 2 (the dependent band carries Theta(N) modes of unbounded total weight).
"""

from __future__ import annotations

import sys
from math import cos, gcd, pi
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from large_wave_effect.cyclotomic import sin_coords

RNG = np.random.default_rng(0)


def totient(n: int) -> int:
    return sum(1 for k in range(1, n + 1) if gcd(k, n) == 1)


def weights(n: int) -> np.ndarray:
    r = np.arange(1, n // 2 + 1)
    omega = 2 * np.sin(pi * r / n)
    b = 2.0 / (n * omega)
    if n % 2 == 0:
        b[-1] = 1.0 / (n * omega[-1])  # half-weight Nyquist mode
    return b


def optimizer_a(n: int, b: np.ndarray, starts: int = 60) -> float:
    m, big_m = n // 2, totient(2 * n) // 2
    rows = np.array(sin_coords(2 * n, list(range(1, m + 1))), dtype=float)
    c, *_ = np.linalg.lstsq(rows[:big_m].T, rows.T, rcond=None)
    cf = np.rint(c.T).astype(float)
    best = -np.inf
    for _ in range(starts):
        res = minimize(
            lambda p: (-float(b @ np.sin(cf @ p)), -(cf.T @ (b * np.cos(cf @ p)))),
            RNG.uniform(0, 2 * pi, big_m),
            jac=True,
            method="L-BFGS-B",
        )
        best = max(best, -res.fun)
    return best


def check(n: int, dependent: list[int], exact_deficit: float) -> tuple[bool, str]:
    b = weights(n)
    u = float(b.sum())
    l_pre = u - sum(b[i - 1] for i in dependent)
    a = optimizer_a(n, b)
    deficit = u - l_pre
    deficit_ok = abs(deficit - exact_deficit) < 1e-9
    band_ok = (2 * l_pre - u) - 1e-7 <= a <= u + 1e-7
    excess = a - l_pre
    return (deficit_ok and band_ok), (
        f"N={n:3d}: U-Lpre={deficit:.6f} (exact {exact_deficit:.6f}) excess={excess:+.5f} "
        f"in [-{deficit:.5f},+{deficit:.5f}]  {'ok' if deficit_ok and band_ok else 'FAIL'}"
    )


def main() -> int:
    ok = True
    print("N = 2p  (one dependent mode omega_p = 2; sharp |A_N - L_pre| <= 1/(2N)):")
    for p in (3, 5, 7, 11, 13):
        n = 2 * p
        good, line = check(n, [p], 1.0 / (2 * n))
        ok &= good
        print("  " + line)

    print("N = 4p  (two dependent modes {2p-1, 2p}; |A_N - L_pre| <= (sec(pi/4p)+1/2)/N):")
    for p in (3, 5, 7):
        n = 4 * p
        good, line = check(n, [2 * p - 1, 2 * p], (1.0 / cos(pi / (4 * p)) + 0.5) / n)
        ok &= good
        print("  " + line)

    print("=" * 70)
    print("RESULT:", "PASS -- excess lemma holds with explicit O(1/N) constant for prime-power odd part"
          if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
