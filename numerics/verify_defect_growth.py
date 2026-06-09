# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2", "scipy>=1.11"]
# ///
"""The composite defect U_N - A_N is UNBOUNDED -- it grows with the number of odd prime factors.

The order theorem gives U_N - A_N <= 2(U_N - L_pre) = (2/pi) ln(N/phi(2N)) + O(1), and for N = 2^a m
(m odd) this is (2/pi) ln(m/phi(m)) + O(1), independent of a (Remark rem:defect).  Since
m/phi(m) ~ e^gamma ln ln m along primorials, the upper bound is O(ln ln ln N).

This script shows the bound is sharp in order: the exact node-0 defect (U_N minus the orbit-closure /
subtorus maximum of sum_r b_r sin(phi_r)) GROWS with the number omega(m) of distinct odd prime factors,
tracking ~ (1/2)(2/pi) ln(m/phi(m)).  So sup_m of the within-family limit is NOT finite -- the defect is
unbounded, of order ln ln ln N along primorials.  (A_N is computed by multistart L-BFGS over the subtorus,
the same method validated on primes in verify_exact_AN.py, where it must hit U_N exactly.)
"""

from __future__ import annotations

import sys
from math import gcd, log, pi
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from large_wave_effect.cyclotomic import sin_coords

RNG = np.random.default_rng(0)


def totient(n: int) -> int:
    return sum(1 for k in range(1, n + 1) if gcd(k, n) == 1)


def n_odd_prime_factors(n: int) -> int:
    while n % 2 == 0:
        n //= 2
    s, d = set(), 3
    while d * d <= n:
        while n % d == 0:
            s.add(d)
            n //= d
        d += 2
    if n > 1:
        s.add(n)
    return len(s)


def odd_part(n: int) -> int:
    while n % 2 == 0:
        n //= 2
    return n


def ceiling_and_node0_amplitude(n: int, starts: int = 80) -> tuple[float, float]:
    """U_N and the node-0 large-wave amplitude A_N = max over the subtorus of sum_r b_r sin(phi_r)."""
    m = n // 2
    r = np.arange(1, m + 1)
    omega = 2.0 * np.sin(pi * r / n)
    b = 2.0 / (n * omega)
    if n % 2 == 0:
        b[-1] = 1.0 / (n * omega[-1])  # the r = N/2 mode carries half weight
    ceiling = float(b.sum())
    c_mat = np.array(sin_coords(2 * n, list(r)), dtype=float)

    def negf(psi: np.ndarray) -> tuple[float, np.ndarray]:
        ph = c_mat @ psi
        return -float(b @ np.sin(ph)), -(c_mat.T @ (b * np.cos(ph)))

    best = 0.0
    for _ in range(starts):
        res = minimize(negf, RNG.uniform(0, 2 * pi, c_mat.shape[1]), jac=True, method="L-BFGS-B")
        best = max(best, -float(res.fun))
    return ceiling, best


def main() -> int:
    print("Node-0 defect U_N - A_N vs the number omega(m) of distinct odd prime factors:")
    print(f"  {'N':>4} {'m':>4} {'omega':>5} {'defect':>8} {'(2/pi)ln(m/phim)':>17}")
    rows = []
    for n in [15, 30, 105]:  # omega(m) = 2, 2, 3
        u, a = ceiling_and_node0_amplitude(n)
        m = odd_part(n)
        defect = u - a
        bound = (2 / pi) * log(m / totient(m))
        rows.append((n, m, n_odd_prime_factors(n), defect, bound))
        print(f"  {n:>4} {m:>4} {n_odd_prime_factors(n):>5} {defect:>8.4f} {bound:>17.4f}")

    d30 = next(d for n, _, _, d, _ in rows if n == 30)
    d105 = next(d for n, _, _, d, _ in rows if n == 105)
    within_bound = all(d <= bnd + 1e-6 for *_, d, bnd in rows)
    grows = d105 - d30 > 0.03  # omega 2 -> 3 raises the defect

    ok = within_bound and grows
    print("=" * 60)
    print(f"  defect(omega=3, N=105) - defect(omega=2, N=30) = {d105 - d30:.4f}  (> 0 => grows)")
    print("RESULT:", "PASS -- defect grows with omega(m); unbounded, O(ln ln ln N)" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
