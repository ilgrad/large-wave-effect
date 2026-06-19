# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2"]
# ///
"""Finite-time large wave: the infinite-time ceiling U_N ~ (1/pi) ln N is approached only over Diophantine
(astronomically long) time scales; the physically reachable amplitude A_N(T) is governed by the horizon T.

The node-0 velocity response is
    u(t) = sum_{r=1}^{N-1} sin(omega_r t) / (2 N sin(pi r/N)),   omega_r = 2 sin(pi r/N),
with |u(t)| <= U_N. Reaching (1-eps) U_N requires sin(omega_r t) ~ 1 for ALL r simultaneously, i.e. a
simultaneous Diophantine approximation of m = floor(N/2) phases to pi/2. By Dirichlet, aligning m phases
to tolerance delta costs time up to ~ delta^{-m}, so the recurrence time

    T_eps(N) = inf { T : A_N(T) >= (1-eps) U_N }

grows at least EXPONENTIALLY in N. We define A_N(T)=max_{t<=T} u(t), measure T_eps(N) for primes, and show
that at any FIXED horizon T the finite-time amplitude saturates well below the (slowly growing) ceiling.
This separates infinite-time amplification (the log law) from finite-time amplification (the physics).
"""

from __future__ import annotations

import numpy as np


def u_ceiling(n: int) -> float:
    r = np.arange(1, n)
    return float(np.sum(1.0 / np.sin(np.pi * r / n)) / (2 * n))


def response_max(n: int, t_grid: np.ndarray) -> np.ndarray:
    """Running max of u(t) over the grid (cumulative), for node 0."""
    r = np.arange(1, n)
    omega = 2.0 * np.sin(np.pi * r / n)
    coeff = 1.0 / (2 * n * np.sin(np.pi * r / n))
    u = (np.sin(np.outer(t_grid, omega)) * coeff).sum(axis=1)
    return np.maximum.accumulate(u)


def t_eps(n: int, eps: float, t_max: float, dt: float = 0.02) -> float:
    """First time the running max reaches (1-eps) U_N; inf (-> t_max+) if not reached within t_max."""
    target = (1.0 - eps) * u_ceiling(n)
    t_grid = np.arange(0.0, t_max, dt)
    rmax = response_max(n, t_grid)
    hit = np.argmax(rmax >= target)
    return float(t_grid[hit]) if rmax[hit] >= target else float("inf")


def main() -> int:
    print("=" * 78)
    print("Finite-time large wave: ceiling U_N ~ (1/pi) ln N reached only on Diophantine time scales")
    print("=" * 78)
    ok = True

    print("\n(A) Recurrence time T_eps(N) to reach (1-eps) U_N grows steeply with N (primes)")
    print(f"    {'N':>4} {'U_N':>7} {'T_0.30':>9} {'T_0.20':>9} {'T_0.10':>9}")
    primes = [5, 7, 11, 13, 17, 19]
    rows = []

    def fmt(x: float) -> str:
        return f"{x:>9.1f}" if np.isfinite(x) else f"{'>cap':>9}"

    for n in primes:
        t30 = t_eps(n, 0.30, 4000.0)
        t20 = t_eps(n, 0.20, 20000.0)
        t10 = t_eps(n, 0.10, 60000.0)
        rows.append((n, t20))
        print(f"    {n:>4} {u_ceiling(n):>7.4f} {fmt(t30)} {fmt(t20)} {fmt(t10)}")
    # fit the asymptotic range only: N=5,7 are pre-asymptotic (T~O(1)) and inflate the slope (~0.63)
    finite = [(n, t) for n, t in rows if np.isfinite(t) and n >= 11]
    if len(finite) >= 3:
        ns = np.array([n for n, _ in finite], float)
        lts = np.log(np.array([t for _, t in finite], float))
        slope = np.polyfit(ns, lts, 1)[0]
        grows = slope > 0.1
        ok &= grows
        print(f"    -> log T_0.20 grows ~ {slope:.2f} per unit N over N>=11 (super-polynomial): {'OK' if grows else 'FAIL'}")
        print("       quantified law T_eps ~ eps^(-a m), a~0.32, over N<=71: verify_finite_time_lower_bound.py")

    print("\n(B) At a FIXED horizon T the finite-time amplitude saturates below the ceiling")
    print(f"    {'N':>4} {'U_N':>7} {'A_N(T=300)':>11} {'A_N(2e4)':>9} {'ratio@300':>10}")
    sat = []
    for n in (5, 11, 17, 23, 29):
        t1 = np.arange(0.0, 300.0, 0.02)
        t2 = np.arange(0.0, 20000.0, 0.05)
        a1 = float(response_max(n, t1)[-1])
        a2 = float(response_max(n, t2)[-1])
        sat.append(a1 / u_ceiling(n))
        print(f"    {n:>4} {u_ceiling(n):>7.4f} {a1:>11.4f} {a2:>9.4f} {a1 / u_ceiling(n):>10.4f}")
    drops = sat[-1] < sat[0]
    ok &= drops
    print(f"    -> at fixed T=300, A_N(T)/U_N DECREASES with N ({sat[0]:.3f} -> {sat[-1]:.3f}):")
    print(f"       more frequencies are harder to align in bounded time: {'OK' if drops else 'FAIL'}")

    print("\n(C) Physical reading")
    print("    Infinite-time amplification:  A_N = U_N ~ (1/pi) ln N (grows without bound).")
    print("    Finite-time amplification:    A_N(T) is horizon-limited; the log ceiling is a")
    print("    mathematical sup attained only after Diophantine-long waiting. A real large wave")
    print("    is the PARTIAL, early alignment -- already a several-fold over-amplification.")

    print("\n" + "=" * 78)
    print("RESULT:", "FINITE-TIME SEPARATION VERIFIED" if ok else "CHECK FAILED")
    print("=" * 78)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
