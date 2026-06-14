# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2", "mpmath"]
# ///
"""Certified finite-time UPPER bound on the recurrence time T_eps(N).

The node-0 velocity response of the ring is
    u(t) = sum_{r=1}^{N-1} sin(omega_r t) / (2 N sin(pi r/N)),   omega_r = 2 sin(pi r/N),   |u| <= U_N,
with the Bohr ceiling U_N = sum_r 1/(2 N sin(pi r/N)). The recurrence time to reach a fraction of the
ceiling is  T_eps(N) = inf{ T : A_N(T) >= (1-eps) U_N },  A_N(T)=max_{t<=T} u(t).

A single witness time gives a rigorous UPPER bound: if u(t*) >= (1-eps) U_N at some explicit t*, then
A_N(t*) >= (1-eps) U_N, hence T_eps(N) <= t*. This script PROVES such a bound for small primes with no
floating-point trust: a witness t* is located numerically, rationalized to t~, and the inequality
u(t~) >= 0.8 U_N is then verified in INTERVAL arithmetic (mpmath.iv) -- every sine, every reciprocal,
and the ceiling itself are rigorous enclosures. The verdict uses lower(u(t~)) >= 0.8 * upper(U_N), so it
is conservative in both directions.

This closes the *existence of an explicit finite reaching time* for each concrete prime N (eps=0.2),
upgrading the numerical recurrence-time data (verify_finite_time.py, Fig. reachability) to a proof. The
FAMILY asymptotic T_eps(N) ~ delta^{-m} stays open: the target pi/2 is inhomogeneous, so a uniform
quantitative Kronecker rate needs a Diophantine lower bound on linear forms in the omega_r (Baker-type),
not the homogeneous Dirichlet pigeonhole.
"""

from __future__ import annotations

from fractions import Fraction

import numpy as np
from mpmath import iv

FRAC = Fraction(4, 5)  # certify A_N(t*) >= 0.8 U_N  (eps = 0.2)
SEARCH_MARGIN = 0.82  # locate the witness where u first exceeds this (cushion absorbs interval width)
TDENOM = 1 << 34  # rationalization grid for the witness time


def find_witness(n: int, t_max: float, dt: float = 0.002) -> tuple[float, float]:
    """Float scan: ceiling U and the first t with u(t) >= SEARCH_MARGIN*U (None if not reached)."""
    r = np.arange(1, n)
    s = np.sin(np.pi * r / n)
    om = 2.0 * s
    b = 1.0 / (2.0 * n * s)
    cap = float(b.sum())
    t = np.arange(0.0, t_max, dt)
    u = (np.sin(np.outer(t, om)) * b).sum(axis=1)
    hit = np.where(u >= SEARCH_MARGIN * cap)[0]
    return cap, (float(t[hit[0]]) if hit.size else float("nan"))


def certify(n: int, t_star: float, dps: int = 60) -> tuple[bool, Fraction]:
    """Interval proof that u(t~) >= FRAC * U_N for a rational t~ near t_star (the verdict is the bool)."""
    iv.dps = dps
    tt_q = Fraction(round(t_star * TDENOM), TDENOM)
    tt = iv.mpf(tt_q.numerator) / iv.mpf(tt_q.denominator)  # exact rational point as a thin interval
    pi = iv.pi
    u_iv = iv.mpf(0)
    cap_iv = iv.mpf(0)
    for r in range(1, n):
        s = iv.sin(pi * r / n)  # sin(pi r / n) > 0
        b = 1 / (2 * n * s)
        u_iv += iv.sin(2 * s * tt) * b
        cap_iv += b
    frac_iv = iv.mpf(FRAC.numerator) / iv.mpf(FRAC.denominator)
    target = frac_iv * cap_iv  # interval enclosing 0.8 U_N
    return bool(u_iv.a >= target.b), tt_q  # conservative: lower(u) >= upper(0.8 U_N)


def u_over_cap(n: int, t: float) -> tuple[float, float]:
    """Float u(t)/U_N and U_N -- for display only, not part of the proof."""
    r = np.arange(1, n)
    s = np.sin(np.pi * r / n)
    b = 1.0 / (2.0 * n * s)
    return float((np.sin(2.0 * s * t) * b).sum() / b.sum()), float(b.sum())


def main() -> int:
    print("=" * 88)
    print("Certified finite-time upper bound:  T_{0.2}(N) <= t*  (witness + interval arithmetic)")
    print("=" * 88)
    # generous horizons (recurrence time grows ~ exp(0.63 N)); keep N small so the witness stays modest
    horizons = {5: 40.0, 7: 120.0, 11: 1500.0, 13: 6000.0}
    print(f"\n{'N':>3} {'U_N':>9} {'witness t~ (<= bound)':>22} {'cert. u(t~)/U_N':>17} {'verdict':>10}")
    ok_all = True
    for n, tmax in horizons.items():
        cap, t_star = find_witness(n, tmax)
        if not np.isfinite(t_star):
            print(f"{n:>3} {cap:>9.4f} {'(0.82 U_N not reached in horizon)':>22}")
            ok_all = False
            continue
        ok, tt_q = certify(n, t_star)
        ok_all &= ok
        ratio, cap_disp = u_over_cap(n, float(tt_q))
        print(f"{n:>3} {cap_disp:>9.4f} {float(tt_q):>22.5f} {ratio:>17.4f} "
              f"{('PROVED' if ok else 'FAIL'):>10}")
    print("\n  Reading: at the listed t~ the amplitude is rigorously >= 0.8 U_N, so the recurrence time")
    print("  T_{0.2}(N) is at most that t~ -- an explicit finite reaching time, proved in interval arithmetic.")
    print("  (The family rate T_eps ~ delta^{-m} stays open: the inhomogeneous target needs a Baker-type")
    print("   Diophantine bound, not the homogeneous Dirichlet pigeonhole.)")
    print("\n" + "=" * 88)
    print("RESULT:", "FINITE-TIME UPPER BOUNDS CERTIFIED" if ok_all else "CHECK FAILED")
    print("=" * 88)
    return 0 if ok_all else 1


if __name__ == "__main__":
    raise SystemExit(main())
