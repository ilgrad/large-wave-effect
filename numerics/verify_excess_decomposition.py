# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2", "scipy>=1.11"]
# ///
"""Verify Lemma lem:triplet (per-relation excess bound) for the composite excess (A1).

PROVED (paper Lemma lem:triplet; this script confirms numerically):
  g(a,b,c) := max_{u,v}[ c sin(u+v) - a(1-sin u) - b(1-sin v) ] <= (c^2/2)(1/a+1/b).
Proof: p=1-sin u, q=1-sin v in [0,2]; sin(u+v) <= sqrt(2p)+sqrt(2q) (|1-p|<=1, sqrt(2p-p^2)<=sqrt(2p));
then max_{p,q>=0}[c(sqrt(2p)+sqrt(2q))-ap-bq] = c^2/(2a)+c^2/(2b).
For a support-3 relation omega_c=omega_a+omega_b the weights b_x=1/(N sin(pi x/N)) satisfy 1/b_a+1/b_b=1/b_c,
so g(b_a,b_b,b_c) <= b_c/2.  This bounds each clean relation IN ISOLATION; the full excess is NOT the sum of
these (comb generators of support>3 + shared prefix phases couple the relations) -- the open core.

CRUX (honest): excess/(U-Lpre) ~ 0.33 along primorials (omega=2,3,4); U-Lpre ~ (1/pi)ln(m/phi(m)) ~ lnlnln N
is ~flat, so O(1) vs Theta(lnlnln N) is numerically UNDECIDABLE -- the lemma's form needs proof, not numerics.
"""
from __future__ import annotations

import sys
import time
from math import pi
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, str(Path(__file__).resolve().parent))
from verify_excess_omega2 import integer_coord_matrix, omega_count, optimize_AN, totient, weights


def g_triplet(a: float, b: float, c: float, starts: int = 20) -> float:
    rng = np.random.default_rng(0)
    best = -np.inf
    for _ in range(starts):
        r = minimize(lambda x: -(c * np.sin(x[0] + x[1]) - a * (1 - np.sin(x[0])) - b * (1 - np.sin(x[1]))),
                     rng.uniform(0, 2 * pi, 2), method="L-BFGS-B")
        best = max(best, -r.fun)
    return best


def main() -> int:
    t0 = time.time()
    ok = True
    print("=" * 88)
    print("Lemma lem:triplet (per-relation excess bound)  +  the excess/(U-Lpre) crux")
    print("=" * 88)

    # (1) the inequality g(a,b,c) <= (c^2/2)(1/a+1/b)
    print("\n[1] g(a,b,c) <= (c^2/2)(1/a+1/b)   (400 random a,b,c>0; max violation must be <= 0)")
    rng = np.random.default_rng(1)
    worst = -np.inf
    for _ in range(400):
        a, b = rng.uniform(0.01, 2.0, 2)
        c = rng.uniform(0.01, min(a, b) + 0.5)
        worst = max(worst, g_triplet(a, b, c, starts=12) - 0.5 * c * c * (1 / a + 1 / b))
    lemma_ok = worst < 1e-6
    ok &= lemma_ok
    print(f"    max(g - bound) = {worst:.2e}   {'OK' if lemma_ok else 'VIOLATED'}")

    # (2) weight identity 1/b_a+1/b_b=1/b_c and the per-triplet bound g(b_a,b_b,b_c) <= b_c/2
    print("\n[2] support-3 relations omega_c=omega_a+omega_b:  1/b_a+1/b_b=1/b_c  and  g <= b_c/2")
    print(f"    {'N':>4} {'#trip':>6} {'reln-id':>8} {'max g_c/(b_c/2)':>16}")
    id_ok = bound_ok = True
    for n in [15, 21, 33, 45]:
        c = integer_coord_matrix(n)
        m, big_m = n // 2, totient(2 * n) // 2
        _u, b, _o = weights(n)
        trip = [(big_m + i, *np.where(row != 0)[0]) for i, row in enumerate(c[big_m:m])
                if (row != 0).sum() == 2 and list(row[row != 0]) == [1, 1]]
        idmax = max((abs((1 / b[a] + 1 / b[bb]) - 1 / b[ci]) for ci, a, bb in trip), default=0.0)
        ratios = [g_triplet(b[a], b[bb], b[ci]) / (b[ci] / 2) for ci, a, bb in trip]
        rmax = max(ratios, default=0.0)
        id_ok &= idmax < 1e-9
        bound_ok &= rmax <= 1.0 + 1e-6
        print(f"    {n:>4} {len(trip):>6} {('OK' if idmax < 1e-9 else 'FAIL'):>8} {rmax:>16.4f}")
    ok &= id_ok and bound_ok

    # (3) the crux: excess/(U-Lpre) along primorials (numerically undecidable O(1) vs Theta(lnlnln N))
    print("\n[3] CRUX  excess/(U-Lpre)  along primorials  (ratio ~0.33; lnlnln-flat => form needs PROOF)")
    print(f"    {'N':>5} {'omega':>5} {'M':>5} {'excess':>8} {'U-Lpre':>8} {'ratio':>7}")
    for n in [15, 105, 1155]:
        starts = 150 if n <= 1155 else 40
        a_n, _sv, bb = optimize_AN(n, seed=0, starts=starts)
        big_m = totient(2 * n) // 2
        u, lp = float(bb.sum()), float(bb[:big_m].sum())
        print(f"    {n:>5} {omega_count(n):>5} {big_m:>5} {a_n - lp:>8.4f} {u - lp:>8.4f} {(a_n - lp) / (u - lp):>7.4f}")

    print("\n" + "=" * 88)
    print("RESULT:", "LEMMA VERIFIED (per-relation bound; full excess open)" if ok else "CHECK FAILED",
          f"[{time.time() - t0:.0f}s]")
    print("=" * 88)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
