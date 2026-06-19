# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2", "scipy>=1.11"]
# ///
"""Verify the chain-split reduction of the excess lemma for 3|N (paper, remark after lem:triplet).

For 3|N the dependent band splits into clean support-3 triplets T (C-row = [1,1], phi_c=psi_a+psi_b) and
higher-support "comb" generators K. The triplets are VERTEX-DISJOINT (no prefix phase in two of them).
Hence for any lambda in [0,1), partitioning the prefix penalty and using lem:triplet
(g((1-lam)b_a,(1-lam)b_b,b_c) <= b_c/(2(1-lam))) with sup-superadditivity gives the RIGOROUS reduction
    A_N - L_pre  <=  (1/(1-lam)) * sum_{c in T} b_c/2  +  F_comb(lam),
    F_comb(lam) = max_psi[ sum_{r in K} b_r sin(phi_r) - lam sum_{s<=M} b_s (1 - sin psi_s) ].
The triplet sum is bounded (for square-free N=3p it is a Riemann sum of (1/2) int_{1/2}^{2/3} csc(pi x) dx
= ln3/(4pi) ~ 0.087; stays < 0.1 across the set). So the excess lemma for 3|N is ISOLATED to the comb
block F_comb = O(1) -- the remaining open core. This script verifies (i) vertex-disjointness, (ii) the
split is a valid upper bound on the true excess, (iii) the triplet bound, on 3|N up to N=315 (M<=240).
"""
from __future__ import annotations

import sys
import time
from math import log, pi
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, str(Path(__file__).resolve().parent))
from verify_excess_omega2 import integer_coord_matrix, optimize_AN, prime_factors, totient, weights


def split_data(n: int):
    m, big_m = n // 2, totient(2 * n) // 2
    c = integer_coord_matrix(n)
    _u, b, _omega = weights(n)
    triplets, combs, used, disjoint = [], [], set(), True
    for i, row in enumerate(c[big_m:m]):
        nz = np.where(row != 0)[0]
        if len(nz) == 2 and sorted(row[nz].tolist()) == [1, 1]:
            a, bb = int(nz[0]), int(nz[1])
            if a in used or bb in used:
                disjoint = False
            used |= {a, bb}
            triplets.append((big_m + i, a, bb))
        else:
            combs.append(big_m + i)
    return m, big_m, c, b, triplets, combs, disjoint


def f_comb(n, lam, combs, c, b, big_m, *, starts=30, seed=2):
    if not combs:
        return 0.0
    cf, bc, bs = c[combs].astype(float), b[combs], b[:big_m]
    rng = np.random.default_rng(seed)
    best = -np.inf
    for _ in range(starts):
        r = minimize(lambda p: (-(float(bc @ np.sin(cf @ p)) - lam * float((bs * (1 - np.sin(p))).sum())),
                                -(cf.T @ (bc * np.cos(cf @ p)) + lam * (bs * np.cos(p)))),
                     rng.uniform(0, 2 * pi, big_m), jac=True, method="L-BFGS-B")
        best = max(best, -r.fun)
    return best


def main() -> int:
    t0, ok = time.time(), True
    print("=" * 92)
    print("Chain-split reduction of the excess lemma for 3|N: excess <= (1/(1-lam)) sum_T b_c/2 + F_comb(lam)")
    print("=" * 92)
    print(f"  {'N':>4} {'#T':>4} {'#K':>4} {'disj':>5} {'cf':>4} {'sum b_c/2':>10} {'splitBound':>11} {'excess':>8} {'bnd>=exc':>9}")
    riem = log(3.0) / (4.0 * pi)  # (1/2) * int_{1/2}^{2/3} csc(pi x) dx = ln3/(4pi)
    for n in [15, 21, 33, 45, 63, 99, 105, 135, 165, 315]:
        if 3 not in prime_factors(n):
            continue
        _m, big_m, c, b, trip, combs, disjoint = split_data(n)
        # proved closed form (identity omega_{N/3+a}=omega_a+omega_{N/3-a}): each clean triplet is (a, N/3-a, N/3+a)
        cf = all((a + 1) + (bb + 1) == n // 3 and (cg + 1) == n // 3 + (a + 1) for cg, a, bb in trip)
        tsum = sum(b[cg] / 2 for cg, _, _ in trip)
        bound = min((1.0 / (1 - lam)) * tsum + f_comb(n, lam, combs, c, b, big_m, starts=20)
                    for lam in np.linspace(0.05, 0.9, 10))
        a_n, _sv, bb = optimize_AN(n, seed=0, starts=80)
        exc = a_n - float(bb[:big_m].sum())
        valid = bound >= exc - 1e-3
        ok &= disjoint and valid and (tsum < 0.1) and cf
        print(f"  {n:>4} {len(trip):>4} {len(combs):>4} {disjoint!s:>5} {cf!s:>4} {tsum:>10.4f} {bound:>11.4f} "
              f"{exc:>8.4f} {valid!s:>9}")
    print(f"\n  (square-free N=3p: sum_T b_c/2 -> ln3/(4pi) = {riem:.4f}; combs carry the open core F_comb)")
    print("=" * 92)
    print("RESULT:", "CHAIN-SPLIT VERIFIED (triplets disjoint; split valid; triplet part <0.1)" if ok
          else "CHECK FAILED", f"[{time.time()-t0:.0f}s]")
    print("=" * 92)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
