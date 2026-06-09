# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2"]
# ///
"""Order theorem on ring products: A = Theta(ln N) on C_N box H for any fixed connected graph H.

The Cartesian product C_N box H has Laplacian eigenvalues lambda_r + mu_h, where lambda_r = 2-2cos(2 pi r/N)
(ring) and mu_h are the eigenvalues of L_H.  Only the gapless branch h=0 (mu_0 = 0, H-constant eigenvector
w_0 = 1/sqrt|H|) has frequencies omega -> 0, so the velocity ceiling U ~ (1/(pi|H|)) ln N comes from it
(the gapped branches add O(1)) -- condition (A) of the order criterion.  That branch's frequencies are
EXACTLY the ring frequencies 2 sin(pi r/N), whose first phi(2N)/2 are Q-independent (the palindromic lemma)
and carry (1-o(1)) U > U/2 -- condition (B).  Hence A = Theta(ln N).

This script confirms, for H = K_2 (prism) and H = P_3 (path), that the alignable-prefix bound
2 L_pre - U is positive and grows like c ln N (=> A = Omega(ln N); with A <= U = O(ln N), A = Theta(ln N)).
"""

from __future__ import annotations

import sys
from math import gcd
from pathlib import Path

import numpy as np
from numpy import pi, sin

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from large_wave_effect.cyclotomic import exact_rank, sin_coords


def totient(n: int) -> int:
    return sum(1 for k in range(1, n + 1) if gcd(k, n) == 1)


def h_spectrum(name: str) -> tuple[np.ndarray, np.ndarray]:
    """(eigenvalues mu_h, first-coordinate eigenvector components w_h(a0=0)) for L_H."""
    if name == "K2":
        lh = np.array([[1.0, -1.0], [-1.0, 1.0]])
    elif name == "P3":
        lh = np.array([[1.0, -1.0, 0.0], [-1.0, 2.0, -1.0], [0.0, -1.0, 1.0]])
    else:
        raise ValueError(name)
    mu, vecs = np.linalg.eigh(lh)
    return mu, vecs[0, :]


def ceiling(n: int, h: str) -> float:
    lam = 2 - 2 * np.cos(2 * pi * np.arange(n) / n)
    mu, w = h_spectrum(h)
    u = 0.0
    for r in range(n):
        for k in range(len(mu)):
            if r == 0 and k == 0:
                continue
            u += (1.0 / n) * w[k] ** 2 / np.sqrt(lam[r] + mu[k])
    return u


def prefix_weight(n: int, h: str) -> float:
    # h=0 ring frequencies, distinct r=1..M (M = phi(2N)/2), each with degeneracy 2 (modes r and N-r)
    big_m = totient(2 * n) // 2
    _, w = h_spectrum(h)
    return sum(2 * (1.0 / n) * w[0] ** 2 / (2 * sin(pi * r / n)) for r in range(1, big_m + 1))


def main() -> int:
    ok = True
    for h, absh in [("K2", 2), ("P3", 3)]:
        print(f"C_N box {h} (|H|={absh}):  prefix bound 2 L_pre - U  (>0 and growing => A=Theta(ln N))")
        print(f"  {'N':>5} {'U':>8} {'L_pre':>8} {'2Lpre-U':>9}")
        vals = []
        for n in [256, 512, 1024, 2048]:
            u, lp = ceiling(n, h), prefix_weight(n, h)
            vals.append(2 * lp - u)
            print(f"  {n:>5} {u:>8.4f} {lp:>8.4f} {2 * lp - u:>9.4f}")
        # positive and increasing -> Omega(ln N)
        if not (min(vals) > 0 and vals[-1] > vals[0] + 0.1):
            ok = False
    # the embedded ring prefix is Q-independent (palindromic lemma), same as the bare ring
    m = totient(2 * 512) // 2
    rank = exact_rank(sin_coords(2 * 512, list(range(1, m + 1))))
    print(f"\nembedded ring prefix Q-independent (N=512): rank {rank} = M={m}: {rank == m}")
    ok &= rank == m
    print("=" * 60)
    print("RESULT:", "PASS -- A=Theta(ln N) on C_N box H (criterion + palindromic lemma)" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
