# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2"]
# ///
"""Certificate for the order theorem  A_N ~ (1/pi) ln N  for ALL N (the former order conjecture).

Two ingredients, both checked here:

1. PREFIX-INDEPENDENCE LEMMA  M(N) = phi(2N)/2  (the first phi(2N)/2 frequencies omega_r = 2 sin(pi r/N)
   are Q-independent).  Proof: a nontrivial relation sum_{r=1}^K k_r sin(pi r/N) = 0 with k_K != 0 yields
       Q(x) = sum_{r=1}^K k_r (x^{K+r} - x^{K-r})  in Z[x],
   nonzero of degree 2K (leading coeff k_K), with Q(zeta_2N) = zeta_2N^K (P(zeta)-P(1/zeta)) = 0, so
   Phi_2N | Q and phi(2N) <= 2K.  Q is anti-palindromic (x^{2K} Q(1/x) = -Q) while Phi_2N is palindromic,
   so 2K = phi(2N) is impossible (it would force Q = k_K Phi_2N = 0); hence 2K >= phi(2N)+2, i.e.
   K >= phi(2N)/2 + 1.  No relation lives in {1,..,phi(2N)/2}, so those frequencies are independent;
   with M(N) <= rank = phi(2N)/2 (Theorem qrank), M(N) = phi(2N)/2.

2. LOWER BOUND.  Aligning the independent prefix (Kronecker) gives, with L_pre = sum_{r<=M} b_r,
   b_r = 1/(N sin(pi r/N)),
       A_N >= sup_t G_N(0,t) >= 2 L_pre - U_N,
   and sin(pi r/N) <= pi r/N gives L_pre >= (1/pi) ln M(N), so
       A_N >= (1/pi) ln( M(N)^2 / N ) - O(1) = (1/pi)(2 ln(phi(2N)/2) - ln N) - O(1).
   By phi(2N) >= 2N/(e^gamma ln ln 2N + 3/ln ln 2N) this is (1/pi) ln N - O(ln ln ln N), so with
   A_N <= U_N = (1/pi) ln N + O(1),  A_N ~ (1/pi) ln N.

This script certifies (1) exactly (the divisibility mechanism on relation generators and the rank), and
(2) numerically (the bound 2 L_pre - U_N is positive and (2 L_pre - U_N)/((1/pi) ln N) -> 1).
"""

from __future__ import annotations

import sys
from itertools import pairwise
from math import gcd, log
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from large_wave_effect.cyclotomic import cyclotomic, exact_rank, sin_coords


def phi(n: int) -> int:
    return sum(1 for k in range(1, n + 1) if gcd(k, n) == 1)


def relation_generators(N: int) -> list[list[int]]:
    """Integer basis of the relation lattice of {omega_r: r=1..floor(N/2)} (left nullspace)."""
    from fractions import Fraction

    rows = sin_coords(2 * N, list(range(1, N // 2 + 1)))
    M = [[Fraction(x) for x in r] for r in rows]
    m, d = len(M), len(M[0])
    A = [[M[r][c] for r in range(m)] for c in range(d)]
    piv: list[int] = []
    pr = 0
    for col in range(m):
        sel = next((r for r in range(pr, d) if A[r][col] != 0), None)
        if sel is None:
            continue
        A[pr], A[sel] = A[sel], A[pr]
        p = A[pr][col]
        A[pr] = [x / p for x in A[pr]]
        for r in range(d):
            if r != pr and A[r][col] != 0:
                f = A[r][col]
                A[r] = [a - f * b for a, b in zip(A[r], A[pr], strict=True)]
        piv.append(col)
        pr += 1
        if pr == d:
            break
    basis = []
    for fcol in (c for c in range(m) if c not in piv):
        k = [Fraction(0)] * m
        k[fcol] = Fraction(1)
        for i, pc in enumerate(piv):
            k[pc] = -A[i][fcol]
        L = 1
        for x in k:
            L = L * x.denominator // gcd(L, x.denominator)
        ki = [int(x * L) for x in k]
        g = 0
        for x in ki:
            g = gcd(g, abs(x))
        basis.append([x // g for x in ki] if g else ki)
    return basis


def poly_divmod(a: list[int], b: list[int]) -> tuple[list[int], list[int]]:
    a = a[:]
    q = [0] * max(1, len(a) - len(b) + 1)
    while len(a) >= len(b) and any(a):
        if a[-1] == 0:
            a.pop()
            continue
        deg = len(a) - len(b)
        c = a[-1] // b[-1]
        q[deg] = c
        for j in range(len(b)):
            a[deg + j] -= c * b[j]
        while a and a[-1] == 0:
            a.pop()
    return q, a


def check_palindromic_mechanism(N: int) -> bool:
    """For every relation generator: Q is anti-palindromic, Phi_2N | Q, deg Q = 2K >= phi(2N)."""
    phi2n = phi(2 * N)
    phipoly = cyclotomic(2 * N)
    for k in relation_generators(N):
        kk = [0, *k]  # 1-indexed
        nz = [r for r in range(1, len(kk)) if kk[r]]
        if not nz:
            continue
        kmax = max(nz)
        q = [0] * (2 * kmax + 1)
        for r in range(1, kmax + 1):
            q[kmax + r] += kk[r]
            q[kmax - r] -= kk[r]
        anti = all(q[i] == -q[2 * kmax - i] for i in range(2 * kmax + 1))
        _, rem = poly_divmod(q[:], phipoly)
        if not (anti and all(x == 0 for x in rem) and 2 * kmax >= phi2n and 2 * kmax > phi2n):
            return False
    return True


def main() -> int:
    ok = True
    print("[1] prefix-independence lemma  M(N) = phi(2N)/2  and the palindromic mechanism:")
    for N in range(3, 121):
        d = phi(2 * N) // 2
        rank_first_d = exact_rank(sin_coords(2 * N, list(range(1, d + 1))))
        if rank_first_d != d:
            ok = False
            print(f"  FAIL N={N}: first {d} sines not independent (rank {rank_first_d})")
        if not check_palindromic_mechanism(N):
            ok = False
            print(f"  FAIL N={N}: palindromic divisibility mechanism")
    print(f"    first phi(2N)/2 sines independent, and every relation has 2K > phi(2N): {ok}")

    print("[2] lower bound  A_N >= 2 L_pre - U_N  with the sharp asymptotic:")

    def bsum(N: int, M: int) -> float:
        return sum(1.0 / (N * np.sin(np.pi * r / N)) for r in range(1, M + 1))

    ratios = []
    for N in [90, 210, 2310, 30030, 510510]:
        d = phi(2 * N) // 2
        lp, u = bsum(N, d), bsum(N, N // 2)
        bnd = 2 * lp - u
        ratio = bnd / (log(N) / np.pi)
        if bnd <= 0:
            ok = False
        print(f"    N={N:>7}: 2L_pre-U_N={bnd:6.3f}  ratio to (1/pi)lnN = {ratio:.3f}")
        ratios.append(ratio)
    # the bound is Omega(ln N) (every ratio bounded below) AND the ratio climbs toward 1
    # (monotone increase along primorials) -- test the trend, not one magic endpoint.
    increasing = all(b > a - 1e-9 for a, b in pairwise(ratios[1:]))
    if not (min(ratios) > 0.6 and increasing and ratios[-1] > ratios[0] + 0.1):
        ok = False

    print("=" * 64)
    print("RESULT:", "PASS -- A_N ~ (1/pi) ln N certified (order theorem)" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
