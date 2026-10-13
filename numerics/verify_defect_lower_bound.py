# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2", "scipy>=1.11"]
# ///
"""Explicit lower bounds for the composite deficit U_N - A_N (Proposition deficit-lower), against the optimizer.

Lemma rel-defect: an integer relation k (sum_r k_r omega_r = 0) with odd coefficient sum forces, at every point
of the orbit closure, sum_r b_r (1 - sin phi_r) >= 1 / (2 sum_r k_r^2 / b_r).  The classes x mod q = 2N/p of an
odd prime p | N give such relations with disjoint supports, so the node-0 deficit is at least the sum of these
terms (>= (s - 1) sin(pi/2p) / (8ps), s = N/p); for 3 | N the vertex-disjoint triplets of Lemma cleancomb with
Lemma triplet give the sharper (1/2) sum_{1 <= a < N/6} b_{N/3+a}, which tends to ln(3)/(4 pi).  Every node
j != 0, N/2 has deficit >= 1/6.  For N = 2p the single relation spans the lattice and
(1 - pi^3/(48p)) pi^2/(8S) <= U_N - A_N <= pi^2/(8S), S = 4p csc(pi/2p), so U_N - A_N ~ (pi^3/16)/N^2.

The bounds are proved in the text and certified in exact arithmetic by exact/pari/defect_lower_bound.gp; this
script checks that they sit below the optimizer's defect U_N - A_N (an upper estimate of the true defect,
since the multistart optimizer can only underestimate A_N), that the N = 2p sandwich holds, and the
node-j bound.
"""

from __future__ import annotations

import math

import numpy as np
from verify_exact_AN import max_over_subtorus, sin_coords


def odd_primes(n: int) -> list[int]:
    out, m, p = [], n, 2
    while p * p <= m:
        if m % p == 0:
            out.append(p)
            while m % p == 0:
                m //= p
        p += 1
    if m > 1:
        out.append(m)
    return [q for q in out if q % 2]


def saturating(n: int) -> bool:
    return not odd_primes(n) or odd_primes(n) == [n]


def weight(n: int, r: int) -> float:
    return 1 / (2 * n) if 2 * r == n else 1 / (n * math.sin(math.pi * r / n))


def orbit_relations(n: int, p: int) -> list[dict[int, int]]:
    s = n // p
    q = 2 * s
    seen: set[int] = set()
    rels = []
    for x in range(1, q):
        if x == s or x in seen:
            continue
        seen |= {x, (-x) % q, (x + s) % q, (s - x) % q}
        k: dict[int, int] = {}
        for j in range(p):
            e = (x + q * j) % (2 * n)
            sg, r = (1, e) if e < n else (-1, e - n)
            r = min(r, n - r)
            k[r] = k.get(r, 0) + sg
        rels.append({r: v for r, v in k.items() if v})
    return rels


def orbit_bound(n: int, p: int) -> float:
    return sum(
        1 / (2 * sum(v * v / weight(n, r) for r, v in k.items())) for k in orbit_relations(n, p)
    )


def triplet_bound(n: int) -> float:
    return 0.5 * sum(weight(n, n // 3 + a) for a in range(1, math.ceil(n / 6)))


def node_defect(n: int, j: int, starts: int = 60) -> tuple[float, float]:
    r = np.arange(1, n)
    omega = 2 * np.sin(np.pi * r / n)
    coeff = (1 / n) * np.cos(2 * np.pi * r * j / n) / omega
    u = float(np.sum(1 / np.sin(np.pi * r / n)) / (2 * n))
    return u, u - max_over_subtorus(coeff, sin_coords(2 * n, r.tolist()), starts)


def main() -> int:
    ok = True
    print(
        "[1] certified lower bounds vs the optimizer defect at node 0 (bound <= defect required):"
    )
    print(f"    {'N':>3} {'U-A0 (opt)':>11} {'orbits':>9} {'triplets':>9}")
    for n in range(6, 61):
        if saturating(n):
            continue
        _, d0 = node_defect(n, 0)
        ob = max(orbit_bound(n, p) for p in odd_primes(n))
        tb = triplet_bound(n) if n % 3 == 0 else 0.0
        good = max(ob, tb) <= d0 + 1e-9
        ok &= good
        if n % 3 == 0 or n in (10, 14, 20, 25, 35, 49):
            print(
                f"    {n:>3} {d0:>11.5f} {ob:>9.5f} {tb:>9.5f}{'' if good else '  <-- VIOLATION'}"
            )

    print(
        "\n[2] every node j != 0, N/2 has deficit >= 1/6 (sum_r b_r |cos(2 pi r j/N)| <= U_N - 1/6):"
    )
    worst = 1.0
    for n in range(6, 201):
        if saturating(n):
            continue
        rr = np.arange(1, n // 2 + 1)
        b = np.array([weight(n, int(r)) for r in rr])
        for j in range(1, n // 2 + (n % 2)):
            if 2 * j == n:
                continue
            worst = min(worst, float(np.sum(b * (1 - np.abs(np.cos(2 * np.pi * rr * j / n))))))
    ok &= worst >= 1 / 6
    print(
        f"    min over 6 <= N <= 200, j != 0, N/2 of the node-j deficit bound: {worst:.5f} (>= 1/6)"
    )

    print(
        "\n[3] N = 2p: (1 - pi^3/48p) pi^2/8S <= U_N - A_N <= pi^2/8S, and N^2 (U_N - A_N) -> pi^3/16:"
    )
    for p in [3, 5, 7, 11, 13, 17, 19, 23, 29, 31]:
        n = 2 * p
        s_val = 4 * p / math.sin(math.pi / (2 * p))
        lo, hi = (1 - math.pi**3 / (48 * p)) * math.pi**2 / (8 * s_val), math.pi**2 / (8 * s_val)
        d = max(node_defect(n, 0)[1], 0.0)
        dd = min(d, node_defect(n, p)[1])
        good = lo - 1e-9 <= dd <= hi + 1e-9
        ok &= good
        print(
            f"    N={n:>3}: {lo:.6f} <= {dd:.6f} <= {hi:.6f}   N^2 * defect = {n * n * dd:.4f}{'' if good else '  <--'}"
        )
    print(f"    pi^3/16 = {math.pi**3 / 16:.4f}")

    print("\n[4] the triplet bound tends to ln(3)/(4 pi):")
    for n in [99, 999, 9999, 99999]:
        print(f"    N={n:>6}: {triplet_bound(n):.6f}")
    lim = math.log(3) / (4 * math.pi)
    ok &= abs(triplet_bound(99999) - lim) < 1e-4
    print(f"    ln(3)/(4 pi) = {lim:.6f}")

    print("=" * 72)
    print(
        "RESULT:", "PASS -- the deficit bounds hold below every optimizer value" if ok else "FAIL"
    )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
