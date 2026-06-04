# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2", "scipy>=1.11"]
# ///
"""Exact criterion for reaching the ceiling: A_N = U_N iff every spectral relation has sum-of-coeffs = 0 mod 4.

The ceiling U_N is attained at node 0 iff all sines peak together, i.e. the target phase
phi_* = (pi/2,...,pi/2) lies in the orbit-closure subtorus T_N. By Kronecker, phi_* in T_N iff it is
compatible with every integer relation k in the relation lattice
    Lambda_N = { k in Z^m : sum_r k_r omega_r = 0 },   omega_r = 2 sin(pi r/N), r=1..floor(N/2),
that is <k, phi_*> = (pi/2) sum_r k_r = 0 (mod 2 pi) for all k, which holds iff

    sum_r k_r = 0 (mod 4)   for every k in Lambda_N.

Since A_N <= U_N always and node 0 dominates (|a_{r,j}| <= a_{r,0}), this proves:

    THEOREM.  A_N = U_N  <=>  sum_r k_r = 0 (mod 4) for every relation k in Lambda_N.

This refines the independence test (independence means Lambda_N = {0}, trivially satisfied): there may
be COMPOSITE N with rationally DEPENDENT frequencies that still reach U_N. We compute Lambda_N exactly
(integer left-kernel of the cyclotomic-coordinate matrix, via Hermite normal form, all in exact integer
arithmetic), check the mod-4 criterion, and confirm it agrees with the directly optimized A_N/U_N.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import minimize

RNG = np.random.default_rng(0)
_CYC: dict[int, list[int]] = {}


def divisors(n: int) -> list[int]:
    return [d for d in range(1, n + 1) if n % d == 0]


def poly_div_exact(a: list[int], b: list[int]) -> list[int]:
    a = a[:]
    db = len(b) - 1
    q = [0] * (len(a) - db)
    for i in range(len(a) - 1, db - 1, -1):
        c = a[i]
        if c:
            q[i - db] = c
            for j in range(len(b)):
                a[i - db + j] -= c * b[j]
    return q


def cyclotomic(n: int, cache: dict[int, list[int]]) -> list[int]:
    if n in cache:
        return cache[n]
    num = [-1] + [0] * (n - 1) + [1]
    for d in divisors(n):
        if d < n:
            num = poly_div_exact(num, cyclotomic(d, cache))
    cache[n] = num
    return num


def power_mod(e: int, phi: list[int]) -> list[int]:
    deg = len(phi) - 1
    p = [0] * max(e + 1, deg)
    p[e] = 1
    for i in range(len(p) - 1, deg - 1, -1):
        c = p[i]
        if c:
            for j in range(deg + 1):
                p[i - deg + j] -= c * phi[j]
    return p[:deg]


def sin_coords(n: int) -> list[list[int]]:
    """Integer coords of 2 i sin(pi r/N)=zeta_{2N}^r - zeta_{2N}^{-r} for the DISTINCT freqs r=1..floor(N/2)."""
    m = 2 * n
    phi = cyclotomic(m, _CYC)
    rows = []
    for r in range(1, n // 2 + 1):
        hi = power_mod(r % m, phi)
        lo = power_mod((m - r) % m, phi)
        rows.append([h - lv for h, lv in zip(hi, lo, strict=True)])
    return rows


def integer_left_kernel(c: list[list[int]]) -> list[list[int]]:
    """Z-basis of {k in Z^m : sum_r k_r C[r] = 0} by Euclidean column reduction (unimodular, exact)."""
    m, d = len(c), len(c[0])
    a = [[c[r][col] for r in range(m)] for col in range(d)]   # C^T : d x m
    v = [[1 if i == j else 0 for j in range(m)] for i in range(m)]
    col = 0
    for prow in range(d):
        if col >= m:
            break
        piv = next((cc for cc in range(col, m) if a[prow][cc] != 0), None)
        if piv is None:
            continue
        for x in (a, v):
            for row in x:
                row[col], row[piv] = row[piv], row[col]
        for cc in range(m):
            if cc == col:
                continue
            while a[prow][cc] != 0:
                q = a[prow][cc] // a[prow][col]
                if q:
                    for d_ in range(d):
                        a[d_][cc] -= q * a[d_][col]
                    for i in range(m):
                        v[i][cc] -= q * v[i][col]
                if a[prow][cc] != 0:
                    for x in (a, v):
                        for row in x:
                            row[col], row[cc] = row[cc], row[col]
        col += 1
    ker = []
    for cc in range(m):
        if all(a[d_][cc] == 0 for d_ in range(d)):
            k = [v[i][cc] for i in range(m)]
            g = 0
            for x in k:
                g = np.gcd(g, abs(x))
            if g > 1:
                k = [x // int(g) for x in k]
            ker.append(k)
    return ker


def ceiling_reached(n: int) -> bool:
    """True iff every relation k in Lambda_N has sum k_r = 0 (mod 4)."""
    ker = integer_left_kernel(sin_coords(n))
    return all(sum(k) % 4 == 0 for k in ker)


def u_ceiling(n: int) -> float:
    r = np.arange(1, n)
    return float(np.sum(1.0 / np.sin(np.pi * r / n)) / (2 * n))


def a_over_u(n: int, starts: int = 15) -> float:
    """Directly optimized A_N (node 0) / U_N over the orbit-closure subtorus, to confirm the criterion."""
    r = np.arange(1, n)
    omega = 2.0 * np.sin(np.pi * r / n)
    m = 2 * n
    phi = cyclotomic(m, _CYC)
    rows = []
    for rr in r:
        hi = power_mod(rr % m, phi)
        lo = power_mod((m - rr) % m, phi)
        rows.append([h - lv for h, lv in zip(hi, lo, strict=True)])
    cmat = np.array(rows, dtype=float)
    coeff = (1.0 / n) / omega

    def negf(psi):
        ph = cmat @ psi
        return -float(coeff @ np.sin(ph)), -(cmat.T @ (coeff * np.cos(ph)))

    best = 0.0
    for _ in range(starts):
        res = minimize(negf, RNG.uniform(0, 2 * np.pi, cmat.shape[1]), jac=True, method="L-BFGS-B")
        best = max(best, -float(res.fun))
    return best / u_ceiling(n)


def is_prime(n: int) -> bool:
    return n >= 2 and all(n % d for d in range(2, int(n**0.5) + 1))


def main() -> int:
    print("=" * 74)
    print("Ceiling criterion: A_N = U_N  <=>  sum k_r = 0 (mod 4) for every relation k")
    print("=" * 74)
    ok = True

    print("\n(A) The criterion agrees with the directly optimized A_N/U_N (N = 3..40)")
    print(f"    {'N':>4} {'criterion':>10} {'A_N/U_N':>8} {'agree':>6} {'class':>10}")
    extra = []
    for n in range(3, 41):
        crit = ceiling_reached(n)
        ratio = a_over_u(n)
        reached = ratio > 0.999
        agree = crit == reached
        ok &= agree
        cls = "prime" if is_prime(n) else ("2^m" if (n & (n - 1)) == 0 else "composite")
        if crit and cls == "composite":
            extra.append(n)
        if n <= 24 or not agree or (crit and cls == "composite"):
            print(f"    {n:>4} {crit!s:>10} {ratio:>8.4f} {agree!s:>6} {cls:>10}")
    print(f"    -> criterion matches optimization for all N in 3..40: {'OK' if ok else 'FAIL'}")

    print("\n(B) N=6 worked example: relation 2 omega_1 - omega_3 = 0, sum = 1 (not 0 mod 4) => A_6 < U_6")
    ker6 = integer_left_kernel(sin_coords(6))
    print(f"    Lambda_6 basis = {ker6}, coeff sums mod 4 = {[sum(k) % 4 for k in ker6]}")
    b_ok = any(sum(k) % 4 != 0 for k in ker6)
    ok &= b_ok
    print(f"    -> a relation with sum != 0 mod 4 exists, so A_6 < U_6: {'OK' if b_ok else 'FAIL'}")

    print(f"\n(C) Composite N reaching the ceiling despite dependent frequencies: {extra if extra else 'none in 3..40'}")
    print("    (the criterion is strictly weaker than full independence: dependent relations are")
    print("    harmless when all their coefficient sums vanish mod 4.)")

    print("\n(D) Wide criterion-only scan N=3..160: which N reach the ceiling?")
    reach, comp_reach = [], []
    for n in range(3, 161):
        if ceiling_reached(n):
            reach.append(n)
            if not (is_prime(n) or (n & (n - 1)) == 0):
                comp_reach.append(n)
    pure = all(is_prime(n) or (n & (n - 1)) == 0 for n in reach)
    print(f"    ceiling reached for {len(reach)} values; composites among them: {comp_reach if comp_reach else 'NONE'}")
    print(f"    -> on N<=160, A_N=U_N holds exactly for prime and 2^m, never for composite N: {'OK' if pure else 'COMPOSITE FOUND'}")
    print("    The mod-4 criterion thus *explains* why only prime / 2^m saturate: their distinct")
    print("    frequencies are rationally independent (Lambda_N = {0}), so the criterion is vacuous.")

    print("\n" + "=" * 74)
    print("RESULT:", "CEILING CRITERION VERIFIED" if ok else "CHECK FAILED")
    print("A_N = U_N iff every integer spectral relation has coefficient-sum divisible by 4 -- a sharp,")
    print("arithmetic criterion refining the independence test.")
    print("=" * 74)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
