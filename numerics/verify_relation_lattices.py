# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2"]
# ///
"""Explicit relation lattices Lambda_N and the mod-4 ceiling criterion, tabulated for small N.

For the distinct frequencies omega_r = 2 sin(pi r/N), r=1..floor(N/2), the relation lattice
Lambda_N = {k in Z^m : sum_r k_r omega_r = 0} is the integer left-kernel of the cyclotomic-coordinate
matrix (2 i sin(pi r/N) = zeta_{2N}^r - zeta_{2N}^{-r}). Theorem (ceiling criterion):
A_N = U_N iff every k in Lambda_N has sum_r k_r divisible by 4. This script prints, for each N, a basis
of Lambda_N and the coefficient-sums mod 4, making the obstruction explicit. (Pure integer arithmetic;
the kernel is exact via Euclidean/Hermite column reduction.)
"""

from __future__ import annotations

import numpy as np

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


def cyclotomic(n: int) -> list[int]:
    if n in _CYC:
        return _CYC[n]
    num = [-1] + [0] * (n - 1) + [1]
    for d in divisors(n):
        if d < n:
            num = poly_div_exact(num, cyclotomic(d))
    _CYC[n] = num
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
    m = 2 * n
    phi = cyclotomic(m)
    rows = []
    for r in range(1, n // 2 + 1):
        hi = power_mod(r % m, phi)
        lo = power_mod((m - r) % m, phi)
        rows.append([h - lv for h, lv in zip(hi, lo, strict=True)])
    return rows


def integer_left_kernel(c: list[list[int]]) -> list[list[int]]:
    m, d = len(c), len(c[0])
    a = [[c[r][col] for r in range(m)] for col in range(d)]
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


def is_prime(n: int) -> bool:
    return n >= 2 and all(n % d for d in range(2, int(n**0.5) + 1))


def main() -> int:
    print("=" * 78)
    print("Relation lattices Lambda_N and the mod-4 ceiling criterion")
    print("=" * 78)
    print(f"\n  {'N':>4} {'class':>10} {'#freq':>6} {'rk Lambda':>10} {'sums mod 4':>14} {'A_N=U_N?':>9}")
    ok = True
    for n in (5, 6, 8, 9, 10, 12, 14, 15, 16, 18, 20, 21, 24):
        ker = integer_left_kernel(sin_coords(n))
        sums = [sum(k) % 4 for k in ker]
        sat = all(s == 0 for s in sums)
        cls = "prime" if is_prime(n) else ("2^m" if (n & (n - 1)) == 0 else "composite")
        # cross-check: prime / 2^m must saturate (empty lattice), others must not (on this range)
        expect = is_prime(n) or (n & (n - 1)) == 0
        ok &= sat == expect
        sums_str = "(none)" if not sums else ",".join(map(str, sums))
        print(f"  {n:>4} {cls:>10} {n // 2:>6} {len(ker):>10} {sums_str:>14} {sat!s:>9}")

    print("\n  Worked examples (a generating relation, read off the lattice basis vector):")

    def relation_str(k: list[int]) -> str:
        terms = [f"{'+' if c > 0 else '-'}{abs(c) if abs(c) != 1 else ''}w{r + 1}"
                 for r, c in enumerate(k) if c]
        s = " ".join(terms).lstrip("+")
        return f"{s} = 0   (coeff-sum {sum(k)} = {sum(k) % 4} mod 4)"

    for n in (6, 12, 15):
        k = integer_left_kernel(sin_coords(n))[0]
        print(f"    N={n:>2}: {relation_str(k)}  ->  A_{n} < U_{n}")

    print("\n" + "=" * 78)
    print("RESULT:", "RELATION LATTICES TABULATED" if ok else "CHECK FAILED")
    print("Prime / 2^m have empty lattice (saturate); composite N carry a relation with sum != 0 mod 4.")
    print("=" * 78)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
