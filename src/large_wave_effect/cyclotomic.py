"""Exact arithmetic in cyclotomic fields: independence ranks of lattice frequencies.

All integer/rational, no floating point in the certificates. Used to certify the linear
independence over Q of the chain frequencies (Theorems for prime / 2^m and the rank map).
"""

from __future__ import annotations

from fractions import Fraction

Poly = list[int]


def divisors(n: int) -> list[int]:
    return [d for d in range(1, n + 1) if n % d == 0]


def _poly_div_exact(a: Poly, b: Poly) -> Poly:
    """Exact quotient a / b for monic b (remainder assumed zero)."""
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


def cyclotomic(n: int, cache: dict[int, Poly] | None = None) -> Poly:
    """Cyclotomic polynomial Phi_n as integer coefficients (low -> high)."""
    cache = {} if cache is None else cache
    if n in cache:
        return cache[n]
    num: Poly = [-1] + [0] * (n - 1) + [1]  # x^n - 1
    for d in divisors(n):
        if d < n:
            num = _poly_div_exact(num, cyclotomic(d, cache))
    cache[n] = num
    return num


def power_mod(e: int, phi: Poly) -> list[int]:
    """Coordinates of x^e mod Phi (monic) in the basis {1, x, ..., x^{deg-1}}."""
    deg = len(phi) - 1
    p = [0] * max(e + 1, deg)
    p[e] = 1
    for i in range(len(p) - 1, deg - 1, -1):
        c = p[i]
        if c:
            for j in range(deg + 1):
                p[i - deg + j] -= c * phi[j]
    return p[:deg]


def sin_coords(modulus: int, exps: list[int]) -> list[list[int]]:
    """Integer coordinates of 2i*sin(pi*e/(modulus/2)) = eta^e - eta^{-e}, eta = zeta_modulus."""
    phi = cyclotomic(modulus)
    rows = []
    for e in exps:
        hi = power_mod(e % modulus, phi)
        lo = power_mod((modulus - e) % modulus, phi)
        rows.append([h - lo_ for h, lo_ in zip(hi, lo, strict=True)])
    return rows


def exact_rank(matrix: list[list[int]]) -> int:
    """Exact rank over Q of an integer matrix (Fraction Gaussian elimination)."""
    m = [[Fraction(x) for x in row] for row in matrix]
    rows, cols = len(m), len(m[0]) if m else 0
    pr = rank = 0
    for col in range(cols):
        sel = next((r for r in range(pr, rows) if m[r][col] != 0), None)
        if sel is None:
            continue
        m[pr], m[sel] = m[sel], m[pr]
        piv = m[pr][col]
        for r in range(rows):
            if r != pr and m[r][col]:
                f = m[r][col] / piv
                m[r] = [a - f * b for a, b in zip(m[r], m[pr], strict=True)]
        pr += 1
        rank += 1
        if pr == rows:
            break
    return rank


def ring_freq_rank(n: int) -> int:
    """Q-rank of {2 sin(pi r/N): r=1..N-1} (= dimension of the orbit-closure subtorus)."""
    return exact_rank(sin_coords(2 * n, list(range(1, n))))
