# /// script
# requires-python = ">=3.11"
# dependencies = []
# ///
"""Exact Q-rank of the ring frequencies {2 sin(pi r/N): r=1..floor(N/2)} for any N.

Lemma 1 says these are linearly independent over Q exactly when the rank equals
floor(N/2). We compute the rank EXACTLY (rational arithmetic, no floats) by expressing
2i sin(pi r/N) = eta^r - eta^{-r} (eta = zeta_{2N}) in the integral basis
{eta^0,...,eta^{phi(2N)-1}} of Q(zeta_{2N}), reducing modulo the cyclotomic polynomial
Phi_{2N}. The rank of the resulting integer matrix is the dimension of the Q-span,
i.e. the number of independent frequencies (= dimension of the Kronecker torus).

  full rank  (= floor(N/2))  <=>  Bohr gives A_N = U_N ~ (1/pi) ln N exactly.
  deficient                  <=>  rational relations among the sines (sub-torus): A_N
                                  may fall below U_N -- "what breaks" on composite N.

Pure standard library: integer polynomial algebra + Fraction Gaussian elimination.
"""

from __future__ import annotations

from fractions import Fraction

Poly = list[int]


def divisors(n: int) -> list[int]:
    return [d for d in range(1, n + 1) if n % d == 0]


def poly_div_exact(a: Poly, b: Poly) -> Poly:
    """Exact quotient a/b for monic b (b[-1] == 1); remainder assumed zero."""
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


def cyclotomic(n: int, cache: dict[int, Poly]) -> Poly:
    """Phi_n(x) as integer coefficients (low->high), via x^n-1 = prod_{d|n} Phi_d."""
    if n in cache:
        return cache[n]
    num = [-1] + [0] * (n - 1) + [1]  # x^n - 1
    for d in divisors(n):
        if d < n:
            num = poly_div_exact(num, cyclotomic(d, cache))
    cache[n] = num
    return num


def power_mod(e: int, phi: Poly) -> list[int]:
    """Coordinates of x^e mod Phi (monic) in {1, x, ..., x^{deg-1}}, integer vector."""
    deg = len(phi) - 1
    p = [0] * max(e + 1, deg)
    p[e] = 1
    for i in range(len(p) - 1, deg - 1, -1):
        c = p[i]
        if c:
            for j in range(deg + 1):
                p[i - deg + j] -= c * phi[j]  # phi monic: cancels p[i]
    return p[:deg]


def sine_matrix(n: int) -> list[list[int]]:
    """Rows = coords of 2i sin(pi r/N) = eta^r - eta^{-r}, eta = zeta_{2N}, r=1..floor(N/2)."""
    cache: dict[int, Poly] = {}
    phi = cyclotomic(2 * n, cache)
    rows = []
    for r in range(1, n // 2 + 1):
        hi = power_mod(r, phi)
        lo = power_mod((2 * n - r) % (2 * n), phi)
        rows.append([h - lv for h, lv in zip(hi, lo, strict=True)])
    return rows


def exact_rank(matrix: list[list[int]]) -> int:
    m = [[Fraction(x) for x in row] for row in matrix]
    rows, cols = len(m), len(m[0])
    pr, rank = 0, 0
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


def report(label: str, ns: list[int]) -> bool:
    print(f"\n[{label}]")
    print(f"  {'N':>4} {'floor(N/2)':>11} {'Q-rank':>7} {'status':>10}")
    all_full = True
    for n in ns:
        half = n // 2
        rank = exact_rank(sine_matrix(n))
        full = rank == half
        all_full &= full
        tag = "full" if full else f"DEFICIENT -{half - rank}"
        print(f"  {n:>4} {half:>11} {rank:>7} {tag:>10}")
    return all_full


def main() -> int:
    print("=" * 64)
    print("Exact Q-rank of ring frequencies {2 sin(pi r/N)} (Lemma 1 map)")
    print("=" * 64)

    primes = [3, 5, 7, 11, 13, 17, 19, 23]
    powers2 = [4, 8, 16, 32]
    composite = [6, 9, 10, 12, 15, 18, 20, 21, 24, 25, 30, 36]

    full_primes = report("odd primes (Theorem 1')", primes)
    full_p2 = report("powers of two (Theorem 1'')", powers2)
    _ = report("composite / prime-power (map)", composite)

    print("\n" + "=" * 64)
    ok = full_primes and full_p2
    print("RESULT:", "PRIMES & 2^m FULL RANK (theorems confirmed)" if ok else "UNEXPECTED")
    print("Composite N show rank deficiency -> rational relations among sin(pi r/N)")
    print("(e.g. N=6: sin(pi/2) = 2 sin(pi/6)) -> sub-torus, A_N may be < U_N. Open case.")
    print("=" * 64)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
