# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2"]
# ///
"""Certificate for the classification theorem A_N = U_N  <=>  N prime or N = 2^m.

The saturation half (prime / 2^m) is Theorems prime, two (Lambda_N = {0}). This script certifies the
other half: every composite N that is not a power of two has A_N < U_N, via the explicit obstruction

    k = the reduction of   sum_{j=0}^{p-1} sin(pi (1 + q j) / N) = 0,   p = odd prime factor, q = 2N/p,

which is a genuine relation in Lambda_N (imaginary part of the vanishing root-of-unity sum
sum_j zeta_{2N}^{1+qj} = zeta_{2N} * (1 + zeta_p + ... + zeta_p^{p-1}) = 0), and whose coefficient sum

    S = sum_r k_r = #{points in (0,N)} - #{points in (N,2N)}

is ODD: none of the p points 1 + q j (mod 2N) lands on {0, N} (since q | 2N, q > 1, and 2s does not
divide p s - 1 for N = p s, s >= 2), so S = p - 0 = p (mod 2) is odd. An odd S gives both
    sum_r k_r       = S            != 0 (mod 4)    [node 0]
    sum_r k_r + 2 T = S + 2 T      != 0 (mod 4)    [antipodal node N/2, even N],
so neither admissible node saturates and A_N < U_N (Theorem ceiling).

This script checks, exactly, every step of that argument for all composite non-2^m N in range.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from large_wave_effect.cyclotomic import sin_coords

MAXN = 600


def odd_prime_factor(n: int) -> int | None:
    """Smallest odd prime factor of n, or None if n is a power of two."""
    while n % 2 == 0:
        n //= 2
    if n == 1:
        return None
    d = 3
    while d * d <= n:
        if n % d == 0:
            return d
        d += 2
    return n


def reduce_index(e: int, n: int) -> tuple[int, int]:
    """sin(pi e/N) = sigma * sin(pi r/N) with r in [1, floor(N/2)], sigma in {-1,0,1}."""
    e %= 2 * n
    if e == 0 or e == n:
        return 0, 0
    base, sigma = (e, 1) if e < n else (e - n, -1)
    return (base if base <= n // 2 else n - base), sigma


def construction(n: int, p: int) -> list[int]:
    """Coefficient vector k (index r = 1..floor(N/2)) of the rotated-root obstruction, a = 1."""
    q = 2 * n // p
    m = n // 2
    k = [0] * (m + 1)
    for j in range(p):
        r, sigma = reduce_index(1 + q * j, n)
        if r:
            k[r] += sigma
    return k[1 : m + 1]


def is_exact_relation(n: int, k: list[int]) -> bool:
    """Exact (integer) check that sum_r k_r * 2 sin(pi r/N) = 0 in Q(zeta_2N)."""
    rows = sin_coords(2 * n, list(range(1, n // 2 + 1)))
    acc = [0] * len(rows[0])
    for r, kr in enumerate(k, start=1):
        if kr:
            row = rows[r - 1]
            acc = [a + kr * c for a, c in zip(acc, row, strict=True)]
    return all(v == 0 for v in acc)


def is_prime(n: int) -> bool:
    return n > 1 and all(n % d for d in range(2, int(n**0.5) + 1))


def is_two_power(n: int) -> bool:
    return n >= 2 and (n & (n - 1)) == 0


def main() -> int:
    ok = True
    checked = 0
    for n in range(4, MAXN + 1):
        if is_prime(n) or is_two_power(n):
            continue
        p = odd_prime_factor(n)
        assert p is not None
        k = construction(n, p)
        s = sum(k)
        t = sum((r + 1) * k[r] for r in range(len(k)))
        rel = is_exact_relation(n, k)
        node0 = s % 4 != 0
        antipodal = (s + 2 * t) % 4 != 0
        s_odd = s % 2 == 1
        if not (rel and s_odd and node0 and antipodal and s != 0):
            ok = False
            print(f"  FAIL N={n}: rel={rel} S={s} node0={node0} antipodal={antipodal}")
        checked += 1

    print(f"composite non-2^m N in [4,{MAXN}]: {checked} cases")
    print("  the rotated-root k is an EXACT relation in Lambda_N,")
    print("  S = sum k_r is ODD, so both the node-0 and antipodal mod-4 branches are violated")
    print("  => A_N < U_N for every such N.")
    print("=" * 60)
    print("RESULT:", "PASS -- classification certified on this range" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
