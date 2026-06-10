# /// script
# requires-python = ">=3.11"
# dependencies = ["mpmath>=1.3", "numpy>=2"]
# ///
"""Degeneracy law of the next-nearest-neighbour ring (Proposition nnn).

NNN ring (g = 1/2): omega_a = 2 sin(pi a/N) sqrt(2 + cos(2 pi a/N)), a = 1..floor(N/2).  With
x = cos(2 pi a/N), omega^2 = f(x) = 2(1-x)(2+x) and f(x_a) - f(x_b) = -2 (x_a - x_b)(1 + x_a + x_b); cos
is injective on (0, pi], so two DISTINCT frequencies collide iff

    cos(2 pi a/N) + cos(2 pi b/N) = -1   <=>   2 + zeta^a + zeta^-a + zeta^b + zeta^-b = 0,

a vanishing length-6 positive sum of N-th roots of unity.  By the Conway-Jones classification of rational
cosine relations the only solution with a != b is {2 pi a/N, 2 pi b/N} = {pi/2, pi}, i.e. (a,b) =
(N/4, N/2) -- so a collision exists iff 4 | N, and then omega_{N/4} = omega_{N/2} = 2.

This script (1) verifies the collision law EXACTLY -- integer polynomial arithmetic modulo Phi_N -- for
all N <= 200, and (2) runs a 100-digit PSLQ search showing no FURTHER rational relation among the NNN
frequencies for 5 <= N <= 40 (dropping the classified duplicate when 4 | N).  On this evidence the NNN
family satisfies condition (B) of the order criterion throughout.
"""

from __future__ import annotations

import sys
from pathlib import Path

import mpmath as mp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from large_wave_effect.cyclotomic import cyclotomic


def poly_rem(a: list[int], b: list[int]) -> list[int]:
    a = a[:]
    while len(a) >= len(b) and any(a):
        if a[-1] == 0:
            a.pop()
            continue
        c = a[-1] // b[-1]
        d = len(a) - len(b)
        for j in range(len(b)):
            a[d + j] -= c * b[j]
        while a and a[-1] == 0:
            a.pop()
    return a


def collision_pairs(n: int, phi_n: list[int], num: int = -1, den: int = 1) -> list[tuple[int, int]]:
    """Exact: pairs a < b <= N/2 with cos(2 pi a/N) + cos(2 pi b/N) = num/den, i.e.
    den (zeta^a + zeta^-a + zeta^b + zeta^-b) - 2 num = 0 in Q(zeta_N).  Default num/den = -1 (g=1/2)."""
    out = []
    for a in range(1, n // 2 + 1):
        for b in range(a + 1, n // 2 + 1):
            p = [0] * (max(n - a, n - b) + 1)
            p[a] += den
            p[(n - a) % n] += den
            p[b] += den
            p[(n - b) % n] += den
            p[0] -= 2 * num
            if not poly_rem(p, phi_n):
                out.append((a, b))
    return out


def nnn_frequencies(n: int) -> list[mp.mpf]:
    return [2 * mp.sin(mp.pi * r / n) * mp.sqrt(2 + mp.cos(2 * mp.pi * r / n)) for r in range(1, n // 2 + 1)]


def main() -> int:
    ok = True

    print("[1] exact collision scan (integer arithmetic mod Phi_N), N <= 200:")
    bad = []
    for n in range(3, 201):
        expect = [(n // 4, n // 2)] if n % 4 == 0 else []
        if collision_pairs(n, cyclotomic(n)) != expect:
            bad.append(n)
    ok &= not bad
    print("    collisions = {(N/4, N/2) : 4 | N} exactly:", "CONFIRMED" if not bad else f"FAIL at {bad[:5]}")

    print("[2] PSLQ (100 digits, maxcoeff 1000): no further relations, 5 <= N <= 40:")
    mp.mp.dps = 100
    extra = []
    for n in range(5, 41):
        om = nnn_frequencies(n)
        if n % 4 == 0:
            del om[n // 2 - 1]  # drop the classified duplicate omega_{N/2}
        if mp.pslq(om, maxcoeff=1000, maxsteps=10**6) is not None:
            extra.append(n)
    ok &= not extra
    print("    full rank apart from the classified collision:", "CONFIRMED" if not extra else f"FAIL at {extra}")

    print("[3] general coupling g (collision <=> cos(2pi a/N)+cos(2pi b/N) = -1/(2g)), exact, N <= 120:")
    # g = 1: the Conway-Jones solutions of cos A + cos B = -1/2 are {pi/3,pi}, {2pi/5,4pi/5}, {pi/2,2pi/3}
    bad_g1 = []
    for n in range(3, 121):
        expect = set()
        if n % 6 == 0:
            expect.add((n // 6, n // 2))
        if n % 5 == 0:
            expect.add((n // 5, 2 * n // 5))
        if n % 12 == 0:
            expect.add(tuple(sorted((n // 4, n // 3))))
        if set(collision_pairs(n, cyclotomic(n), num=-1, den=2)) != expect:
            bad_g1.append(n)
    ok &= not bad_g1
    print("    g=1: families (N/6,N/2), (N/5,2N/5), (N/4,N/3) exactly:",
          "CONFIRMED" if not bad_g1 else f"FAIL at {bad_g1[:5]}")
    # g = 1/4: cos A + cos B = -2 forces A = B = pi -- no distinct-pair collisions at all
    bad_g14 = [n for n in range(3, 121) if collision_pairs(n, cyclotomic(n), num=-2, den=1)]
    ok &= not bad_g14
    print("    g=1/4: no collisions for any N:", "CONFIRMED" if not bad_g14 else f"FAIL at {bad_g14[:5]}")

    print("=" * 64)
    print("RESULT:", "PASS -- NNN degeneracy law (a,b)=(N/4,N/2) iff 4|N; no other relations found"
          if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
