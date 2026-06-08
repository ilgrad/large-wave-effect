# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2"]
# ///
"""Certificate for the Dirichlet-segment rank theorem and saturation classification.

Segment frequencies mu_k = 2 sin(pi k / 2(N+1)) = 2 sin(pi k / N'), k = 1..N, with N' = 2(N+1) -- so
the mu_k are the first N distinct frequencies of the *ring* on N' = 2(N+1) nodes, and the excluded
middle index N+1 = N'/2 carries the rational frequency 2.  Hence, writing M = 4(N+1),

    rank_Q {mu_1,...,mu_N, 2} = phi(2 N')/2 = phi(M)/2          (ring Theorem qrank),

so  rank_Q {mu_k}  =  phi(M)/2 - [2 not in span{mu_k}].  We prove the deficient value by deciding the
bracket:

  * N+1 = 2^m:  the ring on N' = 2^{m+1} is fully Q-independent (Theorem two), so 2 = mu_{N'/2} is
    independent and 2 not in span{mu_k}; thus rank = phi(M)/2 - 1 = N.
  * N+1 composite, odd prime factor p:  the rotated-root relation with start a = N'/2,
        sum_{j=0}^{p-1} zeta_M^{(N+1) + (M/p) j} = 0,
    reduces to a ring-N' relation whose ONLY term hitting the middle index N+1 is j = 0 (because
    (M/p) j = M/2 would need 2 j = p, impossible for p odd), so the middle coefficient is +-1: this
    expresses 2 = mu_{N+1} as a Q-combination of mu_1..mu_N, i.e. 2 in span{mu_k}, and rank = phi(M)/2.
  * N+1 prime:  the mu_k are Q-independent (Proposition 2p), rank = N = phi(M)/2.

Together: rank_Q {mu_k} = min(N, phi(4(N+1))/2) for ALL N, full (= N) iff N+1 is prime or a power of two,
which by Bohr is exactly A_N^Dir = U_N^Dir.  This script certifies every step exactly.
"""

from __future__ import annotations

import sys
from math import gcd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from large_wave_effect.cyclotomic import exact_rank, sin_coords

MAXN = 100


def phi(n: int) -> int:
    return sum(1 for k in range(1, n + 1) if gcd(k, n) == 1)


def is_prime(n: int) -> bool:
    return n > 1 and all(n % d for d in range(2, int(n**0.5) + 1))


def is_two_power(n: int) -> bool:
    return n >= 2 and (n & (n - 1)) == 0


def odd_prime_factor(n: int) -> int | None:
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


def reduce_ring(e: int, np_: int) -> tuple[int, int]:
    e %= 2 * np_
    if e == 0 or e == np_:
        return 0, 0
    base, s = (e, 1) if e < np_ else (e - np_, -1)
    return (base if base <= np_ // 2 else np_ - base), s


def two_in_span_witness(n: int) -> bool:
    """For composite N+1: the a = N'/2 rotated relation has nonzero coeff on the middle index."""
    np_ = 2 * (n + 1)
    p = odd_prime_factor(n + 1)
    if p is None:
        return False  # power of two: no witness (2 is independent)
    q, mid = 2 * np_ // p, n + 1
    k = [0] * (np_ // 2 + 1)
    for j in range(p):
        r, s = reduce_ring(mid + q * j, np_)
        if r:
            k[r] += s
    # exact: is k a genuine ring-N' relation, and is the middle coefficient nonzero?
    rows = sin_coords(2 * np_, list(range(1, np_ // 2 + 1)))
    acc = [0] * len(rows[0])
    for r in range(1, np_ // 2 + 1):
        if k[r]:
            acc = [a + k[r] * c for a, c in zip(acc, rows[r - 1], strict=True)]
    return all(v == 0 for v in acc) and k[n + 1] != 0


def main() -> int:
    ok = True
    for n in range(2, MAXN + 1):
        m_mod = 4 * (n + 1)
        rank = exact_rank(sin_coords(m_mod, list(range(1, n + 1))))
        half_phi = phi(m_mod) // 2
        saturates = is_prime(n + 1) or is_two_power(n + 1)
        if rank != min(n, half_phi):
            ok = False
            print(f"  FAIL rank N={n}: {rank} != min({n},{half_phi})")
        if (rank == n) != saturates:
            ok = False
            print(f"  FAIL class N={n}")
        # the constructive witness for the deficient value (composite N+1)
        composite = not is_two_power(n + 1) and not is_prime(n + 1)
        if composite and not two_in_span_witness(n):
            ok = False
            print(f"  FAIL witness N={n}: a=N'/2 relation does not express 2")

    print(f"Dirichlet segment, N in [2,{MAXN}]:")
    print("  rank_Q {mu_k} = min(N, phi(4(N+1))/2)  (exact, PROVED via the 2-in-span construction)")
    print("  full rank (A^Dir=U^Dir)  <=>  N+1 prime or a power of two")
    print("=" * 64)
    print("RESULT:", "PASS -- segment rank theorem + classification certified" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
