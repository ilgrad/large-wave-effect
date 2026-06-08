# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2"]
# ///
"""Certificate for the Dirichlet-segment classification A_N^Dir = U_N^Dir  <=>  N+1 prime or 2^m.

The simple spectrum mu_k = 2 sin(pi k / 2(N+1)), k = 1..N, makes the segment saturate (A=U, by Bohr)
iff the {mu_k} are Q-independent, i.e. iff rank_Q = N.  Writing 2 i sin(2 pi k/M) = zeta_M^k -
zeta_M^{-k} with M = 4(N+1), the rank is bounded by the minus-eigenspace dimension:

    rank_Q {mu_k : k=1..N}  <=  dim_Q V^-  =  phi(M)/2  =  phi(4(N+1))/2,

and an elementary count gives phi(4(N+1))/2 < N exactly when N+1 is composite and not a power of two
(if N+1 = 2^a m, m odd, then phi/2 = 2^a phi(m) >= N=2^a m -1 iff 2^a (m - phi(m)) <= 1 iff m=1 or
(a=0, m prime)).  So composite non-2^m N+1 forces dependence (A^Dir < U^Dir); for N+1 prime the mu_k are
independent by Proposition 2p, and for N+1 = 2^m by the Theorem-two argument.  This script checks the
exact ranks and the dimension bound over a range.
"""

from __future__ import annotations

import sys
from math import gcd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from large_wave_effect.cyclotomic import exact_rank, sin_coords

MAXN = 60


def phi(n: int) -> int:
    return sum(1 for k in range(1, n + 1) if gcd(k, n) == 1)


def is_prime(n: int) -> bool:
    return n > 1 and all(n % d for d in range(2, int(n**0.5) + 1))


def is_two_power(n: int) -> bool:
    return n >= 2 and (n & (n - 1)) == 0


def main() -> int:
    ok = True
    for n in range(2, MAXN + 1):
        m_mod = 4 * (n + 1)
        rank = exact_rank(sin_coords(m_mod, list(range(1, n + 1))))
        half_phi = phi(m_mod) // 2
        full = rank == n
        saturates = is_prime(n + 1) or is_two_power(n + 1)
        # checks: exact rank formula, the V^- dimension bound, and the classification
        if rank != min(n, half_phi):
            ok = False
            print(f"  FAIL rank N={n}: {rank} != min({n},{half_phi})")
        if rank > half_phi:
            ok = False
            print(f"  FAIL bound N={n}: rank {rank} > phi(M)/2 {half_phi}")
        if full != saturates:
            ok = False
            print(f"  FAIL class N={n}: full={full} but N+1 prime/2^m={saturates}")

    print(f"Dirichlet segment, N in [2,{MAXN}]:")
    print("  rank_Q {mu_k} = min(N, phi(4(N+1))/2), and rank <= phi(4(N+1))/2 (lives in V^-)")
    print("  full rank (A^Dir=U^Dir)  <=>  N+1 prime or a power of two")
    print("=" * 60)
    print("RESULT:", "PASS -- segment classification certified" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
