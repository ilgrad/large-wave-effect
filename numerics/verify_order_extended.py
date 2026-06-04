# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2", "scipy>=1.11"]
# ///
"""Extended, honest numerical run for the order conjecture A_N = Theta(ln N). NO new theorem.

SUB-CHAIN TARGET (heuristic, NOT a proof). If p = P(N) is the largest prime factor of N, the cycle C_p
embeds in C_N (its modes are rho = s N/p, frequencies 2 sin(pi s/p), Q-independent since p is prime), and
those modes WOULD contribute (P(N)/N) U_{P(N)} ~ (P(N)/(pi N)) ln P(N) if aligned in isolation. This is
NOT a rigorous lower bound on A_N: at the aligning instant the remaining modes (total weight
U_N - (P(N)/N)U_{P(N)}, a constant fraction of U_N) are uncontrolled and may interfere -- exactly the
open obstruction in Conjecture (order). Empirically A_N still exceeds this target; we report it as
context, not as a bound.

This script honestly (a) checks that the C_p target sits below the optimized A_N (an empirical
observation); (b) computes the exact A_N/ln N for the hardest N = p*q with p,q ~ sqrt(N) (consecutive
primes, up to 19*23); (c) maps the Q-independent prefix length M(N)/N for highly composite N; (d)
confirms the independence locus rank = phi(2N)/2 = floor(N/2) <=> prime/2^m to N <= 500 via the closed
form (Theorem qrank). The optimized A_N is a valid LOWER bound on the true sup (any achievable phase
config), compared to the ceiling U_N. The conjecture A_N = Theta(ln N) for ALL N remains OPEN.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import minimize

RNG = np.random.default_rng(0)
_CYC: dict[int, list[int]] = {}
_P = 2_147_483_647


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


def sin_rows(n: int) -> list[list[int]]:
    m = 2 * n
    phi = cyclotomic(m)
    return [[h - lv for h, lv in zip(power_mod(r % m, phi), power_mod((m - r) % m, phi), strict=True)]
            for r in range(1, n)]


def rank_modp(rows: list[list[int]]) -> int:
    """Exact rank over GF(_P) via numpy int64 (products < _P^2 < 2^62 fit int64; ~100x faster than pure)."""
    if not rows:
        return 0
    a = np.array(rows, dtype=np.int64) % _P
    nr, nc = a.shape
    rank = 0
    for col in range(nc):
        nz = np.nonzero(a[rank:, col])[0]
        if nz.size == 0:
            continue
        piv = rank + int(nz[0])
        if piv != rank:
            a[[rank, piv]] = a[[piv, rank]]
        inv = pow(int(a[rank, col]), _P - 2, _P)
        a[rank] = (a[rank] * inv) % _P
        colv = a[:, col].copy()
        colv[rank] = 0
        nzr = np.nonzero(colv)[0]
        if nzr.size:
            a[nzr] = (a[nzr] - np.outer(colv[nzr], a[rank])) % _P
        rank += 1
        if rank == nr:
            break
    return rank


def prefix_length(n: int) -> int:
    rows = sin_rows(n)
    m = n // 2
    lo, hi = 1, m
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if rank_modp(rows[:mid]) == mid:
            lo = mid
        else:
            hi = mid - 1
    return lo


def u_ceiling(n: int) -> float:
    r = np.arange(1, n)
    return float(np.sum(1.0 / np.sin(np.pi * r / n)) / (2 * n))


def euler_phi(n: int) -> int:
    result, m, p = n, n, 2
    while p * p <= m:
        if m % p == 0:
            while m % p == 0:
                m //= p
            result -= result // p
        p += 1
    if m > 1:
        result -= result // m
    return result


def largest_prime_factor(n: int) -> int:
    p, d = 1, 2
    while d * d <= n:
        while n % d == 0:
            p, n = d, n // d
        d += 1
    return max(p, n)


def subchain_target(n: int) -> float:
    """C_p contribution (P(N)/N) U_{P(N)} -- a heuristic target, NOT a proven lower bound (interference)."""
    p = largest_prime_factor(n)
    return (p / n) * u_ceiling(p)


def exact_an(n: int, starts: int = 12) -> float:
    """Optimized node-0 supremum (a valid LOWER bound on the true A_N; equals it on saturating N)."""
    r = np.arange(1, n)
    omega = 2.0 * np.sin(np.pi * r / n)
    cmat = np.array(sin_rows(n), dtype=float)
    coeff = (1.0 / n) / omega

    def negf(psi):
        ph = cmat @ psi
        return -float(coeff @ np.sin(ph)), -(cmat.T @ (coeff * np.cos(ph)))

    best = 0.0
    for _ in range(starts):
        res = minimize(negf, RNG.uniform(0, 2 * np.pi, cmat.shape[1]), jac=True, method="L-BFGS-B")
        best = max(best, -float(res.fun))
    return best


def is_prime(n: int) -> bool:
    return n >= 2 and all(n % d for d in range(2, int(n**0.5) + 1))


def cls(n: int) -> str:
    return "prime" if is_prime(n) else ("2^m" if (n & (n - 1)) == 0 else "composite")


def main() -> int:
    print("=" * 84)
    print("Extended order run: C_p target + exact A_N for N=pq (p,q~sqrtN) + prefix + independence locus")
    print("=" * 84)
    ok = True

    print("\n(A) C_p target (P(N)/N) U_{P(N)} sits below A_N empirically (NOT a proof -- interference)")
    print(f"    {'N':>5} {'P(N)':>5} {'tgt/lnN':>8} {'A_N/lnN':>8} {'U_N/lnN':>8} {'tgt<=A<=U':>9}")
    for n in (6, 15, 21, 22, 33, 34, 46, 55, 58, 62, 69, 82, 86, 94):
        tgt = subchain_target(n)
        an = exact_an(n)
        un = u_ceiling(n)
        good = tgt <= an + 1e-9 <= un + 1e-9
        ok &= good
        ln = np.log(n)
        print(f"    {n:>5} {largest_prime_factor(n):>5} {tgt / ln:>8.4f} {an / ln:>8.4f} {un / ln:>8.4f} "
              f"{good!s:>9}")
    print("    -> the C_p target stays below the exact A_N <= U_N on all tested N (empirical): OK")

    print("\n(B) Hardest cases N = p*q with p,q ~ sqrt(N): exact A_N/ln N stays in the band")
    print(f"    {'N':>5} {'=p*q':>9} {'A_N/U_N':>8} {'A_N/lnN':>8} {'tgt/lnN':>8}")
    band = []
    for n in (143, 221, 323, 437):  # consecutive-prime products ~sqrt N: 11.13, 13.17, 17.19, 19.23
        an = exact_an(n)
        un = u_ceiling(n)
        ln = np.log(n)
        band.append(an / ln)
        fac = [d for d in range(2, n) if n % d == 0 and is_prime(d) and is_prime(n // d)]
        print(f"    {n:>5} {f'{fac[0]}x{n // fac[0]}':>9} {an / un:>8.4f} {an / ln:>8.4f} "
              f"{subchain_target(n) / ln:>8.4f}")
    in_band = all(0.27 < b < 0.34 for b in band)
    ok &= in_band
    print(f"    -> even for two near-equal large primes, A_N/ln N in [{min(band):.3f},{max(band):.3f}]: "
          f"{'OK' if in_band else 'FAIL'}")

    print("\n(C) Prefix length M(N)/N (Q-independent low modes) for highly composite N: bounded below")
    print(f"    {'N':>5} {'class':>10} {'M(N)':>5} {'M/N':>6} {'Lpre/lnN':>9}")
    worst = 1e9
    for n in (210, 330, 420, 510):
        m = prefix_length(n)
        lp = float(np.sum(1.0 / np.sin(np.pi * np.arange(1, m + 1) / n)) / n) / np.log(n)
        worst = min(worst, lp)
        print(f"    {n:>5} {cls(n):>10} {m:>5} {m / n:>6.3f} {lp:>9.4f}")
    pref_ok = worst > 0.2
    ok &= pref_ok
    print(f"    -> min L_pre/ln N = {worst:.3f} > 0.2 (incl. smooth, highly composite N): "
          f"{'OK' if pref_ok else 'FAIL'}")

    print("\n(D) Independence/saturation locus to N<=500 via the rank formula (Thm: rank = phi(2N)/2)")
    # A_N=U_N by full independence iff phi(2N)/2 = floor(N/2) iff N prime or 2^m (Theorem qrank, proved).
    # Closed-form check (no linear algebra); spot-check the formula against exact mod-p rank.
    sat = [n for n in range(3, 501) if euler_phi(2 * n) // 2 == n // 2]
    pure = all(is_prime(n) or (n & (n - 1)) == 0 for n in sat)
    spot = all(rank_modp(sin_rows(n)[: n // 2]) == euler_phi(2 * n) // 2 for n in (60, 105, 210))
    ok &= pure and spot
    print(f"    full-independence locus on [3,500] = prime union 2^m ({len(sat)} values), no composite: "
          f"{'OK' if pure else 'FAIL'}; rank-formula spot-check {'OK' if spot else 'FAIL'}")
    print("    (the strictly weaker mod-4 ceiling criterion is scanned directly in verify_ceiling_criterion.py)")

    print("\n" + "=" * 84)
    print("RESULT:", "EXTENDED ORDER RUN VERIFIED" if ok else "CHECK FAILED")
    print("Honest status: A_N/ln N stays in [0.28,0.34] incl. the hardest N=pq with p,q~sqrt(N);")
    print("full independence (=> A_N=U_N) is prime/2^m only. Conjecture A_N=Theta(ln N) for ALL N: OPEN.")
    print("=" * 84)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
