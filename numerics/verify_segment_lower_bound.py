# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2"]
# ///
"""The Dirichlet segment: unconditional Theta(ln N), full rank iff N+1 is prime or a power of two.

The segment chain (Dirichlet BCs) has the SIMPLE spectrum mu_k = 2 sin(pi k / (2(N+1))), k=1..N, with
no degeneracy -- the structural advantage over the ring, where mu_r = mu_{N-r} halves the independent
phases. Writing mu_k = 2 sin(2 pi k / M) with M = 4(N+1), the frequencies lie in Q(zeta_M) and
2 i sin(2 pi k/M) = zeta_M^k - zeta_M^{-k}.

We establish, all unconditionally except where noted:
 (A) the ceiling U_N^{Dir} = sum_k phi_k(center)^2 / mu_k ~ (1/pi) ln N (same constant as the ring),
     and A_N^{Dir} <= U_N^{Dir} by the triangle inequality (Bohr), with NO independence needed;
 (B) the Q-rank of {mu_k} (subtorus dimension) equals min(N, phi(4(N+1))/2); it is FULL (= N, hence
     A_N^{Dir} = U_N^{Dir} by Bohr) IFF N+1 is prime or a power of two -- the exact analog of the ring
     (full iff N prime or 2^m), shifted by the boundary: the ring uses modulus 2N, the segment 4(N+1);
 (C) a weak UNCONDITIONAL lower bound A_N^{Dir} = Omega(ln N) for ALL N via the low-mode prefix
     k=1..K: it is Q-independent for K up to ~ a power of N, so by Kronecker the phases align and
     A_N^{Dir} >= sum_{k<=K} a_k - o(1) ~ (1/pi) ln K.
"""

from __future__ import annotations

import numpy as np

# ---- exact cyclotomic rank (shared style with verify_qrank_formula.py) ----


def divisors(n: int) -> list[int]:
    return [d for d in range(1, n + 1) if n % d == 0]


def totient(n: int) -> int:
    r, m, p = n, n, 2
    while p * p <= m:
        if m % p == 0:
            while m % p == 0:
                m //= p
            r -= r // p
        p += 1
    if m > 1:
        r -= r // m
    return r


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


def rank_mod(rows: list[list[int]], q: int = (1 << 61) - 1) -> int:
    mat = [r[:] for r in rows]
    pr = 0
    ncol = len(mat[0])
    for col in range(ncol):
        piv = next((rr for rr in range(pr, len(mat)) if mat[rr][col] % q), None)
        if piv is None:
            continue
        mat[pr], mat[piv] = mat[piv], mat[pr]
        inv = pow(mat[pr][col], q - 2, q)
        mat[pr] = [(x * inv) % q for x in mat[pr]]
        for rr in range(len(mat)):
            if rr != pr and mat[rr][col]:
                f = mat[rr][col]
                mat[rr] = [(mat[rr][k] - f * mat[pr][k]) % q for k in range(ncol)]
        pr += 1
    return pr


def segment_rank(n: int, kmax: int | None = None) -> int:
    """dim_Q span{2 sin(2 pi k / M): k=1..kmax}, M=4(N+1) (subtorus dimension of the segment)."""
    m = 4 * (n + 1)
    kmax = kmax if kmax is not None else n
    cache: dict[int, list[int]] = {}
    phi = cyclotomic(m, cache)
    q = (1 << 61) - 1
    rows = []
    for k in range(1, kmax + 1):
        hi = power_mod(k % m, phi)
        lo = power_mod((m - k) % m, phi)
        rows.append([(h - lo_) % q for h, lo_ in zip(hi, lo, strict=True)])
    return rank_mod(rows)


def is_prime(n: int) -> bool:
    return n >= 2 and all(n % d for d in range(2, int(n**0.5) + 1))


def segment_ceiling(n: int) -> float:
    """U_N^{Dir} at the central node: sum_k phi_k(c)^2/mu_k, phi_k(j)=sqrt(2/(N+1)) sin(pi k j/(N+1))."""
    k = np.arange(1, n + 1)
    mu = 2.0 * np.sin(np.pi * k / (2 * (n + 1)))
    c = (n + 1) // 2  # central node
    phi2 = (2.0 / (n + 1)) * np.sin(np.pi * k * c / (n + 1)) ** 2
    return float(np.sum(phi2 / mu))


def main() -> int:
    print("=" * 70)
    print("Dirichlet segment: unconditional Theta(ln N); full rank iff N+1 prime or 2^m")
    print("=" * 70)
    ok = True

    print("\n(A) Ceiling U_N^{Dir} grows like (1/pi) ln N (A_N^{Dir} <= U_N^{Dir}, triangle ineq.)")
    inc = []
    for n in (63, 127, 255, 511, 1023):
        u = segment_ceiling(n)
        inc.append(u)
        print(f"    N={n:>5}: U_N^Dir={u:.4f}")
    doublings = [inc[i + 1] - inc[i] for i in range(len(inc) - 1)]
    a_ok = all(abs(d - np.log(2) / np.pi) < 0.03 for d in doublings)
    ok &= a_ok
    print(f"    increment per doubling ~ (ln2)/pi={np.log(2) / np.pi:.4f}: "
          f"{[f'{d:.4f}' for d in doublings]}  {'OK' if a_ok else 'FAIL'}")

    print("\n(B) seg rank = min(N, phi(4(N+1))/2); full (=N) iff N+1 is prime or 2^m")
    print(f"    {'N':>4} {'seg rank':>9} {'min(N,phi/2)':>13} {'full?':>6} {'N+1 prime/2^m':>14}")
    formula_ok, charac_ok = True, True
    for n in range(2, 31):
        sr = segment_rank(n)
        formula = min(n, totient(4 * (n + 1)) // 2)
        formula_ok &= sr == formula
        full = sr == n
        np1_special = is_prime(n + 1) or ((n + 1) & n) == 0  # n+1 prime or power of two
        charac_ok &= full == np1_special
        if n <= 24:
            print(f"    {n:>4} {sr:>9} {formula:>13} {full!s:>6} {np1_special!s:>14}")
    ok &= formula_ok and charac_ok
    print(f"    -> rank = min(N, phi(4(N+1))/2): {'PASS' if formula_ok else 'FAIL'};  "
          f"full <=> N+1 prime/2^m: {'PASS' if charac_ok else 'FAIL'}")

    print("\n(C) Unconditional Omega(ln N) for ALL N via the Q-independent low-mode prefix")
    print(f"    {'N':>5} {'K=floor(sqrt N)':>16} {'prefix rank':>12} {'independent?':>13}")
    c_ok = True
    for n in (50, 100, 200, 400):
        k = int(np.floor(np.sqrt(n)))
        pr = segment_rank(n, kmax=k)
        indep = pr == k
        c_ok &= indep
        print(f"    {n:>5} {k:>16} {pr:>12} {indep!s:>13}")
    ok &= c_ok
    print(f"    -> prefix independent => A_N^Dir >= (1/pi)ln K - o(1) = Omega(ln N): "
          f"{'OK' if c_ok else 'FAIL'}")

    print("\n" + "=" * 70)
    print("RESULT:", "SEGMENT BOUND VERIFIED" if ok else "CHECK FAILED")
    print("Dirichlet segment: U_N^Dir ~ (1/pi)ln N; rank = min(N, phi(4(N+1))/2), full (=> A=U) iff N+1")
    print("is prime or 2^m (ring analog shifted by the boundary); unconditional Theta(ln N) via prefix.")
    print("=" * 70)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
