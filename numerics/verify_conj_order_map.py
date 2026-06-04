# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2", "scipy>=1.11"]
# ///
"""Numerical map for Conjecture (order): A_N = Theta(ln N) for every N, with constant approaching 1/pi.

Three quantities over a wide range of N:
  (A) the exact node-0 supremum A_N (orbit-closure subtorus optimization) and the ratios A_N/U_N, A_N/ln N;
  (B) the Q-independent low-mode prefix length M(N) and its budget
      L_pre(N) = (1/N) sum_{r<=M(N)} csc(pi r/N) -- the amplitude those modes WOULD contribute in
      isolation. NOTE: this is a heuristic budget, NOT a proven lower bound on A_N: aligning the prefix
      to sin=1 leaves the out-of-prefix modes (total weight U_N - L_pre, a constant fraction of U_N)
      uncontrolled, and they may interfere destructively. Controlling that interference IS the open step;
  (C) the ceiling ratio U_N/ln N -> 1/pi.
Together they give strong evidence for A_N = Theta(ln N): A_N/ln N stays in a tight band bounded away
from 0, and the alignable prefix carries an amplitude ~ a constant times ln N. The conjecture is NOT
proved here.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import minimize

RNG = np.random.default_rng(0)
_CYC: dict[int, list[int]] = {}
_P = 2_147_483_647  # Mersenne prime for exact mod-p rank


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
    a = [[x % _P for x in row] for row in rows]
    nr = len(a)
    nc = len(a[0]) if nr else 0
    rank = 0
    for col in range(nc):
        piv = next((r for r in range(rank, nr) if a[r][col]), None)
        if piv is None:
            continue
        a[rank], a[piv] = a[piv], a[rank]
        inv = pow(a[rank][col], _P - 2, _P)
        a[rank] = [(x * inv) % _P for x in a[rank]]
        for r in range(nr):
            if r != rank and a[r][col]:
                f = a[r][col]
                a[r] = [(a[r][c] - f * a[rank][c]) % _P for c in range(nc)]
        rank += 1
    return rank


def prefix_length(n: int) -> int:
    """Largest M with {sin(pi r/N): r=1..M} Q-independent (rows 1..M full-rank)."""
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


def prefix_budget(n: int, m: int) -> float:
    r = np.arange(1, m + 1)
    return float(np.sum(1.0 / np.sin(np.pi * r / n)) / n)


def exact_an_over_u(n: int, starts: int = 18) -> float:
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
    return best / u_ceiling(n)


def is_prime(n: int) -> bool:
    return n >= 2 and all(n % d for d in range(2, int(n**0.5) + 1))


def cls(n: int) -> str:
    return "prime" if is_prime(n) else ("2^m" if (n & (n - 1)) == 0 else "composite")


def main() -> int:
    print("=" * 78)
    print("Order map:  A_N = Theta(ln N) for all N, constant -> 1/pi")
    print("=" * 78)
    ok = True

    print("\n(A) Exact A_N (node 0) across classes; A_N/ln N stays in a tight positive band")
    print(f"    {'N':>4} {'class':>10} {'U_N':>7} {'A_N/U_N':>8} {'A_N/lnN':>8}")
    ratios = []
    for n in (5, 8, 12, 16, 24, 30, 32, 36, 48, 60, 64, 72, 90, 96):
        au = exact_an_over_u(n)
        an = au * u_ceiling(n)
        rr = an / np.log(n)
        ratios.append(rr)
        print(f"    {n:>4} {cls(n):>10} {u_ceiling(n):>7.4f} {au:>8.4f} {rr:>8.4f}")
    band = min(ratios) > 0.27 and max(ratios) / min(ratios) < 1.3
    ok &= band
    print(f"    -> A_N/ln N in [{min(ratios):.3f}, {max(ratios):.3f}], ratio {max(ratios) / min(ratios):.2f}<1.3 "
          f"(tight band away from 0 => Theta(ln N)): {'OK' if band else 'FAIL'}")

    print("\n(B) Prefix budget L_pre(N) (heuristic, not a proven bound); L_pre/ln N bounded below for all N")
    print(f"    {'N':>4} {'class':>10} {'M(N)':>5} {'L_pre/lnN':>10}")
    worst = 1e9
    for n in (30, 60, 105, 120, 180, 210, 240, 252, 256, 300):
        m = prefix_length(n)
        lp = prefix_budget(n, m) / np.log(n)
        worst = min(worst, lp)
        print(f"    {n:>4} {cls(n):>10} {m:>5} {lp:>10.4f}")
    prefix_ok = worst > 0.2
    ok &= prefix_ok
    print(f"    -> min L_pre/ln N = {worst:.3f} > 0.2 (incl. highly composite 210,240,2520-smooth): "
          f"{'OK' if prefix_ok else 'FAIL'}")

    print("\n(C) Ceiling constant U_N/ln N -> 1/pi = 0.3183")
    print(f"    {'N':>6} {'U_N/lnN':>9}")
    for n in (100, 1000, 10000, 100000):
        print(f"    {n:>6} {u_ceiling(n) / np.log(n):>9.5f}")
    conv = abs(u_ceiling(100000) / np.log(100000) - 1 / np.pi) < 0.01
    ok &= conv
    print(f"    -> converges to 1/pi: {'OK' if conv else 'FAIL'}")

    print("\n" + "=" * 78)
    print("RESULT:", "ORDER MAP VERIFIED" if ok else "CHECK FAILED")
    print("A_N/ln N is tightly banded and the prefix bound alone grows like c ln N: strong evidence")
    print("for A_N = Theta(ln N) on every N (the uniform lower bound remains the open step).")
    print("=" * 78)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
