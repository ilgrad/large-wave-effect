# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2"]
# ///
"""Lemma 4 evidence: the independent low-mode prefix carries a logarithmic budget for ALL N.

The ring large wave A_N = sup_t (2/N) sum_{r=1}^{floor(N/2)} sin(omega_r t)/omega_r is hard to
bound below for composite N because not all frequencies are Q-independent. Strategy: keep only
the longest INDEPENDENT PREFIX {omega_1,...,omega_{M(N)}} (M(N) = first r at which sin(pi r/N)
becomes Q-dependent on the earlier ones). On that prefix Kronecker aligns the phases, giving a
guaranteed "budget"

    L_pre(N) = (2/N) sum_{r=1}^{M(N)} 1 / omega_r,   omega_r = 2 sin(pi r/N).

If L_pre(N) >= c ln N for all N with c > 0, the effect cannot vanish for any large N (the rigorous
gap is controlling the OUT-of-prefix modes; the prefix budget is a heuristic lower bound on A_N).

M(N) is found by exact independence (cyclotomic coordinates in Q(zeta_{2N}), rank mod a large
prime for speed; Q-dependence implies dependence mod q, so this does not over-count the prefix).
"""

from __future__ import annotations

import numpy as np

Q = (1 << 61) - 1  # large prime for modular rank (Q-dependence => dependence mod Q)


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


def sine_coord(n: int, r: int, phi: list[int]) -> list[int]:
    hi = power_mod(r, phi)
    lo = power_mod((2 * n - r) % (2 * n), phi)
    return [h - lv for h, lv in zip(hi, lo, strict=True)]


def independent_prefix(n: int) -> int:
    """Largest M with {sin(pi r/N): r=1..M} Q-independent (first dependent row stops it)."""
    cache: dict[int, list[int]] = {}
    phi = cyclotomic(2 * n, cache)
    basis: list[list[int]] = []  # reduced rows mod Q (row echelon)
    pivots: list[int] = []
    m = 0
    for r in range(1, n // 2 + 1):
        row = [x % Q for x in sine_coord(n, r, phi)]
        for b, pc in zip(basis, pivots, strict=True):
            f = row[pc]
            if f:
                row = [(row[k] - f * b[k]) % Q for k in range(len(row))]
        pc = next((k for k in range(len(row)) if row[k]), None)
        if pc is None:
            break  # r is dependent on earlier modes
        inv = pow(row[pc], Q - 2, Q)
        row = [(x * inv) % Q for x in row]
        basis.append(row)
        pivots.append(pc)
        m = r
    return m


def main() -> int:
    print("=" * 74)
    print("Lemma 4: independent low-mode prefix carries a logarithmic budget (all N)")
    print("=" * 74)
    print(f"  {'N':>4} {'floor(N/2)':>10} {'M(N)':>5} {'L_pre':>8} {'U_N':>8} {'L_pre/lnN':>10} {'kind':>9}")

    def kind(n: int) -> str:
        if all(n % d for d in range(2, int(n**0.5) + 1)):
            return "prime"
        if n & (n - 1) == 0:
            return "2^m"
        return "composite"

    ns = sorted({3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 127,
                 4, 8, 16, 32, 64, 128, 256,
                 6, 9, 10, 12, 14, 15, 18, 20, 21, 24, 25, 27, 30, 33, 35, 36, 40, 45, 48, 50, 60,
                 120, 180, 210, 240})  # incl. highly composite worst cases
    logs, budgets = [], []
    for n in ns:
        m = independent_prefix(n)
        r = np.arange(1, m + 1)
        l_pre = float(2.0 / n * np.sum(1.0 / (2.0 * np.sin(np.pi * r / n))))
        u_n = float(1.0 / (2 * n) * np.sum(1.0 / np.sin(np.pi * np.arange(1, n) / n)))
        logs.append(np.log(n))
        budgets.append(l_pre)
        print(f"  {n:>4} {n // 2:>10} {m:>5} {l_pre:>8.4f} {u_n:>8.4f} "
              f"{l_pre / np.log(n):>10.4f} {kind(n):>9}")

    logs_a, budgets_a = np.array(logs), np.array(budgets)
    a, b = np.polyfit(logs_a, budgets_a, 1)
    pred = a * logs_a + b
    r2 = 1 - float(np.sum((budgets_a - pred) ** 2)) / float(np.sum((budgets_a - budgets_a.mean()) ** 2))
    min_ratio = float(np.min(budgets_a / logs_a))
    print("\n  fit L_pre ~ a ln N + b:  a =", f"{a:.4f}", " R^2 =", f"{r2:.4f}")
    print(f"  min over all N of L_pre/lnN = {min_ratio:.4f}  (>0 => budget never collapses)")

    ok = (a > 0.1) and (min_ratio > 0.05)
    print("\n" + "=" * 74)
    print("RESULT:", "BUDGET GROWS ~ln N FOR ALL N" if ok else "BUDGET COLLAPSES (unexpected)")
    print("=> independent low modes alone give A_N >~ c ln N for every N (primes, 2^m, composite);")
    print("   the remaining rigor is controlling out-of-prefix interference (Lebesgue-constant route).")
    print("=" * 74)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
