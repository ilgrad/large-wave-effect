# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2", "scipy>=1.11"]
# ///
"""Numerical investigation of the open problems (EVIDENCE, not proof -- stated honestly).

OP1. Is the sharp constant of A_N ~ c ln N equal to 1/pi for COMPOSITE N too?
     We compute the exact subtorus value A_N and the ceiling U_N and track the defect A_N/U_N. The
     order Theta(ln N) is already proven (prefix). Question: does A_N/U_N -> 1 (universal 1/pi)?
     Focus on the "deficit-1" family N = 2p (p odd prime): rank = p-1, a single rational relation
     (omega_p = 2 = 2 omega_1 only for p=3; for larger p one relation among p frequencies). If the
     defect -> 1 there and for highly composite N, that is strong evidence for a universal 1/pi.

OP2. The constant of B_N = sup_t ||e^{-it L_N}||_{inf->inf} = c_S sqrt(N).
     B_N = sup_t sum_m |K(m,t)|, K(m,t) = (1/N) sum_r e^{-i lambda_r t} e^{2 pi i r m/N}. A t-scan gives
     a LOWER bound on c_S (almost-periodic sup is not scan-reachable); we also probe the Talbot/Gauss
     times t = pi a / N where the kernel is a Gauss sum. We report the apparent constant honestly.

This script PROVES nothing new; it quantifies the open questions so the claims in the paper stay honest.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import minimize

RNG = np.random.default_rng(0)


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


def sin_coords(modulus: int, exps: list[int]) -> np.ndarray:
    cache: dict[int, list[int]] = {}
    phi = cyclotomic(modulus, cache)
    rows = []
    for e in exps:
        hi = power_mod(e % modulus, phi)
        lo = power_mod((modulus - e) % modulus, phi)
        rows.append([h - lv for h, lv in zip(hi, lo, strict=True)])
    return np.array(rows, dtype=float)


def max_over_subtorus(coeff: np.ndarray, c_mat: np.ndarray, starts: int) -> float:
    p = c_mat.shape[1]

    def negf(psi: np.ndarray) -> tuple[float, np.ndarray]:
        ph = c_mat @ psi
        return -float(coeff @ np.sin(ph)), -(c_mat.T @ (coeff * np.cos(ph)))

    best = 0.0
    for _ in range(starts):
        res = minimize(negf, RNG.uniform(0, 2 * np.pi, size=p), jac=True, method="L-BFGS-B")
        best = max(best, -float(res.fun))
    return best


def ring_AN(n: int, starts: int = 30) -> tuple[float, float]:
    r = np.arange(1, n)
    omega = 2.0 * np.sin(np.pi * r / n)
    c_mat = sin_coords(2 * n, list(r))
    u_n = float(np.sum(1.0 / np.sin(np.pi * r / n)) / (2 * n))
    best = 0.0
    for j in range(n // 2 + 1):
        coeff = (1.0 / n) * np.cos(2 * np.pi * r * j / n) / omega
        best = max(best, max_over_subtorus(coeff, c_mat, starts))
    return best, u_n


def is_prime(n: int) -> bool:
    return n >= 2 and all(n % d for d in range(2, int(n**0.5) + 1))


def schrodinger_B(n: int) -> float:
    lam = 4.0 * np.sin(np.pi * np.arange(n) / n) ** 2
    ts = np.concatenate([
        np.linspace(0.01, 6.0 * n, 8000),
        np.pi * np.arange(1, 4 * n) / n,        # Talbot / Gauss-sum times t = pi a / N
    ])
    return max(float(np.sum(np.abs(np.fft.ifft(np.exp(-1j * t * lam))))) for t in ts)


def main() -> int:
    print("=" * 72)
    print("Open problems: numerical EVIDENCE (not proof)")
    print("=" * 72)

    print("\nOP1. Composite defect A_N/U_N -- does the constant approach 1/pi for all N?")
    print(f"    {'N':>4} {'class':>10} {'A_N':>8} {'U_N':>8} {'A_N/U_N':>8} {'A_N/lnN':>8}")
    # deficit-1 family N=2p and some highly composite N
    fam2p = [6, 10, 14, 22, 26, 34, 38, 46]
    composite_other = [12, 24, 30, 36, 48]
    rows = []
    for n in sorted(set(fam2p + composite_other)):
        a, u = ring_AN(n)
        rows.append((n, a, u))
        half = n // 2
        cls = "2p" if (half >= 2 and is_prime(half) and n == 2 * half) else "composite"
        print(f"    {n:>4} {cls:>10} {a:>8.4f} {u:>8.4f} {a / u:>8.4f} {a / np.log(n):>8.4f}")
    # trend of the 2p family defect
    defects_2p = [(n, a / u) for (n, a, u) in rows if n in fam2p]
    trend_up = defects_2p[-1][1] > defects_2p[0][1]
    print(f"\n    deficit-1 family N=2p: A_N/U_N from {defects_2p[0][1]:.4f} (N={defects_2p[0][0]}) "
          f"to {defects_2p[-1][1]:.4f} (N={defects_2p[-1][0]})")
    print(f"    -> defect {'increasing toward 1 (evidence for universal 1/pi)' if trend_up else 'NOT increasing'}; "
          "this is EVIDENCE, not a proof. Sharp constant for composite N remains open.")

    print("\nOP2. Schrodinger constant c_S = lim B_N/sqrt(N)  (t-scan = LOWER bound on c_S)")
    print(f"    {'N':>5} {'B_N':>9} {'B_N/sqrt(N)':>12}")
    cs = []
    for n in (16, 32, 64, 128, 256, 512):
        b = schrodinger_B(n)
        cs.append(b / np.sqrt(n))
        print(f"    {n:>5} {b:>9.4f} {b / np.sqrt(n):>12.4f}")
    print(f"    -> B_N/sqrt(N) clusters near {np.median(cs):.3f} (scan lower bound); the exact constant")
    print("       is an almost-periodic / Gauss-sum supremum -- precise value remains open.")

    print("\n" + "=" * 72)
    print("RESULT: OPEN-PROBLEMS MAP COMPLETE")
    print("Status (kept in sync with the paper): OP1 is now PROVED -- A_N ~ (1/pi) ln N for all N")
    print("(Theorem order), and with U_N - A_N = O(ln ln ln N) the ratio A_N/U_N -> 1; this scan is")
    print("corroboration. OP2: liminf B_N/sqrt(N) >= c_0/sqrt2 is proved and the constant splits by")
    print("parity; the t > N/2 upper bound (which would give limsup = beta_odd) remains open, and any")
    print("t-scan here is a lower bound only.")
    print("=" * 72)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
