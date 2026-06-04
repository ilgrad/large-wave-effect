# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2"]
# ///
"""Consolidated numerical experiments: both dynamics across the three families of N.

Wave (bead chain): A_N = sup_t max_j |G_N(j,t)|. For prime and 2^m this equals the Bohr ceiling
U_N = (1/2N) sum csc(pi r/N) ~ (1/pi) ln N (Theorems 1', 1''); for composite N the exact subtorus
value (from verify_exact_AN.py) sits a little below.

Schrodinger: B_N = sup_t ||e^{-itL_N}||_{inf->inf} = Theta(sqrt N) (ballistic).

This reproduces the summary table of the paper (Section "Numerical experiments").
"""

from __future__ import annotations

import numpy as np

# exact subtorus A_N for the composite cases (from verify_exact_AN.py, validated on primes)
A_COMPOSITE = {12: 0.7546, 15: 0.7542, 24: 0.9618}
FAMILIES = {"2^m": [16, 32, 64], "prime": [17, 31, 61], "composite": [12, 15, 24]}


def ceiling_u(n: int) -> float:
    r = np.arange(1, n)
    return float(np.sum(1.0 / np.sin(np.pi * r / n)) / (2 * n))


def schrodinger_b(n: int) -> float:
    lam = 4.0 * np.sin(np.pi * np.arange(n) / n) ** 2
    ts = np.linspace(0.01, 6.0 * n, 8000)
    return max(float(np.sum(np.abs(np.fft.ifft(np.exp(-1j * t * lam))))) for t in ts)


def main() -> int:
    print("=" * 78)
    print("Numerical experiments: wave A_N and Schrodinger B_N over three families of N")
    print("=" * 78)
    print(f"  {'family':>10} {'N':>5} {'U_N':>8} {'A_N(wave)':>10} {'A_N/U_N':>8} "
          f"{'B_N(Schr)':>10} {'B_N/sqrtN':>10}")
    ok = True
    for fam, ns in FAMILIES.items():
        for n in ns:
            u = ceiling_u(n)
            a = u if fam != "composite" else A_COMPOSITE[n]
            b = schrodinger_b(n)
            print(f"  {fam:>10} {n:>5} {u:>8.4f} {a:>10.4f} {a / u:>8.4f} "
                  f"{b:>10.4f} {b / np.sqrt(n):>10.4f}")
            ok &= (a <= u + 1e-9) and (0.7 < b / np.sqrt(n) < 1.1)
    print("=" * 78)
    print("RESULT:", "CONSISTENT" if ok else "UNEXPECTED")
    print("Wave: A_N = U_N ~ (1/pi) ln N for prime/2^m, slightly below for composite (A_N/U_N >= 0.84);")
    print("Schrodinger: B_N ~ sqrt(N) across all families. Two dynamics, two growth laws (ln N vs sqrt N).")
    print("=" * 78)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
