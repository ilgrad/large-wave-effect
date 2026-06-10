# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2"]
# ///
"""One-sided dynamical self-trapping threshold with the exact constant 4 (Proposition threshold).

For the focusing ring DNLS (eq:dnls)  i u' = -L_N u + gamma |u|^2 u  with single-site data
u_j(0) = sqrt(P) delta_{j j0}, conservation of H = -<L_N u, u> + (gamma/2) sum |u|^4 and of the mass P
give, for ALL t (three lines: <L u0,u0> = 2P, L_N >= 0):

    sum_j |u_j(t)|^4  >=  P^2 - 4P/gamma ,
    max_j |u_j(t)|^2  >=  P - 4/gamma ,          PR(u(t)) <= 1 / (1 - 4/(gamma P)) ,

all nontrivial exactly when gamma P > 4 -- the dynamical threshold constant is the bandwidth 4.  Below
the threshold the bound is vacuous and the excitation may disperse (and numerically does).

This script integrates eq:dnls (mass-conserving Strang splitting) and checks the three bounds along the
flow for gamma P > 4, and exhibits dispersal at gamma P = 3.
"""

from __future__ import annotations

import numpy as np


def evolve_extrema(n: int, p: float, gamma: float, t_final: float, steps: int) -> tuple[float, float, float]:
    """(min_t max_j|u|^2, max_t PR, min_t sum|u|^4) sampled along the flow."""
    lam = 4 * np.sin(np.pi * np.arange(n) / n) ** 2
    u = np.zeros(n, complex)
    u[0] = np.sqrt(p)
    dt = t_final / steps
    half = np.exp(1j * lam * dt / 2)
    min_max2, max_pr, min_q4 = np.inf, 0.0, np.inf
    for s in range(steps):
        u = np.fft.ifft(half * np.fft.fft(u))
        u = np.exp(-1j * gamma * np.abs(u) ** 2 * dt) * u
        u = np.fft.ifft(half * np.fft.fft(u))
        if s % 40 == 0:
            m2 = np.abs(u) ** 2
            q4 = float((m2**2).sum())
            min_max2 = min(min_max2, float(m2.max()))
            max_pr = max(max_pr, float(m2.sum() ** 2 / q4))
            min_q4 = min(min_q4, q4)
    return min_max2, max_pr, min_q4


def main() -> int:
    n, p = 256, 1.0
    ok = True
    print("Single-site launch, eq:dnls; conservation bounds for gamma P > 4 (exact constant 4):")
    print(f"  {'gP':>5} {'P-4/g':>7} {'min max|u|^2':>13} {'PR bound':>9} {'max PR':>7} {'min sum|u|^4':>13} {'>=P^2-4P/g':>11}")
    for gp in [4.5, 6.0, 8.0, 12.0]:
        g = gp / p
        b1, b4 = p - 4 / g, p**2 - 4 * p / g
        prb = 1 / (1 - 4 / gp)
        mm, pr, q4 = evolve_extrema(n, p, g, t_final=200.0, steps=40000)
        good = mm >= b1 - 1e-6 and pr <= prb + 1e-6 and q4 >= b4 - 1e-6
        ok &= good
        print(f"  {gp:>5.1f} {b1:>7.3f} {mm:>13.4f} {prb:>9.2f} {pr:>7.2f} {q4:>13.4f} {b4:>11.3f}  {'ok' if good else 'FAIL'}")
    mm3, pr3, _ = evolve_extrema(n, p, 3.0 / p, t_final=200.0, steps=40000)
    print(f"  below threshold (gP=3): min max|u|^2 = {mm3:.4f} (dispersal; bound vacuous), PR reaches {pr3:.0f}")
    ok &= mm3 < 0.05  # genuinely disperses below the bandwidth
    print("=" * 72)
    print("RESULT:", "PASS -- no dispersal above gamma P = 4 (conservation bounds hold)" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
