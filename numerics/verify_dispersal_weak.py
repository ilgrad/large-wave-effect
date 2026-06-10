# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2", "scipy>=1.11"]
# ///
"""Quantitative dispersal below threshold at weak coupling (Proposition dispersal).

For the focusing ring DNLS (eq:dnls) with single-site data u(0) = sqrt(P) delta_j0, combining
  * Landau's Bessel bound  max_n |J_n(x)| <= c_L x^{-1/3},  c_L = 0.785746...  [Landau, JLMS 61 (2000)],
  * the ring aliasing estimate (copies super-exponentially small for t <= N/5), and
  * the Duhamel persistence bound  ||u - u_lin||_inf <= gamma P^{3/2} t  (Proposition dnls-persist)
gives  max_j |u_j(t)| <= sqrt(P) c_L (2t)^{-1/3} + gamma P^{3/2} t + eps_N.  Minimizing at
t_* = (c_L / (3 * 2^{1/3} gamma P))^{3/4}  yields

    inf_t max_j |u_j(t)|^2  <=  [4 * 3^{-3/4} (c_L 2^{-1/3})^{3/4}]^2 * P sqrt(gamma P) * (1+o(1))
                             =  1.518... * P sqrt(gamma P) * (1+o(1)),

nontrivial (< P) for gamma P < 0.43: the sub-threshold excitation quantitatively disperses, narrowing the
undecided window of the self-trapping threshold to gamma P in [0.43, 4).

This script checks (1) the Landau constant on a wide sample, and (2) the full chain by direct simulation:
max_j |u_j(t_*)|^2 <= bound for gamma in {0.01, 0.05, 0.2}.
"""

from __future__ import annotations

import numpy as np
from scipy.special import jv

C_L = 0.7858  # Landau's constant, rounded up in the safe direction


def evolve_to(n: int, p: float, gamma: float, t_final: float, steps: int) -> np.ndarray:
    lam = 4 * np.sin(np.pi * np.arange(n) / n) ** 2
    u = np.zeros(n, complex)
    u[0] = np.sqrt(p)
    dt = t_final / steps
    half = np.exp(1j * lam * dt / 2)
    for _ in range(steps):
        u = np.fft.ifft(half * np.fft.fft(u))
        u = np.exp(-1j * gamma * np.abs(u) ** 2 * dt) * u
        u = np.fft.ifft(half * np.fft.fft(u))
    return u


def main() -> int:
    ok = True

    print("[1] Landau bound  sup_n |J_n(x)| x^(1/3) <= c_L = 0.7858:")
    worst = 0.0
    for x in [0.5, 1, 2, 5, 10, 30, 100, 300, 1000, 3000]:
        nn = np.arange(0, int(x) + 300)
        worst = max(worst, float(np.abs(jv(nn, x)).max() * x ** (1 / 3)))
    ok &= worst <= C_L
    print(f"    max over sample = {worst:.4f}  ({'ok' if worst <= C_L else 'FAIL'})")

    print("[2] dispersal chain at t_*: max_j |u(t_*)|^2 <= 1.518 P sqrt(gamma P):")
    p = 1.0
    cp = C_L * 2 ** (-1 / 3)
    print(f"    {'gamma':>6} {'t_*':>6} {'N':>5} {'bound':>7} {'max|u(t_*)|^2':>14}")
    for gamma in [0.01, 0.05, 0.2]:
        t_star = (cp / (3 * gamma * p)) ** 0.75
        bound = (4 * 3 ** (-0.75) * cp**0.75 * (gamma * p) ** 0.25) ** 2 * p
        n = max(256, int(np.ceil(5 * t_star)))
        u = evolve_to(n, p, gamma, t_star, max(2000, int(t_star * 400)))
        m2 = float((np.abs(u) ** 2).max())
        good = m2 <= bound + 1e-9
        ok &= good
        print(f"    {gamma:>6.2f} {t_star:>6.2f} {n:>5} {bound:>7.4f} {m2:>14.4f}  {'ok' if good else 'FAIL'}")

    print("=" * 64)
    print("RESULT:", "PASS -- weak-coupling dispersal bound verified (window now [0.43, 4))" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
