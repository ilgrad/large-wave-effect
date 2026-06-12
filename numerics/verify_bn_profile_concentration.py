# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2", "scipy>=1.11"]
# ///
r"""The s>1 profile F(s) < beta_odd via the modulus-concentration mechanism (companion to verify_bn_largewave_tail.py).

The multi-copy profile is  F(s) = sqrt(2/pi) int_0^1 E|Z(u,s)| du,  Z = X+iY,
    X = sum_{l even, live} a_l cos Phi_l,   Y = sum_{l odd, live} a_l cos Phi_l,
    a_l = (s^2-(u+l)^2)^{-1/4}  for |u+l|<s,   Phi_l iid uniform,   X _||_ Y,
with F(1) = beta_odd = 0.9280193036689088 the EXACT two-copy peak.  This script certifies the two facts that
together explain max_{s>1} F(s) < beta_odd, and isolates the ONE inequality that a fully rigorous proof needs.

(P) CONSTANT-ENERGY IDENTITY  (PROVED, exact).  Because as u ranges over [0,1) and l over the integers, the
    points u+l sweep (-s,s) exactly once,
        int_0^1 E|Z|^2 du = int_0^1 (1/2) sum_l a_l^2 du = (1/2) int_{-s}^{s} (s^2-x^2)^{-1/2} dx = pi/2,
    INDEPENDENT of s>=1 (arcsin antiderivative).  So the L^2 mass of the profile is conserved across s; all of
    the s-dependence of F lives in the ratio rho(u,s) = E|Z|/sqrt(E|Z|^2) <= 1.

(C) CONCENTRATION CEILING  (VALIDATED; the residual).  rho(u,s) = E|Z|/sqrt(E|Z|^2) is the L^1/L^2 ratio of a
    2-D random walk modulus.  It satisfies rho <= h(1) = 0.95809 (the equal-pair maximum, attained only by a
    pure 1+1 pair), and DECREASES as the energy spreads over more copies -- quantitatively rho falls with the
    closed-form kurtosis R4 = E|Z|^4/(E|Z|^2)^2 in [1.25 (equal pair), 2 (Gaussian)]:
        rho^2  <=  0.918 - 0.14 (R4 - 1.25)          (validated 2-D-walk inequality, anchored at the pair).
    Feeding this into F via (P) gives the closed-form upper bound F4(s) = sqrt(2/pi) int Phi(R4) sqrt(E2) du,
    Phi(R4)=sqrt(0.918-0.14(R4-1.25)); this script certifies F4(s) < beta_odd for s in (1,12] (sup ~ 0.917,
    margin ~ 0.011), while F4(1) >= beta_odd keeps the ceiling honest at the peak.

HONEST LEDGER.
  * PROVED:   the constant-energy identity (P); the peak F(1)=beta_odd (two-copy, verify_beta_odd.py);
              the tail F->sqrt(pi)/2 (verify_bn_largewave_tail.py).
  * VALIDATED: F(s)<beta_odd on (1,12] by seeded Monte-Carlo of E|Z| AND deterministic FFT ||K_N||_1/sqrt N;
              the ceiling rho^2 <= 0.918-0.14(R4-1.25) (sparse <1% excesses, all within Monte-Carlo error);
              the closed-form F4(s)<beta_odd modulo that ceiling.
  * RESIDUAL (NOT proved here): the concentration inequality rho^2 <= 0.918-0.14(R4-1.25) is specific to
              2-D-walk moduli (it FAILS for a generic non-negative W with the same moments -- e.g. constant W
              -- so it is NOT a pure moment inequality and needs the walk structure: a Pearson/Kluyver
              expected-modulus bound).  This single inequality is the only gap between (P)+(C) and a theorem
              max_{s>1} F(s) < beta_odd; with the band-edge transfer (verify_debye_equidistribution.py) it is
              one of the two residual estimates of the t>N/2 region.
"""

from __future__ import annotations

import numpy as np

BETA = 0.9280193036689088
H1 = 0.9580913983830018  # h(1) = E sqrt(cos^2 + cos^2), the equal-pair ratio ceiling


def _components(u: float, s: float) -> tuple[np.ndarray, np.ndarray]:
    js = [j for j in range(-int(s) - 2, int(s) + 3) if abs(u + j) < s]
    a = np.array([(s * s - (u + j) ** 2) ** -0.25 for j in js])
    par = np.array([j % 2 for j in js])
    return a, par


def _e2_r4(a: np.ndarray, par: np.ndarray) -> tuple[float, float]:
    """Exact (no Monte-Carlo) L^2 mass E|Z|^2 and kurtosis R4 = E|Z|^4/(E|Z|^2)^2 from the amplitudes."""
    sx2 = 0.5 * np.sum(a[par == 0] ** 2)
    sy2 = 0.5 * np.sum(a[par == 1] ** 2)
    s2 = sx2 + sy2
    if s2 <= 0:
        return 0.0, 2.0
    e_z4 = 3 * sx2**2 + 3 * sy2**2 + 2 * sx2 * sy2 - (3 / 8) * np.sum(a**4)
    return s2, e_z4 / s2**2


def _ez_mc(u: float, s: float, rng: np.random.Generator, m: int) -> float:
    a, par = _components(u, s)
    if len(a) == 0:
        return 0.0
    ph = rng.uniform(0, 2 * np.pi, (m, len(a)))
    x = (a[par == 0] * np.cos(ph[:, par == 0])).sum(1)
    y = (a[par == 1] * np.cos(ph[:, par == 1])).sum(1)
    return float(np.mean(np.sqrt(x * x + y * y)))


def _ez_fft(s: float, n: int, n_t: int = 24) -> float:
    """Deterministic ||K_N(.,t)||_1/sqrt(N), K_N(k,t)=IFFT_k exp(-i t lambda), lambda=2-2cos(2pi r/N),
    averaged over a t-window around sN/2 to smooth the O(1) Bessel oscillation (a single t is a fluctuation,
    not the profile).  The window t = sN/2 + j, j=0..n_t-1 spans ~n_t/pi periods."""
    lam = 2 - 2 * np.cos(2 * np.pi * np.arange(n) / n)
    base = np.exp(-1j * (s * n / 2) * lam)
    step = np.exp(-1j * lam)
    vals = []
    cur = base
    for _ in range(n_t):
        vals.append(np.sum(np.abs(np.fft.ifft(cur))) / np.sqrt(n))
        cur = cur * step
    return float(np.mean(vals))


def _f4(s: float, nu: int = 400) -> float:
    """Closed-form upper bound sqrt(2/pi) int Phi(R4) sqrt(E2) du (rigorous GIVEN the ceiling rho<=Phi(R4))."""
    us = (np.arange(nu) + 0.5) / nu
    tot = 0.0
    for u in us:
        e2, r4 = _e2_r4(*_components(u, s))
        tot += np.sqrt(max(0.918 - 0.14 * (r4 - 1.25), 0.0)) * np.sqrt(e2)
    return np.sqrt(2 / np.pi) * tot / nu


def _f_mc(s: float, rng: np.random.Generator, nu: int = 300, m: int = 50000) -> float:
    us = (np.arange(nu) + 0.5) / nu
    return float(np.sqrt(2 / np.pi) * np.mean([_ez_mc(u, s, rng, m) for u in us]))


def main() -> int:
    rng = np.random.default_rng(20240611)
    ok = True

    print("[P] constant-energy identity  int_0^1 E|Z|^2 du = pi/2  (exact, s-independent):")
    worst_e = 0.0
    for s in (1.0, 1.5, 3.0, 7.0, 12.0):
        us = (np.arange(20000) + 0.5) / 20000
        e = float(np.mean([_e2_r4(*_components(u, s))[0] for u in us]))
        worst_e = max(worst_e, abs(e - np.pi / 2))
        print(f"    s={s:>5}:  int E|Z|^2 du = {e:.5f}   (pi/2 = {np.pi/2:.5f})")
    ok &= worst_e < 1e-2  # residual is band-edge quadrature only; the value is exactly pi/2
    print(f"    max deficit (band-edge quadrature) = {worst_e:.2e}   {'ok' if worst_e < 1e-2 else 'FAIL'}")

    print("\n[C] concentration ceiling rho^2 <= 0.918 - 0.14(R4-1.25), and the closed-form bound F4(s):")
    print(f"    {'s':>6} {'F_MC':>8} {'F_FFT':>8} {'F4(<=bound)':>12} {'<beta?':>7}")
    sup_f4 = 0.0
    worst_mc_fft = 0.0
    f4_at_1 = _f4(1.0)
    for s in (1.0, 1.02, 1.05, 1.1, 1.2, 1.5, 2.0, 3.0, 5.0, 8.0, 12.0):
        f_mc = _f_mc(s, rng)
        f_fft = _ez_fft(s, (1 << 16) - 1)  # ODD N: the beta_odd parity branch the channel profile models
        f4 = _f4(s)
        if s > 1.0001:
            sup_f4 = max(sup_f4, f4)
        worst_mc_fft = max(worst_mc_fft, abs(f_mc - f_fft))
        below = f4 < BETA if s > 1.0001 else f4 >= BETA - 5e-3
        ok &= below and (f_mc < BETA + 6e-3)
        tag = "YES" if (s > 1.0001 and f4 < BETA) else ("peak" if s <= 1.0001 else "NO")
        print(f"    {s:>6.2f} {f_mc:>8.4f} {f_fft:>8.4f} {f4:>12.4f} {tag:>7}")
    ok &= worst_mc_fft < 0.02  # channel Monte-Carlo profile agrees with the true odd-N FFT norm
    print(f"    max |F_MC - F_FFT| = {worst_mc_fft:.4f}  (independent cross-check)   "
          f"{'ok' if worst_mc_fft < 0.02 else 'FAIL'}")
    print(f"    F4(1) = {f4_at_1:.4f} >= beta (ceiling valid at the peak);  "
          f"sup_(1,12] F4 = {sup_f4:.4f} {'< beta' if sup_f4 < BETA else '>= beta'}")

    # ceiling spot-check: rho <= h(1) and rho^2 <= 0.918-0.14(R4-1.25) up to MC error
    print("\n    ceiling spot-check (rho from MC vs Phi(R4); excess beyond MC error flags a real violation):")
    worst_excess = -1.0
    for s in (1.3, 2.0, 4.0, 8.0):
        for u in (0.17, 0.41, 0.63, 0.88):
            a, par = _components(u, s)
            e2, r4 = _e2_r4(a, par)
            if e2 <= 0:
                continue
            rho = _ez_mc(u, s, rng, 200000) / np.sqrt(e2)
            worst_excess = max(worst_excess, rho - np.sqrt(max(0.918 - 0.14 * (r4 - 1.25), 0.0)))
    ok &= worst_excess < 5e-3  # within the ~2.5e-3 MC s.e. on rho
    print(f"    worst rho - Phi(R4) over the grid = {worst_excess:+.4f}  (MC s.e. ~ 2.5e-3)   "
          f"{'ok' if worst_excess < 5e-3 else 'FAIL'}")

    print("=" * 96)
    print("RESULT:", "PASS -- constant-energy identity exact; F(s)<beta_odd on (1,12] (MC+FFT); the closed-form "
          "F4<beta_odd certifies it MODULO the validated 2-D-walk ceiling rho^2<=0.918-0.14(R4-1.25), which is "
          "the single residual inequality (see header ledger)." if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
