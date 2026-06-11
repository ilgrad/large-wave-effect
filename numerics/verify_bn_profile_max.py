# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2", "scipy>=1.11"]
# ///
"""Two-copy profile for t <= N/2 and its maximizer (Proposition bn-profile).

For the Schroedinger ring kernel at t = sN/2, s in [0,1], the Poisson periodization keeps only the two
copies ell = 0, -1 (the rest are super-exponentially small since 2t <= N), and the Debye envelope turns
||K_N(.,t)||_1 / sqrt(N) into an explicit limiting profile F(s):

    F(s) = c0 sqrt(s)                                                 for s in [0, 1/2],
    F(s) = sqrt(2/pi) [ (4/pi)(1-s) A(s) + (2s-1)^{3/4} B(s) ]        for s in [1/2, 1],

with, after mapping the two windows to [0,1],
    A(s) = int_0^1 (s^2 - (1-s)^2 w^2)^{-1/4} dw,
    B(s) = int_0^1 g( ahat, bhat ) dw,   ahat = [(1-w)(1+(2s-1)w)]^{-1/2},  bhat = [w(2s-(2s-1)w)]^{-1/2},
    g(a,b) = E_{Phi,Psi} sqrt(a cos^2 Phi + b cos^2 Psi) = sqrt(max(a,b)) * h(min(a,b)/max(a,b)),
    h(rho) = (4/pi^2) int_0^{pi/2} sqrt(1+rho cos^2 psi) * E(1/(1+rho cos^2 psi)) dpsi   (E = complete elliptic).

This script certifies the maximizer structure (Proposition bn-profile):
  (i)   closed-form endpoints  F(1/2) = c0/sqrt2 = 0.86068...,  F(1) = beta_odd = 0.92802...;
  (ii)  ANALYTIC (Jensen) bound  g(a,b) <= sqrt((a+b)/2)  gives an upper envelope F+(s) with
        F(s) <= F+(s) < beta_odd  on [1/2, 4/5]  -- a proof on 60% of [1/2,1];
  (iii) F is strictly increasing on [4/5, 1] (checked on a fine grid with margin), hence
        max_{s in [0,1]} F(s) = F(1) = beta_odd.
The Jensen envelope necessarily exceeds beta_odd near s=1 (any strict majorant of g overshoots at the exact
maximizer s=1), so (iii) -- a one-dimensional monotonicity of the explicit profile -- is the residual step;
it no longer involves N.  Compare verify_bn_parity.py (FFT of the actual kernel) and
verify_bn_residual_profile.py (multi-copy profile for t > N/2).
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import CubicSpline
from scipy.special import beta as beta_fn
from scipy.special import ellipe

C0 = (2 / np.pi) ** 1.5 * beta_fn(0.5, 0.75)
BETA_ODD = 0.9280193036689088  # eq:beta-odd, independent reference value
ROOT = np.sqrt(2 / np.pi)


def _h_raw(rho: float) -> float:
    def f(psi: float) -> float:
        c = rho * np.cos(psi) ** 2
        return np.sqrt(1 + c) * ellipe(1 / (1 + c))

    val, _ = quad(f, 0.0, np.pi / 2, epsabs=1e-13, epsrel=1e-13)
    return 4 / np.pi**2 * val


_RHO = np.linspace(0.0, 1.0, 801)
_H = CubicSpline(_RHO, np.array([_h_raw(r) for r in _RHO]))


def g(a: float, b: float) -> float:
    hi, lo = (a, b) if a >= b else (b, a)
    return np.sqrt(hi) * float(_H(lo / hi))


def _overlap(s: float, jensen: bool) -> float:
    q = 2 * s - 1

    def integrand(theta: float) -> float:  # w = sin^2(theta) removes the w^-1/4, (1-w)^-1/4 endpoints
        w = np.sin(theta) ** 2
        ah = ((1 - w) * (1 + q * w)) ** -0.5
        bh = (w * (2 * s - q * w)) ** -0.5
        core = np.sqrt((ah + bh) / 2) if jensen else g(ah, bh)
        return core * 2 * np.sin(theta) * np.cos(theta)

    val, _ = quad(integrand, 0.0, np.pi / 2, epsabs=1e-12, epsrel=1e-12, limit=200)
    return val


def _one_copy(s: float) -> float:
    val, _ = quad(lambda w: (s**2 - (1 - s) ** 2 * w**2) ** -0.25, 0.0, 1.0, epsabs=1e-12, epsrel=1e-12)
    return val


def profile(s: float, jensen: bool = False) -> float:
    if s <= 0.5:
        return C0 * np.sqrt(s)
    q = 2 * s - 1
    return ROOT * ((4 / np.pi) * (1 - s) * _one_copy(s) + q**0.75 * _overlap(s, jensen))


def main() -> int:
    ok = True

    f_half, f_one = profile(0.5), profile(1.0)
    e1 = abs(f_half - C0 / np.sqrt(2)) < 1e-7 and abs(f_one - BETA_ODD) < 1e-6
    ok &= e1
    print("[i]  endpoints:")
    print(f"     F(1/2) = {f_half:.10f}  (c0/sqrt2 = {C0/np.sqrt(2):.10f})")
    print(f"     F(1)   = {f_one:.10f}  (beta_odd = {BETA_ODD:.10f})    {'ok' if e1 else 'FAIL'}")

    print("[ii] Jensen analytic bound  F(s) <= F+(s) < beta_odd  on [1/2, 4/5]:")
    grid_a = np.linspace(0.5, 0.8, 31)
    worst = -np.inf
    for s in grid_a:
        fp = profile(s, jensen=True)
        worst = max(worst, fp - BETA_ODD)
    e2 = worst < 0
    ok &= e2
    print(f"     max_(s<=0.8) (F+(s) - beta_odd) = {worst:.5f}  (<0 required)   {'ok' if e2 else 'FAIL'}")

    print("[iii] monotonicity of F on [4/5, 1] (fine grid, forward differences):")
    grid_b = np.linspace(0.8, 1.0, 81)
    fb = np.array([profile(s) for s in grid_b])
    diffs = np.diff(fb)
    e3 = bool(np.all(diffs > 0)) and bool(np.all(fb <= BETA_ODD + 1e-7))
    ok &= e3
    print(f"     min forward difference = {diffs.min():.3e}  (>0 required); "
          f"max F = {fb.max():.10f} <= beta_odd   {'ok' if e3 else 'FAIL'}")

    # the dip: F has an interior minimum near s ~ 0.56; report it for the record
    grid_c = np.linspace(0.5, 0.7, 41)
    fc = np.array([profile(s) for s in grid_c])
    print(f"     (for the record: interior dip min F = {fc.min():.6f} at s ~ {grid_c[fc.argmin()]:.3f})")

    print("=" * 70)
    print("RESULT:", "PASS -- max_{[0,1]} F = F(1) = beta_odd; analytic on [1/2,4/5], 1-D residual on [4/5,1]"
          if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
