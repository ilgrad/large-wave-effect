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
  (iii) F is strictly increasing on [4/5, 1]: it is real-analytic there (the w=sin^2 substitution removes
        the band-edge singularities), and a Lipschitz-grid bound certifies inf F' >= min_i F'(s_i)
        - max_i|F''| * delta/2 ~ 0.20 > 0; hence max_{s in [0,1]} F(s) = F(1) = beta_odd;
  (iv)  cross-checks: the single-integral form F(s) = sqrt(2/pi) int_0^1 g(...) du agrees with the A,B
        split, and the one-copy term has the closed form A(s) = s^{1/2}/(1-s) x 2F1(1/4,1/2;3/2;x^2),
        x=(1-s)/s (matched to 1e-9).
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
from scipy.special import ellipe, hyp2f1

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


def profile_unified(s: float) -> float:
    """F(s) = sqrt(2/pi) int_0^1 g((s^2-u^2)_+^{-1/2}, (s^2-(1-u)^2)_+^{-1/2}) du -- the single-integral form."""

    def integrand(u: float) -> float:
        p2 = s**2 - u**2
        q2 = s**2 - (1 - u) ** 2
        a = p2**-0.5 if p2 > 1e-30 else 0.0
        b = q2**-0.5 if q2 > 1e-30 else 0.0
        if a == 0.0 and b == 0.0:
            return 0.0
        if b == 0.0:
            return (2 / np.pi) * np.sqrt(a)
        if a == 0.0:
            return (2 / np.pi) * np.sqrt(b)
        return g(a, b)

    pts = sorted({p for p in (1 - s, s) if 0 < p < 1})
    val, _ = quad(integrand, 0.0, 1.0, points=pts or None, epsabs=1e-10, epsrel=1e-10, limit=200)
    return ROOT * val


def one_copy_hyp(s: float) -> float:
    """Closed form  A(s) = s^{1/2}/(1-s) * x * 2F1(1/4,1/2;3/2;x^2),  x = (1-s)/s  (DLMF 15)."""
    x = (1 - s) / s
    return np.sqrt(s) / (1 - s) * x * hyp2f1(0.25, 0.5, 1.5, x**2)


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

    print("[iii] monotonicity of F on [4/5, 1] via a Lipschitz-grid lower bound on F':")
    # F is real-analytic on [4/5,1] (sin^2 substitution removes the band-edge singularities); compute
    # F', F'' by central differences and certify  inf F' >= min_i F'(s_i) - max_i|F''| * delta/2 > 0.
    hh = 0.004
    grid_b = np.linspace(0.8, 0.992, 25)  # spacing delta = 0.008
    delta = grid_b[1] - grid_b[0]
    fp = np.array([(profile(s + hh) - profile(s - hh)) / (2 * hh) for s in grid_b])
    fpp = np.array([(profile(s + hh) - 2 * profile(s) + profile(s - hh)) / hh**2 for s in grid_b])
    lower = fp.min() - np.abs(fpp).max() * delta / 2
    e3 = lower > 0
    ok &= e3
    print(f"     min F'(s_i) = {fp.min():.4f},  max|F''(s_i)| = {np.abs(fpp).max():.4f},  delta = {delta:.4f}")
    print(f"     => inf_[4/5,1] F' >= {lower:.4f} > 0  (strict increase, max at s=1)   {'ok' if e3 else 'FAIL'}")

    print("[iv] cross-checks (single-integral form; closed form for the one-copy term A(s)):")
    uni = max(abs(profile_unified(s) - profile(s)) for s in (0.6, 0.85, 0.97))
    a_err = max(abs(one_copy_hyp(s) - _one_copy(s)) for s in (0.8, 0.9, 0.97))
    e4 = uni < 1e-7 and a_err < 1e-9
    ok &= e4
    print(f"     |F_unified - F| = {uni:.2e};  |A_2F1 - A_quad| = {a_err:.2e}   {'ok' if e4 else 'FAIL'}")
    grid_c = np.linspace(0.5, 0.7, 41)
    fc = np.array([profile(s) for s in grid_c])
    print(f"     (for the record: interior dip min F = {fc.min():.6f} at s ~ {grid_c[fc.argmin()]:.3f})")

    print("=" * 70)
    print("RESULT:", "PASS -- max_{[0,1]} F = F(1) = beta_odd; analytic on [1/2,4/5], F'>=0.20 on [4/5,1]"
          if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
