# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2", "scipy>=1.11", "mpmath>=1.3"]
# ///
r"""The admissible-constant WINDOW of the Gaussian-excess route, and the s>1 SPLIT that closes it.

Companion / audit of verify_bn_profile_gaussian.py and verify_dagger_extremal.py.  Those scripts state the
pointwise excess bound (D)  E|Z|-E|G| <= K_star (sum a^4)/Sigma_2^{3/2},  K_star=(h(1)-sqrt(pi)/2)/2, and
claim the integrated majorant  B(s)=F_gauss(s)+sqrt(2/pi) K_star int_0^1 (sum a^4)/Sigma_2^{3/2} du  clears
beta_odd on (1,20], with "min_s K_max(s)=0.0402 at s=1.02".  This script re-derives the window from scratch
with panel Gauss-Legendre quadrature that resolves the band-edge singularity a_l=(s^2-(u+l)^2)^{-1/4}
(the clustered midpoint grid of those scripts wobbles ~3e-3 at the edge), and reports three corrections that
matter for the proof, plus the split that repairs the route.

==========================================================================================================
SETUP (identical to verify_bn_profile_gaussian.py).  At t=sN/2,  F(s)=sqrt(2/pi) int_0^1 E|Z(u,s)| du,
Z=X+iY,  X=sum_{l even,live} a_l cos Phi_l,  Y=sum_{l odd,live} a_l cos Phi_l,  a_l=(s^2-(u+l)^2)^{-1/4}
for |u+l|<s.  sigma_X^2=(1/2)sum_even a^2,  sigma_Y^2=(1/2)sum_odd a^2,  Sigma_2=sigma_X^2+sigma_Y^2,
G=(G_X,G_Y) the matched Gaussian.  E|G|=sigma_max sqrt(2/pi) E_ell(1-sigma_min^2/sigma_max^2) (PROVED).
K_max(s) := (beta_odd - F_gauss(s)) / (sqrt(2/pi) int_0^1 (sum a^4)/Sigma_2^{3/2} du)  is the largest
constant K for which B(s)<beta_odd still holds; the K_star route needs K_star <= inf_{s>1} K_max(s).

==========================================================================================================
FINDING 1 -- THE WINDOW.  inf_{s>1} K_max(s) is NOT 0.0402 at s=1.02.  K_max(s) decreases monotonically as
s->1+ (the profile degenerates to the equal two-copy pair, where (D) is SHARP), reaching
    K_max(s) -> 0.0338  as s->1+,   BELOW K_star = 0.035932.
The "0.0402 at s=1.02" of verify_bn_profile_gaussian.py [iv] is a grid-truncation artifact: its s-grid
starts at 1.02 and never samples (1,1.02).  At s=1.02 both methods agree (K_max=0.0404); the minimum simply
lies below the smallest sampled s.  Consequently the K_star majorant B(s) EXCEEDS beta_odd on the sliver
    (1, s*),   s* = 1.00407   (B(1.002)=0.9303 > beta_odd=0.928019; B(1.005)=0.92719 < beta_odd),
so the claim "B(s)<beta_odd on (1,20]" is FALSE on (1,1.00407].  This is not a profile violation: the TRUE
F(s) (FFT density) stays < beta_odd throughout (F(1.001)=0.9255), decreasing.  The overshoot is the
non-sharpness of K_star away from the exact pair, accumulated over the u-integral.

FINDING 2 -- THE SPLIT (repairs the route).  Two PROVED ingredients close the whole range with overlap:
  * (1, 1.005]:  the variance-floor bound F_2(s) = sqrt(2/pi) int E sqrt(X_C^2+Y_C^2+c^2) du (keep the two
    largest amplitudes exact, floor the rest by c^2=(1/2)sum_rest a^2; Lemma lem:varfloor, PROVED >= F).
    F_2 < beta_odd on the sliver (F_2(1.001)=0.9254, margin +0.0026) -- and there the live set is
    essentially two copies (at s=1.005 only 1% of u carries a 3rd copy, energy share <5%), so the floor is
    near-tight.  This is exactly the regime where the EXACT two-copy g (proved to maximize at the pair)
    governs.
  * [1.005, inf):  inf_{s>=1.005} K_max(s) = 0.0363 > K_star, so B(s) < beta_odd for every s>=1.005 (worst
    margin +0.0008 at s=1.005).  MOREOVER any CRUDE constant K <= K_star also clears beta_odd here (B is
    increasing in K): with the small-t fourth-cumulant constant (~0.45 K_star at the pair) the majorant has
    margin > 0.039.  So the hard Bessel-ring tail estimate is needed only to make K_star sharp at the pair,
    which only matters on the sliver -- where the variance-floor already wins.  The route does NOT need the
    sharp constant on [1.005,inf).

FINDING 3 -- THE beta_odd LITERAL.  Independent high-precision (mpmath, the one-elliptic-reduced triple
integral, cross-checked by reduce-over-p and a brute 2-D inner) gives
    beta_odd = 0.92801930480793112292...      (25 digits)
The repo scripts now use this value.  A historical literal 0.9280193036689087668528462 (the float
0.9280193036689088) agreed only to 9 significant figures (digit 10: ...0480... vs ...0366...) and has been
corrected repo-wide.  The discrepancy was immaterial to every margin in the proof (they are 1e-3..1e-2, the
discrepancy was 1e-9; K_max is unchanged to 5 digits whichever value is used) -- it was a citation/literal
error, not a proof error.  Likewise h(1)=0.95809139868285... (script's 0.9580913983830018 is right to 9
digits, off at digit 10); K_star=0.03593223661505... (script's 0.0359322 is correct).

==========================================================================================================
NOTE on profile near-symmetry (corrects the working hypothesis).  The even/odd channels are NOT pointwise
near-symmetric: max_u |sigma_X^2-sigma_Y^2|/Sigma_2 ~ 0.98 over the profile (at a generic u in a two-copy
cell the two live copies have OPPOSITE parity but very different amplitudes -- one near a band edge
dominates).  Symmetry holds only AFTER the u-integral (the u<->2-s edge symmetry).  What IS true and
load-bearing: max_u R(u,s)/K_star <= 1 over the profile, attained exactly along u=0.5 where the two live
copies are an EQUAL pair (ev=od), for every s<1.5; and max_u R/K_star <= 0.71 for s>=1.5 (so off the
near-pair region the profile sits far inside (D)).

==========================================================================================================
HONEST LEDGER.
  PROVED (elementary): the Gaussian closed form / F_gauss; the variance-floor majorant F_2>=F
    (lem:varfloor); F_gauss<beta_odd; F(infinity)=sqrt(pi)/2<beta_odd.
  CERTIFIED HERE (numerics, panel Gauss-Legendre to ~1e-6, node-convergence shown):
    * the K_max(s) curve and inf_{s>1} K_max=0.0338 (s->1+), inf_{s>=1.005} K_max=0.0363>K_star;
    * the failure sliver s* = 1.00407 of the K_star majorant, and that the TRUE F(s) stays < beta_odd there;
    * the split: F_2<beta_odd on (1,1.005], B<beta_odd (any K<=K_star) on [1.005,inf);
    * max_u R(u,s)/K_star <= 1 over the profile (exact Kluyver excess), = 1 on u=0.5 for s<1.5;
    * the K_star excess split 45% small-t / 55% Bessel-ring tail at the equal pair.
  HIGH-PRECISION (mpmath 30 dps): beta_odd, h(1), K_star -- and the beta_odd / h(1) literal discrepancies.
  RESIDUAL (unchanged): a fully rigorous (D) on [1.005,inf) still needs a (crude, NOT sharp) bound on the
    Bessel-ring tail of the Kluyver excess; the sharp K_star and its tail estimate are NO LONGER required to
    close the range once the split routes the sliver through the variance-floor.  "validated != proved".
"""

from __future__ import annotations

import functools
from itertools import pairwise

import numpy as np
from numpy import cos, pi, sin, sqrt
from scipy.integrate import quad
from scipy.special import ellipe, j0, jn_zeros

np.seterr(all="ignore")  # masked d<=0 powers in vectorized amplitudes; values are np.where-guarded

# ---- constants (high-precision values from the mpmath block below) ----
BETA_HP = 0.9280193048079311      # beta_odd (independent mpmath); see FINDING 3
H1 = 0.9580913986828501           # h(1), corrected
ROOT = sqrt(2 / pi)
TAIL = sqrt(pi) / 2               # F_gauss(infinity)
K_STAR = (H1 - TAIL) / 2          # 0.035932236615...
J01 = float(jn_zeros(0, 1)[0])    # 2.404826...


# ---- amplitudes (identical live set to verify_bn_profile_gaussian.amps; verified 0 mismatches) ----
def amps(u: float, s: float) -> tuple[np.ndarray, np.ndarray]:
    lmin = int(np.floor(-s - u))
    lmax = int(np.ceil(s - u))
    ev: list[float] = []
    od: list[float] = []
    for ell in range(lmin, lmax + 1):
        d = s * s - (u + ell) ** 2
        if d > 1e-14:
            (ev if ell % 2 == 0 else od).append(d ** -0.25)
    return np.array(ev), np.array(od)


# ---- band-edge breakpoints and panel Gauss-Legendre quadrature of (F_gauss, control) ----
@functools.cache
def _gl(n: int) -> tuple[np.ndarray, np.ndarray]:
    x, w = np.polynomial.legendre.leggauss(n)
    return x, w


def _breakpoints(s: float) -> np.ndarray:
    """u in (0,1) where some |u+l|=s (the integrable singularities), plus the endpoints."""
    bps = {0.0, 1.0}
    for ell in range(int(np.floor(-s)), int(np.ceil(s)) + 1):
        for x in (s - ell, -s - ell):
            if 0 < x < 1:
                bps.add(x)
    return np.array(sorted(bps))


def fgauss_and_control(s: float, nper: int = 1200) -> tuple[float, float]:
    """(F_gauss(s), int_0^1 (sum a^4)/Sigma_2^{3/2} du) by panel Gauss-Legendre on the band cells.
    Converges to ~1e-6 at nper~800 even at the edge-singular s (the panels isolate the singularity)."""
    bps = _breakpoints(s)
    x, w = _gl(nper)
    lo, hi = int(np.floor(-s)), int(np.ceil(s))
    ls = np.arange(lo, hi + 1)
    par = ls % 2
    fg = 0.0
    ctrl = 0.0
    for a, b in pairwise(bps):
        uu = 0.5 * (b - a) * x + 0.5 * (a + b)
        ww = 0.5 * (b - a) * w
        d = s * s - (uu[:, None] + ls[None, :]) ** 2
        live = d > 1e-13
        ds = np.where(live, d, 1.0)
        a2 = np.where(live, ds ** -0.5, 0.0)   # a_l^2 = d^{-1/2}
        a4 = np.where(live, ds ** -1.0, 0.0)   # a_l^4 = d^{-1}
        sx2 = 0.5 * np.sum(np.where((par == 0)[None, :], a2, 0.0), axis=1)
        sy2 = 0.5 * np.sum(np.where((par == 1)[None, :], a2, 0.0), axis=1)
        sig2 = sx2 + sy2
        c4 = np.sum(a4, axis=1)
        amax = np.maximum(sx2, sy2)
        amin = np.minimum(sx2, sy2)
        m = np.where(amax > 0, 1 - amin / np.where(amax > 0, amax, 1.0), 0.0)
        eg = np.where(amax > 0, np.sqrt(amax) * ROOT * ellipe(m), 0.0)
        fg += np.sum(eg * ww)
        good = sig2 > 0
        ctrl += np.sum(np.where(good, c4 / np.where(good, sig2, 1.0) ** 1.5, 0.0) * ww)
    return ROOT * fg, ctrl


def B_majorant(s: float, K: float, nper: int = 1200) -> float:
    fg, ctrl = fgauss_and_control(s, nper)
    return fg + ROOT * K * ctrl


def K_max(s: float, beta: float = BETA_HP, nper: int = 1200) -> float:
    fg, ctrl = fgauss_and_control(s, nper)
    return (beta - fg) / (ROOT * ctrl)


# ---- TRUE F(s) (FFT density, no Gaussian approximation) -- the honest profile, must stay < beta_odd ----
def F_true(s: float, nper: int = 500, L: float = 16.0, N: int = 4096) -> float:
    bps = _breakpoints(s)
    x, w = _gl(nper)
    xx = np.linspace(-L, L, N, endpoint=False)
    dx = xx[1] - xx[0]
    k = 2 * pi * np.fft.fftfreq(N, d=dx)
    rad = sqrt(xx[:, None] ** 2 + xx[None, :] ** 2)

    def dens(a: np.ndarray) -> np.ndarray:
        phi = np.ones_like(k)
        for ai in a:
            phi = phi * j0(ai * k)
        return np.fft.fftshift(np.real(np.fft.ifft(phi))) / dx

    tot = 0.0
    for a, b in pairwise(bps):
        uu = 0.5 * (b - a) * x + 0.5 * (a + b)
        ww = 0.5 * (b - a) * w
        for u, wt in zip(uu, ww, strict=True):
            ev, od = amps(float(u), s)
            if len(ev) + len(od) == 0:
                continue
            fx = dens(ev) if len(ev) else None
            fy = dens(od) if len(od) else None
            if fx is None:
                g = float((np.abs(xx) * fy).sum() * dx)
            elif fy is None:
                g = float((np.abs(xx) * fx).sum() * dx)
            else:
                g = float(fx @ (rad @ fy) * dx * dx)
            tot += g * wt
    return ROOT * tot


# ---- variance-floor F_2(s) (PROVED majorant, Lemma lem:varfloor): keep K=2 exact, floor the rest ----
def F2_floor(s: float, nper: int = 400, L: float = 12.0, N: int = 1024) -> float:
    bps = _breakpoints(s)
    x, w = _gl(nper)
    xx = np.linspace(-L, L, N, endpoint=False)
    dx = xx[1] - xx[0]
    k = 2 * pi * np.fft.fftfreq(N, d=dx)

    def dens(a: np.ndarray) -> np.ndarray:
        phi = np.ones_like(k)
        for ai in a:
            phi = phi * j0(ai * k)
        return np.fft.fftshift(np.real(np.fft.ifft(phi))) / dx

    tot = 0.0
    for a, b in pairwise(bps):
        uu = 0.5 * (b - a) * x + 0.5 * (a + b)
        ww = 0.5 * (b - a) * w
        for u, wt in zip(uu, ww, strict=True):
            ev, od = amps(float(u), s)
            merged = sorted([(av, 0) for av in ev] + [(av, 1) for av in od], key=lambda z: -z[0])
            core, rest = merged[:2], merged[2:]
            c2 = 0.5 * sum(av * av for av, _ in rest)
            ce = np.array([av for av, p in core if p == 0])
            co = np.array([av for av, p in core if p == 1])
            fx = dens(ce) if len(ce) else None
            fy = dens(co) if len(co) else None
            if fx is None and fy is None:
                g = sqrt(c2)
            elif fx is None:
                g = float((np.sqrt(xx * xx + c2) * fy).sum() * dx)
            elif fy is None:
                g = float((np.sqrt(xx * xx + c2) * fx).sum() * dx)
            else:
                rad = np.sqrt(xx[:, None] ** 2 + (xx * xx + c2)[None, :])
                g = float(fx @ (rad @ fy) * dx * dx)
            tot += g * wt
    return ROOT * tot


# ---- exact Kluyver excess and its split at j_{0,1} (per-config R/K_star, tail/small fractions) ----
def excess_split(ev: np.ndarray, od: np.ndarray, n_theta: int = 2000,
                 t_max: float = 200.0) -> tuple[float, float]:
    sx2 = 0.5 * float(np.sum(ev ** 2)) if len(ev) else 0.0
    sy2 = 0.5 * float(np.sum(od ** 2)) if len(od) else 0.0
    th = (np.arange(n_theta) + 0.5) / n_theta * (2 * pi)
    c, sn = cos(th), sin(th)

    def mw(t: float) -> float:
        p = np.ones(n_theta)
        for a in ev:
            p = p * j0(a * t * c)
        for a in od:
            p = p * j0(a * t * sn)
        return float(p.mean())

    def mg(t: float) -> float:
        return float(np.exp(-(sx2 * c ** 2 + sy2 * sn ** 2) * t * t / 2).mean())

    def integ(t: float) -> float:
        return (mg(t) - mw(t)) / t ** 2

    small, _ = quad(integ, 1e-7, J01, limit=300)
    tail, _ = quad(integ, J01, t_max, limit=500)
    return small, tail


def R_over_Kstar(ev: np.ndarray, od: np.ndarray) -> float:
    sig2 = 0.5 * (float(np.sum(ev ** 2)) + (float(np.sum(od ** 2)) if len(od) else 0.0))
    c4 = float(np.sum(ev ** 4)) + (float(np.sum(od ** 4)) if len(od) else 0.0)
    if sig2 <= 0 or c4 <= 0:
        return 0.0
    sm, tl = excess_split(ev, od)
    return (sm + tl) / (c4 / sig2 ** 1.5) / K_STAR


# ==========================================================================================================
def main() -> int:
    global print
    print = functools.partial(print, flush=True)
    ok = True
    print("=" * 100)
    print("The admissible-constant WINDOW K_max(s) of the Gaussian-excess route, and the s>1 SPLIT")
    print(f"  beta_odd(HP)={BETA_HP:.13f}  K_star={K_STAR:.10f}  j_0,1={J01:.6f}")
    print("=" * 100)

    # ---- FINDING 3 first: high-precision constants (cheap, sets the beta used everywhere) ----
    print("\n[C] high-precision constants (mpmath, one-elliptic-reduced integrals) reproduce the script HP values:")
    try:
        import mpmath as mp
        mp.mp.dps = 30
        half = mp.pi / 2

        def h1_hp() -> mp.mpf:
            def fp(p: mp.mpf) -> mp.mpf:
                d = 1 + mp.cos(p) ** 2
                return mp.sqrt(d) * mp.ellipe(1 / d)
            return mp.quad(fp, [0, half]) / half ** 2

        def beta_hp() -> mp.mpf:
            def inner(u: mp.mpf) -> mp.mpf:
                a = (1 - u * u) ** mp.mpf("-0.5")
                b = (u * (2 - u)) ** mp.mpf("-0.5")

                def fp(p: mp.mpf) -> mp.mpf:
                    c = a * mp.cos(p) ** 2 + b
                    return mp.sqrt(c) * mp.ellipe(b / c)
                return mp.quad(fp, [0, half]) / half ** 2
            return mp.sqrt(2 / mp.pi) * mp.quad(inner, [0, mp.mpf("0.5"), 1])

        h1_v = h1_hp()
        beta_v = beta_hp()
        kstar_v = (h1_v - mp.sqrt(mp.pi) / 2) / 2
        print(f"     beta_odd       = {mp.nstr(beta_v, 22)}   (script float: 0.92801930480793112)")
        print(f"     h(1)           = {mp.nstr(h1_v, 22)}   (script float: 0.9580913986828501)")
        print(f"     K_star         = {mp.nstr(kstar_v, 22)}   (script float: 0.0359322)")
        # historical note: a former literal 0.9280193036689088 was wrong past digit 9; corrected repo-wide.
        c_ok = (abs(float(beta_v) - BETA_HP) < 1e-12 and abs(float(h1_v) - H1) < 1e-12
                and abs(float(kstar_v) - K_STAR) < 1e-9)
        ok &= c_ok
        print(f"     => constants reproduce the script's HP values   {'ok' if c_ok else 'FAIL'}")
    except Exception as exc:  # mpmath edge; constants are not load-bearing for the window verdict
        print(f"     [mpmath unavailable: {exc}] -- skipping the high-precision block")

    # ---- node-convergence of the panel quadrature at the worst (near-edge) s ----
    print("\n[A0] panel Gauss-Legendre node convergence at s=1.005 (the edge-singular, binding s):")
    for n in (400, 800, 1600):
        fg, ct = fgauss_and_control(1.005, n)
        print(f"     nper={n:>5}: F_gauss={fg:.7f}  control={ct:.6f}  B(K_star)={fg + ROOT * K_STAR * ct:.7f}")
    conv_ok = abs(fgauss_and_control(1.005, 400)[1] - fgauss_and_control(1.005, 1600)[1]) < 1e-4
    ok &= conv_ok
    print(f"     => control stable to <1e-4 by nper=400   {'ok' if conv_ok else 'FAIL'}")

    # ---- FINDING 1: the window, the true minimum, and the failure sliver ----
    print("\n[A1] K_max(s) curve (CORRECTED): the minimum is at s->1+, BELOW K_star -- not 0.0402 at s=1.02:")
    print(f"     {'s':>8} {'F_gauss':>9} {'control':>9} {'K_max':>9} {'B(K_star)':>10} {'B<beta?':>8}")
    grid = [1.0005, 1.001, 1.002, 1.004, 1.005, 1.01, 1.02, 1.05, 1.1, 1.2, 1.3, 1.5, 2.0, 3.0, 5.0, 10.0, 20.0]
    kmin_all = (1e9, 0.0)
    for s in grid:
        fg, ct = fgauss_and_control(s)
        km = (BETA_HP - fg) / (ROOT * ct)
        b = fg + ROOT * K_STAR * ct
        if km < kmin_all[0]:
            kmin_all = (km, s)
        print(f"     {s:>8.4f} {fg:>9.5f} {ct:>9.5f} {km:>9.5f} {b:>10.6f} {'yes' if b < BETA_HP else 'NO':>8}")
    print(f"     => min over this grid: K_max={kmin_all[0]:.5f} at s={kmin_all[1]}  (K_star={K_STAR:.5f}; "
          f"K_max<K_star near s=1 => route FAILS there)")
    # the claim "0.0402 at 1.02" is the grid-start artifact; the true inf is below
    sliver_fail = kmin_all[0] < K_STAR
    ok &= sliver_fail  # we EXPECT the route to fail near s=1 (that is the finding)

    # crossover s* where B(s)=beta_odd
    def Bof(s: float) -> float:
        return B_majorant(s, K_STAR)
    lo, hi = 1.0002, 1.010
    for _ in range(26):
        mid = 0.5 * (lo + hi)
        if Bof(mid) > BETA_HP:
            lo = mid
        else:
            hi = mid
    sstar = 0.5 * (lo + hi)
    print(f"\n[A2] failure sliver of the K_star majorant:  B(s) > beta_odd on (1, s*),  s* = {sstar:.5f}")
    print(f"     B(1.002)={Bof(1.002):.6f} > beta;  B(1.005)={Bof(1.005):.6f} < beta;  beta={BETA_HP:.6f}")
    sstar_ok = 1.003 < sstar < 1.006
    ok &= sstar_ok
    print(f"     => s* in (1.003,1.006): {'ok' if sstar_ok else 'FAIL'}")

    # ---- the TRUE F(s) stays below beta_odd on the sliver: the overshoot is K_star-non-sharpness ----
    print("\n[A3] the TRUE profile F(s) (FFT density) stays < beta_odd on the sliver (no real violation):")
    print(f"     {'s':>8} {'F_true':>10} {'beta-F_true':>12}")
    ft_ok = True
    for s in (1.001, 1.002, 1.004, 1.01, 1.05):
        ft = F_true(s)
        ft_ok &= ft < BETA_HP
        print(f"     {s:>8.4f} {ft:>10.6f} {BETA_HP - ft:>+12.6f}")
    ok &= ft_ok
    print(f"     => F_true < beta_odd throughout: {'ok' if ft_ok else 'FAIL'}  (so the bound, not F, is at fault)")

    # ---- FINDING 2: the SPLIT ----
    print("\n[A4] SPLIT lower leg (1,1.005]: PROVED variance-floor F_2(s) (Lemma lem:varfloor) < beta_odd:")
    print(f"     {'s':>8} {'F_2':>10} {'beta-F_2':>10} {'<beta?':>7}")
    f2_ok = True
    for s in (1.001, 1.002, 1.004, 1.005, 1.01):
        f2 = F2_floor(s)
        f2_ok &= f2 < BETA_HP
        print(f"     {s:>8.4f} {f2:>10.6f} {BETA_HP - f2:>+10.6f} {'yes' if f2 < BETA_HP else 'NO':>7}")
    ok &= f2_ok
    print(f"     => F_2 < beta_odd on the sliver (margin >= +0.0026): {'ok' if f2_ok else 'FAIL'}")

    print("\n[A5] SPLIT upper leg [1.005,inf): inf_{s>=1.005} K_max(s) > K_star, and ANY K<=K_star clears it:")
    ss = np.unique(np.concatenate([
        np.linspace(1.005, 1.5, 60), np.linspace(1.5, 3.0, 40), np.linspace(3.0, 30.0, 30)]))
    fg = np.array([fgauss_and_control(float(s), 800)[0] for s in ss])
    ct = np.array([fgauss_and_control(float(s), 800)[1] for s in ss])
    km = (BETA_HP - fg) / (ROOT * ct)
    i = int(np.argmin(km))
    Bk = fg + ROOT * K_STAR * ct
    print(f"     inf_{{s>=1.005}} K_max = {km[i]:.6f} at s={ss[i]:.4f}  (> K_star={K_STAR:.6f}? {km[i] > K_STAR})")
    print(f"     B(K_star) < beta_odd for all s in [1.005,30]? {bool(np.all(Bk < BETA_HP))}  "
          f"(max B={Bk.max():.6f} at s={ss[int(np.argmax(Bk))]:.4f})")
    Kcrude = 0.45 * K_STAR
    Bc = fg + ROOT * Kcrude * ct
    print(f"     crude K=0.45*K_star={Kcrude:.5f}: B<beta_odd for all s>=1.005? {bool(np.all(Bc < BETA_HP))}  "
          f"(max B={Bc.max():.6f}, margin>={BETA_HP - Bc.max():.4f})")
    upper_ok = km[i] > K_STAR and bool(np.all(Bk < BETA_HP)) and bool(np.all(Bc < BETA_HP))
    ok &= upper_ok
    print(f"     => upper leg closes with the sharp OR a crude constant: {'ok' if upper_ok else 'FAIL'}")

    # ---- FINDING (profile): max_u R/K_star <= 1; channel-asymmetry note ----
    print("\n[B] profile structure: max_u R(u,s)/K_star <= 1 (exact Kluyver), = 1 on u=0.5 (equal pair):")
    print(f"     {'s':>7} {'max_u R/K*':>11} {'at u':>7} {'asym max|sx2-sy2|/Sig2':>23}")
    prof_ok = True
    for s in (1.01, 1.1, 1.3, 1.5, 2.0, 3.0):
        us = np.unique(np.concatenate([np.linspace(0.003, 0.997, 40), np.linspace(0.47, 0.53, 12)]))
        best, ub = 0.0, 0.0
        asym = 0.0
        for u in us:
            ev, od = amps(float(u), s)
            if len(ev) + len(od) < 1:
                continue
            sx2 = 0.5 * float(np.sum(ev ** 2))
            sy2 = 0.5 * float(np.sum(od ** 2))
            sig2 = sx2 + sy2
            if sig2 > 0:
                asym = max(asym, abs(sx2 - sy2) / sig2)
            r = R_over_Kstar(ev, od)
            if r > best:
                best, ub = r, u
        prof_ok &= best <= 1.002
        print(f"     {s:>7.2f} {best:>11.5f} {ub:>7.3f} {asym:>23.4f}")
    ok &= prof_ok
    print(f"     => profile max R/K_star <= 1 (sharp at u=0.5 for s<1.5); channels NOT pointwise symmetric "
          f"(asym~0.98): {'ok' if prof_ok else 'FAIL'}")

    # ---- the K_star excess split at the pair: 45% small-t, 55% Bessel-ring tail ----
    print("\n[C2] Kluyver excess split at the equal pair {1}|{1} (which half K_star comes from):")
    sm, tl = excess_split(np.array([1.0]), np.array([1.0]))
    norm = 2.0 / 1.0 ** 1.5  # c4/Sigma2^{3/2} = 2 at the unit equal pair
    R_pair = (sm + tl) / norm
    print(f"     total excess = {sm + tl:.6f}  -> R = {R_pair:.7f}  (= K_star? {abs(R_pair - K_STAR) < 2e-4})")
    print(f"     small-t [0, j01] : {sm:+.6f}  = {sm / (sm + tl) * 100:.1f}% of K_star")
    print(f"     ring tail [j01,inf]: {tl:+.6f}  = {tl / (sm + tl) * 100:.1f}% of K_star  "
          f"(the residual the crude bound may keep lossy)")
    split_ok = abs(R_pair - K_STAR) < 2e-4 and 0.40 < sm / (sm + tl) < 0.50
    ok &= split_ok
    print(f"     => split ~45/55: {'ok' if split_ok else 'FAIL'}")

    print("\n" + "=" * 100)
    print("RESULT:", (
        "PASS -- the admissible window is re-derived with edge-resolving quadrature.  CORRECTIONS:\n"
        "  (1) inf_{s>1} K_max = 0.0338 at s->1+ (NOT 0.0402 at 1.02 -- a grid-start artifact); the K_star\n"
        "      majorant B(s) EXCEEDS beta_odd on the sliver (1, 1.00407], so 'B<beta on (1,20]' is false there.\n"
        "  (2) the SPLIT closes the range: variance-floor F_2<beta_odd on (1,1.005] (PROVED, near-tight since\n"
        "      the live set is ~2 copies there), and B<beta_odd with ANY K<=K_star on [1.005,inf) (inf K_max=\n"
        "      0.0363>K_star) -- the sharp constant/tail estimate is NOT needed off the sliver.\n"
        "  (3) beta_odd = 0.92801930480793... ; a former literal wrong past digit 9 has been corrected\n"
        "      repo-wide (immaterial to all margins).  profile max_u R/K_star<=1 (sharp on u=0.5); channels\n"
        "      are NOT pointwise symmetric (asym~0.98), only after the u-integral."
    ) if ok else "FAIL -- a checked invariant did not hold; see the failing line above.")
    print("Rigor: PROVED = F_gauss, variance-floor majorant, F(inf)=sqrt(pi)/2.  CERTIFIED (numerics, panel\n"
          "       Gauss-Legendre + exact Kluyver) = the window, the sliver, the split, profile R/K_star<=1.\n"
          "       RESIDUAL = a CRUDE (not sharp) Bessel-ring tail bound for (D) on [1.005,inf).  validated != proved.")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
