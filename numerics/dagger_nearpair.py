# /// script
# requires-python = ">=3.11"
# dependencies = ["python-flint>=0.8"]
# ///
r"""RIGOROUS certificate for the near-pair basin of the three-copy (2-vs-1) balanced Gaussian-excess
inequality (dagger).  Together with the singleton band proved by dagger_rigorous.py this covers the
whole three-copy family; configurations with more than three copies are not addressed here.

==========================================================================================================
THE CLAIM.  Two even copies a1, a3 and one odd copy b, X = a1 cos P1 + a3 cos P3, Y = b cos Q (independent
uniform phases), Z = (X, Y), G the matched Gaussian, Sigma_2 = (a1^2 + a3^2 + b^2)/2,
    f = K_star (a1^4 + a3^4 + b^4) / Sigma_2^{3/2} - E|Z| + E|G|      (1-homogeneous).
With u, v, w = a1^2, a3^2, b^2 over Sigma_2 (u + v + w = 2), the band min(u, v, w) in (0.36, 0.7] is
dagger_rigorous.py; here  f >= 0  on the BASIN  min(u, v, w) <= 0.36, with equality only at the pair
(u, v, w) = (1, 0, 1).  The swap a1 <-> a3 leaves f unchanged, so u >= v; then the basin is
    A = {v <= 0.36, w >= 0.36}  in the slice b = 1, coordinates (x, y) = (a1, a3), 0 <= y <= x,
                                 41 y^2 <= 9 (x^2 + 1),  x^2 + y^2 <= 41/9;
    B = {w <= 0.36}             in the slice a1 = 1, coordinates (x, y) = (a3/a1, b/a1),
                                 0 <= x <= 1,  41 y^2 <= 9 (1 + x^2).

==========================================================================================================
THE REPRESENTATION (R1-R5; every identity is exact).
 R1 (odd-copy kernel).  F_b(X) = E_Q sqrt(X^2 + b^2 cos^2 Q) has F_b'' = E[b^2 c^2 (X^2 + b^2 c^2)^{-3/2}]
    (c = cos Q), whose Fourier transform is ghat(b t) with
        ghat(s) = E_Q[2 s |c| K1(s |c|)] = 2 - 2 s I1(s/2) K0(s/2)
    (int e^{itx} (x^2 + beta^2)^{-3/2} dx = 2|t| K1(beta|t|)/beta; d/dz of
    int_0^{pi/2} K0(2 z cos p) dp = (pi/2) I0(z) K0(z); the Wronskian I0 K1 + I1 K0 = 1/z).  x K1(x) decreases
    from 1, so ghat decreases from ghat(0) = 2, and 0 <= 2 - ghat <= 2.  With F_b(0) = 2b/pi:
        F_b(X) = 2b/pi + (1/pi) int_0^inf ghat(b t) (1 - cos tX) / t^2 dt.
 R2 (even copies).  E cos(tX) = J0(a1 t) J0(a3 t) and Tonelli give, in the slice b = 1,
        E|Z| = 2/pi + (1/pi) int_0^inf ghat(t) (1 - J0(a1 t) J0(a3 t)) / t^2 dt.
 R3 (closed form of the subtracted kernel).  For a > c >= 0:
        (1/pi) int_0^inf (2t/(1+t^2)) (1 - J0(at) J0(ct)) / t^2 dt = (2/pi) Lam,
        Lam = ln a + I0(c) K0(a) - ln 2 + gamma.
    (d/da of the left side is (2/pi) int J1(at) J0(ct)/(1+t^2) dt = (2/pi)(1/a - I0(c) K1(a)) by the
    Weber-Schafheitlin value int J1(at) J0(ct) dt = 1/a and d/da of int t J0(at) J0(ct)/(t^2+1) dt =
    I0(c) K0(a); the remaining function of c is constant (Riemann-Lebesgue as a -> inf) and equals
    gamma - ln 2 (a -> 0 at c = 0).)  Scaling t -> t/z and d/dz at z = 1 add the second kernel:
        (1/pi) int_0^inf (t/(1+t^2)^2) (1 - J0 J0) / t^2 dt = (1/pi)(Lam - Lam'/2),
        Lam' = 1 + c I1(c) K0(a) - a I0(c) K1(a).
    With ghat0(t) = 2t/(1+t^2) + 5t/(1+t^2)^2 and ghat1 = ghat - ghat0 = O(t^-5):
        E|Z| = 2/pi + Wc + Wn,  Wc = (1/pi)(7 Lam - (5/2) Lam'),
        Wn = (1/pi) int_0^inf ghat1(t) (1 - J0(a1 t) J0(a3 t)) / t^2 dt.
    Lam is the a > c branch: on boxes meeting the diagonal x = y the formula is an analytic function that
    agrees with f on the closed half x >= y, which is all the certificate uses.
 R4 (Gaussian).  E|G| = pi^{-1/2} int_0^{pi/2} sqrt(sx cos^2 + sy sin^2) with sx, sy the two variances
    times 2; E(m) = (pi/2) 2F1(-1/2, 1/2; 1; m).
 R5 (slice a1 = 1, uniform as b -> 0).  int_0^inf (2 - ghat(s))/s^2 ds = (pi/2) int F_1'' |X| dX = 2, so
        E|Z| = E|X| + (1/pi) int_0^inf (2 - ghat(b t)) J0(t) J0(a3 t) / t^2 dt,
        E|X| = (4/pi^2) (2 E(m) - (1 - m) K(m)),  m = a3^2   (E|cos P1 + c| averaged over c = a3 cos P3).

==========================================================================================================
LEMMAS FOR THE ENCLOSURES.
 L1 (Gauss-Legendre).  As dagger_rigorous L5: |int_{-1}^1 g - GL_m| <= (64/15) M rho^{2-2m}/(rho^2 - 1)
    for g analytic and |g| <= M on the Bernstein ellipse E_rho; times the panel half-width.
 L2 (bounds on the ellipse, Re t >= xm > 0, |Im t| <= ym, |t| <= r).  |J_n(z)| <= e^{|Im z|}; |K_nu(z)| <=
    K_nu(Re z) and |I_nu(z)| <= I_nu(|z|); so |ghat(t)| <= (r/xm) ghat(xm), |2 - ghat(bt)| <=
    2 b r I1(b r/2) K0(b xm/2); |1 + t^2| >= max(xm, 1 - ym)^2; |(1 - J0(x))/x^2| <= (I0(|x|) - 1)/|x|^2.
 L3 (tails).  For t >= T >= 12: |ghat1(t)| <= C_g t^-5 with C_g = (4/pi)(3/8)(4/3)^{5/2} int x^5 K1 + 10
    + T^5 e(T), from ghat(t) = (4/pi) int_0^t x K1(x) (t^2 - x^2)^{-1/2} dx (Taylor remainder of
    (1 - y)^{-1/2} on y <= 1/4; e(T) the exponentially small pieces x >= t/2); |J_nu(x)| <= D/sqrt(x) for
    x >= x0 with D^2 = x0 (J_nu^2 + Y_nu^2)(x0) (nu >= 1) or 2/pi (nu = 0) (Watson 13.74).
 L4 (near t = 0).  On a real ball [0, eps]: ghat in [ghat(eps), 2] (monotone); for R5,
    |int_0^eps (2 - ghat(bt)) J0 J0 / t^2| <= b I0(b eps/2) sqrt(pi b eps) (I1(z) <= z I0(z)/2,
    K0(s) <= sqrt(pi/(2s))).
 L5 (convexity).  If Hess f >= 0 on the convex rectangle U containing p = (1, 0) of the slice b = 1, then
    f >= f(p) + grad f(p).(x - p) = 0 on U: f(p) = 0 defines K_star, f_y(p) = 0 by evenness in a3, and
    f_x(p) = 0 by the swap symmetry f(x, 0) = x f(1/x, 0) of the two-copy line.  Hess f is checked on
    U cap {y >= 0}; evenness in y extends it to the symmetric rectangle.
 L6 (Taylor model).  On a box [xc +- hx] x [yc +- hy]: a point jet at the centre (exact Arb jets of every
    term, the Wn coefficients by L1-L4) and a ball jet one order higher over the box bound each quantity by
    its Taylor polynomial (linear part exact at the centre) plus |coefficients over the box| times h^2.
 L7 (monotonicity).  E|Z| and E|G| are nondecreasing in every amplitude (d E|Z| / d a3 =
    E[cos P3 X/|Z|] >= 0: X/|Z| is nondecreasing in cos P3, which is symmetric), so on a box
    f >= K_lo P_min - E|Z|(upper corner) + E|G|(lower corner).

SCOPE.  Assumed: the classical identities quoted in R1, R3, R5, Watson 13.74, L1, and Arb.  K_star comes
from dagger_rigorous.kstar_enclosure (L10 there).  Coverage of A and B by the certified boxes is checked
in exact rational arithmetic by verify_dagger_nearpair.py.
"""

from __future__ import annotations

import argparse
import functools
import hashlib
import json
import math
import platform
import time
from collections.abc import Callable
from concurrent.futures import FIRST_COMPLETED, Future, ProcessPoolExecutor, wait
from dataclasses import asdict, dataclass
from fractions import Fraction
from math import comb, factorial
from pathlib import Path

import flint
from dagger_rigorous import kstar_enclosure
from flint import acb, arb, ctx

SOURCE = Path(__file__).resolve()
ROOT = SOURCE.parents[1]
OUT = ROOT / "numerics" / "dagger_nearpair.jsonl"
SCHEMA = "beta-odd-dagger-nearpair-v1"


@dataclass(frozen=True)
class Params:
    prec: int = 96
    # region A (slice b = 1): Gauss-Legendre in t
    t_cut: int = 40
    t0: str = "1"
    graded: int = 20  # panels [t0 2^-j-1, t0 2^-j]; [0, t0 2^-graded] by a direct enclosure
    width: str = "1"
    m_main: int = 16
    m_graded: int = 8
    # region B (slice a1 = 1, R5)
    t_cut_b: int = 200
    width_b: str = "1.5"
    m_b: int = 14
    graded_b: int = 40
    # the convexity rectangle U = [u_x0, u_x1] x [0, u_y1] (dyadic) and the lattices
    u_x0: str = "13/16"
    u_x1: str = "3/2"
    u_y1: str = "7/16"
    hess_h: str = "1/64"
    value_side: str = "1/16"
    b_side: str = "1/32"
    max_depth: int = 4


DEFAULT = Params()
RHO = tuple(arb(s) for s in ("1.5", "2", "2.5", "3", "4", "5", "5.5", "6", "8", "10", "12"))


def _up(x: arb) -> arb:
    return arb(x.upper())


def _lo(x: arb) -> arb:
    return arb(x.lower())


def _abs_up(x: arb) -> arb:
    return _up(abs(x))


def _q(x: Fraction) -> arb:
    """Exact arb of a dyadic Fraction."""
    if x.denominator & (x.denominator - 1):
        raise ValueError(f"{x} is not dyadic")
    return arb(x.numerator) / x.denominator


def source_sha256() -> str:
    return hashlib.sha256(SOURCE.read_bytes()).hexdigest()


# ==========================================================================================================
# Bivariate jets: truncated Taylor coefficients c[(i, j)] of a function of (x, y), as Arb balls.
# ==========================================================================================================
class Jet:
    __slots__ = ("N", "c")

    def __init__(self, N: int, c: dict[tuple[int, int], arb]):
        self.N = N
        self.c = {k: v for k, v in c.items() if k[0] + k[1] <= N}

    def get(self, i: int, j: int) -> arb:
        return self.c.get((i, j), arb(0))

    def _lift(self, o: Jet | arb | int | float) -> Jet:
        return o if isinstance(o, Jet) else Jet(self.N, {(0, 0): arb(o)})

    def __add__(self, o: Jet | arb | int | float) -> Jet:
        o = self._lift(o)
        c = dict(self.c)
        for k, v in o.c.items():
            c[k] = c.get(k, arb(0)) + v
        return Jet(self.N, c)

    __radd__ = __add__

    def __neg__(self) -> Jet:
        return Jet(self.N, {k: -v for k, v in self.c.items()})

    def __sub__(self, o: Jet | arb | int | float) -> Jet:
        return self + (-self._lift(o))

    def __rsub__(self, o: Jet | arb | int | float) -> Jet:
        return self._lift(o) - self

    def __mul__(self, o: Jet | arb | int | float) -> Jet:
        if not isinstance(o, Jet):
            return Jet(self.N, {k: v * o for k, v in self.c.items()})
        c: dict[tuple[int, int], arb] = {}
        for (i, j), v in self.c.items():
            for (i2, j2), w in o.c.items():
                if i + j + i2 + j2 <= self.N:
                    c[(i + i2, j + j2)] = c.get((i + i2, j + j2), arb(0)) + v * w
        return Jet(self.N, c)

    __rmul__ = __mul__

    def compose(self, uni: list[arb]) -> Jet:
        """g(self) for g with Taylor coefficients uni at self(0); the constant is dropped EXACTLY (a ball
        minus itself is not zero)."""
        d = Jet(self.N, {k: v for k, v in self.c.items() if k != (0, 0)})
        out = Jet(self.N, {(0, 0): uni[0]})
        p = Jet(self.N, {(0, 0): arb(1)})
        for k in range(1, self.N + 1):
            p = p * d
            out = out + p * uni[k]
        return out


def _var_jets(x: arb, y: arb, N: int) -> tuple[Jet, Jet]:
    return Jet(N, {(0, 0): x, (1, 0): arb(1)}), Jet(N, {(0, 0): y, (0, 1): arb(1)})


def _uni(cs: list[arb], var: int, N: int) -> Jet:
    return Jet(N, {((k, 0) if var == 0 else (0, k)): cs[k] for k in range(N + 1)})


# ==========================================================================================================
# Closed-form terms (R3, R4, P) as jets in the slice b = 1.
# ==========================================================================================================
def _bessel_deriv(vals: list[arb], k: int, sign: int) -> arb:
    """Z0^(k) = (sign/2)^k sum_j C(k,j) Z_|k-2j| for Z = I (sign +1) or K (sign -1)."""
    s = sum((comb(k, j) * vals[abs(k - 2 * j)] for j in range(k + 1)), arb(0))
    return s * (arb(sign) / 2) ** k


def wc_jet(a1: arb, a3: arb, N: int) -> Jet:
    """Wc = (1/pi)(7 Lam - (5/2) Lam') (R3)."""
    bi = [a3.bessel_i(n) for n in range(N + 2)]
    bk = [a1.bessel_k(n) for n in range(N + 2)]
    dI = [_bessel_deriv(bi, k, 1) for k in range(N + 2)]
    dK = [_bessel_deriv(bk, k, -1) for k in range(N + 2)]
    fac = [factorial(k) for k in range(N + 1)]
    r = [dI[k] / fac[k] for k in range(N + 1)]
    p = [dK[k] / fac[k] for k in range(N + 1)]
    s = [(a3 * dI[k + 1] + k * dI[k]) / fac[k] for k in range(N + 1)]  # y I1(y) = y I0'(y)
    q = [-(a1 * dK[k + 1] + k * dK[k]) / fac[k] for k in range(N + 1)]  # x K1(x) = -x K0'(x)
    ln = [a1.log()] + [arb(-1) ** (k + 1) / (k * a1**k) for k in range(1, N + 1)]
    R, P, S, Q, L = (_uni(r, 1, N), _uni(p, 0, N), _uni(s, 1, N), _uni(q, 0, N), _uni(ln, 0, N))
    lam = L + R * P - arb(2).log() + arb.const_euler()
    lamp = 1 + S * P - R * Q
    return (lam * 7 - lamp * arb(2.5)) * (1 / arb.pi())


def _ellip_e_taylor(m0: arb, N: int) -> list[arb]:
    """Taylor coefficients of E(m) = (pi/2) 2F1(-1/2, 1/2; 1; m) at m0 < 1."""
    a, b, c = arb(-0.5), arb(0.5), arb(1)
    out, poch = [], arb(1)
    for k in range(N + 1):
        if k:
            poch *= (a + k - 1) * (b + k - 1) / ((c + k - 1) * k)
        out.append(arb.pi() / 2 * poch * acb(m0).hypgeom_2f1(a + k, b + k, c + k).real)
    return out


def eg_jet(a1: arb, a3: arb, N: int) -> Jet:
    """E|G| = E(1 - a1^2 - a3^2)/sqrt(pi) in the slice b = 1 (R4)."""
    X, Y = _var_jets(a1, a3, N)
    m = 1 - X * X - Y * Y
    return m.compose(_ellip_e_taylor(m.get(0, 0), N)) * (1 / arb.pi().sqrt())


def p_jet(a1: arb, a3: arb, N: int) -> Jet:
    """(a1^4 + a3^4 + 1) / S2^{3/2},  S2 = (a1^2 + a3^2 + 1)/2."""
    X, Y = _var_jets(a1, a3, N)
    xs, ys = X * X, Y * Y
    s2 = (xs + ys + 1) * arb(0.5)
    s0 = s2.get(0, 0)
    pw, bc = [], arb(1)
    for k in range(N + 1):
        if k:
            bc *= (arb(-1.5) - (k - 1)) / k
        pw.append(bc * s0 ** (arb(-1.5) - k))
    return (xs * xs + ys * ys + 1) * s2.compose(pw)


# ==========================================================================================================
# The odd-copy kernel (R1) and the Wn integrals (L1-L4).
# ==========================================================================================================
def k0(z: arb) -> arb:
    """K0(z) for z > 0.  Arb's K0 loses ~2z/ln 2 bits for moderate z, and far more on a ball of nonzero
    radius (relative radius 1e4 at an exact point, 1e9 on a 1e-29 ball, 96 bits, z in [2, 50]).  So
    evaluate at the exact midpoint with 3z + 32 extra bits, and add the radius by the mean-value theorem,
    |K0'| = K1 <= K_{3/2}(x) = sqrt(pi/(2x)) e^{-x} (1 + 1/x) (K_nu increases with nu >= 0)."""
    zm, zl = arb(z.mid()), _lo(z)
    if not (zl > 0):
        raise ValueError("k0 needs z > 0")
    old = ctx.prec
    ctx.prec = old + 3 * math.ceil(float(z.upper())) + 32
    try:
        v = zm.bessel_k(0)
    finally:
        ctx.prec = old
    v = +v
    if z.rad() > 0:
        v += arb(0, _up(z.rad() * (arb.pi() / (2 * zl)).sqrt() * (-zl).exp() * (1 + 1 / zl)))
    return v


def i1k0(z: arb) -> arb:
    return z.bessel_i(1) * k0(z)


def ghat_pt(t: arb) -> arb:
    """ghat(t) = 2 - 2 t I1(t/2) K0(t/2), t > 0."""
    return 2 - 2 * t * i1k0(t / 2)


def ghat_real(t: arb) -> arb:
    """ghat on a real ball t >= 0 that may touch 0 (L4: monotone, ghat(0) = 2)."""
    if t.lower() > 0:
        return ghat_pt(t)
    hi = _up(t)
    return arb(2).union(ghat_pt(hi)) if hi > 0 else arb(2)


def ghat0(t: arb) -> arb:
    d = 1 + t * t
    return 2 * t / d + 5 * t / (d * d)


def _j0_deriv(k: int, J: list[arb]) -> arb:
    """J0^(k) = 2^-k sum_j (-1)^j C(k,j) J_{2j-k},  J_{-n} = (-1)^n J_n."""
    s = arb(0)
    for j in range(k + 1):
        nu = 2 * j - k
        term = comb(k, j) * J[abs(nu)]
        if nu < 0 and nu % 2:
            term = -term
        s += -term if j % 2 else term
    return s / arb(2) ** k


def _phi(x: arb) -> arb:
    """(1 - J0(x))/x^2 = (1/4) 1F2(1; 2, 2; -x^2/4), entire."""
    return acb(-x * x / 4).hypgeom([1], [2, 2]).real / 4


def _core(m: int, n: int, t: arb, J1: list[arb], J3: list[arb], a1: arb, a3: arb, phis) -> arb:
    """d^m_{a1} d^n_{a3} of (1 - J0(a1 t) J0(a3 t))/t^2, written without divisions by t."""
    k = m + n
    if k >= 2:
        return -(t ** (k - 2)) * _j0_deriv(m, J1) * _j0_deriv(n, J3)
    if k == 1:
        if m == 1:
            return a1 * (J1[0] + J1[2]) / 2 * J3[0]
        return a3 * (J3[0] + J3[2]) / 2 * J1[0]
    if phis is not None:
        return a1 * a1 * phis[0] + J1[0] * a3 * a3 * phis[1]
    return (1 - J1[0] * J3[0]) / (t * t)


@functools.cache
def _gl_rule(m: int, prec: int) -> tuple[tuple[arb, arb], ...]:
    old = ctx.prec
    ctx.prec = prec
    try:
        return tuple(arb.legendre_p_root(m, k, weight=True) for k in range(m))
    finally:
        ctx.prec = old


def _panels(t0: arb, graded: int, width: arb, T: int, m_main: int, m_graded: int):
    out = [(t0 / arb(2) ** (j + 1), t0 / arb(2) ** j, m_graded) for j in range(graded)]
    a = t0
    while a < T:
        b = (a + width).min(arb(T))
        out.append((a, b, m_main))
        a = b
    return out


def _ellipse(a: arb, b: arb, rho: arb) -> tuple[arb, arb, arb] | None:
    """(xm, ym, r): min Re, max |Im|, max |t| on the Bernstein ellipse E_rho of [a, b]; None if xm <= 0."""
    c, w = (a + b) / 2, (b - a) / 2
    A, B = w * (rho + 1 / rho) / 2, w * (rho - 1 / rho) / 2
    xm = c - A
    return (xm, B, c + A) if xm > 0 else None


@functools.cache
def _panels_a(prm: Params, prec: int) -> tuple[tuple[tuple, tuple], ...]:
    """Region-A quadrature data independent of (a1, a3): per panel the nodes (t, weight, ghat1(t)) and,
    per admissible ellipse, (L1 factor, xm, ym, r, sup |ghat1|) (L2)."""
    old = ctx.prec
    ctx.prec = prec
    try:
        out = []
        plist = _panels(
            arb(prm.t0), prm.graded, arb(prm.width), prm.t_cut, prm.m_main, prm.m_graded
        )
        for a, b, mm in plist:
            c, w = (a + b) / 2, (b - a) / 2
            nodes = []
            for x, wt in _gl_rule(mm, prec):
                t = c + w * x
                nodes.append((t, w * wt, ghat_pt(t) - ghat0(t)))
            geoms = []
            for rho in RHO:
                g = _ellipse(a, b, rho)
                if g is None:
                    continue
                xm, ym, r = g
                d = xm.max(1 - ym)
                mg = r / xm * ghat_pt(xm) + 2 * r / (d * d) + 5 * r / (d * d * d * d)
                geoms.append(
                    (w * arb(64) / 15 * rho ** (2 - 2 * mm) / (rho * rho - 1), xm, ym, r, mg)
                )
            if not geoms:
                raise ArithmeticError(f"no admissible ellipse for the panel [{a}, {b}]")
            out.append((tuple(nodes), tuple(geoms)))
        return tuple(out)
    finally:
        ctx.prec = old


def _sup_a(
    geom: tuple[arb, arb, arb, arb], a1h: arb, a3h: arb, orders
) -> dict[tuple[int, int], arb]:
    """L2: bounds of |ghat1 * core_mn| on the ellipse (region A); geom = (xm, ym, r, sup |ghat1|)."""
    xm, ym, r, mg = geom
    ex = ((a1h + a3h) * ym).exp()
    out = {}
    for m, n in orders:
        k = m + n
        if k >= 2:
            mc = r ** (k - 2) * ex
        elif k == 1:
            mc = (a1h if m == 1 else a3h) * ex
        else:
            m1 = (1 + ex) / (xm * xm)
            m2 = ((a1h * r).bessel_i(0) - 1 + (a1h * ym).exp() * ((a3h * r).bessel_i(0) - 1)) / (
                r * r
            )
            mc = m1.min(m2)
        out[(m, n)] = mg * mc
    return out


def tail_cg(T: arb) -> arb:
    """L3: C_g with |ghat1(t)| <= C_g t^-5 for all t >= T (T >= 12)."""
    if not (T >= 12):
        raise ValueError("L3 needs T >= 12")
    main = arb(4) / arb.pi() * arb(3) / 8 * (arb(4) / 3) ** arb(2.5) * (arb(45) * arb.pi() / 2)
    e2 = 2 * T / 3 * (T / 2).bessel_k(1)
    e1 = 3 * T / (2 * arb.pi()) * (T / 2).bessel_k(1) + 6 / arb.pi() * T.bessel_k(1) / (
        1 - arb(2.5) / T
    )
    return _up(main + 10 + T**5 * e1.max(e2))


def bessel_decay(x0: arb, nus: set[int]) -> arb:
    """L3: D with |J_nu(x)| <= D / sqrt(x) for every x >= x0 > 0 and nu in nus."""
    out = (2 / arb.pi()).sqrt()
    for nu in nus:
        if nu >= 1:
            out = out.max((x0 * (x0.bessel_j(nu) ** 2 + x0.bessel_y(nu) ** 2)).sqrt())
    return _up(out)


def _tail_a(m: int, n: int, a1l: arb, T: arb, cg: arb) -> arb:
    """L3: |int_T^inf ghat1 core_mn dt| (without the 1/pi of Wn), a1 >= a1l > 0."""
    k = m + n
    if k == 0:
        return _up(cg * 2 / (6 * T**6))
    D1 = bessel_decay(a1l * T, {abs(m - 2 * j) for j in range(m + 1)} | {0, 1, 2})
    if k == 1:
        return _up(cg * D1 / (a1l.sqrt() * arb(5.5) * T ** arb(5.5)))
    return _up(cg * D1 / a1l.sqrt() * T ** (k - arb(6.5)) / (arb(6.5) - k))


def wn_derivs(
    a1: arb, a3: arb, orders: list[tuple[int, int]], prm: Params
) -> dict[tuple[int, int], arb]:
    """{(m, n): enclosure of d^m_{a1} d^n_{a3} Wn} for a1 >= a1l > 0 and a3 >= 0 (point or ball)."""
    numax = max(max(m, n) for m, n in orders) + 2
    a1h, a3h, a1l = _up(abs(a1)), _up(abs(a3)), _lo(a1)
    if not (a1l > 0):
        raise ValueError("Wn needs a1 > 0")
    acc = {o: arb(0) for o in orders}
    err = {o: arb(0) for o in orders}
    for nodes, geoms in _panels_a(prm, ctx.prec):
        for t, wt, g1 in nodes:
            J1 = [(a1 * t).bessel_j(nu) for nu in range(numax + 1)]
            J3 = [(a3 * t).bessel_j(nu) for nu in range(numax + 1)]
            phis = (_phi(a1 * t), _phi(a3 * t)) if (0, 0) in acc and t < 2 else None
            for o in orders:
                acc[o] += wt * g1 * _core(o[0], o[1], t, J1, J3, a1, a3, phis)
        best: dict[tuple[int, int], arb] | None = None
        for fac, *geom in geoms:
            sup = _sup_a(tuple(geom), a1h, a3h, orders)
            cand = {o: _up(fac * sup[o]) for o in orders}
            best = cand if best is None else {o: best[o].min(cand[o]) for o in orders}
        assert best is not None
        for o in orders:
            err[o] += best[o]
    eps = arb(prm.t0) / arb(2) ** prm.graded
    tb = arb(0).union(eps)
    g1 = ghat_real(tb) - ghat0(tb)
    J1 = [(a1 * tb).bessel_j(nu) for nu in range(numax + 1)]
    J3 = [(a3 * tb).bessel_j(nu) for nu in range(numax + 1)]
    phis = (_phi(a1 * tb), _phi(a3 * tb))
    T = arb(prm.t_cut)
    cg = tail_cg(T)
    out = {}
    for o in orders:
        head = eps * g1 * _core(o[0], o[1], tb, J1, J3, a1, a3, phis)
        tot = err[o] + _tail_a(o[0], o[1], a1l, T, cg)
        out[o] = (acc[o] + head) / arb.pi() + arb(0, _up(tot / arb.pi()))
    return out


def f_jet(a1: arb, a3: arb, N: int, K: arb, prm: Params) -> Jet:
    """Jet of f = K P - E|Z| + E|G| (slice b = 1) at a point or over a ball, to order N."""
    orders = [(m, n) for m in range(N + 1) for n in range(N + 1 - m)]
    d = wn_derivs(a1, a3, orders, prm)
    wn = Jet(N, {(m, n): v / (factorial(m) * factorial(n)) for (m, n), v in d.items()})
    ez = wc_jet(a1, a3, N) + wn + 2 / arb.pi()
    return p_jet(a1, a3, N) * K - ez + eg_jet(a1, a3, N)


def ez_a(x: arb, y: arb, prm: Params) -> arb:
    """E|Z| at a point of the slice b = 1 (R3 needs the larger even amplitude first)."""
    a1, a3 = (x, y) if x >= y else (y, x)
    wn = wn_derivs(a1, a3, [(0, 0)], prm)[(0, 0)]
    return 2 / arb.pi() + wc_jet(a1, a3, 0).get(0, 0) + wn


# ==========================================================================================================
# Region B: the slice a1 = 1 with R5.
# ==========================================================================================================
def abs_x(c: arb) -> arb:
    """E|cos P1 + c cos P3| = (4/pi^2)(2 E(m) - (1 - m) K(m)), m = c^2, 0 <= c <= 1."""
    if c == 1:
        return 8 / arb.pi() ** 2
    m = acb(c * c)
    return 4 / arb.pi() ** 2 * (2 * m.elliptic_e() - (1 - m) * m.elliptic_k()).real


@functools.cache
def _nodes_b(prm: Params, prec: int) -> tuple[tuple[arb, arb, arb], ...]:
    """(t, weight, J0(t)) for the region-B quadrature; shared by every point evaluation."""
    old = ctx.prec
    ctx.prec = prec
    try:
        out = []
        for a, b, mm in _panels(
            arb(prm.t0), prm.graded_b, arb(prm.width_b), prm.t_cut_b, prm.m_b, 8
        ):
            c, w = (a + b) / 2, (b - a) / 2
            for x, wt in _gl_rule(mm, prec):
                t = c + w * x
                out.append((t, w * wt, t.bessel_j(0)))
        return tuple(out)
    finally:
        ctx.prec = old


@functools.cache
def _err_b_geoms(
    prm: Params, prec: int
) -> tuple[tuple[arb, int, tuple[tuple[arb, ...], ...]], ...]:
    old = ctx.prec
    ctx.prec = prec
    try:
        out = []
        for a, b, mm in _panels(
            arb(prm.t0), prm.graded_b, arb(prm.width_b), prm.t_cut_b, prm.m_b, 8
        ):
            w = (b - a) / 2
            geoms = []
            for rho in RHO:
                g = _ellipse(a, b, rho)
                if g is not None:
                    geoms.append((rho, *g))
            out.append((w, mm, tuple(geoms)))
        return tuple(out)
    finally:
        ctx.prec = old


def ez_b(c: arb, beta: arb, prm: Params) -> arb:
    """E|Z| at the point (a1, a3, b) = (1, c, beta), 0 <= c <= 1, beta >= 0 (R5, L1, L2, L4)."""
    base = abs_x(c)
    if beta == 0:
        return base
    acc = arb(0)
    for t, wt, j0 in _nodes_b(prm, ctx.prec):
        acc += wt * 2 * beta * i1k0(beta * t / 2) * j0 * (c * t).bessel_j(0) / t
    err = arb(0)
    for w, mm, geoms in _err_b_geoms(prm, ctx.prec):
        best = None
        for rho, xm, ym, r in geoms:
            sup = (
                2
                * beta
                * (beta * r / 2).bessel_i(1)
                * k0(beta * xm / 2)
                * ((1 + c) * ym).exp()
                / xm
            )
            bound = _up(w * arb(64) / 15 * sup * rho ** (2 - 2 * mm) / (rho * rho - 1))
            best = bound if best is None else best.min(bound)
        if best is None:
            raise ArithmeticError("no admissible ellipse in region B")
        err += best
    eps = arb(prm.t0) / arb(2) ** prm.graded_b
    head = beta * (beta * eps / 2).bessel_i(0) * (arb.pi() * beta * eps).sqrt()
    T = arb(prm.t_cut_b)
    D0 = (2 / arb.pi()).sqrt()
    tail = 4 * D0 / (3 * T ** arb(1.5))
    if c > 0:
        tail = tail.min(D0 * D0 / (c.sqrt() * T * T))
    return base + acc / arb.pi() + arb(0, _up((err + head + tail) / arb.pi()))


def eg_b(c: arb, beta: arb) -> arb:
    """E|G| at (1, c, beta): sqrt((1 + c^2)/pi) E(1 - beta^2/(1 + c^2)) (R4)."""
    s = 1 + c * c
    m = acb(1 - beta * beta / s)
    e = arb.pi() / 2 * m.hypgeom_2f1(arb(-0.5), arb(0.5), arb(1)).real
    return (s / arb.pi()).sqrt() * e


def p_b(c: arb, beta: arb) -> arb:
    s2 = (1 + c * c + beta * beta) / 2
    return (1 + c**4 + beta**4) / s2 ** arb(1.5)


def p_a(x: arb, y: arb) -> arb:
    s2 = (1 + x * x + y * y) / 2
    return (1 + x**4 + y**4) / s2 ** arb(1.5)


def _p_min(
    pfun: Callable[[arb, arb], arb], x0: Fraction, x1: Fraction, y0: Fraction, y1: Fraction
) -> arb:
    """Lower bound of P on a box: ball evaluation on a 4 x 4 subdivision."""
    out = None
    for i in range(4):
        for j in range(4):
            xa, xb = x0 + (x1 - x0) * i / 4, x0 + (x1 - x0) * (i + 1) / 4
            ya, yb = y0 + (y1 - y0) * j / 4, y0 + (y1 - y0) * (j + 1) / 4
            v = _lo(pfun(_q(xa).union(_q(xb)), _q(ya).union(_q(yb))))
            out = v if out is None else out.min(v)
    assert out is not None
    return out


# ==========================================================================================================
# Box certificates.
# ==========================================================================================================
def ball(x: arb, digits: int = 30) -> dict[str, str]:
    """Rigorous decimal record: unball(ball(x)) contains x.  The midpoint is printed to `digits` significant
    digits, its rounding error is added to the radius, and the radius is printed rounded up.  (Arb's
    radius=True printing keeps only the digits the radius certifies, which widens a wide ball by up to a
    unit of its last printed digit: enough to lose a certified margin.)"""
    mid = x.mid().str(digits, radius=False)
    r = _up(x.rad() + abs(arb(mid) - x.mid()))
    rec = {"mid": mid, "rad": (r * (1 + arb(2) ** -12)).str(6, radius=False) if r > 0 else "0"}
    if not unball(rec).contains(x):
        raise ArithmeticError(f"the decimal record {rec} does not contain {x}")
    return rec


def unball(rec: dict[str, str]) -> arb:
    return arb(f"[{rec['mid']} +/- {rec['rad']}]")


def hessian_box(xc: Fraction, yc: Fraction, h: Fraction, K: arb, prm: Params) -> dict:
    """L5 + L6: enclosures of f_xx, f_xy, f_yy over [xc +- h] x [yc +- h] (slice b = 1).  The verdict is
    decided from the stored decimal balls, so the checker reproduces it exactly."""
    X, Y, H = _q(xc), _q(yc), _q(h)
    fc = f_jet(X, Y, 3, K, prm)
    fb = f_jet(arb(X, H), arb(Y, H), 4, K, prm)
    c = fc.get

    def b(i: int, j: int) -> arb:
        return _abs_up(fb.get(i, j))

    hess = {"xx": 2 * c(2, 0), "xy": c(1, 1), "yy": 2 * c(0, 2)}
    grad = {
        "xx": (6 * c(3, 0), 2 * c(2, 1)),
        "xy": (2 * c(2, 1), 2 * c(1, 2)),
        "yy": (2 * c(1, 2), 6 * c(0, 3)),
    }
    quad = {
        "xx": (24 * b(4, 0), 6 * b(3, 1), 4 * b(2, 2)),
        "xy": (6 * b(3, 1), 4 * b(2, 2), 6 * b(1, 3)),
        "yy": (4 * b(2, 2), 6 * b(1, 3), 24 * b(0, 4)),
    }
    rec = {}
    for key, v in hess.items():
        g, qd = grad[key], quad[key]
        spread = (_abs_up(g[0]) + _abs_up(g[1])) * H + (qd[0] + 2 * qd[1] + qd[2]) * H * H / 2
        rec[f"f_{key}"] = ball(v + arb(0, _up(spread)))
    xx, xy, yy = (unball(rec[k]) for k in ("f_xx", "f_xy", "f_yy"))
    det_lo = _lo(xx) * _lo(yy) - _abs_up(xy) ** 2
    return {
        "kind": "hessian",
        "box": [str(xc), str(yc), str(h), str(h)],
        **rec,
        "det_lo": ball(det_lo),
        "proves": bool(xx > 0 and yy > 0 and det_lo > 0),
    }


def taylor_box(xc: Fraction, yc: Fraction, h: Fraction, K: arb, prm: Params) -> dict:
    """L6: f(c) - |grad f(c)| h - (1/2) h^T |Hess f(B)| h > 0 on [xc +- h] x [yc +- h] (slice b = 1)."""
    X, Y, H = _q(xc), _q(yc), _q(h)
    fc = f_jet(X, Y, 1, K, prm)
    fb = f_jet(arb(X, H), arb(Y, H), 2, K, prm)
    lin = (_abs_up(fc.get(1, 0)) + _abs_up(fc.get(0, 1))) * H
    quad = (_abs_up(fb.get(2, 0)) + _abs_up(fb.get(1, 1)) + _abs_up(fb.get(0, 2))) * H * H
    rec = {"f_c": ball(fc.get(0, 0)), "lin": ball(lin), "quad": ball(quad)}
    lower = _lo(unball(rec["f_c"])) - _up(unball(rec["lin"])) - _up(unball(rec["quad"]))
    return {
        "kind": "taylor",
        "box": [str(xc), str(yc), str(h), str(h)],
        **rec,
        "lower": ball(lower),
        "proves": bool(lower > 0),
    }


def mono_box(region: str, xc: Fraction, yc: Fraction, h: Fraction, K: arb, prm: Params) -> dict:
    """L7 on the box [xc +- h] x [yc +- h] of region A (slice b = 1) or B (slice a1 = 1)."""
    x0, x1, y0, y1 = xc - h, xc + h, yc - h, yc + h
    x0c, y0c = max(x0, Fraction(0)), max(y0, Fraction(0))
    if region == "A":
        pmin = _p_min(p_a, x0c, x1, y0c, y1)
        ez = ez_a(_q(x1), _q(y1), prm)
        eg = eg_jet(_q(x0c), _q(y0c), 0).get(0, 0)
    else:
        pmin = _p_min(p_b, x0c, min(x1, Fraction(1)), y0c, y1)
        ez = ez_b(_q(min(x1, Fraction(1))), _q(y1), prm)
        eg = eg_b(_q(x0c), _q(y0c))
    rec = {"p_min": ball(pmin), "ez_up": ball(ez), "eg_lo": ball(eg)}
    lower = (
        _lo(K) * _lo(unball(rec["p_min"])) - _up(unball(rec["ez_up"])) + _lo(unball(rec["eg_lo"]))
    )
    return {
        "kind": f"mono_{region}",
        "box": [str(xc), str(yc), str(h), str(h)],
        **rec,
        "lower": ball(lower),
        "proves": bool(lower > 0),
    }


# ==========================================================================================================
# Regions, lattices and the campaign.
# ==========================================================================================================
def meets_a(x0: Fraction, x1: Fraction, y0: Fraction, y1: Fraction) -> bool:
    """Necessary condition for the box to meet A (exact): some corner choice satisfies each constraint."""
    return (
        y0 <= x1
        and 41 * y0 * y0 <= 9 * (x1 * x1 + 1)
        and x0 * x0 + y0 * y0 <= Fraction(41, 9)
        and y1 >= 0
    )


def meets_b(x0: Fraction, x1: Fraction, y0: Fraction, y1: Fraction) -> bool:
    return x0 <= 1 and 41 * y0 * y0 <= 9 * (1 + x1 * x1) and y1 >= 0 and x1 >= 0


def rect_u(prm: Params) -> tuple[Fraction, Fraction, Fraction]:
    return Fraction(prm.u_x0), Fraction(prm.u_x1), Fraction(prm.u_y1)


def inside_u(x0: Fraction, x1: Fraction, y1: Fraction, prm: Params) -> bool:
    ux0, ux1, uy1 = rect_u(prm)
    return x0 >= ux0 and x1 <= ux1 and y1 <= uy1


def relevant(method: str, xc: Fraction, yc: Fraction, h: Fraction, prm: Params) -> bool:
    """Whether a box of this method may be needed: it meets its region (A outside U)."""
    x0, x1, y0, y1 = xc - h, xc + h, yc - h, yc + h
    if method == "A":
        return meets_a(x0, x1, y0, y1) and not inside_u(x0, x1, y1, prm)
    if method == "B":
        return meets_b(x0, x1, y0, y1)
    return True


def base_tasks(prm: Params) -> list[tuple[str, Fraction, Fraction, Fraction]]:
    """(method, xc, yc, h) for the U Hessian boxes, the region-A value boxes and the region-B boxes."""
    tasks: list[tuple[str, Fraction, Fraction, Fraction]] = []
    ux0, ux1, uy1 = rect_u(prm)
    hh = Fraction(prm.hess_h)
    x = ux0
    while x < ux1:
        y = Fraction(0)
        while y < uy1:
            tasks.append(("hessian", x + hh, y + hh, hh))
            y += 2 * hh
        x += 2 * hh
    s = Fraction(prm.value_side)
    xmax = Fraction(9, 4)
    i = 0
    while i * s < xmax:
        x0, x1 = i * s, (i + 1) * s
        j = 0
        while j * s <= x1:
            y0 = j * s
            if relevant("A", x0 + s / 2, y0 + s / 2, s / 2, prm):
                tasks.append(("A", x0 + s / 2, y0 + s / 2, s / 2))
            j += 1
        i += 1
    sb = Fraction(prm.b_side)
    i = 0
    while i * sb < 1:
        x0, x1 = i * sb, (i + 1) * sb
        j = 0
        while j * sb < Fraction(7, 10):
            y0 = j * sb
            if relevant("B", x0 + sb / 2, y0 + sb / 2, sb / 2, prm):
                tasks.append(("B", x0 + sb / 2, y0 + sb / 2, sb / 2))
            j += 1
        i += 1
    return tasks


_K: arb | None = None


def kstar() -> arb:
    """K_star as a ball (dagger_rigorous L10), computed once per process."""
    global _K
    if _K is None:
        k_lo, k_hi = kstar_enclosure()
        _K = k_lo.union(k_hi)
    return _K


def certify(task: tuple[str, str, str, str], prm: Params = DEFAULT) -> dict:
    method, xs, ys, hs = task
    xc, yc, h = Fraction(xs), Fraction(ys), Fraction(hs)
    ctx.prec = prm.prec
    K = kstar()
    t0 = time.time()
    if method == "hessian":
        rec = hessian_box(xc, yc, h, K, prm)
    elif method == "A":
        rec = (
            taylor_box(xc, yc, h, K, prm)
            if xc - h >= Fraction(1, 4)
            else mono_box("A", xc, yc, h, K, prm)
        )
        if not rec["proves"] and rec["kind"] == "taylor":
            alt = mono_box("A", xc, yc, h, K, prm)
            rec = alt if alt["proves"] else rec
    else:
        rec = mono_box("B", xc, yc, h, K, prm)
    rec.update(
        method=method,
        seconds=round(time.time() - t0, 3),
        params=asdict(prm),
        source_sha256=source_sha256(),
        schema=SCHEMA,
        flint=flint.__version__,
        python=platform.python_version(),
    )
    return rec


def run(args: argparse.Namespace) -> int:
    prm = DEFAULT
    tasks = [(m, str(x), str(y), str(h)) for m, x, y, h in base_tasks(prm)]
    if args.limit:
        tasks = tasks[: args.limit]
    out = Path(args.out)
    print(f"near-pair basin: {len(tasks)} base boxes, jobs={args.jobs}, out={out}", flush=True)
    if args.dry_run:
        counts: dict[str, int] = {}
        for t in tasks:
            counts[t[0]] = counts.get(t[0], 0) + 1
        print(counts)
        return 0
    ok = fail = 0
    t0 = time.time()
    with out.open("a", encoding="utf-8") as fh, ProcessPoolExecutor(max_workers=args.jobs) as pool:
        pending: dict[Future, tuple[tuple[str, str, str, str], int]] = {
            pool.submit(certify, t): (t, 0) for t in tasks
        }
        while pending:
            done, _ = wait(pending, return_when=FIRST_COMPLETED)
            for fut in done:
                task, depth = pending.pop(fut)
                rec = fut.result()
                rec["depth"] = depth
                fh.write(json.dumps(rec, sort_keys=True) + "\n")
                fh.flush()
                if rec["proves"]:
                    ok += 1
                    continue
                if depth < prm.max_depth:
                    m, xs, ys, hs = task
                    hh = Fraction(hs) / 2
                    for sx in (-1, 1):
                        for sy in (-1, 1):
                            cx, cy = Fraction(xs) + sx * hh, Fraction(ys) + sy * hh
                            if not relevant(m, cx, cy, hh, prm):
                                continue
                            child = (m, str(cx), str(cy), str(hh))
                            pending[pool.submit(certify, child)] = (child, depth + 1)
                    print(f"  SPLIT {task} -> depth {depth + 1}", flush=True)
                else:
                    fail += 1
                    print(f"  FAIL {task} at max depth", flush=True)
            if (ok + fail) % 50 == 0:
                print(
                    f"  [{ok} ok, {fail} fail, {len(pending)} pending] {time.time() - t0:.0f}s",
                    flush=True,
                )
    print(f"DONE: {ok} proved, {fail} failed [{time.time() - t0:.0f}s] -> {out}")
    return 0 if fail == 0 else 1


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0] if __doc__ else None)
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--out", default=str(OUT))
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument(
        "--box", nargs=4, metavar=("METHOD", "XC", "YC", "H"), help="certify one box and print it"
    )
    args = ap.parse_args()
    if args.box:
        print(json.dumps(certify(tuple(args.box)), indent=1, sort_keys=True))
        return 0
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())
