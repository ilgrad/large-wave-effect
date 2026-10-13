# /// script
# requires-python = ">=3.11"
# dependencies = ["python-flint>=0.8"]
# ///
r"""RIGOROUS per-box certificate for the singleton-channel residual of the balanced Gaussian-excess
inequality (dagger): every bound is an Arb enclosure, there is no finite-difference derivative estimate,
no floating-point grid maximum and no float step on the certified path.  It replaces, for the per-box
claim, the VALIDATED pipeline of verify_dagger_continuum.py (head majorant + theta-ball middle +
theta-Taylor far, lifted to the box by a 25-node finite-difference stencil), which stays available as the
historical reference.

==========================================================================================================
THE CLAIM.  Amplitudes (a1, a3) with b = sqrt(2 - a1^2 - a3^2) (Sigma_2 = 1); the box is
B = [a1c - h, a1c + h] x [a3c - h, a3c + h] with a1c, a3c, h EXACT decimals.  With
    excess(a) = int_0^inf (Mg(t) - M(t)) / t^2 dt                                          (Kluyver (K))
    M(t)  = (2/pi) int_0^{pi/2} J0(a1 t cos th) J0(a3 t cos th) J0(b t sin th) dth,
    Mg(t) = (2/pi) int_0^{pi/2} exp(-(sx2 cos^2 th + sy2 sin^2 th) t^2/2) dth,
    D(a)  = a1^4 + a3^4 + b^4,
certify  excess(a) < K_star D(a)  for EVERY a in B, i.e.  sup_B R < 1  with  R = excess / (K_star D).

COORDINATES.  Work in u = a1^2, v = a3^2, w = b^2 = 2 - u - v.  Then
    J0(a1 y) = 0F1(;1; -u y^2/4)         (entire in u; no square root of b anywhere),
    Mg(t)    = exp(-t^2/4) I0(kappa t^2/4),  kappa = u + v - 1          (closed form, L1),
    D        = u^2 + v^2 + w^2                                          (EXACTLY quadratic in (u, v)).
B maps bijectively onto a rectangle in (u, v); we certify the slightly larger dyadic rectangle
[u_m +- H_u] x [v_m +- H_v] that contains it (box_geometry), so all inputs are exact binary numbers.

==========================================================================================================
THE CERTIFICATE.  excess = G_T + tail_T with G_T = int_0^T F dt, F = (Mg - M)/t^2 (signed, EXACT -- no
head majorant, no region split below T), T = Params.t_cut.  The computable object is the quadrature sum
    Q(u, v) = sum_i om_i [Mg(t_i) - M_N(t_i)] / t_i^2,     M_N(t) = (1/N) sum_j Pi(t, th_j),
Pi = J0 J0 J0 the theta-integrand, th_j = (2j+1) pi/(4N) (midpoint rule), (t_i, om_i) composite
Gauss-Legendre nodes/weights (Arb balls, arb.legendre_p_root).  For every a in the rectangle:
    excess(a) <= P_K(du, dv) + eps_gl + eps_theta + eps_amp + tail_env + tail_mg,
P_K the degree-K Taylor polynomial of Q at (u_m, v_m), computed EXACTLY (Arb, point centre, every
coefficient a ball of radius ~1e-40) by carrying (du, dv) as truncated power series through every J0.
Then  sup_rect [P_K - K_lo D]  is bounded by the exact maximum of a quadratic over a rectangle plus
sum_{i+j>=3} |c_ij| H_u^i H_v^j (with rigorous recentring on sub-rectangles), and
    F_up = that + eps_gl + eps_theta + eps_amp + tail_env + tail_mg.
The box is PROVED iff F_up < 0, and then  sup_B R <= 1 + F_up / (K_hi D_max)  (else 1 + F_up/(K_lo D_min)).

==========================================================================================================
LEMMAS (each used with its hypotheses checked in Arb by the code).
 L1 (Gaussian closed form).  sx2 + sy2 = 1 and sx2 cos^2 + sy2 sin^2 = 1/2 + (kappa/2) cos 2th give
    Mg(t) = exp(-t^2/4) I0(kappa t^2/4) = exp(-t^2/4) 0F1(;1; kappa^2 t^4/64).  |I0(z)| <= I0(|Re z|) <=
    exp(|Re z|) for complex z (I0(z) = (1/pi) int_0^pi e^{z cos phi} dphi).
 L2 (Bessel growth).  |J0(z)| <= exp(|Im z|) for complex z (J0(z) = (1/pi) int_0^pi cos(z sin phi)
    dphi).  For any branch, J0(sqrt(u) y) = 0F1(;1; -u y^2/4), so the bound holds with |Im (sqrt(u) y)|.
 L3 (square-root imaginary part).  For U > 0, |delta| <= R < U:  |Im sqrt(U + delta)| <=
    sqrt(U) rho / (1 + sqrt(1 - rho)),  rho = R/U  (|sqrt(1+z) - 1| = |z| / |sqrt(1+z) + 1| and
    Re sqrt(1+z) >= sqrt(1 - rho)), and |sqrt(U + delta)| <= sqrt(U + R).
 L4 (periodic midpoint rule).  th -> Pi(t, th) is entire, even and pi-periodic.  If |Pi| <= M_eta on the
    strip |Im th| <= eta, then |M(t) - M_N(t)| <= 2 M_eta e^{-4 N eta} / (1 - e^{-4 N eta})  (aliasing:
    the N-point midpoint rule on [0, pi/2] is the 2N-point rule on the period, which integrates
    e^{2ik th} exactly unless 2N | k, and |hat Pi_k| <= M_eta e^{-2|k| eta}).  By L2, for real
    amplitudes and t, M_eta <= exp(t sinh(eta) sqrt((a1 + a3)^2 + b^2)) = exp(t sinh(eta) sqrt(2 + 2 a1 a3))
    (|Im cos(x+iy)| = |sin x| sinh|y|, |Im sin(x+iy)| = |cos x| sinh|y|, Cauchy-Schwarz).
 L5 (Gauss-Legendre, Trefethen ATAP Thm 19.3 for the m-point rule, m >= 2).  If g is analytic in the
    open Bernstein ellipse E_rho (rho > 1) with |g| <= M there, then |int_{-1}^1 g - sum w_k g(x_k)| <=
    (64/15) M rho^{2-2m} / (rho^2 - 1).  (Chebyshev coefficients |a_k| <= 2 M rho^{-k}; the rule is exact
    to degree 2m-1, symmetric, with positive weights summing to 2, and |int T_k| <= 2/15 for even k >= 4.)
    On a panel of half-width w the bound is multiplied by w and M is taken over the affine image.
 L6 (the t-integrand on complex t).  F = (Mg - M)/t^2 is EVEN and ENTIRE in t (Mg, M = 1 - t^2/4 +
    O(t^4)); by L1, L2 for t = x + iy, real amplitudes, |M| <= exp(|y| sqrt(2 + 2 a1 a3)) and
    |Mg| <= exp(max(0, y^2 - x^2)/2), hence |F| <= (exp(max(0, y^2 - x^2)/2) + exp(|y| A)) / x^2 for x > 0,
    and |F(t)| <= (exp((1+|kappa|) R^2/4) + exp(R A)) / R^2 for |t| <= R (maximum modulus principle).
    The same removable singularity holds for the quadrature version (Mg - M_N)/t^2: the midpoint rule
    integrates cos 2th exactly for every N >= 1, so M_N = 1 - t^2/4 + O(t^4) too.
 L7 (bivariate Cauchy estimate).  If Q is analytic on the closed polydisk |du| <= R_u, |dv| <= R_v with
    |Q| <= M_Q there, then for |du| <= H_u, |dv| <= H_v and q = max(H_u/R_u, H_v/R_v) < 1:
    |Q - P_K| <= M_Q sum_{i+j>K} q^{i+j} = M_Q q^{K+1} ((K+2) - (K+1) q) / (1 - q)^2.  M_Q comes from L1-L3,
    L6 applied node by node (positive weights), with Pi bounded by exp(t sqrt((I_u + I_v)^2 + I_w^2)).
 L8 (exact quadratic maximum).  The maximum of a + b.x + x^T C x over a rectangle is attained at a corner,
    on an edge (a 1-D quadratic: endpoints, or its vertex value a - b^2/(4c) when c < 0) or, only if C is
    negative definite, at the critical point.  With Arb coefficients every undecidable branch is resolved
    toward a larger value (convex majorant with c_hi, the global critical value, or the triangle bound).
 L9 (capped envelope tail).  |J0(x)| <= min(1, sqrt(2/(pi x))) for x > 0 (Watson 13.74: x (J0^2 + Y0^2)
    increases to 2/pi), so |M(t)| <= E(t) = (2/pi) int min(1, K2/(t cos th)) min(1, sqrt(2/(pi b t sin th)))
    dth, K2 = 2/(pi sqrt(a1 a3)); E decreases in t, a1, a3, b (so the box minima of the three amplitudes
    bound it), the first factor increases and the second decreases in th, so cell-wise right/left end
    values give an upper Riemann sum, and any cell-wise choice between 1 and the Bessel bound stays an upper
    bound (prefix sums); int_{t_l}^{t_r} E/t^2 <= E(t_l) (1/t_l - 1/t_r); for t >= T_far = Params.tail_far,
    E(t) <= (2/pi) sqrt(2/(pi b)) B0 t^{-1/2}, B0 = int_0^{pi/2} sin^{-1/2} = sqrt(pi) Gamma(1/4)/(2 Gamma(3/4)).
    Mg(t) <= exp(-qmin t^2/2), qmin = min(u+v, w)/2 (L1), and int_T^inf e^{-q t^2/2} t^{-2} dt <=
    e^{-q T^2/2} / (q T^3).
 L10 (the constant).  h(1) = (4/pi^2) int_0^{pi/2} g, g(p) = sqrt(1 + cos^2 p) E(1/(1 + cos^2 p)), g
    decreasing with g(pi/2) = 1; acb.integral encloses int_0^{pi/2 - delta} g (E analytic off [1, inf),
    checked on every ball) and int_{pi/2-delta}^{pi/2} g in [delta, delta g(pi/2 - delta)].  K_star =
    (h(1) - sqrt(pi)/2)/2 in [K_lo, K_hi].

==========================================================================================================
SCOPE.  Assumed, not re-proved here: the identity (K) (verify_dagger_kluyver.py), the two classical facts
quoted in L5 and L9, and the correctness of Arb / python-flint.  That the certified boxes cover the band
lambda in [0.65, 0.82) is checked by verify_dagger_rigorous.py (exact rational interval test); the
near-pair basin lambda >= 0.82 is covered elsewhere (verify_dagger_extremal2.py).  Cost: ~2300 (t, theta)
points per box, ~1.5 s on one core (plus ~6 s once per process for L10), against ~2500 CPU-s per box for
the validated stencil pipeline; on the hardest boxes of the band (a3 ~ 0.61, lambda ~ 0.82) the rigorous
bound is ~0.92 at h = 0.015, where the validated pipeline reported 0.9999 or failed.
"""

from __future__ import annotations

import argparse
import bisect
import functools
import hashlib
import json
import math
import platform
import time
from dataclasses import asdict, dataclass, replace
from fractions import Fraction
from math import comb
from pathlib import Path

import flint
from flint import acb, acb_series, arb, ctx

SOURCE = Path(__file__).resolve()


@dataclass(frozen=True)
class Params:
    """Certificate parameters; every row of the rigorous campaign records them and the checker pins them.
    Exact decimals are strings so that the JSON record reproduces the arithmetic bit for bit."""

    prec: int = 200
    t_cut: str = "24"  # signed integral on [0, t_cut]; capped envelope beyond
    panel: str = "1.5"  # Gauss-Legendre panel width in t
    gl_tol: str = "1e-18"  # per-panel L5 target (chooses the number of nodes)
    gl_max: int = 48
    theta_tol: str = "1e-20"  # per-node L4 target (chooses N)
    theta_min: int = 2
    degree: int = 10  # amplitude Taylor degree K
    q_menu: tuple[str, ...] = ("1/3", "1/4", "1/6", "1/8", "1/10", "1/12", "1/16")
    t_star: str = "2"  # maximum-modulus radius for the small-t nodes in L7
    split_depth: int = 2  # recentring depth of the polynomial maximum (4^depth sub-rectangles)
    tail_segments: tuple[tuple[str, str, int], ...] = (
        ("24", "60", 360),
        ("60", "200", 350),
        ("200", "2000", 450),
        ("2000", "20000", 450),
    )
    tail_cells: int = 16384
    tail_far: str = "20000"


DEFAULT = Params()


def _a(x: str | int | Fraction) -> arb:
    """Exact-as-possible arb from a decimal string, an int or a Fraction (a ball if not dyadic)."""
    if isinstance(x, Fraction):
        return arb(x.numerator) / x.denominator
    return arb(x)


def _pt(x: arb) -> arb:
    """The exact midpoint of a ball, as a radius-0 arb."""
    return arb(x.mid())


def _up(x: arb) -> arb:
    return arb(x.upper())


def _lo(x: arb) -> arb:
    return arb(x.lower())


def _float_up(x: arb) -> float:
    """A float >= every point of x (display only; the certificate uses the Arb ball)."""
    f = float(x.upper())
    return f if arb(f) >= x.upper() else math.nextafter(f, math.inf)


def _max_up(vals: list[arb]) -> arb:
    out = _up(vals[0])
    for v in vals[1:]:
        out = out.max(_up(v))
    return out


def source_sha256() -> str:
    return hashlib.sha256(SOURCE.read_bytes()).hexdigest()


# ==========================================================================================================
# L10: the sharp constant K_star, enclosed.
# ==========================================================================================================
@functools.cache
def kstar_enclosure(prec: int = 128) -> tuple[arb, arb]:
    """(K_lo, K_hi) with K_lo <= K_star <= K_hi, K_star = (h(1) - sqrt(pi)/2)/2 (L10)."""
    old = ctx.prec
    ctx.prec = prec
    try:
        half = arb.pi() / 2
        delta = arb("1e-5")

        def g(p: acb, analytic: bool) -> acb:
            c = p.cos()
            d = 1 + c * c
            m = 1 / d
            if analytic:
                mr, mi = m.real, m.imag
                if mr.upper() >= 1 and mi.lower() <= 0 and mi.upper() >= 0:
                    return acb("nan")
            return d.sqrt(analytic=analytic) * acb.elliptic_e(m)

        body = acb.integral(
            g, 0, half - delta, rel_tol=arb(2) ** -80, eval_limit=10**8, depth_limit=10**5
        ).real
        gend = g(acb(half - delta), False).real
        four_pi2 = 4 / arb.pi() ** 2
        h_lo = (body + delta) * four_pi2
        h_hi = (body + delta * gend) * four_pi2
        root = arb.pi().sqrt() / 2
        k_lo, k_hi = _lo((h_lo - root) / 2), _up((h_hi - root) / 2)
        if not (k_lo.is_finite() and k_hi.is_finite() and k_lo > 0):
            raise ArithmeticError("K_star enclosure failed")
        return k_lo, k_hi
    finally:
        ctx.prec = old


# ==========================================================================================================
# Box geometry: the certified dyadic rectangle in (u, v) and the real ranges the error bounds need.
# ==========================================================================================================
@dataclass(frozen=True)
class Geometry:
    a1c: str
    a3c: str
    h: str
    um: arb  # exact dyadic centre
    vm: arb
    Hu: arb  # exact dyadic half-widths
    Hv: arb

    @property
    def wm(self) -> arb:
        return 2 - self.um - self.vm

    @property
    def u_rng(self) -> tuple[arb, arb]:
        return self.um - self.Hu, self.um + self.Hu

    @property
    def v_rng(self) -> tuple[arb, arb]:
        return self.vm - self.Hv, self.vm + self.Hv

    def amp_max_A(self) -> arb:
        """Upper bound of sqrt(2 + 2 a1 a3) = sqrt((a1 + a3)^2 + b^2) on the rectangle (L4, L6)."""
        return _up((2 + 2 * (self.u_rng[1] * self.v_rng[1]).sqrt()).sqrt())

    def kappa_abs_max(self) -> arb:
        """Upper bound of |kappa| = |u + v - 1| on the rectangle."""
        lo = self.u_rng[0] + self.v_rng[0] - 1
        hi = self.u_rng[1] + self.v_rng[1] - 1
        return _up(lo.abs_upper()).max(_up(hi.abs_upper()))

    def qmin_lower(self) -> arb:
        """Lower bound of qmin = min(u + v, w)/2 on the rectangle (L9)."""
        s_lo = self.u_rng[0] + self.v_rng[0]
        s_hi = self.u_rng[1] + self.v_rng[1]
        return _lo(s_lo / 2).min(_lo((2 - s_hi) / 2))

    def amp_minima(self) -> tuple[arb, arb, arb]:
        """Lower bounds of a1, a3, b on the rectangle (L9 needs the box minimum of each amplitude)."""
        a1 = _lo(self.u_rng[0].sqrt())
        a3 = _lo(self.v_rng[0].sqrt())
        b = _lo((2 - self.u_rng[1] - self.v_rng[1]).sqrt())
        return a1, a3, b

    def lam_range(self) -> tuple[arb, arb]:
        """Enclosure of lambda = 1 - min(u, v, w)/2 over the rectangle."""
        w_lo, w_hi = 2 - self.u_rng[1] - self.v_rng[1], 2 - self.u_rng[0] - self.v_rng[0]
        m_lo = self.u_rng[0].min(self.v_rng[0]).min(w_lo)
        m_hi = self.u_rng[1].min(self.v_rng[1]).min(w_hi)
        return _lo(1 - m_hi / 2), _up(1 - m_lo / 2)


def _dyadic(x: arb, bits: int = 64) -> arb:
    """A dyadic point close to the midpoint of x (exact binary, radius 0)."""
    scale = arb(2) ** bits
    n = (x.mid() * scale).floor().unique_fmpz()
    if n is None:
        raise ArithmeticError(f"the midpoint of {x} times 2^{bits} has no exact floor")
    return arb(int(n)) / scale


def box_geometry(a1c: str, a3c: str, h: str) -> Geometry:
    """Certified rectangle [u_m +- H_u] x [v_m +- H_v] containing {(a1^2, a3^2) : a in B}."""
    A1, A3, Hh = _a(a1c), _a(a3c), _a(h)
    if not (Hh > 0 and A1 - Hh > 0 and A3 - Hh > 0):
        raise ValueError(f"box ({a1c}, {a3c}, {h}) must have positive amplitudes")
    ulo, uhi = _lo((A1 - Hh) ** 2), _up((A1 + Hh) ** 2)
    vlo, vhi = _lo((A3 - Hh) ** 2), _up((A3 + Hh) ** 2)
    um, vm = _dyadic((ulo + uhi) / 2), _dyadic((vlo + vhi) / 2)
    Hu = _up((um - ulo).max(uhi - um))
    Hv = _up((vm - vlo).max(vhi - vm))
    geo = Geometry(a1c, a3c, h, um, vm, Hu, Hv)
    if not (2 - geo.u_rng[1] - geo.v_rng[1] > 0):
        raise ValueError(f"box ({a1c}, {a3c}, {h}) leaves the disc a1^2 + a3^2 < 2")
    return geo


# ==========================================================================================================
# Quadrature plan with rigorous errors (L4, L5, L6).
# ==========================================================================================================
@functools.cache
def _gl_rule(m: int, prec: int) -> tuple[tuple[arb, arb], ...]:
    old = ctx.prec
    ctx.prec = prec
    try:
        return tuple(arb.legendre_p_root(m, k, weight=True) for k in range(m))
    finally:
        ctx.prec = old


_RHO_MENU = ("1.25", "1.5", "2", "2.5", "3", "4", "5", "6", "8", "10", "13", "16", "20", "26", "32")


def _panel_bound(alpha: arb, beta: arb, m: int, A: arb, kap: arb) -> arb:
    """L5 + L6: an upper bound on |int_alpha^beta F - GL_m| valid for every amplitude in the rectangle."""
    w = (beta - alpha) / 2
    c = (alpha + beta) / 2
    bounds: list[arb] = []
    for rs in _RHO_MENU:
        rho = arb(rs)
        ae = w * (rho + 1 / rho) / 2
        be = w * (rho - 1 / rho) / 2
        R = c + ae
        Mdisk = ((1 + kap) * R * R / 4).exp() + (R * A).exp()
        M = Mdisk / (R * R)
        xmin = c - ae
        if xmin > 0:
            gauss = (((be * be - xmin * xmin) / 2).max(arb(0))).exp()
            Mdir = (gauss + (be * A).exp()) / (xmin * xmin)
            M = M.min(Mdir)
        bound = w * arb(64) / 15 * M * rho ** (2 - 2 * m) / (rho * rho - 1)
        bounds.append(_up(bound))
    return functools.reduce(arb.min, bounds)


def _theta_plan(t: arb, A: arb, tol: float, nmin: int) -> tuple[int, arb]:
    """Smallest N >= nmin whose L4 aliasing bound is <= tol, and that bound (valid on the rectangle)."""
    tA = float(_up(t * A))
    n = nmin
    while True:
        r = 4 * n / tA if tA > 0 else math.inf
        if r > 1.0001:
            eta = min(math.acosh(r), 60.0)
            log_bound = math.log(2.0) + tA * math.sinh(eta) - 4 * n * eta
            if log_bound < math.log(tol) - 1.0:
                break
        n += 1
        if n > 4000:
            raise RuntimeError(f"theta plan did not converge at t={float(t.mid())}")
    eta_a = arb(repr(eta))
    e4 = (-4 * n * eta_a).exp()
    bound = 2 * (t * eta_a.sinh() * A).exp() * e4 / (1 - e4)
    return n, _up(bound)


@dataclass(frozen=True)
class Node:
    t: arb
    om: arb
    N: int


def quadrature_plan(geo: Geometry, prm: Params) -> tuple[list[Node], arb, arb]:
    """Composite GL nodes on [0, t_cut], theta counts, and the rigorous quadrature errors
    (eps_gl, eps_theta) valid for every amplitude in the rectangle."""
    A, kap = geo.amp_max_A(), geo.kappa_abs_max()
    T, pw = Fraction(prm.t_cut), Fraction(prm.panel)
    npan = math.ceil(T / pw)
    gl_tol, th_tol = arb(prm.gl_tol), float(prm.theta_tol)
    nodes: list[Node] = []
    eps_gl, eps_th = arb(0), arb(0)
    for k in range(npan):
        alpha, beta = _a(pw * k), _a(min(pw * (k + 1), T))
        m = 2
        while True:
            bnd = _panel_bound(alpha, beta, m, A, kap)
            if bnd < gl_tol:
                break
            m += 1
            if m > prm.gl_max:
                raise RuntimeError(f"GL panel [{alpha}, {beta}] needs more than {prm.gl_max} nodes")
        eps_gl += bnd
        w, c = (beta - alpha) / 2, (alpha + beta) / 2
        for x, wx in _gl_rule(m, prm.prec):
            t, om = c + w * x, w * wx
            N, eb = _theta_plan(t, A, th_tol, prm.theta_min)
            eps_th += om * eb / (t * t)
            nodes.append(Node(t, om, N))
    return nodes, _up(eps_gl), _up(eps_th)


# ==========================================================================================================
# The exact amplitude Taylor polynomial of the quadrature sum Q (point centre, truncated power series).
# ==========================================================================================================
def _j0_lin(z0: arb, z1: arb, n: int) -> list[arb]:
    """Taylor coefficients in x of 0F1(;1; z0 + z1 x) (= J0(sqrt(u) y) with z = -u y^2/4), length n."""
    if ctx.cap < n:  # acb_series silently truncates at ctx.cap; zero-padding would then be wrong
        raise ValueError(f"ctx.cap={ctx.cap} < series length {n}")
    s = acb_series([acb(z0), acb(z1)], prec=n)
    r = acb_series.hypgeom([], [acb(1)], s)  # n = -1: automatic truncation WITH Arb's error bound
    cs = [c.real for c in r.coeffs()]
    return cs + [arb(0)] * (n - len(cs))


def _mg_series(t: arb, kap_m: arb, n: int) -> list[arb]:
    """Taylor coefficients in s = du + dv of Mg = exp(-t^2/4) 0F1(;1; (kap_m + s)^2 t^4/64) (L1)."""
    c = t**4 / 64
    if ctx.cap < n:
        raise ValueError(f"ctx.cap={ctx.cap} < series length {n}")
    s = acb_series([acb(kap_m * kap_m * c), acb(2 * kap_m * c), acb(c)], prec=n)
    r = acb_series.hypgeom([], [acb(1)], s)
    e = (-t * t / 4).exp()
    cs = [e * c.real for c in r.coeffs()]
    return cs + [arb(0)] * (n - len(cs))


def taylor_polynomial(geo: Geometry, nodes: list[Node], K: int) -> dict[tuple[int, int], arb]:
    """c_ij with Q(u_m + du, v_m + dv) = sum_{i+j<=K} c_ij du^i dv^j + O(deg K+1), all exact balls.

    Pi = F1(du) F2(dv) F3(du + dv) with F1 = 0F1(;1; -(u_m + du) y^2/4), F2 likewise in v, F3 =
    0F1(;1; -(w_m - s) y3^2/4), y = t cos th, y3 = t sin th; the trivariate moments
    T[p][q][r] = sum W f1_p f2_q f3_r are accumulated and contracted with binomials at the end."""
    n = K + 1
    if ctx.cap < n:
        ctx.cap = n
    um, vm, wm = geo.um, geo.vm, geo.wm
    kap_m = um + vm - 1
    T = [[[arb(0)] * (n - p - q) for q in range(n - p)] for p in range(n)]
    mg = [arb(0)] * n
    pi = arb.pi()
    trig_cache: dict[int, list[tuple[arb, arb]]] = {}
    for nd in nodes:
        t, inv = nd.t, nd.om / (nd.t * nd.t)
        g = _mg_series(t, kap_m, n)
        for k in range(n):
            mg[k] += inv * g[k]
        W = inv / nd.N
        if nd.N not in trig_cache:
            trig_cache[nd.N] = [
                ((pi * (2 * j + 1) / (4 * nd.N)).cos(), (pi * (2 * j + 1) / (4 * nd.N)).sin())
                for j in range(nd.N)
            ]
        tt4 = t * t / 4
        for cj, sj in trig_cache[nd.N]:
            yc = tt4 * cj * cj
            ys = tt4 * sj * sj
            f1 = _j0_lin(-um * yc, -yc, n)
            f2 = _j0_lin(-vm * yc, -yc, n)
            f3 = _j0_lin(-wm * ys, ys, n)
            f3w = [W * x for x in f3]
            for p in range(n):
                f1p, Tp = f1[p], T[p]
                for q in range(n - p):
                    a = f1p * f2[q]
                    row = Tp[q]
                    for r in range(n - p - q):
                        row[r] += a * f3w[r]
    coef: dict[tuple[int, int], arb] = {}
    for i in range(n):
        for j in range(n - i):
            acc = mg[i + j] * comb(i + j, i)
            for p in range(i + 1):
                for q in range(j + 1):
                    r = (i - p) + (j - q)
                    acc -= T[p][q][r] * comb(r, i - p)
            coef[(i, j)] = acc
    return coef


# ==========================================================================================================
# L7: amplitude Cauchy remainder of Q.
# ==========================================================================================================
def amplitude_remainder(geo: Geometry, nodes: list[Node], prm: Params) -> tuple[arb, str]:
    """min over the q-menu of M_Q q^{K+1}((K+2)-(K+1)q)/(1-q)^2 (L7); returns (bound, chosen q)."""
    K = prm.degree
    um, vm, wm = geo.um, geo.vm, geo.wm
    kap_m = (um + vm - 1).abs_upper()
    ts = _a(prm.t_star)
    best, best_q = None, ""
    for qs in prm.q_menu:
        q = _a(Fraction(qs))
        Ru, Rv = geo.Hu / q, geo.Hv / q
        if not (Ru < um and Rv < vm and Ru + Rv < wm):
            continue

        def im_sqrt(U: arb, R: arb) -> arb:
            rho = R / U
            return U.sqrt() * rho / (1 + (1 - rho).sqrt())

        Iu, Iv, Iw = im_sqrt(um, Ru), im_sqrt(vm, Rv), im_sqrt(wm, Ru + Rv)
        B = ((Iu + Iv) ** 2 + Iw * Iw).sqrt()
        kmax = kap_m + Ru + Rv
        Ac = (((um + Ru).sqrt() + (vm + Rv).sqrt()) ** 2 + (wm + Ru + Rv)).sqrt()
        small = (((1 + kmax) * ts * ts / 4).exp() + (ts * Ac).exp()) / (ts * ts)
        MQ = arb(0)
        for nd in nodes:
            if nd.t < ts:
                MQ += nd.om * small
            else:
                gau = (-(1 - kmax) * nd.t * nd.t / 4).exp()
                MQ += nd.om * (gau + (nd.t * B).exp()) / (nd.t * nd.t)
        tail = q ** (K + 1) * ((K + 2) - (K + 1) * q) / (1 - q) ** 2
        bound = _up(MQ * tail)
        if best is None or bound < best:
            best, best_q = bound, qs
    if best is None:
        raise ValueError("no admissible polydisk radius in q_menu")
    return best, best_q


# ==========================================================================================================
# L9: the tail [t_cut, inf).
# ==========================================================================================================
def tail_bound(geo: Geometry, prm: Params) -> tuple[arb, arb]:
    """(envelope bound on sup_rect int_T^inf |M|/t^2, bound on sup_rect int_T^inf Mg/t^2)  (L9)."""
    a1, a3, b = geo.amp_minima()
    pi = arb.pi()
    K2 = 2 / (pi * (a1 * a3).sqrt())
    c2b = 2 / (pi * b)
    n = prm.tail_cells
    dth = (pi / 2) / n
    # Cell k = [k dth, (k+1) dth]: the first factor is <= min(1, al_k/t), al_k = K2/cos((k+1) dth), and the
    # second <= min(1, sqrt(ga_k/t)), ga_k = c2b/sin(k dth).  The last cell (cos = 0) and the first (sin = 0)
    # only admit the cap 1.  Any choice between 1 and the Bessel bound, cell by cell, is an upper bound, so
    # the float breakpoints below only affect tightness; the sums themselves are Arb prefix sums.
    al = [_up(K2 / (dth * (k + 1)).cos()) for k in range(n - 1)]
    sg = [_up((c2b / (dth * k).sin()).sqrt()) for k in range(1, n)]  # sqrt(ga_k), k = 1..n-1
    al_f = [float(x) for x in al] + [math.inf]  # increasing in k
    neg_ga_f = [-math.inf] + [-(float(x) ** 2) for x in sg]  # -ga_k, increasing in k
    pa, ps, pas = [arb(0)], [arb(0), arb(0)], [arb(0), arb(0)]
    for k in range(n - 1):
        pa.append(pa[-1] + al[k])  # pa[k] = sum_{j<k} al_j
    for k in range(1, n):
        ps.append(ps[-1] + sg[k - 1])  # ps[k] = sum_{1<=j<k} sqrt(ga_j)
        if k < n - 1:
            pas.append(pas[-1] + al[k] * sg[k - 1])  # pas[k] = sum_{1<=j<k} al_j sqrt(ga_j)

    def env(t: arb) -> arb:
        tf = float(t)
        k1 = min(bisect.bisect_left(al_f, tf), n - 1)  # cells >= k1: first factor capped at 1
        k2 = max(bisect.bisect_left(neg_ga_f, -tf), 1)  # cells < k2: second factor capped at 1
        rt = t.sqrt()
        if k1 >= k2:
            acc = (pa[k2] / t) + (pas[k1] - pas[k2]) / (t * rt) + (ps[n] - ps[k1]) / rt
        else:
            acc = (pa[k1] / t) + (k2 - k1) + (ps[n] - ps[k2]) / rt
        return _up(2 / pi * acc * dth)

    total = arb(0)
    t_start = None
    t_prev_end = None
    for lo_s, hi_s, nseg in prm.tail_segments:
        lo, hi = _a(lo_s), _a(hi_s)
        if t_start is None:
            t_start = lo
        elif t_prev_end != lo:
            raise ValueError("tail segments must be contiguous")
        t_prev_end = hi
        dt = (hi - lo) / nseg
        for i in range(nseg):
            tl, tr = lo + dt * i, lo + dt * (i + 1)
            total += env(tl) * (1 / tl - 1 / tr)
    if not (t_start == _a(prm.t_cut) and t_prev_end == _a(prm.tail_far)):
        raise ValueError("tail segments must run from t_cut to tail_far")
    B0 = pi.sqrt() * arb("0.25").gamma() / (2 * arb("0.75").gamma())
    far = (2 / pi) * c2b.sqrt() * B0 * (arb(2) / 3) * _a(prm.tail_far) ** arb("-1.5")
    q = geo.qmin_lower()
    T = _a(prm.t_cut)
    mg = (-q * T * T / 2).exp() / (q * T**3)
    return _up(total + far), _up(mg)


# ==========================================================================================================
# L8: rigorous maximum of the polynomial (minus K D) over the rectangle.
# ==========================================================================================================
def _sup_quad_1d(a: arb, b: arb, c: arb, H: arb) -> arb:
    """Upper bound of sup_{|s|<=H} a + b s + c s^2."""
    chi = _up(c)
    if chi >= 0:
        return _max_up([a + b * H + chi * H * H, a - b * H + chi * H * H])
    vals = [a + b * H + c * H * H, a - b * H + c * H * H]
    s_star = -b / (2 * c)
    if not (s_star > H or s_star < -H):
        vals.append(a - b * b / (4 * c))
    return _max_up(vals)


def _sup_quad_rect(q: dict[tuple[int, int], arb], Hx: arb, Hy: arb) -> arb:
    """Upper bound of sup over [-Hx,Hx]x[-Hy,Hy] of the quadratic part of q (L8)."""
    al, bx, by = q[(0, 0)], q[(1, 0)], q[(0, 1)]
    cxx, cxy, cyy = q[(2, 0)], q[(1, 1)], q[(0, 2)]
    cands = []
    for sx in (Hx, -Hx):
        cands.append(_sup_quad_1d(al + bx * sx + cxx * sx * sx, by + cxy * sx, cyy, Hy))
    for sy in (Hy, -Hy):
        cands.append(_sup_quad_1d(al + by * sy + cyy * sy * sy, bx + cxy * sy, cxx, Hx))
    det = 4 * cxx * cyy - cxy * cxy
    surely_nd = cxx < 0 and det > 0
    surely_not_nd = cxx.lower() >= 0 or cyy.lower() >= 0 or det.upper() <= 0
    if surely_nd:
        xs = (cxy * by - 2 * cyy * bx) / det
        ys = (cxy * bx - 2 * cxx * by) / det
        if not (xs > Hx or xs < -Hx or ys > Hy or ys < -Hy):
            cands.append(al + (bx * xs + by * ys) / 2)
    elif not surely_not_nd:
        cands.append(
            al
            + bx.abs_upper() * Hx
            + by.abs_upper() * Hy
            + cxx.abs_upper() * Hx * Hx
            + cxy.abs_upper() * Hx * Hy
            + cyy.abs_upper() * Hy * Hy
        )
    return _max_up(cands)


def _shift(c: dict[tuple[int, int], arb], x0: arb, y0: arb) -> dict[tuple[int, int], arb]:
    """Exact recentring: coefficients of c(x0 + x, y0 + y)."""
    K = max(i + j for i, j in c)
    out = {}
    for k in range(K + 1):
        for m in range(K + 1 - k):
            acc = arb(0)
            for i in range(k, K + 1):
                for j in range(m, K + 1 - i):
                    if (i, j) in c:
                        acc += c[(i, j)] * comb(i, k) * comb(j, m) * x0 ** (i - k) * y0 ** (j - m)
            out[(k, m)] = acc
    return out


def _higher(c: dict[tuple[int, int], arb], Hx: arb, Hy: arb) -> arb:
    acc = arb(0)
    for (i, j), v in c.items():
        if i + j >= 3:
            acc += v.abs_upper() * Hx**i * Hy**j
    return _up(acc)


def poly_sup(c: dict[tuple[int, int], arb], Hx: arb, Hy: arb, depth: int) -> tuple[arb, arb, arb]:
    """Upper bound of sup over the rectangle of the polynomial c; returns (bound, quadratic-part bound,
    higher-order bound) of the best level (the whole rectangle, or 4^d recentred sub-rectangles)."""
    levels = []
    for d in range(depth + 1):
        k = 2**d
        hx, hy = Hx / k, Hy / k
        cells = []
        for ix in range(k):
            for iy in range(k):
                x0 = -Hx + hx * (2 * ix + 1)
                y0 = -Hy + hy * (2 * iy + 1)
                cc = c if d == 0 else _shift(c, x0, y0)
                qv, hv = _sup_quad_rect(cc, hx, hy), _higher(cc, hx, hy)
                cells.append((_up(qv + hv), qv, hv))
        levels.append(max(cells, key=lambda cell: cell[0]))
    return min(levels, key=lambda level: level[0])


# ==========================================================================================================
# The per-box certificate.
# ==========================================================================================================
def ball(x: arb, digits: int = 25) -> dict[str, str]:
    """Rigorous decimal record of a ball: arb(f'[{mid} +/- {rad}]') contains x."""
    s = x.str(digits, radius=True)
    if s.startswith("["):
        mid, rad = s[1:-1].split(" +/- ")
    else:
        mid, rad = s, "0"
    return {"mid": mid.strip(), "rad": rad.strip()}


def ball_arb(rec: dict[str, str]) -> arb:
    return arb(f"[{rec['mid']} +/- {rec['rad']}]")


ERROR_TERMS = ("eps_gl", "eps_theta", "eps_amp", "tail_envelope", "tail_mg")


def d_quadratic(geo: Geometry) -> dict[tuple[int, int], arb]:
    """D = u^2 + v^2 + w^2 in (du, dv) about the centre: exact, quadratic."""
    um, vm, wm = geo.um, geo.vm, geo.wm
    return {
        (0, 0): um * um + vm * vm + wm * wm,
        (1, 0): 2 * (um - wm),
        (0, 1): 2 * (vm - wm),
        (2, 0): arb(2),
        (1, 1): arb(2),
        (0, 2): arb(2),
    }


def decide(
    coef: dict[tuple[int, int], arb],
    geo: Geometry,
    errs: dict[str, arb],
    kstar: tuple[arb, arb],
    depth: int,
) -> dict[str, arb]:
    """The final step, from the Taylor coefficients and the error terms (shared with the checker):
    F_up >= sup_rect [P_K - K_lo D] + sum(errs) and the ratio bound 1 + F_up / (K D_ref)."""
    k_lo, k_hi = kstar
    Dq = d_quadratic(geo)
    F = dict(coef)
    for key, dv in Dq.items():
        F[key] = F[key] - k_lo * dv
    sup_poly, sup_quad, sup_high = poly_sup(F, geo.Hu, geo.Hv, depth)
    F_up = sup_poly
    for name in ERROR_TERMS:
        F_up = F_up + errs[name]
    F_up = _up(F_up)
    Dmax = _sup_quad_rect(Dq, geo.Hu, geo.Hv)
    Dmin = _lo(-_sup_quad_rect({k: -v for k, v in Dq.items()}, geo.Hu, geo.Hv))
    ratio = 1 + F_up / (k_hi * Dmax) if F_up < 0 else 1 + F_up / (k_lo * Dmin)
    return {
        "F_up": F_up,
        "ratio": _up(ratio),
        "sup_poly_minus_kD": sup_poly,
        "quadratic_part": sup_quad,
        "higher_order": sup_high,
        "Dmin": Dmin,
        "Dmax": Dmax,
    }


def certify_box(a1c: str, a3c: str, h: str, prm: Params = DEFAULT) -> dict:
    """The rigorous certificate for B = [a1c +- h] x [a3c +- h] (exact decimals).  Returns a JSON-ready
    record; 'proves' is True iff the rigorous F_up < 0, i.e. sup_B excess/(K_star D) < 1."""
    t0 = time.time()
    ctx.prec = prm.prec
    geo = box_geometry(a1c, a3c, h)
    nodes, eps_gl, eps_th = quadrature_plan(geo, prm)
    t1 = time.time()
    coef = taylor_polynomial(geo, nodes, prm.degree)
    t2 = time.time()
    eps_amp, q_used = amplitude_remainder(geo, nodes, prm)
    tail_env, tail_mg = tail_bound(geo, prm)
    t3 = time.time()
    kstar = kstar_enclosure()
    ctx.prec = prm.prec
    errs = {
        "eps_gl": eps_gl,
        "eps_theta": eps_th,
        "eps_amp": eps_amp,
        "tail_envelope": tail_env,
        "tail_mg": tail_mg,
    }
    dec = decide(coef, geo, errs, kstar, prm.split_depth)
    # Diagnostic only: the sup of the polynomial model of G_T alone.
    P_sup = poly_sup(coef, geo.Hu, geo.Hv, prm.split_depth)[0]
    lam_lo, lam_hi = geo.lam_range()
    terms = {name: ball(v) for name, v in errs.items()}
    terms.update({k: ball(dec[k]) for k in ("sup_poly_minus_kD", "quadratic_part", "higher_order")})
    terms.update({"Dmin": ball(dec["Dmin"]), "Dmax": ball(dec["Dmax"])})
    terms.update({"K_lo": ball(kstar[0]), "K_hi": ball(kstar[1])})
    terms.update({"G_centre": ball(coef[(0, 0)]), "sup_P": ball(P_sup)})
    return {
        "a1": a1c,
        "a3": a3c,
        "h": h,
        "proves": bool(dec["F_up"] < 0),
        "upper": _float_up(dec["ratio"]),
        "upper_ball": ball(dec["ratio"]),
        "F_up": ball(dec["F_up"]),
        "lam_lo": float(lam_lo),
        "lam_hi": float(lam_hi),
        "terms": terms,
        "rect": {"um": ball(geo.um), "vm": ball(geo.vm), "Hu": ball(geo.Hu), "Hv": ball(geo.Hv)},
        "coef": {f"{i},{j}": ball(v, 30) for (i, j), v in sorted(coef.items())},
        "q": q_used,
        "n_t_nodes": len(nodes),
        "n_points": sum(nd.N for nd in nodes),
        "seconds": {
            "plan": round(t1 - t0, 3),
            "taylor": round(t2 - t1, 3),
            "bounds": round(t3 - t2, 3),
            "total": round(time.time() - t0, 3),
        },
        "params": asdict(prm),
        "code_sha256": source_sha256(),
        "flint": flint.__version__,
        "python": platform.python_version(),
    }


def main() -> int:
    ap = argparse.ArgumentParser(
        description="rigorous per-box certificate (see the module docstring)"
    )
    ap.add_argument("--box", action="append", required=True, help="a1c,a3c,h as exact decimals")
    ap.add_argument("--degree", type=int, default=None, help="override Params.degree (experiments)")
    ap.add_argument("--t-cut", default=None, help="override Params.t_cut (experiments)")
    args = ap.parse_args()
    prm = DEFAULT
    if args.degree is not None:
        prm = replace(prm, degree=args.degree)
    if args.t_cut is not None:
        segs = ((args.t_cut, *prm.tail_segments[0][1:]), *prm.tail_segments[1:])
        prm = replace(prm, t_cut=args.t_cut, tail_segments=segs)
    ok = True
    for spec in args.box:
        a1c, a3c, h = (x.strip() for x in spec.split(","))
        rec = certify_box(a1c, a3c, h, prm)
        ok &= rec["proves"]
        print(json.dumps(rec, indent=1, sort_keys=True), flush=True)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
