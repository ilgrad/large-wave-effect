# /// script
# requires-python = ">=3.11"
# dependencies = ["mpmath>=1.3"]
# ///
"""Validated-numerics proof that max_{[0,1]} F = F(1) = beta_odd for the two-copy profile (Proposition bn-profile).

F(s) = sqrt(2/pi) int_0^1 g(a, b) du,  a = (s^2-u^2)_+^{-1/2},  b = (s^2-(1-u)^2)_+^{-1/2},
g(a, b) = E_{Phi,Psi} sqrt(a cos^2 Phi + b cos^2 Psi),  g(a, 0) = (2/pi) sqrt(a),  beta_odd = F(1).

Writing g(a, b) = sqrt(M) h(m/M), M = max, m = min, the function h(rho) = E sqrt(cos^2 Phi + rho cos^2 Psi)
is increasing and concave on [0, 1] (h' = E[cos^2 Psi / (2 sqrt(.))] > 0, h'' = -E[cos^4 Psi / (4 (.)^{3/2})]
< 0), with h(0) = 2/pi and h'(1) = h(1)/4 (swap Phi and Psi).  Five validated steps cover [0, 1]; every
number below is an interval enclosure (mpmath.iv), never a floating-point estimate.

  (A) h(k/8), k = 1..8, enclosed: h(rho) = (2/pi) int_0^{pi/2} f, f(phi) = (2/pi) sqrt(x + rho) E(rho/(x + rho)),
      x = cos^2 phi, E the complete elliptic integral (computed by the AGM with its error bracket); f decreases
      in phi, so the left- and right-endpoint sums bracket h.
  (B) beta_odd >= beta_lo: the chord interpolant of the concave h through the nine lower bounds lies below h,
      so F(1) >= 2 sqrt(2/pi) int_0^{1/2} sqrt(b) chord(a/b) du (on [0, 1/2], b >= a), enclosed after u = t^4.
  (C) [0, 1/2]: F = c_0 sqrt(s) <= c_0/sqrt(2).  [1/2, 0.575]: sqrt(x^2 + y^2) <= x + y gives
      g(a, b) <= (2/pi)(sqrt(a) + sqrt(b)), hence F(s) <= c_0 sqrt(s) <= c_0 sqrt(0.575); c_0 = 16/Gamma(1/4)^2.
  (D) [0.575, 0.90]: the tangent of h at rho = 1, h(rho) <= h(1)(3 + rho)/4, gives
          F(s) <= 2 sqrt(2/pi) [ (2/pi) int_0^{1-s} (s^2-u^2)^{-1/4} du + (h(1)/4) int_{1-s}^{1/2} (3 sqrt(b) + a/sqrt(b)) du ],
      enclosed for s running over subintervals (adaptive bisection), with v = u-(1-s) = t^4 on the second integral.
  (E) [0.90, 1]: F'(s) > 0.  In the w-form
          F(s) = sqrt(2/pi)[ (4/pi)(1-s) A(s) + (2s-1)^{3/4} B(s) ],   q = 2s-1,
            A(s) = int_0^1 (s^2-(1-s)^2 w^2)^{-1/4} dw,
            B(s) = int_0^1 g(ahat, bhat) dw,   ahat=[(1-w)(1+qw)]^{-1/2}, bhat=[w(2s-qw)]^{-1/2},
          F'(s) = sqrt(2/pi)[ (4/pi)(-A + (1-s)A') + (3/2) q^{-1/4} B + q^{3/4} B' ],
      two elementary inequalities remove the elliptic g:
        * QM-AM:  sqrt(a cos^2 Phi + b cos^2 Psi) >= (a^{1/2}|cos Phi| + b^{1/2}|cos Psi|)/sqrt2, so
          g(a,b) >= (sqrt2/pi)(sqrt(a) + sqrt(b)) and B(s) >= B_lo(s) = int (sqrt2/pi)(sqrt(ahat)+sqrt(bhat)) dw;
        * d_a g = E[cos^2 Phi / (2 sqrt(...))] <= E[|cos Phi|/(2 sqrt a)] = 1/(pi sqrt a), and with
          d_s ahat = -(1-w)w[(1-w)(1+qw)]^{-3/2} (likewise bhat),
              |B'(s)| <= Bp_up(s) = (1/pi) int w(1-w)[ ((1-w)(1+qw))^{-5/4} + (w(2s-qw))^{-5/4} ] dw.
      Since B' < 0, q^{3/4} B' >= -q^{3/4} Bp_up, so
          F'(s) >= F_lo(s) := sqrt(2/pi)[ (4/pi)(-A + (1-s)A') + (3/2) q^{-1/4} B_lo - q^{3/4} Bp_up ],
      all four integrals elementary; the w^{-1/4} endpoint singularities are removed by w = t^4 on [0,1/2] and
      w = 1 - r^4 on [1/2,1], with the t/r powers cancelled by hand so no 0*inf appears.

Every integral is a rigorous interval rectangle sum: the interval evaluation over each box contains the
integrand's range, so box x width contains the contribution, and s itself is an interval in (D) and (E).

RESULT: (C)-(D) give F(s) < beta_lo <= beta_odd on [0, 0.90] and (E) gives F strictly increasing on [0.90, 1],
so max_{[0,1]} F = F(1) = beta_odd, attained only at s = 1.  Companion to verify_bn_profile_max.py (the same
majorants by ordinary quadrature, and the endpoints).
"""

from __future__ import annotations

from itertools import pairwise

import mpmath as mp

iv = mp.iv
iv.prec = 90
mp.mp.prec = 90  # = iv.prec, so converting an interval endpoint to mp.mpf is exact


def _qrt(x):  # x^{1/4} via two sqrt -- the only special op needed
    return iv.sqrt(iv.sqrt(x))


_RT2 = iv.sqrt(iv.mpf([2, 2]))
_PI = iv.pi
_INVPI = 1 / _PI


def _rect(f, a: mp.mpf, b: mp.mpf, n: int):
    """Rigorous interval rectangle-sum enclosure of int_a^b f, f mapping an iv box to an iv value.

    The rounded nodes tile [a, b] exactly (shared endpoints, last node = b), and each box is weighted by an
    interval enclosure of its own width, so rounding in the nodes cannot open a gap or shrink a weight."""
    a, b = mp.mpf(a), mp.mpf(b)
    h = (b - a) / n
    nodes = [a + h * i for i in range(n)] + [b]
    total = iv.mpf([0, 0])
    for lo, hi in pairwise(nodes):
        total = total + f(iv.mpf([lo, hi])) * (iv.mpf([hi, hi]) - iv.mpf([lo, lo]))
    return total


def f_lower(slo: mp.mpf, shi: mp.mpf, n_a: int = 160, n_b: int = 240):
    """Rigorous enclosure of the lower bound F_lo(s) for s in [slo, shi]; if its left end > 0, F' > 0 there."""
    s = iv.mpf([slo, shi])
    q = 2 * s - 1
    one = iv.mpf([1, 1])
    # (1/2)^{1/4}, the substitution endpoint, enclosed: the positive integrands of the lower bound b_lo stop
    # at its lower end (a gap only lowers b_lo), those of the upper bound bp_up run to its upper end (an
    # overlap only raises bp_up), so the junction w = 1/2 is safe either way.
    up_iv = 1 / _qrt(iv.mpf([2, 2]))
    up_lo, up_hi = mp.mpf(up_iv.a), mp.mpf(up_iv.b)

    a_int = _rect(lambda w: 1 / _qrt(s * s - (one - s) ** 2 * w * w), 0, 1, n_a)
    ap_int = _rect(
        lambda w: (
            (iv.mpf([-1, -1]) / 4)
            * (2 * s + 2 * (one - s) * w * w)
            / ((s * s - (one - s) ** 2 * w * w) * _qrt(s * s - (one - s) ** 2 * w * w))
        ),
        0,
        1,
        n_a,
    )

    def blo_t(t):  # w = t^4 on [0,1/2]; the bhat^{-1/4}*4t^3 = 4t^2/qrt(2s-q t^4) (power cancelled)
        return (_RT2 * _INVPI) * (
            4 * t**3 / _qrt((one - t**4) * (one + q * t**4)) + 4 * t * t / _qrt(2 * s - q * t**4)
        )

    def blo_r(r):  # w = 1 - r^4 on [1/2,1]
        return (_RT2 * _INVPI) * (
            4 * r * r / _qrt(one + q * (one - r**4))
            + 4 * r**3 / _qrt((one - r**4) * (2 * s - q * (one - r**4)))
        )

    b_lo = _rect(blo_t, 0, up_lo, n_b) + _rect(blo_r, 0, up_lo, n_b)

    def bpu_t(t):
        p1 = (one - t**4) * (one + q * t**4)
        return _INVPI * (
            4 * t**7 * (one - t**4) / (p1 * _qrt(p1))
            + 4 * t * t * (one - t**4) / ((2 * s - q * t**4) * _qrt(2 * s - q * t**4))
        )

    def bpu_r(r):
        p2 = (one - r**4) * (2 * s - q * (one - r**4))
        return _INVPI * (
            4 * r * r * (one - r**4) / ((one + q * (one - r**4)) * _qrt(one + q * (one - r**4)))
            + 4 * r**7 * (one - r**4) / (p2 * _qrt(p2))
        )

    bp_up = _rect(bpu_t, 0, up_hi, n_b) + _rect(bpu_r, 0, up_hi, n_b)

    pre = iv.sqrt(2 / _PI)
    return pre * (
        (4 / _PI) * (-a_int + (one - s) * ap_int)
        + (iv.mpf([3, 3]) / 2) * (1 / _qrt(q)) * b_lo
        - _qrt(q * q * q) * bp_up
    )


_ONE = iv.mpf([1, 1])


def ellipe_iv(m):
    """Enclosure of E(m) = int_0^{pi/2} sqrt(1 - m sin^2) for a thin interval m in [0, 1), by the AGM.

    a_{n+1} = (a_n + g_n)/2, g_{n+1} = sqrt(a_n g_n), c_{n+1} = c_n^2 / (4 a_{n+1}) (no cancellation), and
    K = pi / (2 AGM) with g_n <= AGM <= a_n, E = K (1 - sum_{n>=0} 2^{n-1} c_n^2).  Once c_N <= g_N, every later
    c shrinks by at least 4, so the neglected tail of the sum is at most 2^{N-1} c_N^2 / 7."""
    a, g, c = _ONE, iv.sqrt(_ONE - m), iv.sqrt(m)
    total = c * c / 2
    p = _ONE
    for _ in range(80):
        a, g = (a + g) / 2, iv.sqrt(a * g)
        c = c * c / (4 * a)
        total = total + p * c * c
        if c.b < mp.mpf("1e-40") and c.b <= g.a:
            break
        p = 2 * p
    else:
        raise RuntimeError(f"AGM did not converge for m = {m}")
    tail = iv.mpf([0, (p * c * c / 7).b])
    k = iv.mpf([(_PI / (2 * a)).a, (_PI / (2 * g)).b])
    return k * (_ONE - total - tail)


def h_bounds(rho: mp.mpf, n: int) -> tuple[mp.mpf, mp.mpf]:
    """Step (A): [lo, hi] containing h(rho) = E sqrt(cos^2 Phi + rho cos^2 Psi), rho in (0, 1] exact."""
    r = iv.mpf([rho, rho])

    def f(i: int):
        x = iv.cos(_PI * i / (2 * n)) ** 2
        return (2 / _PI) * iv.sqrt(x + r) * ellipe_iv(r / (x + r))

    interior = iv.mpf([0, 0])
    for i in range(1, n):
        interior = interior + f(i)
    lo = (interior + (2 / _PI) * iv.sqrt(r)) / n  # f(pi/2) = (2/pi) sqrt(rho) exactly
    hi = (interior + f(0)) / n
    return mp.mpf(lo.a), mp.mpf(hi.b)


def beta_lower(nodes: list[mp.mpf], hlo: list[mp.mpf], n: int) -> mp.mpf:
    """Step (B): a lower bound for beta_odd = F(1) from the chord minorant of the concave h.

    On [0, 1/2], M = b = (u(2-u))^{-1/2} and m = a = (1-u^2)^{-1/2}; u = t^4 turns sqrt(b) du into
    4 t^2 (2-t^4)^{-1/4} dt and a/b into t^2 sqrt((2-t^4)/(1-t^8)), which increases in t."""
    ys = list(hlo)
    for i in range(
        len(ys) - 2, -1, -1
    ):  # a running minimum from the right keeps lower bounds, now monotone
        ys[i] = min(ys[i], ys[i + 1])

    def chord_lower(rho: mp.mpf) -> mp.mpf:
        rho = min(max(rho, mp.mpf(0)), mp.mpf(1))
        for (x0, y0), (x1, y1) in pairwise(zip(nodes, ys, strict=True)):
            if x0 <= rho <= x1:
                x0i, y0i = iv.mpf([x0, x0]), iv.mpf([y0, y0])
                lam = (iv.mpf([rho, rho]) - x0i) / (iv.mpf([x1, x1]) - x0i)
                return mp.mpf((y0i + lam * (iv.mpf([y1, y1]) - y0i)).a)
        raise ValueError(f"rho = {rho} outside [0, 1]")

    def integrand(t):
        t4 = t**4
        rho = t * t * iv.sqrt((2 - t4) / (_ONE - t4 * t4))
        return (
            4 * t * t / _qrt(2 - t4) * chord_lower(mp.mpf(rho.a))
        )  # chord nondecreasing: its value at rho.a

    up = 1 / _qrt(
        iv.mpf([2, 2])
    )  # t at u = 1/2; stopping at its lower end only lowers the positive integral
    return mp.mpf((2 * iv.sqrt(2 / _PI) * _rect(integrand, 0, mp.mpf(up.a), n)).a)


def tangent_upper(slo: mp.mpf, shi: mp.mpf, h1up: mp.mpf, n1: int = 300, n2: int = 1500) -> mp.mpf:
    """Step (D): an upper bound, valid for every s in [slo, shi], of the tangent majorant of F(s).

    Both integrals run over domains enlarged to cover every s in the box (positive integrands), which needs
    shi - slo < slo - 1/2 to keep 2s - 1 - t^4 > 0."""
    s = iv.mpf([slo, shi])
    one_copy = _rect(lambda u: 1 / _qrt(s * s - u * u), 0, mp.mpf((1 - iv.mpf([slo, slo])).b), n1)

    def two_copy(t):  # 3 sqrt(b) + a/sqrt(b) at v = t^4, times dv = 4 t^3 dt
        t4 = t**4
        w = 2 * s - t4
        return 12 * t * t / _qrt(w) + 4 * t4 * _qrt(w) / iv.sqrt((2 * s - 1 - t4) * (1 + t4))

    tmax = _qrt(iv.mpf([shi, shi]) - iv.mpf([0.5, 0.5]))
    two = _rect(two_copy, 0, mp.mpf(tmax.b), n2)
    h1 = iv.mpf([h1up, h1up])
    return mp.mpf((2 * iv.sqrt(2 / _PI) * ((2 / _PI) * one_copy + (h1 / 4) * two)).b)


def certify_tangent(slo: mp.mpf, shi: mp.mpf, h1up: mp.mpf, beta_lo: mp.mpf, margin: mp.mpf):
    """Step (D) driver: bisect [slo, shi] until the majorant clears beta_lo by `margin` on every box."""
    stack, n_boxes, worst = [(slo, shi)], 0, mp.mpf(-1)
    while stack:
        a, b = stack.pop()
        if b - a < (a - mp.mpf(0.5)) / 2:
            up = tangent_upper(a, b, h1up)
            if up < beta_lo - margin:
                n_boxes, worst = n_boxes + 1, max(worst, up)
                continue
        if b - a < mp.mpf("1e-5"):
            print(f"    FAIL: s in [{mp.nstr(a, 8)}, {mp.nstr(b, 8)}] not certified")
            return False, n_boxes, worst
        m = (a + b) / 2
        stack += [(a, m), (m, b)]
    return True, n_boxes, worst


def main() -> int:
    s1, s2 = mp.mpf("0.575"), mp.mpf("0.90")
    print(
        "(A) h(rho) = E sqrt(cos^2 Phi + rho cos^2 Psi) at rho = k/8 (AGM for E(m), monotone sums, 2000 cells):"
    )
    nodes = [mp.mpf(k) / 8 for k in range(9)]
    two_pi = 2 / _PI
    hlo, hhi = [mp.mpf(two_pi.a)], [mp.mpf(two_pi.b)]
    for rho in nodes[1:]:
        lo, hi = h_bounds(rho, 2000)
        hlo.append(lo)
        hhi.append(hi)
        print(f"    h({mp.nstr(rho, 3)}) in [{mp.nstr(lo, 8)}, {mp.nstr(hi, 8)}]")

    beta_lo = beta_lower(nodes, hlo, 12000)
    print(
        f"(B) chord minorant: beta_odd >= {mp.nstr(beta_lo, 8)}   (beta_odd = 0.928019304807931...)"
    )

    c0 = 16 / iv.gamma(iv.mpf([0.25, 0.25])) ** 2
    up_half = mp.mpf((c0 / iv.sqrt(iv.mpf([2, 2]))).b)
    up_tri = mp.mpf((c0 * iv.sqrt(iv.mpf([s1, s1]))).b)
    ok_c = up_half < beta_lo and up_tri < beta_lo
    print(
        f"(C) [0, 1/2]: F <= c_0/sqrt2 <= {mp.nstr(up_half, 8)};  [1/2, {mp.nstr(s1, 4)}]: F <= c_0 sqrt(s) <= "
        f"{mp.nstr(up_tri, 8)}  {'ok' if ok_c else 'FAIL'}"
    )

    ok_d, n_boxes, worst_d = certify_tangent(s1, s2, hhi[-1], beta_lo, mp.mpf("0.001"))
    print(
        f"(D) [{mp.nstr(s1, 4)}, {mp.nstr(s2, 3)}]: tangent majorant <= {mp.nstr(worst_d, 8)} on {n_boxes} boxes, "
        f"each below beta_lo - 0.001  {'ok' if ok_d else 'FAIL'}"
    )

    lo, hi, k = s2, mp.mpf("1.0"), 32
    h = (hi - lo) / k
    worst = mp.mpf("1e9")
    ok_e = True
    print(
        f"(E) validated lower bound on F'(s) over {k} subintervals of [0.90, 1] (left end must be > 0):"
    )
    for i in range(k):
        a, b = lo + h * i, lo + h * (i + 1)
        left = mp.mpf(f_lower(a, b).a)
        worst = min(worst, left)
        ok_e &= left > 0
        if i % 8 == 0 or left <= 0:
            print(
                f"    s in [{mp.nstr(a, 5)}, {mp.nstr(b, 5)}]:  F' >= {mp.nstr(left, 5)}  {'ok' if left > 0 else 'FAIL'}"
            )
    print(
        f"    worst-case F'_lower over [0.90, 1] = {mp.nstr(worst, 6)}  {'ok' if ok_e else 'FAIL'}"
    )

    ok = ok_c and ok_d and ok_e
    print("=" * 72)
    print(
        "RESULT:",
        "PASS -- F < beta_odd on [0, 0.90] and F' > 0 on [0.90, 1], all by validated interval arithmetic: "
        "max_{[0,1]} F = F(1) = beta_odd"
        if ok
        else "FAIL",
    )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
