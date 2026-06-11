# /// script
# requires-python = ">=3.11"
# dependencies = ["mpmath>=1.3"]
# ///
"""Validated-numerics proof that F'(s) > 0 on [0.9366, 1] -- closing the t <= N/2 maximizer (Proposition bn-profile).

This is the residual step of Proposition bn-profile: the tangent-to-concave-h majorant already gives
F(s) < beta_odd analytically on [1/2, 0.9366] (87% of the window), so the only gap was the monotonicity of
the explicit profile F on [0.9366, 1].  Here it is closed by RIGOROUS interval arithmetic.

F(s) = sqrt(2/pi)[ (4/pi)(1-s) A(s) + (2s-1)^{3/4} B(s) ],   q = 2s-1,
  A(s) = int_0^1 (s^2-(1-s)^2 w^2)^{-1/4} dw          (smooth integrand on [0,1]),
  B(s) = int_0^1 g(ahat, bhat) dw,   ahat=[(1-w)(1+qw)]^{-1/2}, bhat=[w(2s-qw)]^{-1/2},
  g(a,b) = E_{Phi,Psi} sqrt(a cos^2 Phi + b cos^2 Psi).

F'(s) = sqrt(2/pi)[ (4/pi)(-A + (1-s)A') + (3/2) q^{-1/4} B + q^{3/4} B' ].  Two ELEMENTARY inequalities
remove the elliptic g (so only sqrt is needed, and the bound is rigorous):

  * QM-AM lower bound on g:   sqrt(a cos^2 Phi + b cos^2 Psi) >= (a^{1/2}|cos Phi| + b^{1/2}|cos Psi|)/sqrt2,
    hence  g(a,b) >= (sqrt2/pi)(sqrt(a) + sqrt(b)),  giving  B(s) >= B_lo(s) = int (sqrt2/pi)(sqrt(ahat)+sqrt(bhat)) dw;
  * derivative bound:   d_a g = E[cos^2 Phi / (2 sqrt(...))] <= E[|cos Phi|/(2 sqrt a)] = 1/(pi sqrt a),
    and with d_s ahat = -(1-w)w[(1-w)(1+qw)]^{-3/2} (likewise bhat),
        |B'(s)| <= Bp_up(s) = (1/pi) int w(1-w)[ ((1-w)(1+qw))^{-5/4} + (w(2s-qw))^{-5/4} ] dw.

Since B' < 0 (ahat, bhat decrease in s), q^{3/4} B' >= -q^{3/4} Bp_up, so

    F'(s) >= F_lo(s) := sqrt(2/pi)[ (4/pi)(-A + (1-s)A') + (3/2) q^{-1/4} B_lo - q^{3/4} Bp_up ].

All four integrals are elementary; F_lo(s) is enclosed by a rigorous interval rectangle sum (mpmath.iv:
the interval evaluation over each box contains the integrand's range, so box x width contains the
contribution -- a guaranteed enclosure).  The w^{-1/4} endpoint singularities are removed by w = t^4 on
[0,1/2] and w = 1 - r^4 on [1/2,1], with the t/r powers cancelled by hand so no 0*inf appears.  Treating s
itself as an interval over subintervals certifies F'(s) > 0 for ALL s in the subinterval.

RESULT: F_lo(s) > 0 on every subinterval of [0.9366, 1] (worst-case lower bound ~ 0.051), so F is strictly
increasing there; with F(1) = beta_odd and the analytic majorant below 0.9366, max_{[0,1]} F = F(1) =
beta_odd is PROVED.  Companion to verify_bn_profile_max.py (the analytic majorant + endpoints).
"""

from __future__ import annotations

import mpmath as mp

iv = mp.iv
iv.prec = 90


def _qrt(x):  # x^{1/4} via two sqrt -- the only special op needed
    return iv.sqrt(iv.sqrt(x))


_RT2 = iv.sqrt(iv.mpf([2, 2]))
_PI = iv.pi
_INVPI = 1 / _PI


def _rect(f, a: mp.mpf, b: mp.mpf, n: int):
    """Rigorous interval rectangle-sum enclosure of int_a^b f, f mapping an iv box to an iv value."""
    a, b = mp.mpf(a), mp.mpf(b)
    h = (b - a) / n
    total = iv.mpf([0, 0])
    for i in range(n):
        total = total + f(iv.mpf([a + h * i, a + h * (i + 1)])) * h
    return total


def f_lower(slo: mp.mpf, shi: mp.mpf, n_a: int = 160, n_b: int = 240):
    """Rigorous enclosure of the lower bound F_lo(s) for s in [slo, shi]; if its left end > 0, F' > 0 there."""
    s = iv.mpf([slo, shi])
    q = 2 * s - 1
    one = iv.mpf([1, 1])
    up = mp.mpf(1) / 2 ** (mp.mpf(1) / 4)  # (1/2)^{1/4}, the substitution endpoint

    a_int = _rect(lambda w: 1 / _qrt(s * s - (one - s) ** 2 * w * w), 0, 1, n_a)
    ap_int = _rect(
        lambda w: (iv.mpf([-1, -1]) / 4)
        * (2 * s + 2 * (one - s) * w * w)
        / ((s * s - (one - s) ** 2 * w * w) * _qrt(s * s - (one - s) ** 2 * w * w)),
        0, 1, n_a,
    )

    def blo_t(t):  # w = t^4 on [0,1/2]; the bhat^{-1/4}*4t^3 = 4t^2/qrt(2s-q t^4) (power cancelled)
        return (_RT2 * _INVPI) * (4 * t**3 / _qrt((one - t**4) * (one + q * t**4)) + 4 * t * t / _qrt(2 * s - q * t**4))

    def blo_r(r):  # w = 1 - r^4 on [1/2,1]
        return (_RT2 * _INVPI) * (
            4 * r * r / _qrt(one + q * (one - r**4)) + 4 * r**3 / _qrt((one - r**4) * (2 * s - q * (one - r**4)))
        )

    b_lo = _rect(blo_t, 0, up, n_b) + _rect(blo_r, 0, up, n_b)

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

    bp_up = _rect(bpu_t, 0, up, n_b) + _rect(bpu_r, 0, up, n_b)

    pre = iv.sqrt(2 / _PI)
    return pre * ((4 / _PI) * (-a_int + (one - s) * ap_int) + (iv.mpf([3, 3]) / 2) * (1 / _qrt(q)) * b_lo - _qrt(q * q * q) * bp_up)


def main() -> int:
    lo, hi, k = mp.mpf("0.9366"), mp.mpf("1.0"), 20
    h = (hi - lo) / k
    worst = mp.mpf("1e9")
    ok = True
    print("Validated interval lower bound on F'(s) over subintervals of [0.9366, 1] (left end must be > 0):")
    for i in range(k):
        a, b = lo + h * i, lo + h * (i + 1)
        left = f_lower(a, b).a
        worst = min(worst, left)
        ok &= left > 0
        if i % 5 == 0 or left <= 0:
            print(f"  s in [{mp.nstr(a, 5)}, {mp.nstr(b, 5)}]:  F' >= {mp.nstr(left, 5)}  {'ok' if left > 0 else 'FAIL'}")
    print("=" * 72)
    print(f"worst-case F'_lower over [0.9366, 1] = {mp.nstr(worst, 6)}")
    print("RESULT:", "PASS -- F'(s) > 0 proved on [0.9366, 1] by validated interval arithmetic; with the "
          "analytic majorant below, max F = F(1) = beta_odd is a theorem" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
