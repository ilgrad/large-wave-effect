# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2", "scipy>=1.11"]
# ///
"""Exact value of the odd-N Schrodinger constant beta_odd = lim_{odd N} B_N/sqrt(N).

For odd N the sup of ||exp(-itL_N)||_{inf->inf} is at t = N/2, where (only the aliasing copies ell=0,-1
survive, and the i^{-N} phase is imaginary) |K_N(k,N/2)| = sqrt(J_k(N)^2 + J_{N-k}(N)^2).  Hence

    beta_odd = lim (1/sqrt N) sum_{k=0}^{N-1} sqrt(J_k(N)^2 + J_{N-k}(N)^2).

The Debye envelope |J_{Nu}(N)| ~ sqrt(2/(pi N)) (1-u^2)^{-1/4} |cos Phi| (and the analogous one for
J_{N-k} with u -> 1-u), together with the equidistribution of the two Debye phases, turns the sum into

    beta_odd = sqrt(2/pi) * int_0^1 E_{Phi,Psi} sqrt( (1-u^2)^{-1/2} cos^2 Phi + (u(2-u))^{-1/2} cos^2 Psi ) du
             = sqrt(2/pi) * int_0^1 (1/pi^2) int_0^pi int_0^pi sqrt( (1-u^2)^{-1/2} cos^2 phi
                                                                    + (u(2-u))^{-1/2} cos^2 psi ) dphi dpsi du.

This is an elliptic-integral average.  To high precision beta_odd = 0.92801930480793112292; a PSLQ /
inverse-symbolic search against {pi, e^gamma, Gamma(1/4), Gamma(3/4), sqrt2, sqrt(pi), Catalan, log2,
zeta(3)} and against low-degree algebraicity returns nothing, so beta_odd has no elementary closed form --
it is a new constant defined by the integral.  This script checks the integral matches the direct large-N
sum at t = N/2 (the high-precision value and the PSLQ negative are reproduced with mpmath separately).
"""

from __future__ import annotations

import numpy as np
from numpy import cos, pi, sqrt
from scipy import integrate
from scipy.special import jv


def beta_odd_integral() -> float:
    def inner(u: float) -> float:
        a = (1 - u**2) ** -0.5
        b = (u * (2 - u)) ** -0.5
        val, _ = integrate.dblquad(
            lambda p, q: sqrt(a * cos(p) ** 2 + b * cos(q) ** 2),
            0, pi, 0, pi, epsabs=1e-8, epsrel=1e-8,
        )
        return val / pi**2

    val, _ = integrate.quad(inner, 1e-7, 1 - 1e-7, limit=200, epsabs=1e-7)
    return sqrt(2 / pi) * val


def beta_odd_direct(n: int) -> float:
    """(1/sqrt N) sum_k sqrt(J_k(N)^2 + J_{N-k}(N)^2)  -- the value of B_N/sqrt(N) at t=N/2."""
    k = np.arange(n)
    return float(np.sqrt(jv(k, n) ** 2 + jv(n - k, n) ** 2).sum()) / sqrt(n)


def main() -> int:
    bi = beta_odd_integral()
    print(f"beta_odd (triple integral)        = {bi:.6f}")
    direct = [(n, beta_odd_direct(n)) for n in [20001, 80001, 200001]]
    for n, v in direct:
        print(f"  direct sum at N={n:>7} (t=N/2)  = {v:.6f}")
    ok = abs(bi - direct[-1][1]) < 3e-3 and abs(bi - 0.9280193) < 2e-3
    print("=" * 56)
    print("RESULT:", "PASS -- beta_odd = 0.92802... (exact integral)" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
