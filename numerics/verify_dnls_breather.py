# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2", "scipy>=1.11"]
# ///
"""Existence of discrete breathers on the ring DNLS, by a finite-dimensional variational argument.

A standing wave z_j(t) = e^{-i Omega t} phi_j (real phi) of the focusing ring DNLS
i z_dot = -L_N z + gamma |z|^2 z solves  (L_N + Omega) phi = gamma |phi|^2 phi.  For Omega > 0 the operator
M = L_N + Omega is positive definite, so minimizing <M phi, phi> over the COMPACT set {||phi||_4^4 = 1}
has a (nonzero) minimizer by continuity + compactness -- a rigorous existence proof on the finite ring,
self-contained (no MacKay-Aubry continuation needed).  The Euler-Lagrange equation of the minimizer is
M phi = 2 mu |phi|^2 phi, i.e. exactly the breather equation with gamma = 2 mu; scaling phi covers all
gamma > 0.  For large Omega the minimizer concentrates on one site (IPR -> 1): the on-site breather.

This script computes the minimizer for several Omega and checks (a) it solves the breather equation
(residual ~ 0) and (b) it localizes (participation IPR -> 1), confirming the existence proposition.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import minimize


def laplacian(phi: np.ndarray) -> np.ndarray:
    return 2 * phi - np.roll(phi, -1) - np.roll(phi, 1)


def breather(n: int, omega: float) -> tuple[float, float, float]:
    def m(p: np.ndarray) -> np.ndarray:
        return laplacian(p) + omega * p

    cons = {"type": "eq", "fun": lambda p: np.sum(p**4) - 1.0}
    p0 = np.zeros(n)
    p0[0], p0[1], p0[-1] = 1.0, 0.3, 0.3
    p0 /= np.sum(p0**4) ** 0.25
    res = minimize(lambda p: p @ m(p), p0, constraints=cons, method="SLSQP",
                   options={"maxiter": 800, "ftol": 1e-13})
    p = res.x
    gamma = float(res.fun)  # ||p||_4^4 = 1 => <Mp,p> = 2 mu => gamma = 2 mu = <Mp,p>
    resid = float(np.linalg.norm(m(p) - gamma * p**3))
    ipr = float(np.sum(p**4) / np.sum(p**2) ** 2)  # 1 = single-site, 1/N = uniform
    return gamma, resid, ipr


def main() -> int:
    n = 64
    print(f"Discrete breather on the ring DNLS (N={n}): variational minimizer solves the breather eqn.")
    print(f"  {'Omega':>6} {'gamma':>8} {'eqn residual':>13} {'IPR':>7}")
    ok = True
    for omega in [1.0, 3.0, 8.0, 20.0, 50.0]:
        gamma, resid, ipr = breather(n, omega)
        if resid > 1e-4:
            ok = False
        print(f"  {omega:>6.1f} {gamma:>8.4f} {resid:>13.2e} {ipr:>7.4f}")
    # localization grows with Omega (on-site breather in the large-Omega / anti-continuum limit)
    if not (breather(n, 50.0)[2] > breather(n, 1.0)[2] + 0.3):
        ok = False
    print("=" * 52)
    print("RESULT:", "PASS -- breather exists (variational); localizes as Omega grows" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
