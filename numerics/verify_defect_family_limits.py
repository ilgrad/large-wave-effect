# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2", "scipy>=1.11"]
# ///
"""Family limits of the composite defect: defect_inf(m) = lim_{a} (U - A)(m 2^a), to high precision.

Within each family N = m 2^a (m odd) the node-0 defect converges geometrically,

    d_a = defect_inf(m) - c1 2^{-a} - c2 4^{-a} + ...,

(the bound (2/pi)ln(m/phi(m)) is a-independent, Remark rem:defect).  This script computes d_a with the
orbit-closure optimizer POLISHED by Newton iterations on the stationarity system (gradient C^T(b cos),
Hessian -C^T diag(b sin) C), so each d_a is accurate to ~1e-12, and fits the two-term geometric model.
Resulting limits (5-point fits, residuals ~6e-9):

    defect_inf(3) = 0.10368161...,   defect_inf(5) = 0.03823846...

A PSLQ / inverse-symbolic search over {1, 1/pi, ln2/pi, ln3/pi, ln5/pi, sqrt3/pi, ln(2+sqrt3)/pi, ...}
finds no small relation: like beta_odd, these family limits appear to be new constants.

The check here uses the first four family members (fit residual and limit reproduced to 5e-5).
"""

from __future__ import annotations

import sys
from math import gcd, pi
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from large_wave_effect.cyclotomic import sin_coords

RNG = np.random.default_rng(0)


def totient(n: int) -> int:
    return sum(1 for k in range(1, n + 1) if gcd(k, n) == 1)


def node0_defect(n: int, starts: int) -> float:
    m, big_m = n // 2, totient(2 * n) // 2
    rows = np.array(sin_coords(2 * n, list(range(1, m + 1))), dtype=float)
    c, *_ = np.linalg.lstsq(rows[:big_m].T, rows.T, rcond=None)
    ci = np.rint(c.T)
    assert np.abs(c.T - ci).max() < 1e-6, "non-integer coordinate matrix"
    r = np.arange(1, m + 1)
    omega = 2 * np.sin(pi * r / n)
    b = 2.0 / (n * omega)
    if n % 2 == 0:
        b[-1] = 1.0 / (n * omega[-1])
    u = float(b.sum())
    cf = ci.astype(float)

    def negf(psi: np.ndarray) -> tuple[float, np.ndarray]:
        ph = cf @ psi
        return -float(b @ np.sin(ph)), -(cf.T @ (b * np.cos(ph)))

    best, bpsi = -np.inf, None
    for _ in range(starts):
        res = minimize(negf, RNG.uniform(0, 2 * pi, big_m), jac=True, method="L-BFGS-B")
        if -res.fun > best:
            best, bpsi = -res.fun, res.x
    psi = bpsi.copy()
    for _ in range(60):  # Newton polish: quadratic convergence to the stationary point
        ph = cf @ psi
        g = cf.T @ (b * np.cos(ph))
        h = -(cf.T * (b * np.sin(ph))) @ cf
        try:
            step = np.linalg.solve(h, -g)
        except np.linalg.LinAlgError:
            break
        psi = psi + step
        if np.linalg.norm(step) < 1e-14:
            break
    return u - float(b @ np.sin(cf @ psi))


def fit_limit(ns: list[int], ds: list[float], m: int) -> tuple[float, float]:
    a0 = int(np.log2(ns[0] / m))
    mat = np.array([[1.0, -(2.0 ** -(a0 + i)), -(4.0 ** -(a0 + i))] for i in range(len(ns))])
    sol, *_ = np.linalg.lstsq(mat, np.array(ds), rcond=None)
    return float(sol[0]), float(np.abs(mat @ sol - ds).max())


def main() -> int:
    targets = {3: 0.1036816, 5: 0.0382385}
    ok = True
    print("Family defect limits (node 0), geometric model d_a = L - c1 2^-a - c2 4^-a:")
    for m, ns, starts in [(3, [24, 48, 96, 192], 150), (5, [20, 40, 80, 160], 150)]:
        ds = [node0_defect(n, starts) for n in ns]
        lim, resid = fit_limit(ns, ds, m)
        good = abs(lim - targets[m]) < 5e-5 and resid < 1e-5
        ok &= good
        seq = ", ".join(f"{d:.7f}" for d in ds)
        print(f"  m={m}: d = [{seq}]  ->  L = {lim:.7f} (target {targets[m]}), fit residual {resid:.1e}  "
              f"{'ok' if good else 'FAIL'}")
    print("=" * 64)
    print("RESULT:", "PASS -- family limits reproduced; new constants (no closed form found)" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
