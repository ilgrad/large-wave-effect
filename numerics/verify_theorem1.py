# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2"]
# ///
"""Numerical verification of Theorem 1 (docs/bead_chain_formulas_claude.md, §8).

Claim: for the periodic bead chain (m = kappa = 1) the maximal deflection of a
single bead under fixed energy E0 on the zero-mean subspace is

    M_N(E0) = sqrt( 2 * E0 * (L_N^+)_jj ) = sqrt( E0 * (N^2 - 1) / (6N) ),

independent of time t and of node j, because (L_N^+)_jj = (N^2 - 1) / (12N).

Checks:
  1. Identities: sum csc^2(pi r/N) = (N^2-1)/3, sum 1/lambda_r = (N^2-1)/12,
     (L_N^+)_jj = (N^2-1)/(12N) -- spectral form AND dense pinv agree.
  2. Invariance: a(t)^T H^+ a(t) = (L_N^+)_jj for all sampled t, j.
  3. Reachability: build the optimal IC, evolve exactly, confirm the focusing
     bead hits M_N, energy is conserved, and no node ever exceeds M_N.
  4. Contrast: the diploma step IC stays O(1), far below the sqrt(N) ceiling.

All operators f(L_N) act diagonally in the DFT basis (circulant), applied via FFT.
"""

from __future__ import annotations

import math

import numpy as np

RTOL = 1e-9


def operators(n: int):
    """Return DFT-domain multipliers and an apply() for the chain of size n."""
    r = np.arange(n)
    theta = np.pi * r / n
    lam = 4.0 * np.sin(theta) ** 2          # eigenvalues of L_N, lam[0] = 0
    omega = 2.0 * np.abs(np.sin(theta))     # sqrt(lam)
    nz = lam > 1e-12                         # non-zero (drop the r=0 mode)

    lpinv = np.zeros(n)
    lpinv[nz] = 1.0 / lam[nz]                # L_N^+ multiplier

    def apply(mult: np.ndarray, x: np.ndarray) -> np.ndarray:
        return np.real(np.fft.ifft(mult * np.fft.fft(x)))

    def m_cos(t: float) -> np.ndarray:       # cos(t sqrt(L))
        return np.cos(omega * t)

    def m_lm12_sin(t: float) -> np.ndarray:  # L^{-1/2} sin(t sqrt(L)), 0 at r=0
        out = np.zeros(n)
        out[nz] = np.sin(omega[nz] * t) / omega[nz]
        return out

    def m_sqrtl_sin(t: float) -> np.ndarray:  # sqrt(L) sin(t sqrt(L)) for velocity
        return omega * np.sin(omega * t)

    return {
        "n": n, "lam": lam, "omega": omega, "lpinv": lpinv, "apply": apply,
        "m_cos": m_cos, "m_lm12_sin": m_lm12_sin, "m_sqrtl_sin": m_sqrtl_sin,
    }


def diag_lpinv_spectral(op) -> float:
    """(L_N^+)_jj = (1/N) sum_r lpinv_r  (same for every j by translation symmetry)."""
    return float(np.mean(op["lpinv"]))


def diag_lpinv_dense(n: int) -> np.ndarray:
    """Diagonal of the Moore-Penrose pseudoinverse of the dense circulant L_N."""
    lap = 2.0 * np.eye(n) - np.eye(n, k=1) - np.eye(n, k=-1)
    lap[0, -1] = lap[-1, 0] = -1.0
    return np.diag(np.linalg.pinv(lap))


def energy(op, u: np.ndarray, udot: np.ndarray) -> float:
    """E = 1/2 udot^T udot + 1/2 u^T L u  (m = kappa = 1)."""
    lu = op["apply"](op["lam"], u)
    return float(0.5 * udot @ udot + 0.5 * u @ lu)


def check(label: str, got: float, want: float, rtol: float = RTOL, atol: float = 0.0) -> bool:
    # Relative tolerance when comparing to a non-zero target; absolute tolerance
    # (atol) when the target is exactly zero (std, energy drift, overshoot).
    diff = abs(got - want)
    ok = diff <= atol + rtol * abs(want)
    metric = f"abs={diff:.2e}" if abs(want) < 1e-30 else f"rel={diff / abs(want):.2e}"
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}: got={got:.12g} want={want:.12g} {metric}")
    return ok


def check_identities(n: int) -> bool:
    print(f"(1) Identities, N={n}")
    op = operators(n)
    r = np.arange(1, n)
    csc2 = np.sum(1.0 / np.sin(np.pi * r / n) ** 2)
    inv_lam = np.sum(1.0 / (4.0 * np.sin(np.pi * r / n) ** 2))
    diag_spec = diag_lpinv_spectral(op)
    diag_dense = diag_lpinv_dense(n)
    ok = True
    ok &= check("sum csc^2 = (N^2-1)/3", csc2, (n * n - 1) / 3)
    ok &= check("sum 1/lambda = (N^2-1)/12", inv_lam, (n * n - 1) / 12)
    ok &= check("(L+)_jj spectral = (N^2-1)/(12N)", diag_spec, (n * n - 1) / (12 * n))
    ok &= check("dense pinv diag is constant", float(diag_dense.std()), 0.0, rtol=0.0, atol=1e-6)
    ok &= check("dense pinv diag = spectral", float(diag_dense.mean()), diag_spec, rtol=1e-7)
    return ok


def check_invariance(n: int) -> bool:
    print(f"(2) Invariance a(t)^T H^+ a(t) = (L+)_jj for all t, j, N={n}")
    op = operators(n)
    apply, lpinv = op["apply"], op["lpinv"]
    target = diag_lpinv_spectral(op)
    ok = True
    rng = np.random.default_rng(0)
    js = [0, 1, n // 3, n // 2, n - 1]
    ts = [0.0, 0.37, 1.0, 3.1415, 12.5, 99.0, float(rng.uniform(0, 100))]
    worst = 0.0
    for j in js:
        e = np.zeros(n)
        e[j] = 1.0
        for t in ts:
            c = apply(op["m_cos"](t), e)            # cos(t sqrt L) e_j
            s = apply(op["m_lm12_sin"](t), e)       # L^{-1/2} sin(t sqrt L) e_j
            val = c @ apply(lpinv, c) + s @ s       # a^T H^+ a
            worst = max(worst, abs(val - target) / target)
    ok &= check(f"max rel dev over {len(js)}x{len(ts)} (t,j) samples", worst, 0.0, rtol=0.0, atol=1e-9)
    return ok


def check_reachability(n: int, e0: float = 1.0, tstar: float = 7.3, j0: int = 0) -> bool:
    print(f"(3) Reachability: optimal IC focuses to M_N, N={n}, E0={e0}, t*={tstar}, j0={j0}")
    op = operators(n)
    apply, lpinv = op["apply"], op["lpinv"]
    diag = diag_lpinv_spectral(op)
    m_n = math.sqrt(2 * e0 * diag)

    e = np.zeros(n)
    e[j0] = 1.0
    cstar = apply(op["m_cos"](tstar), e)
    sstar = apply(op["m_lm12_sin"](tstar), e)
    a_h_a = cstar @ apply(lpinv, cstar) + sstar @ sstar
    scale = math.sqrt(2 * e0 / a_h_a)
    u0 = scale * apply(lpinv, cstar)     # u0* = c L^+ cos(t* sqrt L) e_j0
    v0 = scale * sstar                   # v0* = c L^{-1/2} sin(t* sqrt L) e_j0

    def state(t: float):
        u = apply(op["m_cos"](t), u0) + apply(op["m_lm12_sin"](t), v0)
        udot = -apply(op["m_sqrtl_sin"](t), u0) + apply(op["m_cos"](t), v0)
        return u, udot

    ok = True
    u_star, _ = state(tstar)
    ok &= check("u_{j0}(t*) = M_N", float(u_star[j0]), m_n)
    ok &= check("M_N = sqrt(E0(N^2-1)/6N)", m_n, math.sqrt(e0 * (n * n - 1) / (6 * n)))

    tgrid = np.linspace(0.0, 3.0 * tstar, 1500)
    e_dev, global_max = 0.0, 0.0
    for t in tgrid:
        u, udot = state(t)
        e_dev = max(e_dev, abs(energy(op, u, udot) - e0) / e0)
        global_max = max(global_max, float(np.max(np.abs(u))))
    ok &= check("energy conserved on trajectory", e_dev, 0.0, rtol=0.0, atol=1e-8)
    # ceiling: no node, no time may exceed M_N (allow tiny grid overshoot tol)
    overshoot = max(0.0, (global_max - m_n) / m_n)
    ok &= check("global max over (j,t) does not exceed M_N", overshoot, 0.0, rtol=0.0, atol=1e-6)
    print(f"       global max/M_N = {global_max / m_n:.6f}  (=1 at the focusing point)")
    return ok


def contrast_step(ns: list[int]) -> None:
    print("(4) Contrast: diploma step IC vs the sqrt(N) ceiling (illustrative)")
    print(f"  {'N':>6} {'E0':>10} {'M_N(ceil)':>12} {'peak(grid)':>12} {'peak/M_N':>10} {'peak/sqrtN':>11}")
    for n in ns:
        op = operators(n)
        diag = diag_lpinv_spectral(op)
        alpha = np.ones(n)
        alpha[0] = 0.0
        alpha = alpha - alpha.mean()                 # zero-mean displacement
        e0 = 0.5 * alpha @ op["apply"](op["lam"], alpha)
        m_n = math.sqrt(2 * e0 * diag)
        tgrid = np.linspace(0.0, 4.0 * n, 4000)      # accessible-time grid scan
        peak = 0.0
        for t in tgrid:
            u = op["apply"](op["m_cos"](t), alpha)   # v0 = 0
            peak = max(peak, float(np.max(np.abs(u))))
        print(f"  {n:>6} {e0:>10.4f} {m_n:>12.4f} {peak:>12.4f} "
              f"{peak / m_n:>10.4f} {peak / math.sqrt(n):>11.5f}")
    print("  -> observed peak is O(1) (lower bound on sup_t); ceiling grows ~sqrt(N):")
    print("     simple ICs stay far from the variational optimum (cf. log N growth, §10).")


def main() -> int:
    print("=" * 72)
    print("Theorem 1 numerical verification  (periodic bead chain, m=kappa=1)")
    print("=" * 72)
    all_ok = True
    for n in (8, 16, 32, 64, 128, 257, 1024):
        all_ok &= check_identities(n)
    for n in (16, 64, 257):
        all_ok &= check_invariance(n)
    for n in (16, 64, 256):
        all_ok &= check_reachability(n)
    print()
    contrast_step([16, 32, 64, 128, 256, 512, 1024])
    print()
    print("=" * 72)
    print("RESULT:", "ALL CHECKS PASSED" if all_ok else "SOME CHECKS FAILED")
    print("=" * 72)
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
