# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2"]
# ///
"""Persistence of the linear large wave under weak nonlinearity (Proposition dnls-persist).

DNLS on the ring:  i u_t = L_N u + gamma |u|^2 u,  L_N the graph Laplacian (eigenvalues 4 sin^2(pi r/N)).
Duhamel + mass conservation give the rigorous bound

    || u(t) - e^{-i t L_N} u(0) ||_2  <=  |gamma| t || u(0) ||_2^3 .

We integrate the DNLS with a mass-conserving Strang split-step scheme (linear part exact in Fourier,
nonlinear part an exact pointwise phase rotation) and verify (i) the power ||u||_2^2 is conserved, and
(ii) the deviation from the linear flow stays under |gamma| t ||u0||_2^3 -- so the linear large wave
survives for small gamma.
"""

from __future__ import annotations

import numpy as np


def lin_step(u: np.ndarray, lam: np.ndarray, dt: float) -> np.ndarray:
    return np.fft.ifft(np.exp(-1j * lam * dt) * np.fft.fft(u))


def integrate(u0: np.ndarray, lam: np.ndarray, gamma: float, t: float, steps: int) -> np.ndarray:
    dt = t / steps
    u = u0.copy()
    for _ in range(steps):
        u = lin_step(u, lam, dt / 2)
        u = np.exp(-1j * gamma * np.abs(u) ** 2 * dt) * u  # nonlinear phase (exact, mass-preserving)
        u = lin_step(u, lam, dt / 2)
    return u


def main() -> int:
    rng = np.random.default_rng(0)
    n = 64
    lam = 4 * np.sin(np.pi * np.arange(n) / n) ** 2
    ok = True
    print(f"{'gamma':>8} {'t':>6} {'||u||2 drift':>13} {'||u-u_lin||2':>13} {'bound |g|t||u0||2^3':>19}")
    for gamma in [0.01, 0.05, 0.2]:
        u0 = rng.standard_normal(n) + 1j * rng.standard_normal(n)
        m0 = np.linalg.norm(u0)
        t = 8.0
        u = integrate(u0, lam, gamma, t, 4000)
        u_lin = lin_step(u0, lam, t)
        dev = np.linalg.norm(u - u_lin)
        bound = abs(gamma) * t * m0**3
        mass_drift = abs(np.linalg.norm(u) - m0)
        held = dev <= bound + 1e-9 and mass_drift < 1e-9
        ok &= held
        print(f"{gamma:>8.2f} {t:>6.1f} {mass_drift:>13.2e} {dev:>13.4f} {bound:>19.4f}  {'ok' if held else 'FAIL'}")
    print("=" * 64)
    print("RESULT:", "PASS -- Duhamel bound holds, mass conserved (persistence)" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
