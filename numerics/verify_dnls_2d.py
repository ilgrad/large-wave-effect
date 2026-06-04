# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2"]
# ///
"""Two-dimensional DNLS: a discrete breather self-traps in 2D, where the LINEAR large wave cannot.

Part I showed the linear large wave is one-dimensional: on the d-torus the ceiling U_N^{(d)} is O(1)
for d >= 2, so a localized linear disturbance stays bounded and spreads. Nonlinearity changes this. On
the 2-torus (Z_N)^2 the DNLS
    i psi_dot_{j,k} = -(L_N^{(2)} psi)_{j,k} + gamma |psi_{j,k}|^2 psi_{j,k},
with the BCCB Laplacian L_N^{(2)} (eigenvalues lambda_{r,s} = 4 sin^2(pi r/N) + 4 sin^2(pi s/N), band
[0,8]), still self-traps a single-site excitation into a stationary, localized 2D discrete breather once
the on-site nonlinear frequency exceeds the 2D band top, gamma|psi|^2 > 8.

A norm-conserving 2D split-step integrator (2D FFT linear step, pointwise nonlinear step) gives, as the
correctness test, exact norm conservation. We verify:
 (A) norm conserved to ~1e-12;
 (B) the 2D participation ratio P = (sum|psi|^2)^2 / sum|psi|^4 collapses from ~N^2 (dispersed) to O(1)
     (trapped) as gamma crosses the 2D band top ~8 -- a self-trapping transition in two dimensions;
 (C) the trapped profile is a stationary 2D discrete breather (a localized hump steady in time).

So nonlinearity restores a large, localized amplitude in 2D, by a mechanism (self-trapping) absent from
the linear theory -- the nonlinear counterpart of the dimensional threshold of Part I.
"""

from __future__ import annotations

import numpy as np


def step2d(psi: np.ndarray, dt: float, gamma: float, lam: np.ndarray) -> np.ndarray:
    psi = psi * np.exp(-1j * gamma * np.abs(psi) ** 2 * dt)
    return np.fft.ifft2(np.exp(1j * lam * dt) * np.fft.fft2(psi))


def participation(psi: np.ndarray) -> float:
    a2 = np.abs(psi) ** 2
    return float(a2.sum() ** 2 / np.sum(a2**2))


def lam2d(n: int) -> np.ndarray:
    d = 4.0 * np.sin(np.pi * np.arange(n) / n) ** 2
    return d[:, None] + d[None, :]


def evolve(n: int, gamma: float, amp: float, dt: float, steps: int) -> dict:
    lam = lam2d(n)
    psi = np.zeros((n, n), dtype=complex)
    psi[0, 0] = amp
    norm0 = float(np.sum(np.abs(psi) ** 2))
    ndrift = 0.0
    profiles = []
    for s in range(steps):
        psi = step2d(psi, dt, gamma, lam)
        if s % (steps // 10) == 0:
            ndrift = max(ndrift, abs(np.sum(np.abs(psi) ** 2) - norm0) / norm0)
        if s >= steps - 4 * (steps // 10) and s % (steps // 10) == 0:
            profiles.append(np.abs(psi).copy())
    return {"P": participation(psi), "ndrift": ndrift, "profiles": profiles, "final": np.abs(psi)}


def main() -> int:
    print("=" * 70)
    print("2D DNLS: discrete breather self-traps where the linear large wave cannot")
    print("=" * 70)
    ok = True
    n, dt, steps = 64, 0.01, 20000

    r = evolve(n, gamma=16.0, amp=1.0, dt=dt, steps=steps)
    print(f"\n(A) Norm conservation (gamma=16): drift = {r['ndrift']:.2e}")
    a_ok = r["ndrift"] < 1e-10
    ok &= a_ok
    print(f"    2D split-step is norm-exact: {'OK' if a_ok else 'FAIL'}")

    print("\n(B) 2D self-trapping transition: participation ratio P vs gamma (N^2 = 4096 sites)")
    print(f"    {'gamma':>6} {'P_final':>9} {'regime':>10}")
    pf = {}
    for g in (2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 16.0):
        pf[g] = evolve(n, gamma=g, amp=1.0, dt=dt, steps=steps)["P"]
        regime = "trapped" if pf[g] < 8 else ("partial" if pf[g] < 60 else "dispersed")
        print(f"    {g:>6.1f} {pf[g]:>9.1f} {regime:>10}")
    b_ok = pf[2.0] > 60 and pf[16.0] < 8 and pf[16.0] < pf[2.0]
    ok &= b_ok
    print(f"    -> transition near the 2D band top gamma~8: dispersed (P~{pf[2.0]:.0f}) -> trapped "
          f"(P~{pf[16.0]:.1f}): {'OK' if b_ok else 'FAIL'}")

    r2 = evolve(n, gamma=16.0, amp=1.0, dt=dt, steps=steps)
    profs = np.array(r2["profiles"])
    var = float(np.max(np.std(profs, axis=0)))
    c_ok = var < 0.05
    ok &= c_ok
    print(f"\n(C) 2D discrete breather (gamma=16): max temporal std of |psi| = {var:.4f} "
          f"(stationary if <0.05) {'OK' if c_ok else 'FAIL'}")

    print("\n" + "=" * 70)
    print("RESULT:", "2D DNLS BREATHER VERIFIED" if ok else "CHECK FAILED")
    print("Nonlinearity self-traps a 2D discrete breather above the 2D band top (gamma~8); the linear")
    print("large wave, by contrast, is bounded for d>=2 (Part I). Nonlinearity restores 2D localization.")
    print("=" * 70)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
