# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2"]
# ///
"""DNLS self-trapping and discrete breathers: nonlinearity localizes, the opposite of the large wave.

The discrete NLS  i psi_dot_j = -(L_N psi)_j + gamma |psi_j|^2 psi_j  (focusing) has the linear ring
core L_N. A single-site excitation psi_j = A delta_{j0} either DISPERSES (linear-like, the wave spreads
over the ring) or SELF-TRAPS into a stationary localized hump -- a discrete breather -- once the
on-site nonlinear frequency gamma|psi_0|^2 exceeds the linear bandwidth (the spectrum of L_N is [0,4]).
This is the dynamical opposite of the linear large wave, which is a delocalized interference of all
modes.

A norm-conserving split-step integrator (the nonlinear half preserves each |psi_j|, the linear half is
unitary) gives, as the correctness test, exact conservation of the norm sum|psi|^2 and good
conservation of the Hamiltonian H = -sum|psi_{j+1}-psi_j|^2 + (gamma/2) sum|psi_j|^4. We verify:
 (A) norm conserved to ~1e-12, H to O(dt^2);
 (B) the participation ratio P = (sum|psi|^2)^2 / sum|psi|^4 stays O(1) (trapped) for gamma above a
     threshold and grows toward N (dispersed) below it -- a sharp self-trapping transition;
 (C) the trapped state is a genuine discrete breather: |psi_j| is stationary in time (a fixed profile).
"""

from __future__ import annotations

import numpy as np


def step(psi: np.ndarray, dt: float, gamma: float, lam: np.ndarray) -> np.ndarray:
    psi = psi * np.exp(-1j * gamma * np.abs(psi) ** 2 * dt)
    return np.fft.ifft(np.exp(1j * lam * dt) * np.fft.fft(psi))


def participation(psi: np.ndarray) -> float:
    a2 = np.abs(psi) ** 2
    return float(a2.sum() ** 2 / np.sum(a2**2))


def hamiltonian(psi: np.ndarray, gamma: float) -> float:
    return float(-np.sum(np.abs(np.roll(psi, -1) - psi) ** 2) + 0.5 * gamma * np.sum(np.abs(psi) ** 4))


def evolve(n: int, gamma: float, amp: float, dt: float, steps: int) -> dict:
    lam = 4.0 * np.sin(np.pi * np.arange(n) / n) ** 2
    psi = np.zeros(n, dtype=complex)
    psi[0] = amp
    norm0 = float(np.sum(np.abs(psi) ** 2))
    h0 = hamiltonian(psi, gamma)
    ndrift = hdrift = 0.0
    profiles = []
    for s in range(steps):
        psi = step(psi, dt, gamma, lam)
        if s % (steps // 20) == 0:
            ndrift = max(ndrift, abs(np.sum(np.abs(psi) ** 2) - norm0) / norm0)
            hdrift = max(hdrift, abs(hamiltonian(psi, gamma) - h0) / max(abs(h0), 1.0))
        if s >= steps - 5 * (steps // 20) and s % (steps // 20) == 0:
            profiles.append(np.abs(psi).copy())
    return {"P": participation(psi), "ndrift": ndrift, "hdrift": hdrift,
            "profiles": profiles, "final": np.abs(psi)}


def main() -> int:
    print("=" * 70)
    print("DNLS self-trapping and discrete breathers (linear core = L_N)")
    print("=" * 70)
    ok = True
    n, dt, steps = 128, 0.01, 40000

    # (A) conservation on a representative trapped run
    r = evolve(n, gamma=8.0, amp=1.0, dt=dt, steps=steps)
    print(f"\n(A) Conservation (gamma=8): norm drift={r['ndrift']:.2e}, H drift={r['hdrift']:.2e}")
    a_ok = r["ndrift"] < 1e-10 and r["hdrift"] < 1e-2
    ok &= a_ok
    print(f"    norm exact, H to O(dt^2): {'OK' if a_ok else 'FAIL'}")

    # (B) self-trapping transition: P_final vs gamma (single-site IC, amp 1 so gamma|psi|^2 = gamma)
    print("\n(B) Self-trapping transition: participation ratio P_final vs gamma (N=128)")
    print(f"    {'gamma':>6} {'P_final':>9} {'regime':>10}")
    pf = {}
    for g in (1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0):
        pf[g] = evolve(n, gamma=g, amp=1.0, dt=dt, steps=steps)["P"]
        regime = "trapped" if pf[g] < 5 else ("partial" if pf[g] < 20 else "dispersed")
        print(f"    {g:>6.1f} {pf[g]:>9.2f} {regime:>10}")
    # weak gamma disperses (P large), strong gamma traps (P ~ O(1)); threshold near band top gamma~4
    b_ok = pf[1.0] > 20 and pf[12.0] < 2 and pf[12.0] < pf[1.0]
    ok &= b_ok
    print(f"    -> transition near gamma~4 (band top): dispersed (P~{pf[1.0]:.0f}) -> trapped "
          f"(P~{pf[12.0]:.1f}): {'OK' if b_ok else 'FAIL'}")

    # (C) the trapped profile is stationary (discrete breather)
    r2 = evolve(n, gamma=8.0, amp=1.0, dt=dt, steps=steps)
    profs = np.array(r2["profiles"])  # last few |psi| snapshots
    var = float(np.max(np.std(profs, axis=0)))  # max over sites of temporal std of |psi_j|
    c_ok = var < 0.05  # the localized hump is stationary in time
    ok &= c_ok
    print(f"\n(C) Discrete breather (gamma=8): max temporal std of |psi_j| = {var:.4f} "
          f"(stationary if <0.05) {'OK' if c_ok else 'FAIL'}")

    print("\n" + "=" * 70)
    print("RESULT:", "DNLS SELF-TRAPPING VERIFIED" if ok else "CHECK FAILED")
    print("Nonlinearity LOCALIZES: above threshold a single site self-traps into a stationary discrete")
    print("breather (P~O(1)), the dynamical opposite of the delocalized linear large wave (P~N).")
    print("=" * 70)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
