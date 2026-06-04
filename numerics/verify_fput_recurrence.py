# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2"]
# ///
"""FPUT-beta recurrence on the ring: energy returns to the seeded mode (the FPUT paradox).

The Fermi-Pasta-Ulam-Tsingou beta lattice on the periodic ring is
    u_ddot_j = (u_{j+1} - 2u_j + u_{j-1}) + beta[(u_{j+1}-u_j)^3 - (u_j-u_{j-1})^3],
with Hamiltonian H = sum p_j^2/2 + sum [ (u_{j+1}-u_j)^2/2 + beta (u_{j+1}-u_j)^4/4 ]. Its linear part
is the same ring Laplacian L_N (omega_r = 2 sin(pi r/N)) as the linear theory; the cubic term couples
the normal modes. Seeding all energy in mode 1, the energy does NOT thermalize on short times: it
flows to a few low modes and returns nearly to mode 1 (the FPUT recurrence).

We verify, with a symplectic velocity-Verlet integrator (energy conservation is the correctness test):
 (A) H is conserved to O(dt^2) and stays flat over the whole run (no secular drift);
 (B) the mode-1 harmonic energy E_1(t) drops and then RECURS to a large fraction of its initial value
     -- a genuine return, not mere persistence;
 (C) on these timescales the dynamics is quasi-periodic among a few modes (the participation in mode
     space stays small), the near-integrable face of the nonlinear lattice.
"""

from __future__ import annotations

import numpy as np


def force(u: np.ndarray, beta: float) -> np.ndarray:
    dr = np.roll(u, -1) - u  # u_{j+1} - u_j
    # F_j = (u_{j+1}-2u_j+u_{j-1}) + beta[(u_{j+1}-u_j)^3 - (u_j-u_{j-1})^3] = (dr - dr_{j-1}) + beta(dr^3 - dr^3_{j-1})
    return (dr - np.roll(dr, 1)) + beta * (dr**3 - np.roll(dr**3, 1))


def hamiltonian(u: np.ndarray, p: np.ndarray, beta: float) -> float:
    dr = np.roll(u, -1) - u
    return float(np.sum(p**2) / 2 + np.sum(dr**2) / 2 + beta * np.sum(dr**4) / 4)


def mode_energies(u: np.ndarray, p: np.ndarray, n: int) -> np.ndarray:
    omega2 = 4.0 * np.sin(np.pi * np.arange(n) / n) ** 2
    uh = np.fft.fft(u) / np.sqrt(n)
    ph = np.fft.fft(p) / np.sqrt(n)
    return 0.5 * (np.abs(ph) ** 2 + omega2 * np.abs(uh) ** 2)


def run(n: int, beta: float, amp: float, dt: float, steps: int) -> dict:
    j = np.arange(n)
    u = amp * np.sin(2 * np.pi * j / n)  # seed mode r=1
    p = np.zeros(n)
    h0 = hamiltonian(u, p, beta)
    e1_0 = mode_energies(u, p, n)[1] + mode_energies(u, p, n)[n - 1]  # r=1 and its degenerate r=N-1
    hmax = 0.0
    e1_series, t_series = [], []
    f = force(u, beta)
    for s in range(steps):
        p += 0.5 * dt * f
        u += dt * p
        f = force(u, beta)
        p += 0.5 * dt * f
        if s % 20 == 0:
            me = mode_energies(u, p, n)
            e1_series.append(float(me[1] + me[n - 1]))
            t_series.append(s * dt)
            hmax = max(hmax, abs(hamiltonian(u, p, beta) - h0) / abs(h0))
    return {
        "h0": h0, "e1_0": e1_0, "hdrift": hmax,
        "t": np.array(t_series), "e1": np.array(e1_series),
    }


def main() -> int:
    print("=" * 70)
    print("FPUT-beta recurrence on the ring (symplectic; linear core = L_N)")
    print("=" * 70)
    ok = True
    n, beta, amp, dt, steps = 32, 0.25, 4.0, 0.05, 600000
    r = run(n, beta, amp, dt, steps)

    print(f"\n(A) Energy conservation: max |H(t)-H_0|/|H_0| = {r['hdrift']:.2e} over T={steps * dt:.0f}")
    a_ok = r["hdrift"] < 1e-3
    ok &= a_ok
    print(f"    symplectic, no secular drift: {'OK' if a_ok else 'FAIL'}")

    e1 = r["e1"] / r["e1_0"]  # mode-1 energy as fraction of initial
    e1_min = float(e1.min())
    # recurrence: after dropping below 0.6, E_1 returns above 0.85
    dropped = np.where(e1 < 0.6)[0]
    recur = 0.0
    if len(dropped):
        after = e1[dropped[0]:]
        recur = float(after.max())
    print(f"\n(B) Mode-1 energy E_1(t)/E_1(0): min={e1_min:.3f}, recurs to {recur:.3f} after dropping")
    b_ok = e1_min < 0.6 and recur > 0.85
    ok &= b_ok
    print(f"    genuine recurrence (drops below 0.6, returns above 0.85): {'OK' if b_ok else 'FAIL'}")

    # (C) only a few low modes ever carry appreciable energy (quasi-periodic, not thermalized)
    j = np.arange(n)
    u = amp * np.sin(2 * np.pi * j / n)
    p = np.zeros(n)
    f = force(u, beta)
    peak_share = np.zeros(n)
    for s in range(200000):
        p += 0.5 * dt * f
        u += dt * p
        f = force(u, beta)
        p += 0.5 * dt * f
        if s % 200 == 0:
            me = mode_energies(u, p, n)
            peak_share = np.maximum(peak_share, me / me.sum())
    active = int(np.sum(peak_share[: n // 2] > 0.02))  # modes ever above 2% of energy
    print(f"\n(C) Active low modes (ever >2% of energy) = {active} of {n // 2} "
          f"(quasi-periodic, not thermalized)")
    c_ok = active <= 8
    ok &= c_ok
    print(f"    energy confined to a few modes: {'OK' if c_ok else 'FAIL'}")

    print("\n" + "=" * 70)
    print("RESULT:", "FPUT RECURRENCE VERIFIED" if ok else "CHECK FAILED")
    print("The nonlinear lattice (linear core L_N) is near-integrable at this energy: energy returns")
    print("to the seeded mode (FPUT recurrence) instead of thermalizing -- the regular regime.")
    print("=" * 70)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
