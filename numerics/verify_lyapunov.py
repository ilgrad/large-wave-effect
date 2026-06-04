# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2"]
# ///
"""Route to chaos in the FPUT-beta ring: the largest Lyapunov exponent turns on with energy.

The linear large wave lives on integrable tori (zero Lyapunov exponent). Adding the cubic coupling, the
FPUT-beta lattice u_ddot_j = (u_{j+1}-2u_j+u_{j-1}) + beta[(u_{j+1}-u_j)^3 - (u_j-u_{j-1})^3] keeps its
linear core L_N but acquires genuine chaos ABOVE an energy threshold. We measure the maximal Lyapunov
exponent lambda_max by the standard Benettin method: a tangent vector (delta u, delta p) is co-evolved
with the variational equations and periodically renormalized,
    lambda_max = (1 / T) sum log( ||tangent|| / d0 ).

We verify:
 (A) at low energy lambda_max ~ 0 (regular / quasi-periodic, the FPUT-recurrence regime);
 (B) lambda_max grows with energy and is clearly positive at high energy (chaos);
 (C) the underlying trajectory conserves H (symplectic), so the positive exponent is dynamical, not a
     numerical artefact.
"""

from __future__ import annotations

import numpy as np

RNG = np.random.default_rng(0)


def force(u: np.ndarray, beta: float) -> np.ndarray:
    dr = np.roll(u, -1) - u
    return (dr - np.roll(dr, 1)) + beta * (dr**3 - np.roll(dr**3, 1))


def varforce(u: np.ndarray, du: np.ndarray, beta: float) -> np.ndarray:
    dr = np.roll(u, -1) - u
    ddr = np.roll(du, -1) - du
    g = dr**2 * ddr
    return (ddr - np.roll(ddr, 1)) + 3.0 * beta * (g - np.roll(g, 1))


def hamiltonian(u: np.ndarray, p: np.ndarray, beta: float) -> float:
    dr = np.roll(u, -1) - u
    return float(np.sum(p**2) / 2 + np.sum(dr**2) / 2 + beta * np.sum(dr**4) / 4)


def lyapunov(n: int, beta: float, amp: float, dt: float, t_steps: int) -> tuple[float, float, float]:
    j = np.arange(n)
    u = amp * np.sin(2 * np.pi * j / n)
    p = np.zeros(n)
    h0 = hamiltonian(u, p, beta)
    d0 = 1e-7
    du = RNG.standard_normal(n)
    dp = RNG.standard_normal(n)
    nrm = np.sqrt(np.sum(du**2) + np.sum(dp**2))
    du *= d0 / nrm
    dp *= d0 / nrm
    f, vf = force(u, beta), varforce(u, du, beta)
    renorm_every = max(1, round(1.0 / dt))
    lsum, ncount, hmax = 0.0, 0, 0.0
    for s in range(1, t_steps + 1):
        p += 0.5 * dt * f
        dp += 0.5 * dt * vf
        u += dt * p
        du += dt * dp
        f, vf = force(u, beta), varforce(u, du, beta)
        p += 0.5 * dt * f
        dp += 0.5 * dt * vf
        if s % renorm_every == 0:
            d = np.sqrt(np.sum(du**2) + np.sum(dp**2))
            lsum += np.log(d / d0)
            du *= d0 / d
            dp *= d0 / d
            ncount += 1
            hmax = max(hmax, abs(hamiltonian(u, p, beta) - h0) / abs(h0))
    lam = lsum / (ncount * renorm_every * dt)
    return lam, h0, hmax


def main() -> int:
    print("=" * 70)
    print("Lyapunov exponent and the route to chaos in the FPUT-beta ring")
    print("=" * 70)
    n, beta, dt, t_steps = 32, 0.25, 0.02, 250000
    print(f"\n(A,B) lambda_max vs energy (N={n}, beta={beta}, T={t_steps * dt:.0f})")
    print(f"    {'amp':>5} {'energy H':>10} {'lambda_max':>11} {'Hdrift':>9} {'regime':>10}")
    rows = []
    for amp in (1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0):
        lam, h0, hdrift = lyapunov(n, beta, amp, dt, t_steps)
        regime = "chaotic" if lam > 0.02 else ("weak" if lam > 0.005 else "regular")
        rows.append((amp, h0, lam, hdrift))
        print(f"    {amp:>5.1f} {h0:>10.2f} {lam:>11.4f} {hdrift:>9.1e} {regime:>10}")

    lam_lo = rows[0][2]   # amp=1 (low energy)
    lam_hi = rows[-1][2]  # amp=10 (high energy)
    hdrift_max = max(r[3] for r in rows)
    a_ok = lam_lo < 0.01
    b_ok = lam_hi > 0.02 and lam_hi > 8 * max(lam_lo, 1e-4)
    c_ok = hdrift_max < 1e-2
    ok = a_ok and b_ok and c_ok
    print(f"\n  (A) low energy lambda_max={lam_lo:.4f} ~ 0 (regular): {'OK' if a_ok else 'FAIL'}")
    print(f"  (B) high energy lambda_max={lam_hi:.4f} > 0 (chaos), grows with E: {'OK' if b_ok else 'FAIL'}")
    print(f"  (C) trajectory H conserved (max drift {hdrift_max:.1e}): {'OK' if c_ok else 'FAIL'}")

    print("\n" + "=" * 70)
    print("RESULT:", "LYAPUNOV TRANSITION VERIFIED" if ok else "CHECK FAILED")
    print("The same lattice (linear core L_N) is integrable-like at low energy (lambda~0, FPUT")
    print("recurrence) and genuinely chaotic at high energy (lambda>0) -- the regular-to-chaotic route.")
    print("=" * 70)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
