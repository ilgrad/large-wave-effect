# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2"]
# ///
"""Stability of the ring-DNLS breathers: site-centred stable, bond-centred unstable.

For the standing wave z = e^{-i Omega t} phi, (L_N + Omega) phi = gamma phi^3, the linearization
z = e^{-i Omega t}(phi + u + i v) obeys  u_dot = -L_- v,  v_dot = L_+ u  with
    L_+ = L_N + Omega - 3 gamma phi^2 ,   L_- = L_N + Omega - gamma phi^2 .
Spectral stability <=> the operator JL = [[0, -L_-], [L_+, 0]] has purely imaginary spectrum.

Two checks, matching the classical DNLS picture (and Vakhitov-Kolokolov / Grillakis-Shatah-Strauss):
  * SITE-centred branch: P(Omega) = ||phi||^2 has dP/dOmega > 0 (P ~ Omega/gamma + O(1)), the VK slope
    condition, and max |Re lambda| ~ 0 -> spectrally & orbitally STABLE.
  * BOND-centred branch: a real eigenvalue pair, max Re lambda >> 0 -> UNSTABLE.
"""

from __future__ import annotations

import numpy as np

N, GAMMA = 64, 1.0


def lap() -> np.ndarray:
    m = np.zeros((N, N))
    for j in range(N):
        m[j, j] = 2.0
        m[j, (j + 1) % N] -= 1.0
        m[j, (j - 1) % N] -= 1.0
    return m


LN = lap()


def breather(omega: float, seed: np.ndarray) -> np.ndarray:
    phi = seed.copy()
    moi = LN + omega * np.eye(N)
    for _ in range(200):
        f = moi @ phi - GAMMA * phi**3
        jac = moi - 3 * GAMMA * np.diag(phi**2)
        d = np.linalg.solve(jac, f)
        phi = phi - d
        if np.linalg.norm(d) < 1e-13:
            break
    return phi


def max_re_lambda(phi: np.ndarray, omega: float) -> float:
    moi = LN + omega * np.eye(N)
    lp = moi - 3 * GAMMA * np.diag(phi**2)
    lm = moi - GAMMA * np.diag(phi**2)
    jl = np.block([[np.zeros((N, N)), -lm], [lp, np.zeros((N, N))]])
    return float(np.max(np.abs(np.linalg.eigvals(jl).real)))


def main() -> int:
    onsite = np.zeros(N)
    onsite[0] = 1.0
    bond = np.zeros(N)
    bond[0] = bond[1] = 1.0

    oms = np.array([2.0, 4.0, 6.0, 8.0, 12.0, 16.0, 20.0])
    powers = np.array([breather(o, onsite * np.sqrt(o + 2)) @ breather(o, onsite * np.sqrt(o + 2)) for o in oms])
    dP = np.gradient(powers, oms)

    print("SITE-centred breather (VK slope dP/dOmega and linearization spectrum):")
    print(f"  {'Omega':>6} {'P':>9} {'dP/dOmega':>10} {'max|Re lambda|':>14}")
    site_ok = True
    for i, o in enumerate(oms):
        phi = breather(o, onsite * np.sqrt(o + 2))
        mr = max_re_lambda(phi, o)
        if not (dP[i] > 0 and mr < 1e-4):
            site_ok = False
        print(f"  {o:>6.1f} {powers[i]:>9.4f} {dP[i]:>10.4f} {mr:>14.2e}")

    print("BOND-centred breather (expect a real unstable eigenvalue):")
    print(f"  {'Omega':>6} {'P':>9} {'max|Re lambda|':>14}")
    bond_unstable = True
    for o in [4.0, 8.0, 12.0, 20.0]:
        phi = breather(o, bond * np.sqrt((o + 2) / 2))
        mr = max_re_lambda(phi, o)
        if mr < 0.1:
            bond_unstable = False
        print(f"  {o:>6.1f} {phi @ phi:>9.4f} {mr:>14.3f}")

    ok = site_ok and bond_unstable
    print("=" * 56)
    print("RESULT:", "PASS -- site-centred STABLE (VK), bond-centred UNSTABLE" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
