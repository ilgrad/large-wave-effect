# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2"]
# ///
"""Modulational instability of the ring DNLS, tied to the L_N spectrum (Proposition dnls-mi).

For the ring DNLS  i z_dot = -L_N z + gamma |z|^2 z  (paper convention), the uniform plane wave
z_j(t) = A e^{-i gamma A^2 t} is an exact solution.  Linearizing z = (A + b) e^{-i gamma A^2 t} with
b = x + i y gives the closed real system

    x_dot = -L_N y ,        y_dot = (L_N - 2 gamma A^2) x ,

so in each Fourier mode Q (eigenvalue lam_Q = 2 - 2 cos(2 pi Q / N) of L_N) the 2x2 block has eigenvalues
mu with  mu^2 = lam_Q (2 gamma A^2 - lam_Q).  Hence for gamma > 0 the mode grows like e^{sigma_Q t} with

    sigma_Q = sqrt( lam_Q (2 gamma A^2 - lam_Q) )   on the band 0 < lam_Q < 2 gamma A^2,

and is neutrally stable otherwise; the fastest growth is sigma_max = gamma A^2 at lam_Q = gamma A^2.
This script builds the exact 2N x 2N linearized operator and checks its spectrum equals these sigma_Q.
"""

from __future__ import annotations

import numpy as np
from numpy import cos, pi, sqrt


def main() -> int:
    n, gamma, a = 128, 1.0, 1.0
    lam = 2 - 2 * cos(2 * pi * np.arange(n) / n)
    lmat = np.zeros((n, n))
    for j in range(n):
        lmat[j, j] = 2.0
        lmat[j, (j + 1) % n] -= 1.0
        lmat[j, (j - 1) % n] -= 1.0
    # [x_dot; y_dot] = [[0, -L],[L - 2 g A^2 I, 0]] [x; y]
    op = np.block([[np.zeros((n, n)), -lmat], [lmat - 2 * gamma * a**2 * np.eye(n), np.zeros((n, n))]])
    ev = np.linalg.eigvals(op)
    max_growth = float(np.max(ev.real))

    sigma = np.array([sqrt(lam[q] * (2 * gamma * a**2 - lam[q])) if 0 < lam[q] < 2 * gamma * a**2 else 0.0
                      for q in range(n)])
    band = [q for q in range(n) if 0 < lam[q] < 2 * gamma * a**2]

    print(f"DNLS modulational instability (N={n}, gamma={gamma}, A={a}):")
    print(f"  fastest growth: operator = {max_growth:.5f},  discrete max sigma_Q = {sigma.max():.5f}"
          f"  (continuum bound gamma A^2 = {gamma * a**2:.2f})")
    print(f"  unstable band |{{Q: lam_Q < 2 gamma A^2}}| = {len(band)} modes")
    print(f"  {'Q':>3} {'lam_Q':>7} {'sigma_Q':>9} {'in operator spectrum?':>22}")
    ok = abs(max_growth - sigma.max()) < 1e-6  # operator spectrum matches the formula exactly
    for q in [1, 3, 8, 12, 20, 30]:
        present = bool(np.any(np.abs(ev.real - sigma[q]) < 1e-6)) if sigma[q] > 0 else True
        ok &= present
        print(f"  {q:>3} {lam[q]:>7.4f} {sigma[q]:>9.5f} {present!s:>22}")

    print("=" * 60)
    print("RESULT:", "PASS -- MI growth rate sigma_Q = sqrt(lam_Q(2 g A^2 - lam_Q))" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
