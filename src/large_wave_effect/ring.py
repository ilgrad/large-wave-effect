"""Periodic bead chain (ring): spectrum, velocity Green's function, ceilings.

Normalization m = kappa = 1, so the ring Laplacian L_N = circ(2,-1,0,...,0,-1) has eigenvalues
lambda_r = 4 sin^2(pi r/N) and frequencies omega_r = 2 sin(pi r/N). Operators f(L_N) are applied
in the Fourier basis via the FFT (O(N log N)).
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

FloatArray = npt.NDArray[np.float64]


def frequencies(n: int) -> tuple[FloatArray, FloatArray]:
    """Return (lambda_r, omega_r) for r = 0..N-1."""
    theta = np.pi * np.arange(n) / n
    lam = 4.0 * np.sin(theta) ** 2
    return lam, 2.0 * np.abs(np.sin(theta))


def diag_lpinv(n: int) -> float:
    """(L_N^+)_jj = (1/N) sum_{r!=0} 1/lambda_r = (N^2-1)/(12 N), independent of j."""
    r = np.arange(1, n)
    return float(np.sum(1.0 / (4.0 * np.sin(np.pi * r / n) ** 2)) / n)


def variational_max(n: int, e0: float = 1.0) -> float:
    """Theorem 1: max deflection M_N = sqrt(2 E0 (L_N^+)_jj) = sqrt(E0 (N^2-1)/6N)."""
    return float(np.sqrt(2.0 * e0 * diag_lpinv(n)))


def ceiling_U(n: int) -> float:
    """Bohr ceiling U_N = (1/2N) sum_{r=1}^{N-1} csc(pi r/N) ~ (1/pi) ln N + (gamma+ln(2/pi))/pi."""
    r = np.arange(1, n)
    return float(np.sum(1.0 / np.sin(np.pi * r / n)) / (2 * n))


def green_velocity(n: int, t: float) -> FloatArray:
    """Velocity Green's function u_j(t) = G_N(j,t) for the zero-mean impulse v0 = e_0 - 1/N."""
    _, omega = frequencies(n)
    nz = omega > 1e-12
    coeff = np.zeros(n)
    coeff[nz] = np.sin(omega[nz] * t) / omega[nz]
    return np.real(np.fft.ifft(coeff))  # ifft carries 1/N -> equals G_N


def prefix_independent_length(n: int) -> int:
    """M(N): largest M with {sin(pi r/N): r=1..M} Q-independent (modular incremental rank)."""
    from .cyclotomic import cyclotomic, power_mod

    q = (1 << 61) - 1
    phi = cyclotomic(2 * n)
    basis: list[list[int]] = []
    pivots: list[int] = []
    m = 0
    for r in range(1, n // 2 + 1):
        hi = power_mod(r % (2 * n), phi)
        lo = power_mod((2 * n - r) % (2 * n), phi)
        row = [(h - lo_) % q for h, lo_ in zip(hi, lo, strict=True)]
        for b, pc in zip(basis, pivots, strict=True):
            f = row[pc]
            if f:
                row = [(row[k] - f * b[k]) % q for k in range(len(row))]
        pc = next((k for k in range(len(row)) if row[k]), None)
        if pc is None:
            break
        inv = pow(row[pc], q - 2, q)
        basis.append([(x * inv) % q for x in row])
        pivots.append(pc)
        m = r
    return m
