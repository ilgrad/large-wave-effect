"""Discrete Schrodinger dynamics on the ring: i dz/dt = L_N z, z(t) = exp(-i t L_N) z(0).

Phases are governed by lambda_r = 4 sin^2(pi r/N) = 2 - 2 cos(2 pi r/N).
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

ComplexArray = npt.NDArray[np.complex128]


def kernel(n: int, t: float) -> ComplexArray:
    """Propagator kernel K_N(k,t) (= z_k(t) for the delta initial state z(0)=e_0)."""
    lam = 4.0 * np.sin(np.pi * np.arange(n) / n) ** 2
    return np.fft.ifft(np.exp(-1j * t * lam))  # ifft carries 1/N -> equals K_N


def ell_inf_norm(n: int, t: float) -> float:
    """||exp(-i t L_N)||_{inf->inf} = sum_k |K_N(k,t)| (circulant row sum)."""
    return float(np.sum(np.abs(kernel(n, t))))


def amplification_B(n: int, t_max_mult: float = 6.0, samples: int = 6000) -> float:
    """B_N = sup_t ||exp(-i t L_N)||_{inf->inf} (grid lower bound; ballistic ~ sqrt(N))."""
    ts = np.linspace(0.01, t_max_mult * n, samples)
    return max(ell_inf_norm(n, float(t)) for t in ts)


def max_local_amplitude(n: int, t_max_mult: float = 8.0, samples: int = 6000) -> float:
    """max_j sup_t |z_j(t)| for the delta state; <= 1 by unitarity (no large wave from a point)."""
    ts = np.linspace(0.0, t_max_mult * n, samples)
    return max(float(np.max(np.abs(kernel(n, float(t))))) for t in ts)
