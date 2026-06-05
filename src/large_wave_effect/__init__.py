"""Large-wave amplification on the periodic discrete chain: reusable core routines.

Experiments and certificates live in ``numerics/``; this package holds the tested core they build on.
"""

from __future__ import annotations

from .cyclotomic import cyclotomic, exact_rank, power_mod, ring_freq_rank, sin_coords
from .ring import (
    ceiling_U,
    diag_lpinv,
    frequencies,
    green_velocity,
    prefix_independent_length,
    variational_max,
)
from .schrodinger import amplification_B, ell_inf_norm, kernel, max_local_amplitude

__all__ = [
    "amplification_B",
    "ceiling_U",
    "cyclotomic",
    "diag_lpinv",
    "ell_inf_norm",
    "exact_rank",
    "frequencies",
    "green_velocity",
    "kernel",
    "max_local_amplitude",
    "power_mod",
    "prefix_independent_length",
    "ring_freq_rank",
    "sin_coords",
    "variational_max",
]
