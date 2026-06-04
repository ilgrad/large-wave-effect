"""Discrete Schrodinger on the ring: unitarity ceiling and ballistic ell^inf amplification."""

from __future__ import annotations

import numpy as np
import pytest

from large_wave_effect import amplification_B, kernel, max_local_amplitude


@pytest.mark.parametrize("n", [16, 64, 128])
def test_delta_state_not_amplified(n: int) -> None:
    # Unitarity: a localized state never exceeds amplitude 1 (no large wave from a point).
    assert max_local_amplitude(n) <= 1.0 + 1e-9


@pytest.mark.parametrize("n", [16, 64, 256])
def test_kernel_is_unitary_in_l2(n: int) -> None:
    k = kernel(n, 5.3)
    assert float(np.sum(np.abs(k) ** 2)) == pytest.approx(1.0, rel=1e-10)


def test_ell_inf_amplification_grows_like_sqrt_n() -> None:
    # B_N is increasing and consistent with Theta(sqrt N): B_N / sqrt(N) stays in a tight band.
    ns = [16, 64, 256]
    b = [amplification_B(n) for n in ns]
    assert b[0] < b[1] < b[2]
    ratios = [bi / np.sqrt(ni) for bi, ni in zip(b, ns, strict=True)]
    assert max(ratios) / min(ratios) < 1.2  # sqrt(N) scaling, not ln N or N
