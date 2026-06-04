"""Bead chain: variational ceiling (Theorem 1), Bohr ceiling U_N, Green's-function invariants."""

from __future__ import annotations

import math

import numpy as np
import pytest

from large_wave_effect import (
    ceiling_U,
    diag_lpinv,
    green_velocity,
    prefix_independent_length,
    variational_max,
)


@pytest.mark.parametrize("n", [8, 16, 31, 64, 257, 1024])
def test_diag_lpinv_exact_identity(n: int) -> None:
    # (L_N^+)_jj = (N^2 - 1) / (12 N), the core of Theorem 1.
    assert diag_lpinv(n) == pytest.approx((n * n - 1) / (12 * n), rel=1e-12)


@pytest.mark.parametrize("n", [8, 16, 64, 256])
def test_variational_max_formula(n: int) -> None:
    # M_N = sqrt(E0 (N^2-1)/6N) at E0 = 1.
    assert variational_max(n, 1.0) == pytest.approx(math.sqrt((n * n - 1) / (6 * n)), rel=1e-12)


@pytest.mark.parametrize("n", [256, 1024, 4096, 16384])
def test_ceiling_asymptotic_constant(n: int) -> None:
    # U_N - (1/pi) ln N -> (gamma + ln(2/pi)) / pi ~ 0.03999.
    c = (np.euler_gamma + math.log(2 / math.pi)) / math.pi
    assert ceiling_U(n) - math.log(n) / math.pi == pytest.approx(c, abs=2e-3)


@pytest.mark.parametrize("n", [16, 32, 64])
def test_green_velocity_real_and_zero_mean(n: int) -> None:
    u = green_velocity(n, 3.7)
    assert u.shape == (n,)
    assert float(np.sum(u)) == pytest.approx(0.0, abs=1e-10)  # zero-mean impulse


@pytest.mark.parametrize(("n", "expect_full"), [(7, True), (16, True), (6, False), (12, False)])
def test_prefix_length(n: int, expect_full: bool) -> None:
    m = prefix_independent_length(n)
    assert m >= 1
    assert (m == n // 2) is expect_full  # full prefix iff prime/2^m here
