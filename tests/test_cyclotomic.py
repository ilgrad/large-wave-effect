"""Exact Q-rank of the ring frequencies: independence theorems (prime, 2^m) and composite defect."""

from __future__ import annotations

import pytest

from large_wave_effect import ring_freq_rank


@pytest.mark.parametrize(("p", "expected"), [(5, 2), (7, 3), (11, 5), (13, 6), (17, 8), (23, 11)])
def test_prime_frequencies_independent(p: int, expected: int) -> None:
    # Theorem 1': {sin(pi r/p)} are Q-independent, so rank = (p-1)/2 (full half).
    assert ring_freq_rank(p) == expected == (p - 1) // 2


@pytest.mark.parametrize(("n", "expected"), [(4, 2), (8, 4), (16, 8), (32, 16)])
def test_power_of_two_independent(n: int, expected: int) -> None:
    # Theorem 1'': rank = N/2 for N = 2^m.
    assert ring_freq_rank(n) == expected == n // 2


@pytest.mark.parametrize(("n", "rank", "half"), [(6, 2, 3), (12, 4, 6), (15, 4, 7), (9, 3, 4)])
def test_composite_is_deficient(n: int, rank: int, half: int) -> None:
    # Composite N: rational relations among sines -> rank strictly below floor(N/2).
    r = ring_freq_rank(n)
    assert r == rank
    assert r < half
