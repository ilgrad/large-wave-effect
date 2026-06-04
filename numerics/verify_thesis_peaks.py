# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2"]
# ///
"""Validation against the 2015 diploma's numerical peaks (independent reproduction).

The diploma solved the free wave equation on the ring with the step displacement IC
alpha_j = 1 for 1 <= j <= N-1, alpha_0 = 0, beta = 0, and reported max_{j,t} y_j(t):

    N    32     64    128    256    512   1024
    peak 1.617  1.562 1.505  1.445  1.397 1.397

Here we recompute the peak by the EXACT Bohr supremum (no time scan). For N = 2^m the
frequencies omega_r = 2 sin(pi r/N) are Q-independent (Theorem 1''), so for u(t) = cos(t sqrt L) alpha,

    max_t u_j(t) = abar + sum_r |amplitude_{r,j}|,   abar = (1/N) sum_j alpha_j,

and the large wave is the maximum over nodes j. The exact supremum turns out to be ~2 for all N.

FINDING: this EXCEEDS the diploma's 1.617 -> 1.397, and explains them: those values are finite-time
scan lower bounds (the diploma scanned up to t ~ 2e8). The exact sup lives at near-recurrence times
that the scan does not reach, and the undershoot worsens with N -- exactly the diploma's decreasing
trend. So the present analysis reproduces the EFFECT (large wave, peak > 1) and corrects the value
(true sup ~ 2), without any time scan.
"""

from __future__ import annotations

import numpy as np

THESIS = {32: 1.617, 64: 1.562, 128: 1.505, 256: 1.445, 512: 1.397, 1024: 1.397}


def bohr_peak_step(n: int) -> float:
    """Exact sup_{t} max_j cos(t sqrt L_N) applied to the step displacement IC (N = 2^m)."""
    alpha = np.ones(n)
    alpha[0] = 0.0
    h = np.fft.fft(alpha)  # modal coefficients
    abar = h[0].real / n
    j = np.arange(n)
    # complex per-mode contribution g[r,j] = (1/N) h_r e^{2 pi i r j / N}
    g = (h[:, None] / n) * np.exp(2j * np.pi * np.outer(np.arange(n), j) / n)
    amp = np.zeros((n, n))
    amp[1 : n // 2] = 2.0 * np.real(g[1 : n // 2])  # paired modes r and N-r
    amp[n // 2] = np.real(g[n // 2])  # middle mode (N even)
    sup_per_node = abar + np.sum(np.abs(amp[1 : n // 2 + 1]), axis=0)
    return float(np.max(sup_per_node))


def main() -> int:
    print("=" * 66)
    print("Diploma 2015 step-IC peak: exact Bohr sup vs the reported scan value")
    print("=" * 66)
    print(f"  {'N':>5} {'exact sup':>10} {'diploma scan':>13} {'sup>scan?':>10}")
    exceeds = True
    near_two = True
    for n, ref in THESIS.items():
        p = bohr_peak_step(n)
        exceeds &= p > ref
        near_two &= 1.9 <= p <= 2.001
        print(f"  {n:>5} {p:>10.4f} {ref:>13.3f} {('yes' if p > ref else 'NO'):>10}")
    ok = exceeds and near_two
    print("=" * 66)
    print("RESULT:", "CONSISTENT (sup ~2 exceeds the finite-scan diploma values)" if ok else "UNEXPECTED")
    print("Interpretation: the true large wave is ~2 for all N (exact, no scan). The diploma's")
    print("1.617 -> 1.397 are time-scan lower bounds whose decay with N is the expected worsening")
    print("of scan undershoot. The effect (peak > 1) is reproduced; the magnitude is corrected.")
    print("=" * 66)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
