# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2"]
# ///
"""Statistics of the ring: density of states (arcsine law) and value distribution (Bohr-Jessen).

Density of states. The eigenvalues lambda_r = 4 sin^2(pi r/N) are the pushforward of the uniform law
on r/N by lambda(x)=4 sin^2(pi x). As N -> infinity their empirical distribution converges to the
ARCSINE law on [0,4],
    rho(lambda) = 1 / (pi sqrt(lambda (4 - lambda))),   CDF F(lambda) = (2/pi) arcsin(sqrt(lambda)/2),
the classical density of states of the 1D discrete Laplacian. Likewise the frequencies
omega_r = 2 sin(pi r/N) have density rho(omega) = 2/(pi sqrt(4 - omega^2)) on [0,2].

Value distribution. For prime N the impulse response u_0(t)=(1/N) sum_r sin(omega_r t)/omega_r is an
almost-periodic function with Q-independent frequencies; by the Bohr-Jessen theorem its values have a
limiting distribution. We verify it is symmetric with mean 0 and variance = 1/12 (the exact RMS^2 from
the ring degeneracy), and NON-Gaussian -- platykurtic (kurtosis < 3), because the lowest mode
dominates (weight 1/omega_1 ~ N) and a single sinusoid is arcsine-distributed (kurtosis 1.5). The
large wave is the rare near-alignment supremum (growing like ln N), not a heavy tail of the typical value.

All three checks are quantitative (Kolmogorov-Smirnov distances and moments).
"""

from __future__ import annotations

import numpy as np


def ks_distance(samples: np.ndarray, cdf) -> float:
    xs = np.sort(samples)
    n = len(xs)
    emp = np.arange(1, n + 1) / n
    return float(np.max(np.abs(emp - cdf(xs))))


def main() -> int:
    print("=" * 70)
    print("Statistics: density of states (arcsine) and value distribution (Bohr-Jessen)")
    print("=" * 70)
    ok = True

    # (A) DOS of lambda_r -> arcsine on [0,4]
    print("\n(A) Density of states: lambda_r ~ arcsine,  F(l)=(2/pi) arcsin(sqrt(l)/2)")
    for n in (256, 1024, 4096):
        lam = 4.0 * np.sin(np.pi * np.arange(1, n) / n) ** 2
        ks = ks_distance(lam, lambda x: (2 / np.pi) * np.arcsin(np.sqrt(np.clip(x, 0, 4)) / 2))
        ok &= ks < 2.0 / np.sqrt(n)  # KS -> 0 like 1/sqrt(N)
        print(f"  N={n:>5}: KS(lambda, arcsine) = {ks:.4f}  (< 2/sqrt(N)={2 / np.sqrt(n):.4f})  "
              f"{'OK' if ks < 2 / np.sqrt(n) else 'FAIL'}")

    # (B) frequency density omega_r -> 2/(pi sqrt(4-omega^2)) on [0,2]
    print("\n(B) Frequency density: omega_r,  F(w)=(2/pi) arcsin(w/2)")
    for n in (1024, 4096):
        om = 2.0 * np.sin(np.pi * np.arange(1, n) / n)
        ks = ks_distance(om, lambda w: (2 / np.pi) * np.arcsin(np.clip(w, 0, 2) / 2))
        ok &= ks < 2.0 / np.sqrt(n)
        print(f"  N={n:>5}: KS(omega, law) = {ks:.4f}  {'OK' if ks < 2 / np.sqrt(n) else 'FAIL'}")

    # (C) value distribution of u_0(t) (Bohr-Jessen): mean 0, var 1/12, symmetric, heavy tails
    print("\n(C) Value distribution of u_0(t) (prime N): mean 0, var 1/12, symmetric, platykurtic")
    n = 61
    om = 2.0 * np.sin(np.pi * np.arange(1, n) / n)
    ts = np.linspace(0.0, 200000.0, 400000)
    b = 1.0 / n / om  # u_0(t) = (1/N) sum_r sin(omega_r t)/omega_r
    vals = (np.sin(np.outer(ts, om)) * b).sum(axis=1)
    mean, var = float(vals.mean()), float(vals.var())
    skew = float(((vals - mean) ** 3).mean() / var**1.5)
    kurt = float(((vals - mean) ** 4).mean() / var**2)  # Gaussian = 3, arcsine = 1.5
    cok = abs(mean) < 5e-3 and abs(var - 1 / 12) < 5e-3 and abs(skew) < 0.05 and 1.5 <= kurt < 3.0
    ok &= cok
    print(f"  N={n}: mean={mean:+.4f} (0)  var={var:.4f} (1/12={1 / 12:.4f})  "
          f"skew={skew:+.3f} (0)  kurtosis={kurt:.2f} (Gauss 3, arcsine 1.5)")
    print(f"  -> symmetric, var=1/12, platykurtic (lowest-mode dominated): {'OK' if cok else 'FAIL'}")

    print("\n" + "=" * 70)
    print("RESULT:", "STATISTICS VERIFIED" if ok else "CHECK FAILED")
    print("Spectrum follows the arcsine density of states; the large-wave amplitude has a Bohr-Jessen")
    print("value law (mean 0, var 1/12, platykurtic, lowest-mode dominated). The large wave is the rare")
    print("near-alignment supremum (~ ln N), not a heavy tail of the typical value.")
    print("=" * 70)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
