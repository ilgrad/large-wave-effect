# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2", "scipy>=1.11"]
# ///
"""Bessel propagation kernel for the discrete Schrodinger equation, and Poisson summation.

On the infinite chain Z, i dz/dt = L z with z(0)=delta_0 has the EXACT kernel
    K_inf(n,t) = (1/2pi) int e^{i n theta} e^{-i t (2 - 2 cos theta)} dtheta = e^{-2 i t} i^n J_n(2t),
so |K_inf(n,t)| = |J_n(2t)|: ballistic light cone |n| ~ 2t, Airy front max_n |J_n(2t)| ~ 0.6748 (2t)^{-1/3}.

On the ring the kernel is the PERIODIZATION (Poisson summation):
    K_N(j,t) = sum_{l in Z} K_inf(j + l N, t).

Checks (scipy.special.jv for J_n):
 (A) the finite-N ring kernel matches e^{-2it} i^n J_n(2t) while the front stays inside the ring (2t < N/2);
 (B) Poisson summation reproduces the ring kernel once several images l are included;
 (C) the Airy front constant max_n |J_n(2t)| (2t)^{1/3} -> 0.6748.
"""

from __future__ import annotations

import numpy as np
from scipy.special import jv


def ring_kernel(n: int, t: float) -> np.ndarray:
    lam = 4.0 * np.sin(np.pi * np.arange(n) / n) ** 2
    return np.fft.ifft(np.exp(-1j * t * lam))  # K_N(k,t), k=0..N-1


def k_inf(n: np.ndarray, t: float) -> np.ndarray:
    return np.exp(-2j * t) * (1j ** n) * jv(n, 2 * t)


def main() -> int:
    print("=" * 70)
    print("Bessel kernel for discrete Schrodinger + Poisson summation")
    print("=" * 70)
    ok = True

    # (A) ring kernel == e^{-2it} i^n J_n(2t) while the front 2t is well inside the ring
    big, t = 2048, 180.0  # 2t = 360 < N/2 = 1024
    kr = ring_kernel(big, t)
    sites = np.arange(0, 60)
    err = float(np.max(np.abs(kr[sites] - k_inf(sites, t))))
    print(f"\n(A) infinite-chain identity K(n,t)=e^(-2it) i^n J_n(2t),  N={big}, t={t}")
    print(f"    max |K_N(n,t) - e^(-2it) i^n J_n(2t)|, n=0..59 : {err:.2e}")
    print(f"    -> {'PASS' if err < 1e-9 else 'FAIL'}")
    ok &= err < 1e-9

    # (B) Poisson summation: ring = periodization of the infinite kernel
    n_small, t2 = 16, 22.0
    kr2 = ring_kernel(n_small, t2)
    images = np.arange(-60, 61)
    poisson = np.array([np.sum(k_inf(j + images * n_small, t2)) for j in range(n_small)])
    perr = float(np.max(np.abs(kr2 - poisson)))
    print(f"\n(B) Poisson sum K_N(j,t)=sum_l K_inf(j+lN,t),  N={n_small}, t={t2}")
    print(f"    max |ring - periodized| over j : {perr:.2e}")
    print(f"    -> {'PASS' if perr < 1e-9 else 'FAIL'}")
    ok &= perr < 1e-9

    # (C) Airy front constant
    print("\n(C) Airy front:  max_n |J_n(2t)| (2t)^{1/3} -> 0.6748 (Airy)")
    print(f"    {'2t':>7} {'max|J|':>10} {'x^(1/3) max':>12}")
    const = None
    for z in (50.0, 200.0, 800.0, 3200.0):
        nv = np.arange(0, int(z * 1.3) + 5)
        mx = float(np.max(np.abs(jv(nv, z))))
        const = z ** (1 / 3) * mx
        print(f"    {z:>7.0f} {mx:>10.5f} {const:>12.5f}")
    airy_ok = abs(const - 0.6748) < 5e-3
    print(f"    -> matches 0.6748: {'PASS' if airy_ok else 'FAIL'}")
    ok &= airy_ok

    print("\n" + "=" * 70)
    print("RESULT:", "BESSEL KERNEL VERIFIED" if ok else "CHECK FAILED")
    print("The ring Schrodinger kernel is the periodized Bessel kernel; the ballistic front and its")
    print("(2t)^{-1/3} decay underlie B_N ~ sqrt(N) and the vanishing of the effect on Z.")
    print("=" * 70)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
