# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2", "scipy>=1.11"]
# ///
"""Spectral zeta function of the ring Laplacian: a unifying thread for the large wave.

Define the spectral zeta of L_N (zero mode removed):
    zeta_{L_N}(s) = sum_{r=1}^{N-1} lambda_r^{-s},   lambda_r = 4 sin^2(pi r/N).

The large-wave constants are special values:
    * Bohr ceiling      U_N           = zeta_{L_N}(1/2) / N         (large wave A_N ~ U_N ~ (1/pi) ln N)
    * variational M_N   M_N^2         = 2 E0 zeta_{L_N}(1) / N      (Theorem 1; (L_N^+)_jj = zeta(1)/N)

Closed forms at integer s come from cosecant power sums (Bernoulli numbers):
    zeta(1) = (N^2-1)/12,   zeta(2) = (N^2-1)(N^2+11)/720.

The N -> infinity density is the symbol integral (Szego):
    (1/N) zeta_{L_N}(s) -> I(s) = (1/2pi) int_0^{2pi} (2-2cos th)^{-s} dth
                                = 4^{-s} Gamma(1/2 - s) / (sqrt(pi) Gamma(1 - s)),   Re s < 1/2,
which diverges as s -> 1/2 (the logarithmic edge that produces U_N ~ (1/pi) ln N) and is invalid at
s=1 (where zeta(1) ~ N^2/12 gives the sqrt(N) variational ceiling). Everything below is verified.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import quad
from scipy.special import gamma


def zeta_LN(n: int, s: float) -> float:
    r = np.arange(1, n)
    lam = 4.0 * np.sin(np.pi * r / n) ** 2
    return float(np.sum(lam ** (-s)))


def U_N(n: int) -> float:
    r = np.arange(1, n)
    return float(np.sum(1.0 / np.sin(np.pi * r / n)) / (2 * n))


def I_gamma(s: float) -> float:
    return float(4.0 ** (-s) * gamma(0.5 - s) / (np.sqrt(np.pi) * gamma(1.0 - s)))


def I_quad(s: float) -> float:
    val, _ = quad(lambda th: (2 - 2 * np.cos(th)) ** (-s), 0, np.pi, points=[0], limit=200)
    return val / np.pi  # (1/2pi) int_0^{2pi} = (1/pi) int_0^pi by symmetry


def main() -> int:
    print("=" * 74)
    print("Spectral zeta function of the ring Laplacian")
    print("=" * 74)
    ok = True

    # (A) exact integer values
    print("\n(A) Closed forms at integer s")
    for n in (16, 64, 256, 1024):
        z1, z2 = zeta_LN(n, 1), zeta_LN(n, 2)
        e1, e2 = (n * n - 1) / 12, (n * n - 1) * (n * n + 11) / 720
        a = abs(z1 - e1) / e1 < 1e-10 and abs(z2 - e2) / e2 < 1e-10
        ok &= a
        print(f"  N={n:>5}: zeta(1)={z1:.4f} [exp {e1:.4f}], zeta(2)={z2:.6g} [exp {e2:.6g}]  {'OK' if a else 'FAIL'}")

    # (B) large-wave constants as special values
    print("\n(B) U_N = zeta(1/2)/N  and  M_N^2 = 2 zeta(1)/N  (E0=1)")
    for n in (32, 128, 512):
        lhs_u, rhs_u = U_N(n), zeta_LN(n, 0.5) / n
        m2, z1n = (n * n - 1) / (6 * n), 2 * zeta_LN(n, 1) / n
        bu = abs(lhs_u - rhs_u) < 1e-10
        bm = abs(m2 - z1n) < 1e-9
        ok &= bu and bm
        print(f"  N={n:>4}: U_N={lhs_u:.6f}=zeta(1/2)/N={rhs_u:.6f} {'OK' if bu else 'FAIL'} | "
              f"M_N^2={m2:.5f}=2zeta(1)/N={z1n:.5f} {'OK' if bm else 'FAIL'}")

    # (C) symbol-integral density and the gamma closed form
    print("\n(C) (1/N) zeta_{L_N}(s) -> I(s) = 4^-s Gamma(1/2-s)/(sqrt(pi)Gamma(1-s)), Re s<1/2")
    print(f"  {'s':>5} {'I_gamma':>10} {'I_quad':>10} {'(1/N)zeta, N=2^14':>18}")
    for s in (0.15, 0.25, 0.35):
        ig, iq = I_gamma(s), I_quad(s)
        zn = zeta_LN(16384, s) / 16384
        gq = abs(ig - iq) / abs(ig) < 1e-3
        ok &= gq
        print(f"  {s:>5.2f} {ig:>10.5f} {iq:>10.5f} {zn:>18.5f}   gamma=quad:{'OK' if gq else 'FAIL'}")
    # convergence (1/N)zeta(0.25) -> I(0.25) as N grows
    s0 = 0.25
    errs = [abs(zeta_LN(n, s0) / n - I_gamma(s0)) for n in (256, 1024, 4096, 16384)]
    conv = errs[-1] < errs[0] / 3 and errs[-1] < 0.02
    ok &= conv
    print(f"  convergence at s=0.25: errors {[f'{e:.4f}' for e in errs]}  -> {'OK' if conv else 'FAIL'}")

    # (D) s=1/2 edge: zeta(1/2)/N - (1/pi)lnN -> (gamma+ln(2/pi))/pi (the U_N log)
    c = (np.euler_gamma + np.log(2 / np.pi)) / np.pi
    edge = zeta_LN(65536, 0.5) / 65536 - np.log(65536) / np.pi
    ed = abs(edge - c) < 5e-3
    ok &= ed
    print(f"\n(D) s=1/2 edge: zeta(1/2)/N - lnN/pi = {edge:.5f} -> {c:.5f}  {'OK' if ed else 'FAIL'}")

    print("\n" + "=" * 74)
    print("RESULT:", "SPECTRAL ZETA VERIFIED" if ok else "CHECK FAILED")
    print("The large-wave growth laws are the s=1/2 (log, Bohr ceiling) and s=1 (sqrt N, variational)")
    print("behaviour of one object: the spectral zeta of the discrete Laplacian. Series <-> spectrum.")
    print("=" * 74)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
