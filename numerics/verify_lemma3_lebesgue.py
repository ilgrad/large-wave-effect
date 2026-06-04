# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2"]
# ///
"""Lemma 3 (exact Lebesgue-type constant) + the continuum / Bessel picture.

(A) RIGOROUS ceiling constant. Euler-Maclaurin on the cosecant sum gives
        U_N = (1/2N) sum_{r=1}^{N-1} csc(pi r/N) = (1/pi)(ln N + gamma + ln(2/pi)) + o(1),
    i.e. U_N - (1/pi) ln N -> C := (gamma + ln(2/pi))/pi ~ 0.03997. This is the ring's
    Lebesgue-type constant; it is unconditional (pure asymptotics of a sum, no number theory).

(B) Does the continuum/Lebesgue route close the LOWER bound for composite N? Test a fixed
    "continuum" instant t = N/4. There omega_r * t ~ pi r/2, so sin(omega_r t) ~ sin(pi r/2)
    alternates in sign and the coherent log collapses to the convergent series
        u_0(N/4) ~ (1/2pi) sum_r sin(pi r/2)/r = (1/2pi)(pi/4) = 1/8 = O(1).
    We confirm u_0 and max_j stay O(1) (do NOT grow with N) -> a fixed instant gives no log;
    the logarithm genuinely needs phase coherence (Kronecker), which is the open arithmetic part.

(C) Bessel / infinite-chain limit. The single-site propagator on Z is the Bessel kernel
    J_n(2t); its peak over sites decays, max_n |J_n(2t)| ~ c (2t)^{-1/3} (Airy front). So the
    large wave VANISHES on the infinite chain: it is a finite-size resonance of the ring's
    discrete spectrum (no large wave on Z). (J_n via its integral.)
"""

from __future__ import annotations

import numpy as np

GAMMA = 0.5772156649015329


def ring_U(n: int) -> float:
    r = np.arange(1, n)
    return float(np.sum(1.0 / np.sin(np.pi * r / n)) / (2 * n))


def u_field(n: int, t: float) -> np.ndarray:
    """u_j(t) = G_N^{(s)}(j,t) for the unit impulse at node 0, all nodes j."""
    r = np.arange(n)
    omega = 2.0 * np.sin(np.pi * r / n)
    coeff = np.zeros(n)
    nz = omega > 1e-12
    coeff[nz] = np.sin(omega[nz] * t) / omega[nz]
    return np.fft.ifft(coeff).real  # ifft carries the 1/N -> equals G^{(s)}


def bessel_J(n_vals: np.ndarray, z: float, m: int = 4000) -> np.ndarray:
    """J_n(z) = (1/pi) int_0^pi cos(n theta - z sin theta) dtheta, vectorised over n."""
    theta = np.linspace(0.0, np.pi, m)
    integrand = np.cos(n_vals[:, None] * theta[None, :] - z * np.sin(theta)[None, :])
    return np.trapezoid(integrand, theta, axis=1) / np.pi


def main() -> int:
    print("=" * 74)
    print("Lemma 3 (Lebesgue constant) + continuum/Bessel picture")
    print("=" * 74)
    ok = True
    inv_pi = 1.0 / np.pi
    C = (GAMMA + np.log(2.0 / np.pi)) * inv_pi

    # (A) exact asymptotic constant
    print(f"\n(A) U_N - (1/pi) ln N -> C = (gamma+ln(2/pi))/pi = {C:.5f}")
    print(f"  {'N':>7} {'U_N':>10} {'U_N-(lnN)/pi':>14}")
    last = None
    for n in (64, 256, 1024, 4096, 16384, 65536):
        u = ring_U(n)
        diff = u - inv_pi * np.log(n)
        print(f"  {n:>7} {u:>10.5f} {diff:>14.5f}")
        last = diff
    okA = abs(last - C) < 5e-3
    print(f"  -> converges to C: {'PASS' if okA else 'FAIL'}")
    ok &= okA

    # (B) fixed continuum instant t = N/4 gives O(1), not log
    print("\n(B) Fixed instant t=N/4: u_0 and max_j stay O(1) (coherence is required for log)")
    print(f"  {'N':>6} {'u_0(N/4)':>10} {'max_j':>9} {'U_N(ceiling)':>13} {'kind':>10}")
    grow = []
    for n in (24, 48, 96, 192, 384, 768):  # composite (rich) N
        u = u_field(n, n / 4.0)
        kind = "2^m*3" if n % 3 == 0 else "composite"
        print(f"  {n:>6} {u[0]:>10.4f} {np.max(np.abs(u)):>9.4f} {ring_U(n):>13.4f} {kind:>10}")
        grow.append(np.max(np.abs(u)))
    bounded = max(grow) < 1.0  # stays O(1) while U_N ~ (1/pi) ln N grows
    print(f"  -> max_j at t=N/4 bounded (O(1)) as N grows: {'PASS' if bounded else 'FAIL'}")
    print("     (analytic: u_0(N/4) -> (1/2pi) sum sin(pi r/2)/r = 1/8; alternating, no log)")
    ok &= bounded

    # (C) Bessel peak decays -> effect vanishes on the infinite chain
    print("\n(C) Infinite-chain Bessel kernel: max_n |J_n(2t)| ~ c (2t)^{-1/3} (decays)")
    print(f"  {'2t':>7} {'max_n|J_n|':>11} {'x^{1/3}*max':>12}")
    zs = np.array([20.0, 60.0, 200.0, 600.0, 2000.0])
    peaks = []
    for z in zs:
        nv = np.arange(0, int(z * 1.3) + 5)
        mx = float(np.max(np.abs(bessel_J(nv, z))))
        peaks.append(mx)
        print(f"  {z:>7.0f} {mx:>11.5f} {z ** (1 / 3) * mx:>12.5f}")
    slope = float(np.polyfit(np.log(zs), np.log(peaks), 1)[0])
    print(f"  -> slope d ln(max)/d ln(2t) = {slope:.3f} (Airy front: -1/3 = {-1/3:.3f})")
    okC = abs(slope + 1 / 3) < 0.05
    print(f"  -> matches -1/3 (peak decays, no large wave on Z): {'PASS' if okC else 'FAIL'}")
    ok &= okC

    print("\n" + "=" * 74)
    print("RESULT:", "ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED")
    print("Honest status: Lemma 3 (ceiling U_N = (1/pi)ln N + C, C=0.0400) is RIGOROUS and")
    print("unconditional. But the continuum/Lebesgue route does NOT by itself close the lower")
    print("bound for composite N: a fixed instant gives O(1), so the log needs phase coherence")
    print("(Kronecker) -- which holds for prime/2^m (Thms 1',1'') and is the open part otherwise.")
    print("The Bessel decay shows the effect is a finite-size resonance (vanishes on Z).")
    print("=" * 74)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
