# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2"]
# ///
"""Staggered-momentum identity and the interaction-Morawetz reduction of the threshold window (problem 3).

The virial convexity (eq:virial) is governed by the staggered momentum Re sum_n conj(a_n) a_{n+2} = P<cos 2k>
(since <L2 a,a> = 2P - 2 Re sum conj(a_n) a_{n+2} = 2P(1 - <cos 2k>)).  The linear flow preserves every
momentum moment |hat a_k|^2, so its evolution is PURELY NONLINEAR and -- unlike the virial -- carries NO
position weight:

    d/dt Re sum_n conj(a_n) a_{n+2}  =  - gamma sum_n (|a_n|^2 - |a_{n+2}|^2) Im(conj(a_n) a_{n+2}),
    |d/dt Re sum_n conj(a_n) a_{n+2}|  <=  2 gamma ||a||_4^4      (Cauchy-Schwarz).

The identity is proven SYMBOLICALLY (Maxima and sympy, ring of N=7) and verified here to machine precision
along the flow by evaluating d/dt directly from the ODE.  Integrating the bound controls the drift by the
L^4 spacetime norm  gamma * int_0^T ||a||_4^4 dt,  but that norm diverges even on the FREE flow:
||a(t)||_4^4 = P^2 sum_n J_n(2t)^4 ~ P^2 (ln 2t + c)/(2 pi^2 t), c ~ 4.05 (PARI), so the integral grows like
P^2 (ln 2T + c)^2/(4 pi^2) for every gamma >= 0 -- and linearly in T above gamma P = 4, where
||a||_4^4 >= P^2 - 4P/gamma.  This norm therefore cannot close the dispersal argument.  (The integrable
Ablowitz-Ladik cousin conserves a staggered-momentum functional that motivates the identity, but it drifts
under the non-integrable flow -- the O(1) AL gap.)  Compare verify_virial_window.py.
"""

from __future__ import annotations

import math
from itertools import pairwise

import numpy as np


def _lap(z: np.ndarray) -> np.ndarray:
    return np.roll(z, -1) - 2 * z + np.roll(z, 1)


def _adot(z: np.ndarray, gamma: float) -> np.ndarray:
    return -1j * (_lap(z) + gamma * np.abs(z) ** 2 * z)


def _stagger(z: np.ndarray) -> float:
    return float(np.real(np.sum(np.conj(z) * np.roll(z, -2))))


def _dstagger_ode(z: np.ndarray, gamma: float) -> float:
    zd = _adot(z, gamma)
    return float(np.real(np.sum(np.conj(z) * np.roll(zd, -2) + np.conj(zd) * np.roll(z, -2))))


def _dstagger_formula(z: np.ndarray, gamma: float) -> float:
    return -gamma * float(
        np.sum(
            (np.abs(z) ** 2 - np.abs(np.roll(z, -2)) ** 2) * np.imag(np.conj(z) * np.roll(z, -2))
        )
    )


def _l4(z: np.ndarray) -> float:
    return float(np.sum(np.abs(z) ** 4))


def _l4_integrals(n: int, gamma: float, times: tuple[float, ...], dt: float) -> list[float]:
    """int_0^T ||a||_4^4 dt at each T in `times` (single-site seed, P = 1, Strang splitting)."""
    lam = 2 - 2 * np.cos(2 * np.pi * np.arange(n) / n)
    z = np.zeros(n, complex)
    z[0] = 1.0
    half = np.exp(1j * lam * dt / 2)
    integral, out = 0.0, []
    for k in range(1, round(times[-1] / dt) + 1):
        z = np.fft.ifft(half * np.fft.fft(z))
        z = np.exp(-1j * gamma * np.abs(z) ** 2 * dt) * z
        z = np.fft.ifft(half * np.fft.fft(z))
        integral += _l4(z) * dt
        if any(abs(k * dt - t) < dt / 2 for t in times):
            out.append(integral)
    return out


def main() -> int:
    rng = np.random.default_rng(1)
    ok = True

    print(
        "[1] staggered-momentum identity  d/dt Re sum a_n* a_{n+2} = -g sum(|a_n|^2-|a_{n+2}|^2)Im(a_n* a_{n+2})"
    )
    print("    (ODE-level, machine precision; symbolically proven in Maxima/sympy):")
    worst_id, worst_bound = 0.0, True
    for _ in range(4):
        nn = 40
        z = rng.standard_normal(nn) + 1j * rng.standard_normal(nn)
        z /= np.sqrt(np.sum(np.abs(z) ** 2))
        g = float(rng.uniform(1, 4))
        diff = abs(_dstagger_ode(z, g) - _dstagger_formula(z, g))
        worst_id = max(worst_id, diff)
        worst_bound &= abs(_dstagger_formula(z, g)) <= 2 * g * _l4(z) + 1e-12
    ok &= worst_id < 1e-10 and worst_bound
    print(
        f"    max |ODE - formula| = {worst_id:.1e};  bound |rate| <= 2g||a||_4^4 holds: {worst_bound}   "
        f"{'ok' if worst_id < 1e-10 and worst_bound else 'FAIL'}"
    )

    # N = 2048 keeps the front |n| <= 2T = 320 far from the wrap-around at N/2 = 1024
    times = (40.0, 80.0, 160.0)
    print(
        "[2] the L^4 spacetime norm int_0^T ||a||_4^4 dt diverges, even on the free flow (P=1, N=2048):"
    )
    print(
        f"    {'gP':>4} "
        + " ".join(f"{'T=' + str(int(t)):>8}" for t in times)
        + "   increments per doubling"
    )
    free_pred = [(math.log(2 * t) + 4.05) ** 2 / (4 * math.pi**2) for t in times]
    grows = True
    for g in (0.0, 1.0, 3.0, 4.5):
        vals = _l4_integrals(2048, g, times, 0.01)
        inc = [b - a for a, b in pairwise(vals)]
        grows &= inc[1] > inc[0] > 0
        print(
            f"    {g:>4.1f} "
            + " ".join(f"{v:>8.3f}" for v in vals)
            + "   "
            + "  ".join(f"{d:.3f}" for d in inc)
        )
        if g == 0.0:
            pred_inc = [b - a for a, b in pairwise(free_pred)]
            law = all(abs(d - q) < 0.01 for d, q in zip(inc, pred_inc, strict=True))
            grows &= law
            print(
                f"    free increments vs (ln 2T + 4.05)^2/(4 pi^2): {'  '.join(f'{q:.3f}' for q in pred_inc)}"
                f"  {'ok' if law else 'FAIL'}"
            )
    ok &= grows
    print(
        f"    increments grow for every gP (no bounded L^4 norm below 4): {'ok' if grows else 'FAIL'}"
    )

    print("=" * 76)
    print(
        "RESULT:",
        "PASS -- staggered-momentum identity verified; its L^4 spacetime control diverges like "
        "(ln T)^2 even on the free flow, so it cannot close the dispersal argument"
        if ok
        else "FAIL",
    )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
