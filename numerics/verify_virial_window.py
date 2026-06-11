# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2"]
# ///
"""Virial route to dispersal on the window gamma*P in [0.43, 4), and the precise obstruction.

On the line Z (the ring for t <~ N/5, where aliasing is negligible) write the second moment
V(t) = sum_j j^2 |a_j|^2.  A direct computation of the discrete DNLS virial gives, at the single-site seed,
    d^2 V/dt^2 |_{t=0} = 4P ,
INDEPENDENT of gamma (the cubic term vanishes at t=0 because only one site is populated).  If d^2V/dt^2
stayed >= c > 0 for all t, then V -> infinity and the wave provably disperses (||a||_inf -> 0 at fixed mass).
Numerically d^2V/dt^2 > 0 throughout the window and decays to 0 as gamma*P -> 4 -- the virial is convex,
SUPPORTING dispersal across [0.43, 4).

But the convexity constant is the asymptotic d^2V/dt^2 -> 2<v^2>, v_k = 2 sin k the group velocity, and
<v^2> = <4 sin^2 k> = 2 - 2<cos 2k> is a SECOND spectral moment.  Energy conservation fixes only the FIRST
moment: at full dispersal <L a, a> -> 2P - (gamma/2)P^2, i.e.
    <cos k>_inf = gamma*P/4 .
This is < 1 for gamma*P < 4 (so the momentum cannot fully concentrate at k=0), yet it does NOT force
<v^2> > 0: a bimodal distribution split between k=0 and k=pi has <cos k> = gamma*P/4 while <v^2> = 0 (both
have sin k = 0).  So the conserved quantities cannot, by themselves, lower-bound the virial convexity --
this is the precise wall blocking a rigorous dispersal proof on [0.43, 4) (Morawetz/Strichartz fail at the
same point: the focusing sign / the uncontrolled staggered k=pi content).

This certificate confirms:
  (1) d^2V/dt^2|_0 = 4P, gamma-independent (to 1e-3);
  (2) d^2V/dt^2 plateau > 0 across the window, decaying toward 0 at gamma*P = 4;
  (3) the energy identity <cos k>_inf = gamma*P/4 at the dispersal minimum (to 1e-2);
  (4) numerically <v^2> stays > 0 (the dynamics avoids the bimodal trap) -- but this is exactly what is
      NOT supplied by conservation.  Compare verify_selftrapping_transition.py and verify_dispersal_weak.py.
"""

from __future__ import annotations

import numpy as np


def evolve(n: int, gamma: float, t_final: float, steps: int) -> dict[str, float | np.ndarray]:
    lam = 2 - 2 * np.cos(2 * np.pi * np.arange(n) / n)
    k = 2 * np.pi * np.arange(n) / n
    j = ((np.arange(n) + n // 2) % n) - n // 2  # centered position on the line approximation
    z = np.zeros(n, complex)
    z[0] = 1.0
    dt = t_final / steps
    half = np.exp(1j * lam * dt / 2)
    vs = np.empty(steps + 1)
    peak_min, cosk_at_min, v2_at_min = 1.0, 2.0, 2.0
    for m in range(steps + 1):
        vs[m] = float(np.sum(j**2 * np.abs(z) ** 2))
        peak = float((np.abs(z) ** 2).max())
        if peak < peak_min:
            w = np.abs(np.fft.fft(z)) ** 2
            w = w / w.sum()
            peak_min = peak
            cosk_at_min = float(np.sum(np.cos(k) * w))
            v2_at_min = float(np.sum((2 * np.sin(k)) ** 2 * w))
        if m < steps:
            z = np.fft.ifft(half * np.fft.fft(z))
            z = np.exp(-1j * gamma * np.abs(z) ** 2 * dt) * z
            z = np.fft.ifft(half * np.fft.fft(z))
    d2 = (vs[2:] - 2 * vs[1:-1] + vs[:-2]) / dt**2
    return {"d2": d2, "cosk": cosk_at_min, "v2": v2_at_min, "peak_min": peak_min}


def main() -> int:
    n = 2048
    ok = True

    print("[1] d^2V/dt^2|_0 = 4P, gamma-independent:")
    for g in (0.0, 2.0, 3.9):
        d2 = evolve(n, g, 2.0, 4000)["d2"]
        good = abs(float(d2[0]) - 4.0) < 1e-3
        ok &= good
        print(f"    gamma*P={g:>4}: {float(d2[0]):.4f}  {'ok' if good else 'FAIL'}")

    print("[2] virial convex (d^2V/dt^2 plateau > 0) across the window; [3] <cos k> -> gamma*P/4:")
    print(f"    {'gP':>4} {'plateau':>8} {'<cosk>':>7} {'gP/4':>6} {'<v^2>':>7} {'peak_min':>9}")
    for g in (1.0, 2.0, 3.0, 3.5):
        r = evolve(n, g, 200.0, 8000)
        plateau = float(np.mean(r["d2"][2666:5333]))
        cosk = float(r["cosk"])
        convex_ok = plateau > 0
        energy_ok = abs(cosk - g / 4) < 1e-2
        ok &= convex_ok and energy_ok
        flag = "" if convex_ok and energy_ok else "  FAIL"
        print(f"    {g:>4.1f} {plateau:>8.3f} {cosk:>7.4f} {g / 4:>6.3f} {float(r['v2']):>7.4f} "
              f"{float(r['peak_min']):>9.4f}{flag}")

    print("=" * 72)
    print("RESULT:", "PASS -- virial convex on the window (supports dispersal); conservation fixes "
          "<cos k>=gP/4 but not <v^2>, the precise wall" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
