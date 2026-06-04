# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2", "scipy>=1.11"]
# ///
"""The discrete-Schrodinger wavefront is a fold caustic: the Airy function, scaling as t^{-1/3}.

Part I's free discrete Schrodinger kernel on Z is K(n,t) = e^{-2it} i^n J_n(2t). The ballistic front
sits at n = 2t, where the lattice group velocity v_g(k) = 2 sin k is extremal (v_g' = 0 at k = pi/2): a
degenerate stationary-phase point, i.e. a FOLD CAUSTIC. Catastrophe theory predicts the universal Airy
profile with width ~ t^{1/3} and height ~ t^{-1/3}. The uniform Bessel asymptotics give precisely
    J_n(2t) ~ t^{-1/3} Ai( 2^{1/3} (n - 2t) / (2t)^{1/3} ),
so the discrete large-wave front IS the Airy caustic. This is the lattice instance of the same
energy-concentration-by-refraction that focuses ocean swell into rogue hot spots: rays (here, the band
of group velocities) pile up on a caustic, where the amplitude concentrates with the universal -1/3
exponent of a fold.

We verify, with scipy's Bessel and Airy functions:
 (A) the front peak max_n |J_n(2t)| scales as t^{-1/3} (the fold-caustic exponent);
 (B) under the caustic scaling s = 2^{1/3}(n-2t)/(2t)^{1/3} the rescaled front t^{1/3} J_n(2t) collapses
     onto the Airy function Ai(s), with the error vanishing as t grows.
"""

from __future__ import annotations

import numpy as np
from scipy.special import airy, jv


def front_window(t: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n = np.arange(int(2 * t) - 60, int(2 * t) + 18)
    j = jv(n, 2 * t)
    s = 2.0 ** (1 / 3) * (n - 2 * t) / (2 * t) ** (1 / 3)
    return n, j, s


def main() -> int:
    print("=" * 70)
    print("Discrete-Schrodinger front as a fold caustic: Airy profile, t^{-1/3}")
    print("=" * 70)
    ok = True

    # (A) peak scaling ~ t^{-1/3}
    print("\n(A) front peak max_n |J_n(2t)| ~ C t^{-1/3} (fold-caustic exponent -1/3)")
    ts = np.array([50.0, 100.0, 200.0, 400.0, 800.0])
    peaks = np.array([np.abs(front_window(t)[1]).max() for t in ts])
    slope = np.polyfit(np.log(ts), np.log(peaks), 1)[0]
    print(f"    {'t':>6} {'peak':>10} {'peak*t^(1/3)':>13}")
    for t, p in zip(ts, peaks, strict=True):
        print(f"    {t:>6.0f} {p:>10.5f} {p * t ** (1 / 3):>13.5f}")
    a_ok = abs(slope + 1 / 3) < 0.02
    ok &= a_ok
    print(f"    log-log slope = {slope:.4f}  (Airy fold = -1/3 = {-1 / 3:.4f}): {'OK' if a_ok else 'FAIL'}")
    print(f"    peak*t^(1/3) -> {peaks[-1] * ts[-1] ** (1 / 3):.4f}  (Airy max |Ai| = 0.5357)")

    # (B) Airy collapse: t^{1/3} J_n(2t) -> Ai(s)
    print("\n(B) caustic-scaled front t^{1/3} J_n(2t) collapses onto Ai(s)")
    print(f"    {'t':>6} {'max|scaled - Ai(s)|':>20}")
    errs = []
    for t in ts:
        _, j, s = front_window(t)
        ai = airy(s)[0]
        scaled = t ** (1 / 3) * j
        win = s > -6  # the caustic/illuminated side and the first oscillations
        err = float(np.max(np.abs(scaled[win] - ai[win])))
        errs.append(err)
        print(f"    {t:>6.0f} {err:>20.4f}")
    b_ok = errs[-1] < 0.02 and errs[-1] < errs[0]
    ok &= b_ok
    print(f"    -> rescaled front converges to the Airy function (err -> 0): {'OK' if b_ok else 'FAIL'}")

    print("\n" + "=" * 70)
    print("RESULT:", "AIRY CAUSTIC VERIFIED" if ok else "CHECK FAILED")
    print("The discrete large-wave front is a fold caustic: J_n(2t) ~ t^{-1/3} Ai(...). Energy")
    print("concentrates on the caustic with the universal -1/3 exponent -- the lattice form of the")
    print("refractive focusing that makes ocean rogue hot spots.")
    print("=" * 70)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
