# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2"]
# ///
"""A linear caustic SEEDS nonlinear self-focusing: the two energy-concentration routes cooperate.

Ocean rogue waves combine the two mechanisms of this paper: a current/bathymetry lens concentrates
energy by linear refraction (a caustic), raising the local steepness until it crosses the modulational
threshold, whereupon nonlinear self-focusing takes over -- rogue waves are "the nonlinear stage of the
modulational instability" (Zakharov & Shamin). We model this on the focusing NLS
    i psi_t + (1/2) psi_xx + |psi|^2 psi = 0,
seeding a Gaussian beam with a converging (lens) phase, psi(x,0) = A exp(-x^2/2w^2) exp(-i b x^2/2),
which focuses at a finite distance under the LINEAR flow. We compare three runs at the same weak
amplitude A (below the bare modulational threshold over this time):
  (i)   plane wave + nonlinearity     -> stays flat (a clean uniform wave has no seed for MI);
  (ii)  lens + linear flow (gamma=0)  -> diffractive focusing concentrates energy ~7x (the caustic);
  (iii) lens + nonlinearity           -> self-focusing sharpens the focal peak on top of the caustic.

So the caustic does the bulk of the concentration and the nonlinearity adds the self-focusing stage
at the focus; the lens is what localizes the energy in the first place (without it, the flat wave never
spikes). In 1D the nonlinear enhancement is modest (~10-15%, the NLS being L^2-subcritical, no
collapse); the dramatic collapsing self-focusing is a 2D effect (Zakharov wave collapse), the natural
extension.

We verify (norm-conserving split-step): (A) norm conserved; (B) the lens focuses linearly (peak rises
~7x at a focal time z_f); (C) the lens+nonlinear peak exceeds the linear-caustic peak (self-focusing
adds), while the lensless plane wave stays flat (the lens is required to localize).
"""

from __future__ import annotations

import numpy as np


def step(psi: np.ndarray, dt: float, gamma: float, k2: np.ndarray) -> np.ndarray:
    psi = psi * np.exp(1j * gamma * np.abs(psi) ** 2 * dt / 2)
    psi = np.fft.ifft(np.exp(-1j * k2 * dt / 2) * np.fft.fft(psi))
    return psi * np.exp(1j * gamma * np.abs(psi) ** 2 * dt / 2)


def run(x: np.ndarray, psi0: np.ndarray, gamma: float, dt: float, steps: int) -> dict:
    n = len(x)
    k2 = (2 * np.pi * np.fft.fftfreq(n, d=(x[1] - x[0]))) ** 2
    psi = psi0.astype(complex).copy()
    norm0 = float(np.sum(np.abs(psi) ** 2))
    peak, t_peak, ndrift = float(np.abs(psi).max()), 0.0, 0.0
    for s in range(steps):
        psi = step(psi, dt, gamma, k2)
        a = float(np.abs(psi).max())
        if a > peak:
            peak, t_peak = a, s * dt
        if s % (steps // 10) == 0:
            ndrift = max(ndrift, abs(np.sum(np.abs(psi) ** 2) - norm0) / norm0)
    return {"peak": peak, "t_peak": t_peak, "ndrift": ndrift, "final": np.abs(psi)}


def main() -> int:
    print("=" * 72)
    print("A linear caustic seeds nonlinear self-focusing (lens + NLS)")
    print("=" * 72)
    ok = True
    n = 4096
    x = np.linspace(-80, 80, n, endpoint=False)
    amp, w, b = 0.42, 14.0, 0.28  # weak amplitude, beam width, lens curvature
    dt, steps = 0.01, 2500
    lens = amp * np.exp(-x**2 / (2 * w**2)) * np.exp(-1j * b * x**2 / 2)
    plane = amp * np.ones(n, dtype=complex)

    r_lin = run(x, lens, gamma=0.0, dt=dt, steps=steps)    # caustic only
    r_nl = run(x, lens, gamma=1.0, dt=dt, steps=steps)     # caustic + self-focusing
    r_pw = run(x, plane, gamma=1.0, dt=dt, steps=steps)    # plane wave + nonlinearity

    print(f"\n(A) Norm conservation (lens+NLS): drift = {r_nl['ndrift']:.2e}")
    a_ok = r_nl["ndrift"] < 1e-9
    ok &= a_ok
    print(f"    norm-conserving split-step: {'OK' if a_ok else 'FAIL'}")

    print("\n(B) Linear lens focuses: caustic peak vs the flat input amplitude")
    b_ok = r_lin["peak"] > 1.5 * amp and r_lin["t_peak"] > 0
    ok &= b_ok
    print(f"    input |psi|={amp:.2f} -> linear caustic peak={r_lin['peak']:.3f} at z_f={r_lin['t_peak']:.1f}  "
          f"{'OK' if b_ok else 'FAIL'}")

    print("\n(C) Caustic SEEDS self-focusing: lens+NLS peak exceeds both references")
    print(f"    {'run':>26} {'peak':>8} {'/input':>8}")
    print(f"    {'(i) plane wave + NLS':>26} {r_pw['peak']:>8.3f} {r_pw['peak'] / amp:>8.2f}")
    print(f"    {'(ii) lens, linear (caustic)':>26} {r_lin['peak']:>8.3f} {r_lin['peak'] / amp:>8.2f}")
    print(f"    {'(iii) lens + NLS':>26} {r_nl['peak']:>8.3f} {r_nl['peak'] / amp:>8.2f}")
    enh = r_nl["peak"] / r_lin["peak"] - 1
    c_ok = r_nl["peak"] > 1.08 * r_lin["peak"] and r_pw["peak"] < 1.3 * amp
    ok &= c_ok
    print(f"    -> self-focusing adds {enh:.0%} on top of the caustic; the lensless plane wave stays "
          f"flat ({r_pw['peak'] / amp:.2f}x): {'OK' if c_ok else 'FAIL'}")

    print("\n" + "=" * 72)
    print("RESULT:", "CAUSTIC-SEEDED SELF-FOCUSING VERIFIED" if ok else "CHECK FAILED")
    print("The linear lens (caustic) concentrates energy past the modulational threshold, where")
    print("self-focusing takes over: rogue waves as the nonlinear stage of the instability.")
    print("=" * 72)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
