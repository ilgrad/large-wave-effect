# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2"]
# ///
"""Branched flow: a wave through a weak random current focuses into caustic hot spots (rogue spots).

Ocean rogue "hot spots" (the Agulhas current off SE Africa being the classic example) arise when swell
crosses a spatially varying current: the current refracts the wave rays, neighbouring rays cross, and
energy concentrates on CAUSTICS. A weak random current field produces a cascade of such caustics --
"branched flow" -- whose first focal distance and heavy-tailed intensity statistics are the LINEAR
(refractive) route to large local amplitude, complementary to the nonlinear self-focusing of Part II.

We model this with the paraxial wave (Schrodinger) equation in one transverse dimension,
    i d_z psi = -(1/2) d_xx psi + V(x) psi,
with V a weak, smooth, static Gaussian random field (the frozen current pattern), propagated by a
unitary split-step integrator. We verify:
 (A) the propagator is unitary (norm conserved to ~1e-12);
 (B) starting from a uniform wave, the intensity develops caustic HOT SPOTS: the scintillation index
     max_x|psi|^2 / mean rises from 1 to >> 1, peaking at a finite focal distance z_f;
 (C) at the caustic the intensity is HEAVY-TAILED -- its exceedance lies far above the exponential law
     of a structureless random-phase field (intensity kurtosis >> 2) -- the rogue-hot-spot signature.
"""

from __future__ import annotations

import numpy as np

RNG = np.random.default_rng(0)


def random_potential(n: int, vrms: float, corr: float) -> np.ndarray:
    """Smooth static Gaussian random field on the ring: white noise filtered to correlation length corr."""
    k = 2 * np.pi * np.fft.fftfreq(n)
    filt = np.exp(-0.5 * (k * corr) ** 2)
    w = np.fft.fft(RNG.standard_normal(n))
    v = np.real(np.fft.ifft(w * filt))
    return vrms * v / np.std(v)


def propagate(n: int, vrms: float, corr: float, dz: float, zsteps: int) -> dict:
    v = random_potential(n, vrms, corr)
    k2 = (2 * np.pi * np.fft.fftfreq(n)) ** 2
    kin = np.exp(-1j * k2 * dz / 2)
    half = np.exp(-1j * v * dz / 2)
    psi = np.ones(n, dtype=complex)  # uniform incident wave
    norm0 = float(np.sum(np.abs(psi) ** 2))
    ndrift, sci, zs, fields = 0.0, [], [], []
    best_sci, best_field = 0.0, np.abs(psi) ** 2
    for s in range(zsteps):
        psi = half * psi
        psi = np.fft.ifft(kin * np.fft.fft(psi))
        psi = half * psi
        inten = np.abs(psi) ** 2
        rr = float(inten.max() / inten.mean())
        sci.append(rr)
        zs.append(s * dz)
        if rr > best_sci:
            best_sci, best_field = rr, inten.copy()
        if s % 6 == 0:
            fields.append(inten.copy())
        if s % (zsteps // 10) == 0:
            ndrift = max(ndrift, abs(np.sum(np.abs(psi) ** 2) - norm0) / norm0)
    return {"v": v, "sci": np.array(sci), "z": np.array(zs), "ndrift": ndrift,
            "field": np.array(fields), "best_field": best_field, "final": np.abs(psi) ** 2}


def main() -> int:
    print("=" * 70)
    print("Branched flow: caustic hot spots of a wave in a random current")
    print("=" * 70)
    ok = True
    n, vrms, corr, dz, zsteps = 1024, 0.12, 15.0, 0.05, 3000
    r = propagate(n, vrms, corr, dz, zsteps)

    print(f"\n(A) Unitary paraxial split-step: norm drift = {r['ndrift']:.2e}")
    a_ok = r["ndrift"] < 1e-10
    ok &= a_ok
    print(f"    norm conserved: {'OK' if a_ok else 'FAIL'}")

    sci, zz = r["sci"], r["z"]
    sci_max = float(sci.max())
    z_f = float(zz[int(sci.argmax())])
    print("\n(B) Scintillation index max|psi|^2/mean vs propagation distance z")
    print(f"    {'z':>7} {'max/mean':>9}")
    for zt in (0.0, z_f / 2, z_f, 2 * z_f):
        i = int(np.argmin(np.abs(zz - zt)))
        print(f"    {zz[i]:>7.1f} {sci[i]:>9.2f}")
    b_ok = sci[0] < 1.2 and sci_max > 8 and 0 < z_f < zz[-1]
    ok &= b_ok
    print(f"    -> hot spots form: index 1 -> {sci_max:.1f} at focal distance z_f={z_f:.1f} "
          f"{'OK' if b_ok else 'FAIL'}")

    # (C) heavy-tailed intensity at the focal (max-scintillation) plane
    inten = r["best_field"] / r["best_field"].mean()
    kurt = float((inten**2).mean() / inten.mean() ** 2)  # <I^2>/<I>^2; =2 for exponential (Rayleigh)
    frac = float(np.mean(inten > 8))  # P(I > 8<I>): exponential gives e^-8=3.4e-4
    c_ok = kurt > 2.3 and frac > 5 * np.exp(-8.0)
    ok &= c_ok
    print(f"\n(C) Intensity at the focal plane: <I^2>/<I>^2 = {kurt:.2f} (exponential field = 2.0), "
          f"P(I>8<I>) = {frac:.2e}")
    print(f"    -> heavy-tailed: variance above exponential, tail {frac / np.exp(-8.0):.0f}x the "
          f"exponential law e^-8={np.exp(-8.0):.1e}: {'OK' if c_ok else 'FAIL'}")

    print("\n" + "=" * 70)
    print("RESULT:", "BRANCHED FLOW VERIFIED" if ok else "CHECK FAILED")
    print("A weak random current refracts a uniform wave into caustic hot spots with heavy-tailed")
    print("intensity -- the linear refractive route to rogue waves (the Agulhas-type lens), distinct")
    print("from and complementary to the nonlinear self-focusing of Part II.")
    print("=" * 70)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
