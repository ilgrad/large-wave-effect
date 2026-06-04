# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2"]
# ///
"""Ray tracing in a current jet: parallel rays bend and cross on a caustic (geometric optics).

Geometric optics is the skeleton of the refractive focusing of Section ``Caustics''. A wave packet
follows rays x_dot = d omega/d k, k_dot = -d omega/d x. For surface waves crossing a current U(x), the
Doppler term k.U(x) refracts the rays toward the slow (opposing-current) region; in the paraxial limit
this is the ray equation x'' = -K x exp(-x^2/2 sigma^2) of a focusing lens. A bundle of initially
parallel rays therefore converges and, where neighbouring rays cross, forms a CAUSTIC: the ray map
x0 -> x(z) folds (its Jacobian J = dx/dx0 passes through zero), the ray density 1/|J| diverges, and the
linear amplitude concentrates. This is the same fold that the discrete-Schrodinger front realizes as
the Airy function (verify_caustic_airy.py).

We verify, with a symplectic (leapfrog) ray integrator:
 (A) the Jacobian J(z) of the ray map starts at 1 and first vanishes at a finite focal distance z_c
     (the caustic), the same focal distance at which the ray bundle width is minimal;
 (B) past z_c the ray map is multivalued -- up to three rays reach a transverse point -- the fold
     signature of a caustic;
 (C) outside the lens (K=0) rays stay parallel and no caustic forms (J == 1).
"""

from __future__ import annotations

import numpy as np


def trace(x0: np.ndarray, k: float, sigma: float, dz: float, zsteps: int) -> np.ndarray:
    """Leapfrog ray positions x(z); returns array (zsteps+1, n_rays). x'' = -k x exp(-x^2/2 sigma^2)."""
    x = x0.copy()
    v = np.zeros_like(x0)
    acc = -k * x * np.exp(-x**2 / (2 * sigma**2))
    out = [x.copy()]
    for _ in range(zsteps):
        v += 0.5 * dz * acc
        x = x + dz * v
        acc = -k * x * np.exp(-x**2 / (2 * sigma**2))
        v += 0.5 * dz * acc
        out.append(x.copy())
    return np.array(out)


def main() -> int:
    print("=" * 70)
    print("Ray tracing in a current jet: parallel rays cross on a caustic")
    print("=" * 70)
    ok = True
    k, sigma, dz, zsteps = 0.04, 12.0, 0.1, 1200
    x0 = np.linspace(-30, 30, 1201)
    xz = trace(x0, k, sigma, dz, zsteps)  # (Z+1, n)
    z = np.arange(zsteps + 1) * dz

    # (A) Jacobian J = dx/dx0 (finite differences); first zero crossing = caustic
    jac = np.gradient(xz, x0, axis=1)
    jmin = jac.min(axis=1)  # most-folded ray at each z
    caustic_idx = int(np.argmax(jmin < 0))  # first z where some J<0 (rays crossed)
    z_c = z[caustic_idx]
    central = np.abs(x0) < sigma
    spread = xz[:, central].std(axis=1)
    z_focus = z[int(spread.argmin())]
    a_ok = 0 < z_c < z[-1] and abs(z_c - z_focus) < 0.2 * z_c + 5
    ok &= a_ok
    print(f"\n(A) Caustic: Jacobian first vanishes at z_c = {z_c:.1f}; bundle focal z = {z_focus:.1f}  "
          f"{'OK' if a_ok else 'FAIL'}")

    # (B) fold: at a plane past z_c, count rays reaching a target transverse point (1 -> 3)
    zp = int((caustic_idx + zsteps) / 2)
    xp = xz[zp]
    target = 0.0
    crossings = int(np.sum(np.abs(np.diff(np.sign(xp - target))) > 0))
    b_ok = crossings >= 3
    ok &= b_ok
    print(f"(B) Past the caustic (z={z[zp]:.0f}): {crossings} rays reach x=0 (fold gives >=3): "
          f"{'OK' if b_ok else 'FAIL'}")

    # (C) no lens (K=0): rays stay parallel, no caustic
    xz0 = trace(x0, 0.0, sigma, dz, zsteps)
    j0 = np.gradient(xz0, x0, axis=1)
    c_ok = bool(np.all(j0 > 0.99))
    ok &= c_ok
    print(f"(C) No lens (K=0): min Jacobian = {j0.min():.3f} (==1, rays parallel): "
          f"{'OK' if c_ok else 'FAIL'}")

    print("\n" + "=" * 70)
    print("RESULT:", "RAY CAUSTIC VERIFIED" if ok else "CHECK FAILED")
    print("A current-jet lens bends parallel rays into a caustic (fold), where the ray density and the")
    print("linear wave amplitude concentrate -- the geometric-optics skeleton of refractive focusing.")
    print("=" * 70)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
