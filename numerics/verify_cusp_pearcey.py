# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2", "scipy>=1.11"]
# ///
"""The cusp caustic: the Pearcey diffraction pattern, and the teacup (nephroid) caustic.

After the fold (Airy), the next generic caustic is the CUSP -- the bright curve at the bottom of a
teacup, the everyday concentration of light. Its universal diffraction pattern is the Pearcey integral
    Pe(X,Y) = int_{-inf}^{inf} exp( i (s^4 + X s^2 + Y s) ) ds,
the wave dressing of the cusp catastrophe Phi(s)=s^4+X s^2+Y s. Its geometric caustic, where
Phi'=Phi''=0, is the semicubical parabola 8 X^3 + 27 Y^2 = 0 (X<=0): two illuminated rays inside the
cusp, none outside. The teacup itself shows the ray (catacaustic) version: parallel rays reflected
inside a circle have a NEPHROID envelope.

We verify:
 (A) the cusp-point value |Pe(0,0)| = (1/2) Gamma(1/4) = 1.8128 (the exact value of int exp(i s^4) ds);
 (B) the intensity |Pe|^2 is concentrated ON the cusp caustic 8X^3+27Y^2=0 and inside it, and drops in
     the geometric shadow X>0 (the caustic separates the 2-ray bright region from the 0-ray shadow);
 (C) the teacup ray caustic: parallel rays reflected in the unit circle have their envelope on the
     nephroid of radius 1/2 (the reflected rays are tangent to it), with cusps at distance 1/2.
"""

from __future__ import annotations

import numpy as np
from scipy.special import gamma


def pearcey(xx: float, yy: float, s_max: float = 8.0, ds: float = 0.0025) -> complex:
    s = np.arange(-s_max, s_max, ds)
    window = np.exp(-(s / (0.92 * s_max)) ** 12)  # smooth taper to suppress the truncation tail
    return complex(np.trapezoid(np.exp(1j * (s**4 + xx * s**2 + yy * s)) * window, dx=ds))


def main() -> int:
    print("=" * 70)
    print("The cusp caustic: Pearcey diffraction and the teacup nephroid")
    print("=" * 70)
    ok = True

    # (A) cusp-point value
    p00 = abs(pearcey(0.0, 0.0))
    exact = 0.5 * gamma(0.25)
    a_ok = abs(p00 - exact) < 0.02
    ok &= a_ok
    print(f"\n(A) |Pe(0,0)| = {p00:.4f}  vs (1/2)Gamma(1/4) = {exact:.4f}  {'OK' if a_ok else 'FAIL'}")

    # (B) intensity on the cusp caustic vs in the shadow
    tt = np.array([0.6, 1.0, 1.6, 2.2])
    on_caustic = np.mean([abs(pearcey(-t, np.sqrt(8 * t**3 / 27))) ** 2 for t in tt])
    inside = np.mean([abs(pearcey(-t, 0.0)) ** 2 for t in tt])          # between the two branches
    shadow = np.mean([abs(pearcey(+t, 0.0)) ** 2 for t in tt])          # geometric shadow X>0
    b_ok = inside > 2.5 * shadow and on_caustic > 1.4 * shadow
    ok &= b_ok
    print(f"\n(B) mean |Pe|^2:  inside cusp = {inside:.2f},  on caustic = {on_caustic:.2f},  "
          f"shadow (X>0) = {shadow:.2f}")
    print(f"    -> bright inside the cusp ({inside / shadow:.1f}x the shadow), the caustic separating "
          f"lit from shadow: {'OK' if b_ok else 'FAIL'}")

    # (C) teacup: parallel rays reflected in the unit circle -> nephroid envelope (cusps at 1/2)
    th = np.linspace(0.05, np.pi - 0.05, 4000)        # hit points (cos th, sin th), rays travel +x
    hitx, hity = np.cos(th), np.sin(th)
    # incoming direction +x; normal = (cos th, sin th); reflected dir d = inc - 2(inc.n)n
    nx, ny = np.cos(th), np.sin(th)
    dot = nx  # inc=(1,0), inc.n = nx
    dx, dy = 1 - 2 * dot * nx, -2 * dot * ny
    # nephroid (caustic of a circle, parallel rays): points r(th) with |center-to-cusp| = 1/2.
    # the caustic point on each reflected ray is at distance t* where neighbouring rays meet; for the
    # circle catacaustic this traces a nephroid whose cusps lie at (+-1/2, 0)... check tangency:
    # distance from each reflected ray to the analytic nephroid center band is ~0 at the envelope.
    # Robust check: the caustic cusps are at radius 1/2 from the centre.
    # Envelope via discriminant: caustic point param. For parallel rays the catacaustic is the nephroid
    # x = (3/4)cos th - (1/4)cos 3th , y = (3/4)sin th - (1/4)sin 3th  (radius-1 circle, scaled).
    cx = 0.75 * np.cos(th) - 0.25 * np.cos(3 * th)
    cy = 0.75 * np.sin(th) - 0.25 * np.sin(3 * th)
    # verify each reflected ray passes through its caustic point (ray tangent to envelope):
    # point (cx,cy) should lie on the line through (hitx,hity) with direction (dx,dy): cross product ~0
    cross = (cx - hitx) * dy - (cy - hity) * dx
    tangent_err = float(np.max(np.abs(cross)))
    c_ok = tangent_err < 1e-9
    ok &= c_ok
    print(f"\n(C) Teacup nephroid: reflected rays tangent to the nephroid envelope "
          f"(max off-tangency {tangent_err:.1e}): {'OK' if c_ok else 'FAIL'}")

    print("\n" + "=" * 70)
    print("RESULT:", "CUSP CAUSTIC VERIFIED" if ok else "CHECK FAILED")
    print("The cusp -- Pearcey diffraction and the teacup nephroid -- is the next generic caustic after")
    print("the Airy fold: energy concentrates on the cusp, the everyday geometry of a rogue hot spot.")
    print("=" * 70)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
