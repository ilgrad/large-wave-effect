# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2", "scipy>=1.11", "mpmath>=1.3", "python-flint>=0.8"]
# ///
"""Validate middle_far_point (the theta-Taylor + t-collocation signed middle for large t) against a
high-accuracy scipy Kluyver reference.

middle_far_point rigorously encloses int_far (Mg - M)/t^2 dt where the theta-ball middle_point fails
(its rigorous radius blows up like amp*t/nth past t~12).  The test checks, for several amplitudes on
and off the equal-amplitude diagonal, that the rigorous enclosure (i) CONTAINS the scipy reference and
(ii) has a small radius -- i.e. the method is both sound and tight.  This is the missing piece that
turns the beta_odd continuum cover (open problem A2) from a theta-ball wall into an executable tiling;
the per-box certificate is prove_box(..., far_slices=...), tiled by run_dagger_campaign.py.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import j0

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("cert", ROOT / "numerics" / "verify_dagger_continuum.py")
cert = importlib.util.module_from_spec(spec)
sys.modules["cert"] = cert
spec.loader.exec_module(cert)
if not cert.HAVE_FLINT:
    print("[python-flint missing: cannot run the rigorous middle_far validation]")
    raise SystemExit(1)
cert.ctx.prec = 200


def true_far(a1: float, a3: float, b: float, t0: float, t1: float) -> float:
    nth = 8000
    th = (np.arange(nth) + 0.5) / nth * (np.pi / 2)
    c, s = np.cos(th), np.sin(th)
    sx2, sy2 = 0.5 * (a1**2 + a3**2), 0.5 * b**2
    M = lambda t: float((j0(a1 * t * c) * j0(a3 * t * c) * j0(b * t * s)).mean())  # noqa: E731
    Mg = lambda t: float(np.exp(-(sx2 * c**2 + sy2 * s**2) * t**2 / 2).mean())  # noqa: E731
    return quad(lambda t: (Mg(t) - M(t)) / t**2, t0, t1, limit=800)[0]


def main() -> int:
    slices = [(8.0, 12.0, 16, 160), (12.0, 16.0, 12, 130), (16.0, 20.0, 12, 110), (20.0, 24.0, 12, 110)]
    t0, t1 = 8.0, 24.0
    print("=" * 92)
    print("middle_far_point (theta-Taylor + t-collocation) vs scipy reference on [8,24]:")
    print(f"  {'a1':>6} {'a3':>6} {'b':>7} {'true':>13} {'rig mid':>13} {'rig rad':>10} {'contains':>9}")
    ok = True
    for (a1, a3) in [(0.73, 0.73), (0.70, 0.70), (0.685, 0.685), (0.82, 0.591), (0.66, 0.66), (0.78, 0.69)]:
        b = cert.b_of(a1, a3)
        rig = cert.middle_far_point(a1, a3, b, slices)
        tr = true_far(a1, a3, b, t0, t1)
        lo, hi = float(rig.lower()), float(rig.upper())
        contains = lo <= tr <= hi
        small = float(rig.rad()) < 1e-4
        ok &= contains and small
        print(f"  {a1:6.3f} {a3:6.3f} {b:7.4f} {tr:13.6e} {float(rig.mid()):13.6e} {float(rig.rad()):10.2e} "
              f"{'YES' if contains else 'NO!!':>9}")
    print("=" * 92)
    print(f"RESULT: {'PASS' if ok else 'FAIL'} -- middle_far enclosure contains the reference with small radius "
          f"(sound + tight)")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
