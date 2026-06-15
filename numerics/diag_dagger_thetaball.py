# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2", "scipy>=1.11", "mpmath>=1.3", "python-flint>=0.8"]
# ///
"""Diagnostic: why the signed (theta-ball) middle cannot extend past t~12.

`middle_point` encloses the theta-integral with theta-BALLS of half-width dth/2.  For t>~12 the
integrand J0(amp*t*cos theta)... oscillates fast in theta (frequency ~ amp*t), so the ball radius
of the enclosure scales like amp*t/nth -- it is the rigorous RADIUS, not the midpoint, that blows up.
This probe evaluates the signed middle on a single large-t slice [12,20] for the worst (largest-b)
amplitude and shows:
  * the midpoint is correct and stable across precision (the 0F1 J0-series midpoint is fine);
  * the RADIUS is precision-INDEPENDENT and shrinks only like 1/nth -> ~1e7 theta-points would be
    needed to reach a usable ~1e-2 radius.  Hence t=12 is the theta-ball viability limit, and going
    further needs a theta-TAYLOR model (see middle_far_point).

    uv run --script numerics/diag_dagger_thetaball.py
"""

import importlib.util
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CERTP = ROOT / "numerics" / "verify_dagger_continuum.py"
A = 0.685  # large-b worst node near the residual diagonal


def _load(prec):
    spec = importlib.util.spec_from_file_location("cert", CERTP)
    m = importlib.util.module_from_spec(spec)
    sys.modules["cert"] = m
    spec.loader.exec_module(m)
    m.ctx.prec = prec
    return m


def run(args):
    prec, nth = args
    cert = _load(prec)
    b = cert.b_of(A, A)
    t0 = time.time()
    v = cert.middle_point(A, A, b, [(16.0, 20.0, 16, nth)])
    return (prec, nth, float(v.mid()), float(v.rad()), time.time() - t0)


if __name__ == "__main__":
    cert = _load(200)
    b = cert.b_of(A, A)
    print(f"signed middle on [16,20], a1=a3={A}, b={b:.4f} (z=b*20={b * 20:.1f}); midpoint vs RADIUS:")
    print(f"{'prec':>6} {'nth':>7} {'mid':>14} {'radius':>12} {'sec':>6}")
    jobs = [(200, 500), (200, 1000), (200, 2000), (400, 1000), (800, 1000)]
    with ProcessPoolExecutor(max_workers=min(8, len(jobs))) as pool:
        for (prec, nth, mid, rad, dt) in pool.map(run, jobs):
            print(f"{prec:>6} {nth:>7} {mid:14.6e} {rad:12.3e} {dt:6.1f}")
    print("\nradius ~ 1/nth and ~independent of prec => theta-ball limit; use theta-Taylor beyond t~12.")
