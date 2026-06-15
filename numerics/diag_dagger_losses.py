# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2", "scipy>=1.11", "mpmath>=1.3", "python-flint>=0.8"]
# ///
"""Diagnostic: the TRUE margin and the rigorous-enclosure LOSS budget of the dagger continuum
certificate, at a chosen equal-amplitude box.

It prints, in R/K* units (so the numbers are additive and directly comparable to 1):
  * the TRUE (scipy Kluyver) head/middle/tail split and total R/K* -- the actual value, which
    is what the certificate must come in under;
  * the RIGOROUS middle enclosure at a couple of theta resolutions (the O(1/nth) theta-disc loss);
  * the RIGOROUS capped-|M| tail vs the true (signed) tail (the cancellation the monotone envelope
    cannot see).
This is the measurement behind the finding that the residual is enclosure-limited (true margin ~0.12),
not near-tight, and that the binding losses are the capped tail and the middle theta-discretization.

    uv run --script numerics/diag_dagger_losses.py            # default a1=a3=0.73
    uv run --script numerics/diag_dagger_losses.py 0.70 0.70
"""

import importlib.util
import sys
import time
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import j0

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("cert", ROOT / "numerics" / "verify_dagger_continuum.py")
cert = importlib.util.module_from_spec(spec)
sys.modules["cert"] = cert
spec.loader.exec_module(cert)
cert.ctx.prec = 200

a1 = float(sys.argv[1]) if len(sys.argv) > 1 else 0.73
a3 = float(sys.argv[2]) if len(sys.argv) > 2 else 0.73
b = cert.b_of(a1, a3)
D = a1**4 + a3**4 + b**4
KD = cert.K_STAR * D
print(f"point a1={a1} a3={a3} b={b:.5f} lambda={cert.lam_conc(a1, a3, b):.4f}  K*={cert.K_STAR:.8f}")

nth = 8000
th = (np.arange(nth) + 0.5) / nth * (np.pi / 2)
c, s = np.cos(th), np.sin(th)
sx2, sy2 = 0.5 * (a1**2 + a3**2), 0.5 * b**2
M = lambda t: float((j0(a1 * t * c) * j0(a3 * t * c) * j0(b * t * s)).mean())  # noqa: E731
Mg = lambda t: float(np.exp(-(sx2 * c**2 + sy2 * s**2) * t**2 / 2).mean())  # noqa: E731
f = lambda t: (Mg(t) - M(t)) / t**2  # noqa: E731
th_true = quad(f, 1e-7, 1.5, limit=400)[0] / KD
md_true = quad(f, 1.5, 12.0, limit=800)[0] / KD
tl_true = quad(f, 12.0, 200.0, limit=800)[0] / KD
print(f"\nTRUE (scipy)  head[0,1.5]={th_true:+.4f}  mid[1.5,12]={md_true:+.4f}  tail[12,inf]={tl_true:+.4f}"
      f"   total R/K*={th_true + md_true + tl_true:.4f}   (margin to 1 = {1 - (th_true + md_true + tl_true):.4f})")

full_mid = [(1.5, 2.4048, 8, 700), (2.4048, 4.0, 12, 1100), (4.0, 6.0, 12, 1900),
            (6.0, 8.0, 12, 3000), (8.0, 10.0, 14, 5500), (10.0, 12.0, 16, 11000)]
for mult in (1, 2):
    sl = [(lo, hi, nb, nt * mult) for (lo, hi, nb, nt) in full_mid]
    t0 = time.time()
    rig = float(cert.middle_point(a1, a3, b, sl).upper()) / KD
    print(f"RIGOROUS mid[1.5,12] nth x{mult}: {rig:+.4f}   theta-loss vs true = {rig - md_true:+.4f}  "
          f"({time.time() - t0:.0f}s)")

full_tail = [(12.0, 24.0, 20, 1200), (24.0, 60.0, 16, 2400), (60.0, 200.0, 12, 6000)]
tail = cert.tail_modulus_box(a1, a3, b, 12.0, full_tail, 200.0)
mgt = cert.mg_tail_box(min(sx2, sy2), 12.0)
rig_tail = float((tail + mgt).upper()) / KD
print(f"RIGOROUS tail[12,inf] capped-|M|: {rig_tail:+.4f}   loss vs true signed tail = {rig_tail - tl_true:+.4f}")
