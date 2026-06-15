# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2", "scipy>=1.11", "mpmath>=1.3", "python-flint>=0.8"]
# ///
"""Diagnostic: per-node middle values across the 25-node amplitude stencil of a box, separating the
[1.5,12] part from a [12,20] extension, and showing what build_stencil does with the node MIDPOINTS
vs the node RADII.

It demonstrates that the extension's failure is NOT a finite-difference blow-up of the midpoints
(build_stencil on the smooth midpoints is clean) but the node ENCLOSURE RADIUS (the theta-ball loss
of the large-t middle, see diag_dagger_thetaball.py).  Reduced theta is used for speed.

    uv run --script numerics/diag_dagger_nodes.py            # default box 0.70,0.70,0.0075
    uv run --script numerics/diag_dagger_nodes.py 0.73 0.73 0.0075
"""

import importlib.util
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CERTP = ROOT / "numerics" / "verify_dagger_continuum.py"
BASE = [(1.5, 2.4048, 8, 200), (2.4048, 4.0, 12, 300), (4.0, 6.0, 12, 500),
        (6.0, 8.0, 12, 700), (8.0, 10.0, 14, 1100), (10.0, 12.0, 16, 2000)]
EXT = [(12.0, 16.0, 16, 1200), (16.0, 20.0, 16, 1000)]
_CERT = None


def _load():
    spec = importlib.util.spec_from_file_location("cert", CERTP)
    m = importlib.util.module_from_spec(spec)
    sys.modules["cert"] = m
    spec.loader.exec_module(m)
    m.ctx.prec = 200
    return m


def _init():
    global _CERT
    _CERT = _load()


def node(coord):
    a1, a3 = coord
    cert = _CERT
    b = cert.b_of(a1, a3)
    m12 = cert.middle_point(a1, a3, b, BASE)
    mx = cert.middle_point(a1, a3, b, EXT)
    return (a1, a3, b, float(m12.mid()), float(mx.mid()), float(mx.rad()))


if __name__ == "__main__":
    a1c = float(sys.argv[1]) if len(sys.argv) > 1 else 0.70
    a3c = float(sys.argv[2]) if len(sys.argv) > 2 else 0.70
    h = float(sys.argv[3]) if len(sys.argv) > 3 else 0.0075
    cert = _load()
    from flint import arb
    coords = [(a1c + p * h, a3c + q * h) for p in range(-2, 3) for q in range(-2, 3)]
    with ProcessPoolExecutor(max_workers=14, initializer=_init) as pool:
        res = list(pool.map(node, coords))
    print(f"box ({a1c},{a3c},h={h})  reduced-theta per-node middle:")
    print(f"{'a1':>7} {'a3':>7} {'b':>7} {'mid[1.5,12]':>13} {'mid[12,20]':>13} {'rad[12,20]':>11}")
    m12d, mxd = {}, {}
    for (a1, a3, b, m12, mx, rx) in res:
        flag = "  <-- HUGE RADIUS" if rx > 1e-2 else ""
        print(f"{a1:7.4f} {a3:7.4f} {b:7.4f} {m12:13.6e} {mx:13.6e} {rx:11.2e}{flag}")
        m12d[(round(a1, 9), round(a3, 9))] = (m12, 0.0)
        mxd[(round(a1, 9), round(a3, 9))] = (m12 + mx, 0.0)
    KD = cert.K_STAR * (a1c**4 + a3c**4 + cert.b_of(a1c, a3c)**4)
    for label, dat in [("[1.5,12] mids", m12d), ("[1.5,20] mids", mxd)]:
        look = lambda a1, a3: arb(*dat[(round(a1, 9), round(a3, 9))])  # noqa: E731,B023
        sm = cert.build_stencil(look, a1c, a3c, h, h, kspan=2)
        print(f"build_stencil on {label:16}: sup/K*D = {sm.sup() / KD:+.5f}  rem={sm.rem:.2e}  "
              f"(clean => the failure is the node RADIUS, not the midpoints)")
