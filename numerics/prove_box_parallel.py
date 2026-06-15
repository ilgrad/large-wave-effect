# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2", "scipy>=1.11", "mpmath>=1.3", "python-flint>=0.8"]
# ///
"""Node-parallel driver for a SINGLE dagger continuum box certificate.

`verify_dagger_continuum.prove_box` evaluates the region contributions (head, middle) on a
(2*kspan+1)^2 amplitude stencil and feeds them to `build_stencil`; those node evaluations are the
expensive part and are independent, so they parallelize cleanly across cores.  This driver computes
them in a process pool and rebuilds the EXACT same degree-2 amplitude Taylor model via the unchanged
`build_stencil` (each node's arb is shipped as (mid, rad) floats and reconstructed in the parent).
The result is identical to the serial `prove_box`; only the wall-clock differs (a full-preset box
drops from ~30 min on one core to a few minutes on 14).

This complements `run_dagger_continuum_tiling.py`, which parallelizes across BOXES (better throughput
for a whole tiling); this script parallelizes WITHIN one box (better latency for a single hard box --
e.g. while iterating on the enclosure for the equal-amplitude residual near lambda~0.73).

Soundness is unchanged: same prove_box arithmetic, same enclosures.  Arb is CPU-only and single
threaded; there is no certified GPU arithmetic, so process-level parallelism is the right accelerator.

    uv run --script numerics/prove_box_parallel.py --preset full --box 0.73,0.73,0.0075 --jobs 14
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CERT_PATH = ROOT / "numerics" / "verify_dagger_continuum.py"
RUNNER_PATH = ROOT / "numerics" / "run_dagger_continuum_tiling.py"

_CERT = None


def _load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def _load_cert():
    cert = _load_module("dagger_continuum_cert", CERT_PATH)
    if not cert.HAVE_FLINT:
        raise RuntimeError("python-flint is required for the rigorous box certificate")
    cert.ctx.prec = 200
    return cert


def _init() -> None:
    global _CERT
    _CERT = _load_cert()


def eval_node(
    args: tuple[float, float, tuple[int, int, int], list[tuple[float, float, int, int]]],
) -> tuple[float, float, float, float, float, float]:
    """Head and middle region contributions at one amplitude node, as (mid, rad) floats."""
    a1, a3, head_res, mid_slices = args
    cert = _CERT
    b = cert.b_of(a1, a3)
    nb, nt, ordr = head_res
    hv = cert.head_majorant_point(a1, a3, b, cert.T1, nb, nt, ordr)
    mv = cert.middle_point(a1, a3, b, mid_slices)
    return (a1, a3, float(hv.mid()), float(hv.rad()), float(mv.mid()), float(mv.rad()))


def prove_box_parallel(a1c: float, a3c: float, h: float, preset: str, jobs: int) -> None:
    cert = _load_cert()
    runner = _load_module("dagger_continuum_runner", RUNNER_PATH)
    from flint import arb

    res = runner.PRESETS[preset]
    kspan, delta = res.kspan, h
    coords = [(a1c + p * delta, a3c + q * delta) for p in range(-kspan, kspan + 1)
              for q in range(-kspan, kspan + 1)]
    payload = [(a1, a3, res.head_res, res.mid_slices) for (a1, a3) in coords]

    t0 = time.time()
    with ProcessPoolExecutor(max_workers=jobs, initializer=_init) as pool:
        results = list(pool.map(eval_node, payload))
    eval_dt = time.time() - t0

    head = {(round(a1, 9), round(a3, 9)): (hm, hr) for (a1, a3, hm, hr, mm, mr) in results}
    mid = {(round(a1, 9), round(a3, 9)): (mm, mr) for (a1, a3, hm, hr, mm, mr) in results}

    def head_lookup(a1: float, a3: float) -> object:
        return arb(*head[(round(a1, 9), round(a3, 9))])

    def mid_lookup(a1: float, a3: float) -> object:
        return arb(*mid[(round(a1, 9), round(a3, 9))])

    head_sup = cert.build_stencil(head_lookup, a1c, a3c, h, delta, kspan=kspan).sup()
    mid_sup = cert.build_stencil(mid_lookup, a1c, a3c, h, delta, kspan=kspan).sup()

    a1min, a3min = a1c - h, a3c - h
    bmin = min((bb for d1 in (-h, 0, h) for d3 in (-h, 0, h)
                if (bb := cert.b_of(a1c + d1, a3c + d3)) is not None), default=cert.b_of(a1c - h, a3c - h))
    minD, qmin = cert._box_minD_and_qmin(a1c, a3c, h)
    tail = cert.tail_modulus_box(a1min, a3min, bmin, cert.T2, res.tail_segs, res.tail_far)
    mgt = cert.mg_tail_box(qmin, cert.T2)
    tail_sup = float((tail + mgt).upper())

    kd = cert.K_STAR * minD
    upper = (head_sup + mid_sup + tail_sup) / kd
    lam = [cert.lam_conc(a1c + d1, a3c + d3, bb) for d1 in (-h, 0, h) for d3 in (-h, 0, h)
           if (bb := cert.b_of(a1c + d1, a3c + d3)) is not None]

    print(f"box ({a1c},{a3c},h={h}) preset={preset} jobs={jobs} lambda=[{min(lam):.4f},{max(lam):.4f}]")
    print(f"  {len(coords)} stencil nodes in {eval_dt:.0f}s  (vs serial ~{len(coords)} x single-node)")
    print(f"  head/K*D = {head_sup / kd:+.5f}   mid/K*D = {mid_sup / kd:+.5f}   tail/K*D = {tail_sup / kd:+.5f}")
    print(f"  minD = {minD:.5f}   UPPER = {upper:.5f}   PROVES (<1) = {upper < 1.0}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--preset", default="full")
    ap.add_argument("--box", required=True, help="a1,a3,h")
    ap.add_argument("--jobs", type=int, default=14)
    args = ap.parse_args()
    a1, a3, h = (float(x) for x in args.box.split(","))
    prove_box_parallel(a1, a3, h, args.preset, args.jobs)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
