# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2", "scipy>=1.11", "mpmath>=1.3", "python-flint>=0.8"]
# ///
"""Adaptive, all-cores tiling campaign for the beta_odd continuum cover (open problem A2), using the
theta-Taylor middle_far so the equal-amplitude singleton residual near lambda~0.73 actually closes.

Per box the certificate is the split signed middle  mid[1.5,8] (theta-ball middle_point) +
far[8,24] (theta-Taylor middle_far_point) + capped-|M| tail from t=24 -- exactly prove_box with
far_slices, but the 25 stencil-node evaluations run in a process pool (node-level parallelism, all
cores).  Boxes that do not certify upper<1 are split into 4 (h/2) and retried to a depth cap; every
result is checkpointed to JSONL so the multi-hour sweep is resumable.

The worst box (0.70,0.70,h=0.005) certifies upper=0.948 (validated); most boxes are easier, so the
adaptive grid closes the band.  Soundness is the unmodified prove_box arithmetic (middle_far validated
by numerics/verify_middle_far.py); Arb is CPU-only, so process parallelism is the accelerator.

    uv run --script numerics/run_dagger_campaign.py --dry-run
    uv run --script numerics/run_dagger_campaign.py --h 0.0075 --step 0.015 --jobs 14
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CERTP = ROOT / "numerics" / "verify_dagger_continuum.py"
OUT = ROOT / "numerics" / "dagger_campaign.jsonl"

HEAD_RES = (30, 600, 6)
MID_SLICES = [(1.5, 2.4048, 8, 1400), (2.4048, 4.0, 12, 2200), (4.0, 6.0, 12, 3800), (6.0, 8.0, 12, 6000)]
FAR_SLICES = [(8.0, 12.0, 16, 160), (12.0, 16.0, 12, 130), (16.0, 20.0, 12, 110), (20.0, 24.0, 12, 110)]
TAIL_SEGS = [(24.0, 60.0, 16, 2400), (60.0, 200.0, 12, 6000)]
T_TAIL, TAIL_FAR = 24.0, 200.0
_CERT = None


def _load():
    spec = importlib.util.spec_from_file_location("dagger_cert", CERTP)
    m = importlib.util.module_from_spec(spec)
    sys.modules["dagger_cert"] = m
    spec.loader.exec_module(m)
    if not m.HAVE_FLINT:
        raise RuntimeError("python-flint is required")
    m.ctx.prec = 200
    return m


def _init():
    global _CERT
    _CERT = _load()


def eval_node(coord):
    a1, a3 = coord
    cert = _CERT
    b = cert.b_of(a1, a3)
    nb, nt, ordr = HEAD_RES
    hv = cert.head_majorant_point(a1, a3, b, cert.T1, nb, nt, ordr)
    mv = cert.middle_point(a1, a3, b, MID_SLICES)
    fv = cert.middle_far_point(a1, a3, b, FAR_SLICES)
    return (a1, a3, float(hv.mid()), float(hv.rad()), float(mv.mid()), float(mv.rad()),
            float(fv.mid()), float(fv.rad()))


def prove(cert, pool, a1c, a3c, h):
    from flint import arb
    coords = [(a1c + p * h, a3c + q * h) for p in range(-2, 3) for q in range(-2, 3)]
    res = list(pool.map(eval_node, coords))
    d = {tag: {(round(a1, 9), round(a3, 9)): (v[i], v[i + 1]) for v in res for a1, a3 in [(v[0], v[1])]}
         for tag, i in [("h", 2), ("m", 4), ("f", 6)]}
    look = lambda t: (lambda a1, a3: arb(*d[t][(round(a1, 9), round(a3, 9))]))  # noqa: E731
    hs = cert.build_stencil(look("h"), a1c, a3c, h, h, kspan=2).sup()
    ms = cert.build_stencil(look("m"), a1c, a3c, h, h, kspan=2).sup()
    fs = cert.build_stencil(look("f"), a1c, a3c, h, h, kspan=2).sup()
    a1min, a3min = a1c - h, a3c - h
    bmin = min((bb for d1 in (-h, 0, h) for d3 in (-h, 0, h)
                if (bb := cert.b_of(a1c + d1, a3c + d3)) is not None), default=cert.b_of(a1c - h, a3c - h))
    minD, qmin = cert._box_minD_and_qmin(a1c, a3c, h)
    tail = cert.tail_modulus_box(a1min, a3min, bmin, T_TAIL, TAIL_SEGS, TAIL_FAR)
    mgt = cert.mg_tail_box(qmin, T_TAIL)
    kd = cert.K_STAR * minD
    upper = (hs + ms + fs + float((tail + mgt).upper())) / kd
    lam = [cert.lam_conc(a1c + d1, a3c + d3, bb) for d1 in (-h, 0, h) for d3 in (-h, 0, h)
           if (bb := cert.b_of(a1c + d1, a3c + d3)) is not None]
    return {"a1": float(a1c), "a3": float(a3c), "h": float(h), "upper": float(upper),
            "proves": bool(upper < 1.0), "lam_lo": float(min(lam)), "lam_hi": float(max(lam))}


def base_grid(cert, h, step):
    n = round((1.08 - 0.55) / step) + 1
    vals = [round(0.55 + i * step, 6) for i in range(n)]
    out = []
    for a1 in vals:
        for a3 in vals:
            if a1 + 1e-12 < a3 or cert.b_of(a1, a3) is None:
                continue
            lr = [cert.lam_conc(a1 + d1, a3 + d3, bb) for d1 in (-h, 0, h) for d3 in (-h, 0, h)
                  if (bb := cert.b_of(a1 + d1, a3 + d3)) is not None]
            if not lr or max(lr) < cert.LAM_LO or min(lr) >= cert.LAM_HI:
                continue
            out.append((a1, a3, h))
    out.sort(key=lambda b: abs(b[0] - b[1]))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h", type=float, default=0.015)
    ap.add_argument("--step", type=float, default=0.03)
    ap.add_argument("--jobs", type=int, default=14)
    ap.add_argument("--max-depth", type=int, default=2)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--box", action="append", help="target box a1,a3,h (repeatable); skips the grid")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    cert = _load()
    if args.box:
        boxes = [tuple(float(x) for x in spec.split(",")) for spec in args.box]
        with ProcessPoolExecutor(max_workers=args.jobs, initializer=_init) as pool:
            for (a1, a3, h) in boxes:
                t0 = time.time()
                row = prove(cert, pool, a1, a3, h)
                print(f"box ({a1},{a3},h={h}) UPPER={row['upper']:.5f} PROVES={row['proves']} "
                      f"lambda=[{row['lam_lo']:.4f},{row['lam_hi']:.4f}] ({time.time() - t0:.0f}s)")
        return 0
    queue = [(b, 0) for b in base_grid(cert, args.h, args.step)]
    if args.limit:
        queue = queue[: args.limit]
    print(f"base boxes={len(queue)} h={args.h} step={args.step} jobs={args.jobs} max_depth={args.max_depth}")
    if args.dry_run:
        for (b, _d) in queue[:12]:
            print(f"  {b}")
        return 0
    ok = fail = 0
    t0 = time.time()
    with OUT.open("a", encoding="utf-8") as fh, ProcessPoolExecutor(max_workers=args.jobs, initializer=_init) as pool:
        while queue:
            (a1, a3, h), depth = queue.pop(0)
            row = prove(cert, pool, a1, a3, h)
            row["depth"] = depth
            fh.write(json.dumps(row, sort_keys=True) + "\n")
            fh.flush()
            if row["proves"]:
                ok += 1
            elif depth < args.max_depth:
                hh = h / 2
                queue[:0] = [((round(a1 + s1 * hh, 7), round(a3 + s3 * hh, 7), hh), depth + 1)
                             for s1 in (-1, 1) for s3 in (-1, 1)]
                print(f"  SPLIT ({a1},{a3},{h}) upper={row['upper']:.4f} -> depth {depth + 1}", flush=True)
            else:
                fail += 1
                print(f"  FAIL ({a1},{a3},{h}) upper={row['upper']:.4f} at max depth", flush=True)
            print(f"  [{ok} ok, {fail} fail, {len(queue)} queued] {time.time() - t0:.0f}s "
                  f"last=({a1:.3f},{a3:.3f},h={h}) upper={row['upper']:.4f}", flush=True)
    print(f"\nDONE: {ok} certified, {fail} failed at max depth [{time.time() - t0:.0f}s]  ->  {OUT}")
    return 0 if fail == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
