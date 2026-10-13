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

The (0.70,0.70) box certifies upper=0.9915; the hardest boxes are the high-lambda edge (a3~0.61), each
closed by one h-refinement (worst certified leaf 0.9999), so the adaptive grid closes the band.  This
validates the residual >=3-copy inequality across the band at the validated standard (the amplitude
remainder uses finite-difference third derivatives, see build_stencil) -- it is not an unconditional proof
of limsup B_N/sqrt(N)=beta_odd.  Soundness is the unmodified prove_box arithmetic (middle_far validated
by numerics/verify_middle_far.py); Arb is CPU-only, so process parallelism is the accelerator.  Node
values are cached, so adjacent boxes share their common stencil nodes.  --coverage checks that the
certified boxes in the JSONL cover the band.

--rigorous replaces the per-box certificate by numerics/dagger_rigorous.py (every bound an Arb enclosure:
exact signed Kluyver integral on [0, 24] by Gauss-Legendre x periodic midpoint rule with analytic error
bounds, an exact degree-10 amplitude Taylor polynomial with a Cauchy remainder, an exact quadratic maximum
over the box, the capped-envelope tail).  Boxes are exact decimals; the base grid is every lattice box of
half-width h whose interval enclosure meets the band, so a run in which every box proves covers the band by
construction; one box costs ~1.5 s, so boxes (not nodes) are the parallel unit.  Rows go to
dagger_campaign_rigorous.jsonl with the Arb balls of every term, the parameters and the SHA-256 of the
certificate source; numerics/verify_dagger_rigorous.py re-checks them and the band coverage.

    uv run --script numerics/run_dagger_campaign.py --dry-run
    uv run --script numerics/run_dagger_campaign.py --h 0.0075 --step 0.015 --jobs 14
    uv run --script numerics/run_dagger_campaign.py --coverage
    uv run --script numerics/run_dagger_campaign.py --rigorous --jobs 8
"""

from __future__ import annotations

import argparse
import functools
import importlib.util
import json
import math
import sys
import time
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from decimal import Decimal, localcontext
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CERTP = ROOT / "numerics" / "verify_dagger_continuum.py"
OUT = ROOT / "numerics" / "dagger_campaign.jsonl"
CERT_RIG = ROOT / "numerics" / "dagger_rigorous.py"
OUT_RIG = ROOT / "numerics" / "dagger_campaign_rigorous.jsonl"
ORIGIN = Fraction("0.535")  # box edges sit on ORIGIN + k * (2h)
M_BAND = (Fraction("0.36"), Fraction("0.7"))  # lambda = 1 - min(a1^2, a3^2, b^2)/2 in [0.65, 0.82)

HEAD_RES = (30, 600, 6)
MID_SLICES = [
    (1.5, 2.4048, 8, 1400),
    (2.4048, 4.0, 12, 2200),
    (4.0, 6.0, 12, 3800),
    (6.0, 8.0, 12, 6000),
]
FAR_SLICES = [
    (8.0, 12.0, 16, 160),
    (12.0, 16.0, 12, 130),
    (16.0, 20.0, 12, 110),
    (20.0, 24.0, 12, 110),
]
TAIL_SEGS = [(24.0, 60.0, 16, 2400), (60.0, 200.0, 12, 6000)]
T_TAIL, TAIL_FAR = 24.0, 200.0
_NODES: dict[tuple[float, float], tuple] = {}


def _load():
    spec = importlib.util.spec_from_file_location("dagger_cert", CERTP)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {CERTP}")
    m = importlib.util.module_from_spec(spec)
    sys.modules["dagger_cert"] = m
    spec.loader.exec_module(m)
    if not m.HAVE_FLINT:
        raise RuntimeError("python-flint is required")
    m.ctx.prec = 200
    return m


@functools.cache
def _worker_cert():
    return _load()


def _init():
    _worker_cert()


def eval_node(coord):
    a1, a3 = coord
    cert = _worker_cert()
    b = cert.b_of(a1, a3)
    nb, nt, ordr = HEAD_RES
    hv = cert.head_majorant_point(a1, a3, b, cert.T1, nb, nt, ordr)
    mv = cert.middle_point(a1, a3, b, MID_SLICES)
    fv = cert.middle_far_point(a1, a3, b, FAR_SLICES)
    return (
        a1,
        a3,
        float(hv.mid()),
        float(hv.rad()),
        float(mv.mid()),
        float(mv.rad()),
        float(fv.mid()),
        float(fv.rad()),
    )


def _key(a1, a3):
    return round(a1, 9), round(a3, 9)


def prove(cert, pool, a1c, a3c, h):
    from flint import arb

    coords = [(a1c + p * h, a3c + q * h) for p in range(-2, 3) for q in range(-2, 3)]
    todo = [c for c in coords if _key(*c) not in _NODES]
    _NODES.update((_key(v[0], v[1]), v) for v in pool.map(eval_node, todo))
    res = [_NODES[_key(*c)] for c in coords]
    d = {
        tag: {
            (round(a1, 9), round(a3, 9)): (v[i], v[i + 1]) for v in res for a1, a3 in [(v[0], v[1])]
        }
        for tag, i in [("h", 2), ("m", 4), ("f", 6)]
    }
    look = lambda t: lambda a1, a3: arb(*d[t][(round(a1, 9), round(a3, 9))])  # noqa: E731
    hs = cert.build_stencil(look("h"), a1c, a3c, h, h, kspan=2).sup()
    ms = cert.build_stencil(look("m"), a1c, a3c, h, h, kspan=2).sup()
    fs = cert.build_stencil(look("f"), a1c, a3c, h, h, kspan=2).sup()
    a1min, a3min = a1c - h, a3c - h
    bmin = min(
        (
            bb
            for d1 in (-h, 0, h)
            for d3 in (-h, 0, h)
            if (bb := cert.b_of(a1c + d1, a3c + d3)) is not None
        ),
        default=cert.b_of(a1c - h, a3c - h),
    )
    minD, qmin = cert._box_minD_and_qmin(a1c, a3c, h)
    tail = cert.tail_modulus_box(a1min, a3min, bmin, T_TAIL, TAIL_SEGS, TAIL_FAR)
    mgt = cert.mg_tail_box(qmin, T_TAIL)
    kd = cert.K_STAR * minD
    upper = (hs + ms + fs + float((tail + mgt).upper())) / kd
    lam = [
        cert.lam_conc(a1c + d1, a3c + d3, bb)
        for d1 in (-h, 0, h)
        for d3 in (-h, 0, h)
        if (bb := cert.b_of(a1c + d1, a3c + d3)) is not None
    ]
    return {
        "a1": float(a1c),
        "a3": float(a3c),
        "h": float(h),
        "upper": float(upper),
        "proves": bool(upper < 1.0),
        "lam_lo": float(min(lam)),
        "lam_hi": float(max(lam)),
    }


def base_grid(cert, h, step):
    # lambda < LAM_HI needs every square > 2 (1 - LAM_HI), so the band reaches this a1
    a1_max = math.sqrt(2 - 4 * (1 - cert.LAM_HI))
    n = math.ceil((a1_max - h - 0.55) / step) + 1
    vals = [round(0.55 + i * step, 6) for i in range(n)]
    out = []
    for a1 in vals:
        for a3 in vals:
            if a1 + 1e-12 < a3 or cert.b_of(a1, a3) is None:
                continue
            lr = [
                cert.lam_conc(a1 + d1, a3 + d3, bb)
                for d1 in (-h, 0, h)
                for d3 in (-h, 0, h)
                if (bb := cert.b_of(a1 + d1, a3 + d3)) is not None
            ]
            if not lr or max(lr) < cert.LAM_LO or min(lr) >= cert.LAM_HI:
                continue
            out.append((a1, a3, h))
    out.sort(key=lambda b: abs(b[0] - b[1]))
    return out


def uncovered_cells(cert, rows, origin=0.535):
    """Cells of the box-edge lattice that meet the band (a1 >= a3, lambda in [LAM_LO, LAM_HI)) but lie in
    no certified box.  Every box edge sits on the lattice origin + k*cell, cell = the smallest certified h, so
    the certified union is a union of cells; each cell is tested against the band by an interval enclosure
    of the smallest squared amplitude (lambda = 1 - min a^2 / 2), which can only over-report."""
    boxes = [(r["a1"], r["a3"], r["h"]) for r in rows if r["proves"]]
    cell = min(h for _, _, h in boxes)
    covered = set()
    for a1, a3, h in boxes:
        edges = [a1 - h, a1 + h, a3 - h, a3 + h]
        k = [round((e - origin) / cell) for e in edges]
        if any(abs(origin + ki * cell - e) > 1e-9 for ki, e in zip(k, edges, strict=True)):
            raise ValueError(f"box ({a1}, {a3}, h={h}) is off the lattice {origin} + k*{cell}")
        covered.update((i, j) for i in range(k[0], k[1]) for j in range(k[2], k[3]))
    m_lo_band, m_hi_band = 2 * (1 - cert.LAM_HI), 2 * (1 - cert.LAM_LO)
    out = []
    for i in range(math.ceil((math.sqrt(2) - origin) / cell)):
        x0, x1 = origin + i * cell, origin + (i + 1) * cell
        for j in range(i + 1):
            y0, y1 = origin + j * cell, origin + (j + 1) * cell
            m_lo = min(x0 * x0, y0 * y0, 2 - x1 * x1 - y1 * y1)
            m_hi = min(x1 * x1, y1 * y1, 2 - x0 * x0 - y0 * y0)
            if m_hi > m_lo_band - 1e-12 and m_lo <= m_hi_band + 1e-12 and (i, j) not in covered:
                out.append((x0, y0, cell))
    return out


def _load_rigorous():
    spec = importlib.util.spec_from_file_location("dagger_rigorous", CERT_RIG)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {CERT_RIG}")
    m = importlib.util.module_from_spec(spec)
    sys.modules["dagger_rigorous"] = m
    spec.loader.exec_module(m)
    return m


@functools.cache
def _worker_rigorous():
    m = _load_rigorous()
    m.kstar_enclosure()
    return m


def _init_rigorous():
    _worker_rigorous()


def certify_rigorous(task):
    (a1, a3, h), depth = task
    rec = _worker_rigorous().certify_box(a1, a3, h)
    rec["depth"] = depth
    return rec


def dec(x: Fraction) -> str:
    """Exact decimal string of a Fraction with denominator 2^i 5^j."""
    d = x.denominator
    for p in (2, 5):
        while d % p == 0:
            d //= p
    if d != 1:
        raise ValueError(f"{x} has no finite decimal expansion")
    with localcontext() as c:
        c.prec = 80
        s = format((Decimal(x.numerator) / Decimal(x.denominator)).normalize(), "f")
    if Fraction(s) != x:
        raise ArithmeticError(f"decimal conversion of {x} is inexact")
    return s


def meets_band(x0: Fraction, x1: Fraction, y0: Fraction, y1: Fraction) -> bool:
    """Exact interval test: can min(a1^2, a3^2, b^2) on [x0,x1]x[y0,y1] fall in the band (0.36, 0.7]?"""
    m_lo = min(x0 * x0, y0 * y0, 2 - x1 * x1 - y1 * y1)
    m_hi = min(x1 * x1, y1 * y1, 2 - x0 * x0 - y0 * y0)
    return m_hi > M_BAND[0] and m_lo <= M_BAND[1]


def band_boxes(h: Fraction) -> list[tuple[str, str, str]]:
    """Every box [a1c +- h] x [a3c +- h] of the lattice ORIGIN + 2hk with a1c >= a3c, inside the disc,
    whose interval enclosure meets the band (exact rational arithmetic)."""
    out = []
    k = 0
    while ORIGIN + 2 * h * k < Fraction(3, 2):
        x0 = ORIGIN + 2 * h * k
        for j in range(k + 1):
            y0 = ORIGIN + 2 * h * j
            x1, y1 = x0 + 2 * h, y0 + 2 * h
            if x1 * x1 + y1 * y1 < 2 and meets_band(x0, x1, y0, y1):
                out.append((dec(x0 + h), dec(y0 + h), dec(h)))
        k += 1
    out.sort(key=lambda b: (abs(Fraction(b[0]) - Fraction(b[1])), b))
    return out


def run_rigorous(args) -> int:
    out = Path(args.out) if args.out else OUT_RIG
    if args.box:
        base = [tuple(dec(Fraction(x.strip())) for x in spec.split(",")) for spec in args.box]
    else:
        base = band_boxes(Fraction(str(args.h)))
    tasks = [(b, 0) for b in base]
    if args.limit:
        tasks = tasks[: args.limit]
    print(
        f"rigorous: base boxes={len(tasks)} h={args.h} jobs={args.jobs} max_depth={args.max_depth}"
    )
    print(f"  output {out}")
    if args.dry_run:
        for b, _d in tasks[:12]:
            print(f"  {b}")
        return 0
    ok = fail = 0
    t0 = time.time()
    with (
        out.open("a", encoding="utf-8") as fh,
        ProcessPoolExecutor(max_workers=args.jobs, initializer=_init_rigorous) as pool,
    ):
        pending = {pool.submit(certify_rigorous, t): t for t in tasks}
        while pending:
            done, _ = wait(pending, return_when=FIRST_COMPLETED)
            for fut in done:
                (a1, a3, h), depth = pending.pop(fut)
                rec = fut.result()
                fh.write(json.dumps(rec, sort_keys=True) + "\n")
                fh.flush()
                if rec["proves"]:
                    ok += 1
                elif depth < args.max_depth:
                    hh = Fraction(h) / 2
                    for s1 in (-1, 1):
                        for s3 in (-1, 1):
                            child = (
                                dec(Fraction(a1) + s1 * hh),
                                dec(Fraction(a3) + s3 * hh),
                                dec(hh),
                            )
                            pending[pool.submit(certify_rigorous, (child, depth + 1))] = (
                                child,
                                depth + 1,
                            )
                    print(f"  SPLIT ({a1},{a3},{h}) upper={rec['upper']:.5f} -> depth {depth + 1}")
                else:
                    fail += 1
                    print(f"  FAIL ({a1},{a3},{h}) upper={rec['upper']:.5f} at max depth")
                print(
                    f"  [{ok} ok, {fail} fail, {len(pending)} pending] {time.time() - t0:.0f}s "
                    f"({a1},{a3},h={h}) upper={rec['upper']:.5f} ({rec['seconds']['total']:.1f}s)",
                    flush=True,
                )
    print(
        f"\nDONE: {ok} certified, {fail} failed at max depth [{time.time() - t0:.0f}s]  ->  {out}"
    )
    return 0 if fail == 0 else 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h", type=float, default=0.015)
    ap.add_argument("--step", type=float, default=0.03)
    ap.add_argument("--jobs", type=int, default=14)
    ap.add_argument("--max-depth", type=int, default=2)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument(
        "--box", action="append", help="target box a1,a3,h (repeatable); replaces the base grid"
    )
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument(
        "--coverage", action="store_true", help="check that the certified boxes cover the band"
    )
    ap.add_argument(
        "--rigorous",
        action="store_true",
        help="use the fully rigorous per-box certificate numerics/dagger_rigorous.py",
    )
    ap.add_argument("--out", default=None, help="JSONL output (rigorous mode)")
    args = ap.parse_args()
    if args.rigorous:
        if args.coverage:
            print("rigorous rows: run numerics/verify_dagger_rigorous.py (exact coverage + checks)")
            return 2
        return run_rigorous(args)
    cert = _load()
    if args.coverage:
        rows = [json.loads(s) for s in OUT.read_text(encoding="utf-8").splitlines() if s.strip()]
        miss = uncovered_cells(cert, rows)
        print(
            f"{len(miss)} band cells outside the {sum(r['proves'] for r in rows)} certified boxes of {OUT.name}"
        )
        for x0, y0, c in miss[:20]:
            print(f"  a1 in [{x0:.5f}, {x0 + c:.5f}], a3 in [{y0:.5f}, {y0 + c:.5f}]")
        return 0 if not miss else 1
    if args.box:
        base = [tuple(float(x) for x in spec.split(",")) for spec in args.box]
    else:
        base = base_grid(cert, args.h, args.step)
    queue = [(b, 0) for b in base]
    if args.limit:
        queue = queue[: args.limit]
    print(
        f"base boxes={len(queue)} h={args.h} step={args.step} jobs={args.jobs} max_depth={args.max_depth}"
    )
    if args.dry_run:
        for b, _d in queue[:12]:
            print(f"  {b}")
        return 0
    ok = fail = 0
    t0 = time.time()
    with (
        OUT.open("a", encoding="utf-8") as fh,
        ProcessPoolExecutor(max_workers=args.jobs, initializer=_init) as pool,
    ):
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
                queue[:0] = [
                    ((round(a1 + s1 * hh, 7), round(a3 + s3 * hh, 7), hh), depth + 1)
                    for s1 in (-1, 1)
                    for s3 in (-1, 1)
                ]
                print(
                    f"  SPLIT ({a1},{a3},{h}) upper={row['upper']:.4f} -> depth {depth + 1}",
                    flush=True,
                )
            else:
                fail += 1
                print(f"  FAIL ({a1},{a3},{h}) upper={row['upper']:.4f} at max depth", flush=True)
            print(
                f"  [{ok} ok, {fail} fail, {len(queue)} queued] {time.time() - t0:.0f}s "
                f"last=({a1:.3f},{a3:.3f},h={h}) upper={row['upper']:.4f}",
                flush=True,
            )
    print(
        f"\nDONE: {ok} certified, {fail} failed at max depth [{time.time() - t0:.0f}s]  ->  {OUT}"
    )
    return 0 if fail == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
