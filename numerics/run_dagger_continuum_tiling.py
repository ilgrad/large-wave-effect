# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2", "scipy>=1.11", "mpmath>=1.3", "python-flint>=0.8"]
# ///
"""Offline runner for the singleton-continuum (dagger) box certificates.

This is intentionally a runner, not a new theorem.  `verify_dagger_continuum.py`
contains the actual Arb certificate machinery (`prove_box`).  This script adds
the missing operational layer needed for the expensive part:

* reproducible resolution presets (`demo`, `moderate`, `full`);
* deterministic box-grid generation over the singleton band;
* JSONL checkpointing after every box, so multi-hour tilings can be resumed;
* a single-box mode for targeted full-resolution jobs.

Examples:
    uv run --script numerics/run_dagger_continuum_tiling.py --dry-run --limit 20
    uv run --script numerics/run_dagger_continuum_tiling.py --preset demo --box 0.88,0.88,0.015
    uv run --script numerics/run_dagger_continuum_tiling.py --preset full --limit 16 --jobs 4
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CERT_PATH = ROOT / "numerics" / "verify_dagger_continuum.py"
DEFAULT_OUT = ROOT / "numerics" / "dagger_continuum_tiling.jsonl"


def load_certificate_module():
    spec = importlib.util.spec_from_file_location("dagger_continuum_cert", CERT_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {CERT_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    if not module.HAVE_FLINT:
        raise RuntimeError("python-flint is required for the rigorous box certificates")
    module.ctx.prec = 200
    return module


@dataclass(frozen=True)
class Resolution:
    head_res: tuple[int, int, int]
    mid_slices: list[tuple[float, float, int, int]]
    tail_segs: list[tuple[float, float, int, int]]
    tail_far: float
    kspan: int


PRESETS: dict[str, Resolution] = {
    "demo": Resolution(
        head_res=(16, 250, 6),
        mid_slices=[
            (1.5, 2.4048, 5, 120),
            (2.4048, 4.0, 6, 180),
            (4.0, 6.0, 6, 300),
            (6.0, 8.0, 6, 450),
            (8.0, 10.0, 8, 700),
            (10.0, 12.0, 8, 1300),
        ],
        tail_segs=[(12.0, 24.0, 12, 500), (24.0, 60.0, 10, 1000), (60.0, 200.0, 8, 2500)],
        tail_far=200.0,
        kspan=1,
    ),
    "moderate": Resolution(
        head_res=(30, 600, 6),
        mid_slices=[
            (1.5, 2.4048, 8, 400),
            (2.4048, 4.0, 10, 650),
            (4.0, 6.0, 10, 1100),
            (6.0, 8.0, 10, 1700),
            (8.0, 10.0, 12, 2600),
            (10.0, 12.0, 14, 5000),
        ],
        tail_segs=[(12.0, 24.0, 16, 900), (24.0, 60.0, 12, 1800), (60.0, 200.0, 10, 4500)],
        tail_far=200.0,
        kspan=2,
    ),
    "full": Resolution(
        head_res=(30, 600, 6),
        mid_slices=[
            (1.5, 2.4048, 8, 700),
            (2.4048, 4.0, 12, 1100),
            (4.0, 6.0, 12, 1900),
            (6.0, 8.0, 12, 3000),
            (8.0, 10.0, 14, 5500),
            (10.0, 12.0, 16, 11000),
        ],
        tail_segs=[(12.0, 24.0, 20, 1200), (24.0, 60.0, 16, 2400), (60.0, 200.0, 12, 6000)],
        tail_far=200.0,
        kspan=2,
    ),
}


@dataclass(frozen=True)
class Box:
    a1: float
    a3: float
    h: float

    @property
    def key(self) -> str:
        return f"{self.a1:.6f},{self.a3:.6f},{self.h:.6f}"


def parse_box(raw: str) -> Box:
    parts = [float(x.strip()) for x in raw.split(",")]
    if len(parts) != 3:
        raise argparse.ArgumentTypeError("--box must be a1,a3,h")
    return Box(*parts)


def lambda_range(cert, box: Box) -> tuple[float, float] | None:
    vals = []
    for d1 in (-box.h, 0.0, box.h):
        for d3 in (-box.h, 0.0, box.h):
            b = cert.b_of(box.a1 + d1, box.a3 + d3)
            if b is not None:
                vals.append(cert.lam_conc(box.a1 + d1, box.a3 + d3, b))
    return (min(vals), max(vals)) if vals else None


def grid_boxes(cert, h: float, step: float, skip_handoff: bool) -> list[Box]:
    boxes: list[Box] = []
    # The singleton family is symmetric under a1 <-> a3; keep a1 >= a3 to avoid duplicated work.
    n = round((1.08 - 0.55) / step) + 1
    vals = [0.55 + i * step for i in range(n)]
    for a1 in vals:
        for a3 in vals:
            if a1 + 1e-12 < a3:
                continue
            box = Box(round(a1, 6), round(a3, 6), h)
            if cert.b_of(box.a1, box.a3) is None:
                continue
            lr = lambda_range(cert, box)
            if lr is None:
                continue
            lam_lo, lam_hi = lr
            if lam_hi < cert.LAM_LO or lam_lo >= cert.LAM_HI:
                continue
            if skip_handoff and lam_hi >= cert.LAM_HI:
                continue
            boxes.append(box)
    # Prioritize low reference-risk boxes: away from the singleton/handoff boundary and with a1~a3 first.
    boxes.sort(key=lambda b: (abs(b.a1 - b.a3), lambda_range(cert, b)[1]))
    return boxes


def load_done(path: Path) -> set[str]:
    done: set[str] = set()
    if not path.exists():
        return done
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        try:
            row = json.loads(line)
        except json.JSONDecodeError:
            continue
        if row.get("status") == "ok" and "key" in row:
            done.add(row["key"])
    return done


def run_box(cert, box: Box, preset_name: str, res: Resolution) -> dict[str, object]:
    started = time.time()
    result = cert.prove_box(
        box.a1,
        box.a3,
        box.h,
        head_res=res.head_res,
        mid_slices=res.mid_slices,
        tail_segs=res.tail_segs,
        tail_far=res.tail_far,
        delta=box.h,
        kspan=res.kspan,
    )
    row = asdict(result)
    row.update(
        {
            "status": "ok",
            "key": box.key,
            "preset": preset_name,
            "seconds": time.time() - started,
        }
    )
    return row


def run_box_worker(payload: tuple[Box, str, Resolution]) -> dict[str, object]:
    box, preset_name, res = payload
    cert = load_certificate_module()
    return run_box(cert, box, preset_name, res)


def json_ready(value: object) -> object:
    """Convert NumPy scalars from certificate arithmetic into plain JSON values."""
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    if isinstance(value, dict):
        return {str(k): json_ready(v) for k, v in value.items()}
    if isinstance(value, list | tuple):
        return [json_ready(v) for v in value]
    return value


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preset", choices=sorted(PRESETS), default="demo")
    parser.add_argument("--box", action="append", type=parse_box, help="target box a1,a3,h; repeatable")
    parser.add_argument("--h", type=float, default=0.015, help="grid box half-width")
    parser.add_argument("--step", type=float, default=0.03, help="grid center spacing")
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--jobs", type=int, default=1, help="parallel worker processes")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--include-handoff", action="store_true", help="include boxes with lambda_hi >= 0.82")
    args = parser.parse_args()

    if args.jobs < 1:
        parser.error("--jobs must be >= 1")

    cert = load_certificate_module()
    res = PRESETS[args.preset]
    boxes = args.box if args.box else grid_boxes(cert, args.h, args.step, skip_handoff=not args.include_handoff)
    done = load_done(args.out)
    boxes = [box for box in boxes if box.key not in done]
    if args.limit is not None:
        boxes = boxes[: args.limit]

    print(f"preset={args.preset} boxes={len(boxes)} jobs={args.jobs} out={args.out}")
    for box in boxes[:10]:
        lam = lambda_range(cert, box)
        print(f"  {box.key} lambda=[{lam[0]:.3f},{lam[1]:.3f}]")
    if args.dry_run:
        return 0

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("a", encoding="utf-8") as fh:
        def write_row(row: dict[str, object]) -> None:
            fh.write(json.dumps(json_ready(row), sort_keys=True) + "\n")
            fh.flush()
            if row["status"] == "ok":
                print(
                    f"    upper={row['upper']:.4f} proves={row['proves']} "
                    f"lambda=[{row['lam_lo']:.3f},{row['lam_hi']:.3f}] seconds={row['seconds']:.1f}",
                    flush=True,
                )
            else:
                print(f"    ERROR {row['error']}", flush=True)

        if args.jobs > 1:
            with ProcessPoolExecutor(max_workers=args.jobs) as pool:
                futures = {pool.submit(run_box_worker, (box, args.preset, res)): box for box in boxes}
                for i, fut in enumerate(as_completed(futures), start=1):
                    box = futures[fut]
                    print(f"[{i}/{len(boxes)}] {box.key}", flush=True)
                    try:
                        row = fut.result()
                    except Exception as exc:
                        row = {
                            "status": "error",
                            "key": box.key,
                            "preset": args.preset,
                            "error": repr(exc),
                        }
                    write_row(row)
            return 0

        for i, box in enumerate(boxes, start=1):
            print(f"[{i}/{len(boxes)}] {box.key}", flush=True)
            try:
                row = run_box(cert, box, args.preset, res)
            except Exception as exc:  # keep a long tiling run resumable even if one box fails
                row = {
                    "status": "error",
                    "key": box.key,
                    "preset": args.preset,
                    "error": repr(exc),
                }
            write_row(row)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
