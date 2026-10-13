# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2", "scipy>=1.11", "python-flint>=0.8"]
# ///
"""Checker and manifest for the RIGOROUS beta_odd continuum cover (the role of a verify_all.py).

It re-reads numerics/dagger_campaign_rigorous.jsonl (written by run_dagger_campaign.py --rigorous) and
trusts neither the driver nor the stored verdicts:
  1. every row was produced by the current certificate source (SHA-256 of numerics/dagger_rigorous.py)
     with the pinned parameters dagger_rigorous.DEFAULT;
  2. every leaf is RE-DECIDED from its stored Arb balls: the certified rectangle is recomputed from the
     decimal box and must lie inside the stored balls, K_star is re-enclosed, and dagger_rigorous.decide()
     rebuilds  F_up >= sup_rect [P_K - K_lo D] + eps_gl + eps_theta + eps_amp + tail_env + tail_mg  from the
     stored Taylor coefficients and error terms; the leaf counts only if F_up < 0 and the ratio ball < 1
     (Arb comparisons, no float);
  3. the re-decided leaves cover the band lambda in [0.65, 0.82), a1 >= a3, by an exact rational interval
     test on the box-edge lattice (run_dagger_campaign.meets_band);
  4. --recompute N re-runs certify_box on the N worst leaves and requires bit-identical balls;
  5. --crosscheck compares, on the two worst leaves, the stored Taylor polynomial with an independent scipy
     evaluation of the signed integral at the box corners, and the scipy ratio on a grid with the bound.
--write-manifest then records numerics/dagger_rigorous_manifest.json (source and data SHA-256, parameters,
versions, precision, leaf count, worst bound, coverage); without it an existing manifest is verified.

    uv run --script numerics/verify_dagger_rigorous.py [--recompute 3] [--crosscheck] [--write-manifest]
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import platform
import sys
from dataclasses import asdict
from fractions import Fraction
from pathlib import Path

import flint
from flint import arb, ctx

ROOT = Path(__file__).resolve().parents[1]
NUM = ROOT / "numerics"
DATA = NUM / "dagger_campaign_rigorous.jsonl"
MANIFEST = NUM / "dagger_rigorous_manifest.json"
SCHEMA = "beta-odd-dagger-rigorous-cover-v1"
SOURCES = ("dagger_rigorous.py", "run_dagger_campaign.py", "verify_dagger_rigorous.py")


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def redecide(dr, row: dict, kstar: tuple[arb, arb]) -> tuple[bool, arb, arb, list[str]]:
    """Re-derive the verdict of one row from its stored balls; returns (proved, F_up, ratio, problems)."""
    problems = []
    geo = dr.box_geometry(row["a1"], row["a3"], row["h"])
    for key, val in (("um", geo.um), ("vm", geo.vm), ("Hu", geo.Hu), ("Hv", geo.Hv)):
        if not dr.ball_arb(row["rect"][key]).contains(val):
            problems.append(f"stored rect.{key} does not contain the recomputed value")
    coef = {}
    for key, rec in row["coef"].items():
        i, j = (int(x) for x in key.split(","))
        coef[(i, j)] = dr.ball_arb(rec)
    K = row["params"]["degree"]
    if set(coef) != {(i, j) for i in range(K + 1) for j in range(K + 1 - i)}:
        problems.append("coefficient set is not the full degree-K triangle")
        return False, arb("nan"), arb("nan"), problems
    errs = {name: dr.ball_arb(row["terms"][name]) for name in dr.ERROR_TERMS}
    if not all(e >= 0 for e in errs.values()):
        problems.append("an error term is not a nonnegative ball")
    dec = dr.decide(coef, geo, errs, kstar, row["params"]["split_depth"])
    proved = bool(dec["F_up"] < 0 and dec["ratio"] < 1)
    # Re-deriving from the 25-30 digit decimal balls lands ~1e-27 above the stored value (more
    # conservative); a larger gap means a corrupted or hand-edited row.
    if abs(float(dr.ball_arb(row["upper_ball"])) - float(dec["ratio"])) > 1e-12:
        problems.append("stored upper_ball disagrees with the re-derived ratio")
    return proved, dec["F_up"], dec["ratio"], problems


def uncovered(rc, leaves: list[tuple[Fraction, Fraction, Fraction]]) -> list[tuple[Fraction, ...]]:
    """Band cells (a1 >= a3) of the lattice ORIGIN + k*cell, cell = the smallest leaf h, lying in no leaf."""
    cell = min(h for _, _, h in leaves)
    covered = set()
    for a1, a3, h in leaves:
        edges = [a1 - h, a1 + h, a3 - h, a3 + h]
        ks = [(e - rc.ORIGIN) / cell for e in edges]
        if any(k.denominator != 1 for k in ks):
            raise ValueError(f"leaf ({a1}, {a3}, {h}) is off the lattice {rc.ORIGIN} + k*{cell}")
        k = [int(x) for x in ks]
        covered.update((i, j) for i in range(k[0], k[1]) for j in range(k[2], k[3]))
    out = []
    i = 0
    while rc.ORIGIN + i * cell < Fraction(3, 2):
        x0, x1 = rc.ORIGIN + i * cell, rc.ORIGIN + (i + 1) * cell
        for j in range(i + 1):
            y0, y1 = rc.ORIGIN + j * cell, rc.ORIGIN + (j + 1) * cell
            if rc.meets_band(x0, x1, y0, y1) and (i, j) not in covered:
                out.append((x0, y0, cell))
        i += 1
    return out


def crosscheck(dr, rows: list[dict]) -> list[str]:
    """Independent scipy check of the stored polynomial and of the bound (not part of the proof)."""
    import numpy as np
    from scipy.integrate import quad
    from scipy.special import i0e, j0

    def F(t: float, u: float, v: float) -> float:
        th = (np.arange(6000) + 0.5) / 6000 * (np.pi / 2)
        c, s = np.cos(th), np.sin(th)
        m = (
            j0(np.sqrt(u) * t * c) * j0(np.sqrt(v) * t * c) * j0(np.sqrt(2 - u - v) * t * s)
        ).mean()
        x = abs(u + v - 1) * t * t / 4
        return (np.exp(-t * t / 4 + x) * i0e(x) - m) / t**2

    def integral(u: float, v: float, a: float, b: float) -> float:
        return quad(lambda t: F(t, u, v), a, b, limit=2000, epsabs=1e-13)[0]

    k_star = float(dr.kstar_enclosure()[0])
    lines = []
    for row in rows:
        um, vm = float(dr.ball_arb(row["rect"]["um"])), float(dr.ball_arb(row["rect"]["vm"]))
        Hu, Hv = float(dr.ball_arb(row["rect"]["Hu"])), float(dr.ball_arb(row["rect"]["Hv"]))
        T = float(row["params"]["t_cut"])
        coef = {
            tuple(int(x) for x in k.split(",")): float(dr.ball_arb(v))
            for k, v in row["coef"].items()
        }
        worst = 0.0
        for su, sv in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
            du, dv = su * Hu, sv * Hv
            p = sum(c * du**i * dv**j for (i, j), c in coef.items())
            worst = max(worst, abs(p - integral(um + du, vm + dv, 1e-6, T)))
        a1c, a3c, h = (float(row[k]) for k in ("a1", "a3", "h"))
        rmax = 0.0
        for a1 in np.linspace(a1c - h, a1c + h, 3):
            for a3 in np.linspace(a3c - h, a3c + h, 3):
                u, v = a1 * a1, a3 * a3
                ex = integral(u, v, 1e-6, T) + integral(u, v, T, 400.0)
                rmax = max(rmax, ex / (k_star * (u * u + v * v + (2 - u - v) ** 2)))
        ok = worst < 1e-9 and rmax <= row["upper"]
        lines.append(
            f"  ({row['a1']},{row['a3']},{row['h']}): |P_K - scipy| at corners <= {worst:.1e}; "
            f"scipy max R on 3x3 = {rmax:.5f} <= upper {row['upper']:.5f}: {'ok' if ok else 'FAIL'}"
        )
        if not ok:
            lines.append("CROSSCHECK FAILED")
    return lines


def main() -> int:
    ap = argparse.ArgumentParser(description="checker for the rigorous beta_odd continuum cover")
    ap.add_argument("--jsonl", default=str(DATA))
    ap.add_argument(
        "--recompute", type=int, default=0, help="re-run certify_box on the N worst leaves"
    )
    ap.add_argument("--crosscheck", action="store_true", help="independent scipy cross-check")
    ap.add_argument("--write-manifest", action="store_true")
    args = ap.parse_args()
    if not __debug__:
        raise RuntimeError("run without -O")
    dr = _load("dagger_rigorous", NUM / "dagger_rigorous.py")
    rc = _load("run_dagger_campaign", NUM / "run_dagger_campaign.py")
    ctx.prec = dr.DEFAULT.prec
    data = Path(args.jsonl)
    rows = [json.loads(s) for s in data.read_text(encoding="utf-8").splitlines() if s.strip()]
    code = sha256(NUM / "dagger_rigorous.py")
    params = json.loads(json.dumps(asdict(dr.DEFAULT)))
    problems: list[str] = []
    kstar = dr.kstar_enclosure()
    ctx.prec = dr.DEFAULT.prec
    leaves: dict[tuple[Fraction, Fraction, Fraction], tuple[arb, dict]] = {}
    versions = set()
    for n, row in enumerate(rows):
        tag = f"row {n} ({row['a1']},{row['a3']},{row['h']})"
        if row["code_sha256"] != code:
            problems.append(f"{tag}: produced by another certificate source")
            continue
        if row["params"] != params:
            problems.append(f"{tag}: parameters differ from dagger_rigorous.DEFAULT")
            continue
        versions.add(row["flint"])
        proved, _F_up, ratio, probs = redecide(dr, row, kstar)
        problems += [f"{tag}: {p}" for p in probs]
        if row["proves"] and not proved:
            problems.append(
                f"{tag}: stored verdict 'proves' is not reproduced from the stored balls"
            )
        if proved:
            key = (Fraction(row["a1"]), Fraction(row["a3"]), Fraction(row["h"]))
            if key not in leaves or ratio < leaves[key][0]:
                leaves[key] = (ratio, row)
    if not leaves:
        problems.append("no certified leaf")
        miss = []
    else:
        miss = uncovered(rc, list(leaves))
    for x0, y0, c in miss[:20]:
        problems.append(
            f"uncovered band cell a1 in [{rc.dec(x0)}, {rc.dec(x0 + c)}], "
            f"a3 in [{rc.dec(y0)}, {rc.dec(y0 + c)}]"
        )
    worst_key = max(leaves, key=lambda k: float(leaves[k][0])) if leaves else None
    if args.recompute and leaves:
        order = sorted(leaves, key=lambda k: -float(leaves[k][0]))[: args.recompute]
        for key in order:
            row = leaves[key][1]
            fresh = dr.certify_box(row["a1"], row["a3"], row["h"])
            same = fresh["upper_ball"] == row["upper_ball"] and fresh["F_up"] == row["F_up"]
            print(f"  recompute ({row['a1']},{row['a3']},{row['h']}): bit-identical={same}")
            if not same:
                problems.append(f"recompute of {key} differs from the stored row")
    if args.crosscheck and leaves:
        order = sorted(leaves, key=lambda k: -float(leaves[k][0]))[:2]
        lines = crosscheck(dr, [leaves[k][1] for k in order])
        print("\n".join(lines))
        if "CROSSCHECK FAILED" in lines:
            problems.append("scipy cross-check failed")
    summary = {
        "schema": SCHEMA,
        "status": "pass" if not problems else "fail",
        "data": data.name,
        "data_sha256": sha256(data),
        "rows": len(rows),
        "leaves": len(leaves),
        "uncovered_band_cells": len(miss),
        "worst_leaf": None
        if worst_key is None
        else {
            "box": [rc.dec(x) for x in worst_key],
            "upper_ball": dr.ball(leaves[worst_key][0]),
        },
        "K_star": [dr.ball(kstar[0]), dr.ball(kstar[1])],
        "params": params,
        "source_sha256": {name: sha256(NUM / name) for name in SOURCES},
        "python_flint_version_rows": sorted(versions),
        "python_flint_version": flint.__version__,
        "python_version": platform.python_version(),
        "precision_bits": dr.DEFAULT.prec,
    }
    if args.write_manifest:
        if problems:
            print("refusing to write the manifest: the checks failed")
        else:
            MANIFEST.write_text(
                json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
            )
            print(f"wrote {MANIFEST.name}")
    elif MANIFEST.exists() and data.resolve() == DATA.resolve():
        man = json.loads(MANIFEST.read_text(encoding="utf-8"))
        for key in ("schema", "data_sha256", "source_sha256", "params", "leaves", "worst_leaf"):
            if man.get(key) != summary[key]:
                problems.append(f"manifest field '{key}' does not match the current files")
        summary["status"] = "pass" if not problems else "fail"
    for p in problems[:40]:
        print("PROBLEM:", p)
    print(
        json.dumps(
            {
                k: summary[k]
                for k in ("status", "rows", "leaves", "uncovered_band_cells", "worst_leaf")
            },
            indent=1,
        )
    )
    return 0 if not problems else 1


if __name__ == "__main__":
    raise SystemExit(main())
