# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2", "scipy>=1.11", "python-flint>=0.8"]
# ///
"""Checker and manifest for the RIGOROUS near-pair basin certificate (dagger_nearpair.py).

It re-reads numerics/dagger_nearpair.jsonl and trusts neither the driver nor the stored verdicts:
  1. every row was produced by the current certificate source (SHA-256 of numerics/dagger_nearpair.py)
     with the pinned parameters dagger_nearpair.DEFAULT;
  2. every row is RE-DECIDED from its stored Arb balls (Arb comparisons, no float): a Hessian box needs
     f_xx > 0, f_yy > 0 and lo(f_xx) lo(f_yy) - up|f_xy|^2 > 0; a Taylor box lo(f_c) - lin - quad > 0 with
     lin, quad >= 0; a monotone box K_lo P_min - E|Z|_up + E|G|_lo > 0, K_star re-enclosed by
     dagger_rigorous.kstar_enclosure;
  3. coverage, in exact rationals on the finest lattice of each family: the proved Hessian boxes lie in
     and tile U = [u_x0, u_x1] x [0, u_y1], which contains p = (1, 0) and lies in {x > y}; every lattice cell
     that may meet A lies in U or in a proved region-A box, and every cell that may meet B lies in a proved
     region-B box.  "May meet" is the necessary corner test (cells have nonnegative corners):
         A = {0 <= y <= x, 41 y^2 <= 9 (x^2 + 1), x^2 + y^2 <= 41/9}:  y0 <= x1, 41 y0^2 <= 9 (x1^2 + 1),
             x0^2 + y0^2 <= 41/9;  A lies in [0, 9/4]^2 since 41/9 < (9/4)^2;
         B = {0 <= x <= 1, 41 y^2 <= 9 (1 + x^2)}:  x0 <= 1, 41 y0^2 <= 9 (1 + x1^2);  B lies in
             [0, 1] x [0, 7/10] since 18/41 < (7/10)^2;
  4. --recompute N re-runs dagger_nearpair.certify on the N tightest rows of each kind and requires
     bit-identical balls;
  5. --crosscheck recomputes E|Z|, E|G| and f by an independent route (scipy: the odd copy through the
     complete elliptic integral E_Q sqrt(X^2 + b^2 cos^2 Q) = (2/pi) sqrt(X^2 + b^2) E(b^2/(X^2 + b^2)),
     adaptive quadrature over the even phases, finite differences for the Hessian) at the pair and at the
     tightest boxes, and requires agreement with the stored balls (not part of the proof).
--write-manifest records numerics/dagger_nearpair_manifest.json (source and data SHA-256, parameters,
versions, row counts, tightest margins, coverage); without it an existing manifest is verified.

    uv run --script numerics/verify_dagger_nearpair.py [--recompute 3] [--crosscheck] [--write-manifest]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
from dataclasses import asdict
from fractions import Fraction
from pathlib import Path

import dagger_nearpair as dn
import dagger_rigorous as dr
import flint
from flint import arb, ctx

ROOT = Path(__file__).resolve().parents[1]
NUM = ROOT / "numerics"
DATA = NUM / "dagger_nearpair.jsonl"
MANIFEST = NUM / "dagger_nearpair_manifest.json"
SCHEMA = "beta-odd-dagger-nearpair-cover-v1"
SOURCES = ("dagger_nearpair.py", "verify_dagger_nearpair.py", "dagger_rigorous.py")
FAMILY = {"hessian": "hessian", "taylor": "A", "mono_A": "A", "mono_B": "B"}
BALLS = {
    "hessian": ("f_xx", "f_xy", "f_yy", "det_lo"),
    "taylor": ("f_c", "lin", "quad", "lower"),
    "mono_A": ("p_min", "ez_up", "eg_lo", "lower"),
    "mono_B": ("p_min", "ez_up", "eg_lo", "lower"),
}

Box = tuple[Fraction, Fraction, Fraction]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _lo(x: arb) -> arb:
    return arb(x.lower())


def _up(x: arb) -> arb:
    return arb(x.upper())


def redecide(row: dict, k_lo: arb) -> tuple[bool, arb, list[str]]:
    """(proved, margin, problems) re-derived from the stored balls of one row."""
    b = {key: dr.ball_arb(row[key]) for key in BALLS[row["kind"]]}
    problems = []
    if row["kind"] == "hessian":
        margin = _lo(b["f_xx"]) * _lo(b["f_yy"]) - _up(abs(b["f_xy"])) ** 2
        return bool(b["f_xx"] > 0 and b["f_yy"] > 0 and margin > 0), margin, problems
    if row["kind"] == "taylor":
        if not (b["lin"] >= 0 and b["quad"] >= 0):
            problems.append("a Taylor remainder term is not a nonnegative ball")
        margin = _lo(b["f_c"]) - _up(b["lin"]) - _up(b["quad"])
        return bool(margin > 0) and not problems, margin, problems
    if not (b["p_min"] >= 0):
        problems.append("P_min is not a nonnegative ball")
    margin = k_lo * _lo(b["p_min"]) - _up(b["ez_up"]) + _lo(b["eg_lo"])
    return bool(margin > 0) and not problems, margin, problems


def covered_cells(
    boxes: list[Box], ox: Fraction, oy: Fraction, cell: Fraction
) -> set[tuple[int, int]]:
    out: set[tuple[int, int]] = set()
    for xc, yc, h in boxes:
        ks = [
            (xc - h - ox) / cell,
            (xc + h - ox) / cell,
            (yc - h - oy) / cell,
            (yc + h - oy) / cell,
        ]
        if any(k.denominator != 1 for k in ks):
            raise ValueError(f"box ({xc}, {yc}, {h}) is off the lattice ({ox}, {oy}) + k*{cell}")
        i0, i1, j0, j1 = (int(k) for k in ks)
        out.update((i, j) for i in range(i0, i1) for j in range(j0, j1))
    return out


def may_meet_a(x0: Fraction, x1: Fraction, y0: Fraction) -> bool:
    return y0 <= x1 and 41 * y0 * y0 <= 9 * (x1 * x1 + 1) and x0 * x0 + y0 * y0 <= Fraction(41, 9)


def may_meet_b(x0: Fraction, x1: Fraction, y0: Fraction) -> bool:
    return x0 <= 1 and 41 * y0 * y0 <= 9 * (1 + x1 * x1)


def coverage(proved: dict[str, list[Box]], prm: dn.Params) -> list[str]:
    """Exact-rational coverage of U, A and B (item 3); returns the problems found."""
    problems = []
    ux0, ux1, uy1 = Fraction(prm.u_x0), Fraction(prm.u_x1), Fraction(prm.u_y1)
    if not (ux0 <= 1 <= ux1 and uy1 < ux0):
        problems.append("U must contain p = (1, 0) and lie in {x > y}")
    hess = proved["hessian"]
    if not hess:
        return [*problems, "no proved Hessian box"]
    for xc, yc, h in hess:
        if not (xc - h >= ux0 and xc + h <= ux1 and yc - h >= 0 and yc + h <= uy1):
            problems.append(f"Hessian box ({xc}, {yc}, {h}) is not inside U")
    cell = min(2 * h for _, _, h in hess)
    if ((ux1 - ux0) / cell).denominator != 1 or (uy1 / cell).denominator != 1:
        problems.append("U is not a union of Hessian lattice cells")
        return problems
    got = covered_cells(hess, ux0, Fraction(0), cell)
    miss = [
        (i, j)
        for i in range(int((ux1 - ux0) / cell))
        for j in range(int(uy1 / cell))
        if (i, j) not in got
    ]
    problems += [
        f"U cell {ux0 + i * cell}, {j * cell} (side {cell}) not tiled" for i, j in miss[:10]
    ]

    for fam, meets, xmax, ymax in (
        ("A", may_meet_a, Fraction(9, 4), Fraction(9, 4)),
        ("B", may_meet_b, Fraction(1), Fraction(7, 10)),
    ):
        boxes = proved[fam]
        if not boxes:
            problems.append(f"no proved region-{fam} box")
            continue
        cell = min(2 * h for _, _, h in boxes)
        if fam == "A" and any((v / cell).denominator != 1 for v in (ux0, ux1, uy1)):
            problems.append("U is not a union of region-A lattice cells")
            continue
        got = covered_cells(boxes, Fraction(0), Fraction(0), cell)
        miss = []
        i = 0
        while i * cell < xmax:
            x0, x1 = i * cell, (i + 1) * cell
            j = 0
            while j * cell < ymax:
                y0, y1 = j * cell, (j + 1) * cell
                in_u = fam == "A" and x0 >= ux0 and x1 <= ux1 and y1 <= uy1
                if meets(x0, x1, y0) and not in_u and (i, j) not in got:
                    miss.append((x0, y0))
                j += 1
            i += 1
        problems += [f"region-{fam} cell ({x}, {y}) side {cell} uncovered" for x, y in miss[:10]]
    return problems


def crosscheck(rows: dict[str, list[tuple[arb, dict]]], K: float) -> list[str]:
    """Independent scipy evaluation (item 5); returns report lines, 'CROSSCHECK FAILED' on a mismatch."""
    import numpy as np
    from scipy import integrate, special

    def odd(X: float, b: float) -> float:
        s = X * X + b * b
        return 2 / np.pi * np.sqrt(s) * special.ellipe(b * b / s) if s > 0 else 0.0

    def ez(a1: float, a3: float, b: float) -> float:
        val, _ = integrate.dblquad(
            lambda p1, p3: odd(a1 * np.cos(p1) + a3 * np.cos(p3), b),
            0,
            np.pi,
            0,
            np.pi,
            epsabs=1e-13,
            epsrel=1e-13,
        )
        return val / np.pi**2

    def eg(a1: float, a3: float, b: float) -> float:
        sx, sy = (a1 * a1 + a3 * a3) / 2, b * b / 2
        val, _ = integrate.quad(
            lambda th: np.sqrt(sx * np.cos(th) ** 2 + sy * np.sin(th) ** 2),
            0,
            np.pi / 2,
            epsabs=1e-15,
        )
        return np.sqrt(np.pi / 2) * (2 / np.pi) * val

    def f(a1: float, a3: float, b: float = 1.0) -> float:
        s2 = (a1 * a1 + a3 * a3 + b * b) / 2
        return K * (a1**4 + a3**4 + b**4) / s2**1.5 - ez(a1, a3, b) + eg(a1, a3, b)

    def hess(x: float, y: float) -> tuple[float, float, float]:
        def fd(h: float) -> tuple[float, float, float]:
            c = f(x, y)
            xx = (f(x + h, y) - 2 * c + f(x - h, y)) / (h * h)
            yy = (f(x, y + h) - 2 * c + f(x, y - h)) / (h * h)
            xy = (f(x + h, y + h) - f(x + h, y - h) - f(x - h, y + h) + f(x - h, y - h)) / (
                4 * h * h
            )
            return xx, xy, yy

        a, b = fd(0.02), fd(0.01)
        return (4 * b[0] - a[0]) / 3, (4 * b[1] - a[1]) / 3, (4 * b[2] - a[2]) / 3

    def near(val: float, rec: dict, tol: float) -> bool:
        x = dr.ball_arb(rec)
        return float(x.lower()) - tol <= val <= float(x.upper()) + tol

    lines, ok = [], True
    h1 = 0.9580913986828502
    zp = ez(1.0, 0.0, 1.0)
    fxx, _, fyy = hess(1.0, 0.0)
    good = abs(zp - h1) < 1e-10 and abs(f(1.0, 0.0)) < 1e-10
    ok &= good
    lines.append(
        f"  pair: E|Z| = {zp:.13f} (h(1) = {h1}), f = {f(1.0, 0.0):.2e}, "
        f"f_xx = {fxx:.7f}, f_yy = {fyy:.7f}  {'ok' if good else 'MISMATCH'}"
    )
    for kind in ("hessian", "taylor", "mono_A", "mono_B"):
        for _, row in rows.get(kind, [])[:2]:
            xc, yc, h = (Fraction(s) for s in row["box"][:3])
            x, y = float(xc), float(yc)
            if kind == "hessian":
                hxx, hxy, hyy = hess(x, y)
                good = all(
                    near(v, row[k], 1e-6) for v, k in ((hxx, "f_xx"), (hxy, "f_xy"), (hyy, "f_yy"))
                )
                msg = f"Hessian at the centre ({hxx:.6f}, {hxy:.6f}, {hyy:.6f}) in the box balls"
            elif kind == "taylor":
                fc = f(x, y)
                good = near(fc, row["f_c"], 1e-9)
                msg = f"f(centre) = {fc:.12f} vs stored f_c {row['f_c']['mid'][:16]}"
            else:
                x1, y1 = float(xc + h), float(yc + h)
                x0, y0 = max(float(xc - h), 0.0), max(float(yc - h), 0.0)
                if kind == "mono_A":
                    zu, gl = ez(x1, y1, 1.0), eg(x0, y0, 1.0)
                else:
                    zu, gl = ez(1.0, min(x1, 1.0), y1), eg(1.0, x0, y0)
                good = near(zu, row["ez_up"], 1e-9) and near(gl, row["eg_lo"], 1e-12)
                msg = f"E|Z|(upper) = {zu:.10f}, E|G|(lower) = {gl:.12f} in the stored balls"
            ok &= good
            lines.append(f"  {kind} {row['box'][:3]}: {msg}  {'ok' if good else 'MISMATCH'}")
    if not ok:
        lines.append("CROSSCHECK FAILED")
    return lines


def main() -> int:
    ap = argparse.ArgumentParser(description="checker for the rigorous near-pair basin certificate")
    ap.add_argument("--jsonl", default=str(DATA))
    ap.add_argument("--recompute", type=int, default=0, help="re-run the N tightest rows per kind")
    ap.add_argument("--crosscheck", action="store_true", help="independent scipy cross-check")
    ap.add_argument("--write-manifest", action="store_true")
    args = ap.parse_args()
    if not __debug__:
        raise RuntimeError("run without -O")
    prm = dn.DEFAULT
    ctx.prec = prm.prec
    data = Path(args.jsonl)
    rows = [json.loads(s) for s in data.read_text(encoding="utf-8").splitlines() if s.strip()]
    code = sha256(NUM / "dagger_nearpair.py")
    params = json.loads(json.dumps(asdict(prm)))
    k_lo, k_hi = dr.kstar_enclosure()
    ctx.prec = prm.prec
    problems: list[str] = []
    versions = set()
    best: dict[tuple[str, Box], tuple[arb, dict]] = {}
    for n, row in enumerate(rows):
        tag = f"row {n} {row.get('kind')} {row.get('box')}"
        if row.get("schema") != dn.SCHEMA or row.get("source_sha256") != code:
            problems.append(f"{tag}: produced by another certificate source or schema")
            continue
        if row["params"] != params:
            problems.append(f"{tag}: parameters differ from dagger_nearpair.DEFAULT")
            continue
        kind, (xs, ys, hx, hy) = row["kind"], row["box"]
        if kind not in FAMILY or hx != hy or FAMILY[kind] != row["method"]:
            problems.append(f"{tag}: inconsistent kind, method or box")
            continue
        versions.add(row["flint"])
        proved, margin, probs = redecide(row, k_lo)
        problems += [f"{tag}: {p}" for p in probs]
        if row["proves"] != proved:
            problems.append(f"{tag}: stored verdict is not reproduced from the stored balls")
        if proved:
            key = (kind, (Fraction(xs), Fraction(ys), Fraction(hx)))
            if key not in best or margin > best[key][0]:
                best[key] = (margin, row)
    proved_boxes: dict[str, list[Box]] = {"hessian": [], "A": [], "B": []}
    by_kind: dict[str, list[tuple[arb, dict]]] = {}
    for (kind, box), (margin, row) in best.items():
        proved_boxes[FAMILY[kind]].append(box)
        by_kind.setdefault(kind, []).append((margin, row))
    for lst in by_kind.values():
        lst.sort(key=lambda mr: float(mr[0]))
    problems += coverage(proved_boxes, prm)
    if args.recompute:
        for kind, lst in sorted(by_kind.items()):
            for _, row in lst[: args.recompute]:
                xc, yc, h = row["box"][:3]
                fresh = dn.certify((row["method"], xc, yc, h))
                same = all(fresh[k] == row[k] for k in (*BALLS[kind], "kind", "proves"))
                print(f"  recompute {kind} {row['box'][:3]}: bit-identical={same}")
                if not same:
                    problems.append(
                        f"recompute of {kind} {row['box'][:3]} differs from the stored row"
                    )
    if args.crosscheck:
        lines = crosscheck(by_kind, float(k_lo.union(k_hi).mid()))
        print("\n".join(lines))
        if "CROSSCHECK FAILED" in lines:
            problems.append("scipy cross-check failed")
    tight = {
        kind: {"box": lst[0][1]["box"][:3], "margin": dn.ball(lst[0][0])}
        for kind, lst in sorted(by_kind.items())
    }
    summary = {
        "schema": SCHEMA,
        "status": "pass" if not problems else "fail",
        "data": data.name,
        "data_sha256": sha256(data),
        "rows": len(rows),
        "proved_boxes": {kind: len(lst) for kind, lst in sorted(by_kind.items())},
        "tightest": tight,
        "U": [prm.u_x0, prm.u_x1, prm.u_y1],
        "K_star": [dn.ball(k_lo), dn.ball(k_hi)],
        "params": params,
        "source_sha256": {name: sha256(NUM / name) for name in SOURCES},
        "python_flint_version_rows": sorted(versions),
        "python_flint_version": flint.__version__,
        "python_version": platform.python_version(),
        "precision_bits": prm.prec,
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
        for key in ("schema", "data_sha256", "source_sha256", "params", "proved_boxes", "tightest"):
            if man.get(key) != summary[key]:
                problems.append(f"manifest field '{key}' does not match the current files")
        summary["status"] = "pass" if not problems else "fail"
    for p in problems[:40]:
        print("PROBLEM:", p)
    print(
        json.dumps(
            {k: summary[k] for k in ("status", "rows", "proved_boxes", "tightest")}, indent=1
        )
    )
    return 0 if not problems else 1


if __name__ == "__main__":
    raise SystemExit(main())
