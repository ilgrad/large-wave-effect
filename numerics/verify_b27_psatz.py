# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2", "python-flint>=0.8"]
# ///
r"""Exact check of the certificate for B_27 < sqrt 27 (paper, Theorem bn-saturation (4)) written by
numerics/b27_psatz_certificate.py.  Nothing is taken from the generator but the integers in the certificate.

 (a) c(r) are the coordinates of 2 cos(2 pi r/27) in the basis 1, mu_1, ..., mu_8 (floating point here, exactly
     in exact/pari/b27_psatz_check.gp), and R_l(X) = sum_k X^{c(k+l) - c(k)} is built from them;
 (b) Q is expanded from one value per G-orbit of index pairs of B_3, the monomials of word length <= 3;
 (c) sigma = m^* Q m, sum_l g_l R_l and the residual r = sigma + sum_l g_l R_l + c are computed in exact integer
     arithmetic, coefficient by coefficient;
 (d) lambda_min(Q) >= lambda: for a floating Cholesky factor L of Q - tau I (tau = lambda + 1/2), rounded to
     L~ = round(2^30 L), the matrix E = Q - tau I - 2^-60 L~ L~^T is computed exactly, so Q >= (tau - ||E||_F) I
     since L~ L~^T >= 0.  The product L~ L~^T is formed twice, by FLINT and by exact float splitting in numpy;
 (e) lambda |B_3| + c > ||r||_1, which rules out a common zero of the R_l on T^9 (see the generator);
 (f) a perturbed Gram value, a perturbed multiplier and an inflated lambda must each be rejected.

    uv run --script numerics/verify_b27_psatz.py [CERT]          (a few minutes)
"""

from __future__ import annotations

import argparse
import gzip
import sys
import time
from fractions import Fraction
from math import isqrt
from pathlib import Path

import flint
import numpy as np

N = 27
DEFAULT_CERT = Path(__file__).resolve().parent.parent / "exact" / "pari" / "b27_psatz_cert.gp.gz"
T0 = time.time()


def log(msg: str) -> None:
    print(f"{msg}  ({time.time() - T0:.1f}s)", flush=True)


def mu_coordinates() -> np.ndarray:
    c = np.zeros((N, 9), dtype=np.int64)
    c[0, 0] = 2
    for r in range(1, 9):
        c[r, r] = 1
    c[9, 0] = -1
    for i in range(1, 5):
        c[9 + i, i] = c[9 + i, 9 - i] = -1
    for r in range(14, N):
        c[r] = c[N - r]
    mu = 2 * np.cos(2 * np.pi * np.arange(N) / N)
    assert np.abs(c @ np.concatenate(([1.0], mu[1:9])) - mu).max() < 1e-12
    return c


C = mu_coordinates()


def group() -> list[np.ndarray]:
    """X^e -> X^{M e}: the multipliers 2^a (M c(r) = c(2^a r)), each with and without inversion."""
    m2 = np.zeros((9, 9), dtype=np.int64)
    m2[0, 0] = 1
    for b in range(1, 9):
        m2[:, b] = C[(2 * b) % N]
    assert all(np.array_equal(m2 @ C[r], C[(2 * r) % N]) for r in range(N))
    out, m = [], np.eye(9, dtype=np.int64)
    for _ in range(9):
        out += [m.copy(), -m]
        m = m2 @ m
    return out


G = group()


def fold(ell: int) -> int:
    ell %= N
    return min(ell, N - ell)


def autocorrelations() -> dict[int, tuple[np.ndarray, np.ndarray]]:
    out = {}
    for ell in range(1, 14):
        terms: dict[tuple[int, ...], int] = {}
        for k in range(N):
            e = tuple(int(v) for v in C[(k + ell) % N] - C[k])
            terms[e] = terms.get(e, 0) + 1
        out[ell] = (
            np.array(list(terms), dtype=np.int64),
            np.array(list(terms.values()), dtype=np.int64),
        )
    return out


def b3() -> np.ndarray:
    gens = [C[r] - C[0] for r in range(1, 14)]
    gens += [-g for g in gens]
    seen = {(0,) * 9}
    layer = [np.zeros(9, dtype=np.int64)]
    for _ in range(3):
        nxt = []
        for e in layer:
            for g in gens:
                f = e + g
                if tuple(f) not in seen:
                    seen.add(tuple(f))
                    nxt.append(f)
        layer = nxt
    return np.array(sorted(seen), dtype=np.int64)


def code(a: np.ndarray) -> np.ndarray:
    """Injective integer key of exponent vectors with entries in (-32, 32)."""
    assert np.abs(a).max() < 32
    out = np.zeros(a.shape[:-1], dtype=np.int64)
    for i in range(9):
        out = out * 64 + (a[..., i] + 32)
    return out


def read_certificate(path: Path) -> dict:
    header, values, rows = None, [], []
    with gzip.open(path, "rt", encoding="ascii") as f:
        for line in f:
            if line.startswith("\\\\") or not line.strip():
                continue
            v = [int(t) for t in line.strip().strip("[]").split(",")]
            if header is None:
                header = v
            elif v[0] == 1:
                values += v[1:]
            elif v[0] == 2:
                rows.append(v[1:])
            else:
                raise ValueError(f"unknown record {v[:3]}")
    assert header is not None
    n, degree, scale, c, lam, n_values, n_rows = header
    assert (n, degree) == (N, 3) and len(values) == n_values and len(rows) == n_rows
    return dict(
        scale=scale,
        c=c,
        lam=Fraction(lam),
        values=np.array(values, dtype=np.int64),
        rows=np.array(rows, dtype=np.int64),
    )


def expand_gram(B3: np.ndarray, values: np.ndarray) -> np.ndarray:
    n = len(B3)
    order = np.argsort(code(B3))
    perms = []
    for m in G:
        img = code(B3 @ m.T)
        pos = order[np.searchsorted(code(B3)[order], img)]
        assert np.array_equal(code(B3)[pos], img), "B_3 is not G-invariant"
        perms.append(pos)
    ii, jj = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
    key = np.full((n, n), n * n, dtype=np.int64)
    for p in perms:
        a, b = p[ii], p[jj]
        key = np.minimum(key, np.minimum(a, b) * n + np.maximum(a, b))
    keys = np.unique(key)
    assert len(keys) == len(values), (len(keys), len(values))
    return values[np.searchsorted(keys, key)]


def residual(B3: np.ndarray, QI: np.ndarray, rows: np.ndarray, c: int, RL: dict) -> np.ndarray:
    """Exact integer coefficients of sigma + sum_l g_l R_l + c (all scaled by 2^scale)."""
    # code is affine in the exponents, so code(beta - alpha) = code(beta) - code(alpha) + code(0)
    assert np.abs(B3).max() < 16
    cb = code(B3)
    diff = (cb[None, :] - cb[:, None] + code(np.zeros(9, dtype=np.int64))).reshape(-1)
    keys_sigma = [diff]
    vals_sigma = [QI.reshape(-1)]
    keys_g, vals_g = [], []
    for row in rows:
        j, e, y = int(row[0]), row[1:10], int(row[10])
        for a, m in enumerate(G):
            T, cf = RL[fold(2 ** (a // 2) * j)]
            keys_g.append(code((m @ e)[None, :] + T))
            vals_g.append(cf * y)
    keys = np.concatenate(keys_sigma + keys_g + [code(np.zeros((1, 9), dtype=np.int64))])
    vals = np.concatenate(vals_sigma + vals_g + [np.array([c], dtype=np.int64)])
    uniq, inv = np.unique(keys, return_inverse=True)
    out = np.zeros(len(uniq), dtype=np.int64)
    np.add.at(out, inv, vals)
    return out


def gram_product_flint(Lt: np.ndarray) -> list:
    M = flint.fmpz_mat(Lt.shape[0], Lt.shape[1], Lt.reshape(-1).tolist())
    return (M * M.transpose()).entries()


def gram_product_split(Lt: np.ndarray) -> np.ndarray:
    """L~ L~^T exactly: 13-bit limbs, each limb product exact in float64 (sums < 2^53)."""
    assert np.abs(Lt).max() < 2**38
    lo = Lt % 2**13
    mid = (Lt // 2**13) % 2**13
    hi = Lt // 2**26
    limbs = [(lo, 0), (mid, 13), (hi, 26)]
    total = np.zeros(Lt.shape[0] ** 2, dtype=object)
    for x, sx in limbs:
        for z, sz in limbs:
            P = x.astype(np.float64) @ z.astype(np.float64).T
            assert np.abs(P).max() < 2**52
            total = total + P.astype(np.int64).reshape(-1).astype(object) * (1 << (sx + sz))
    return total


def psd_bound(QI: np.ndarray, scale: int, lam: Fraction, both: bool) -> Fraction | None:
    """A rigorous lower bound >= lam for lambda_min(QI / 2^scale), or None if the Cholesky step fails."""
    tau = lam + Fraction(1, 2)
    T = tau * 2**scale
    assert T.denominator == 1
    A = QI - int(T) * np.eye(len(QI), dtype=np.int64)
    try:
        L = np.linalg.cholesky(A / 2**scale)
    except np.linalg.LinAlgError:
        return None
    s = 30
    Lt = np.round(L * 2**s).astype(np.int64)
    shift = 1 << (2 * s - scale)
    target = [int(v) * shift for v in A.reshape(-1).tolist()]
    prod = gram_product_flint(Lt)
    ss = sum((t - int(p)) ** 2 for t, p in zip(target, prod, strict=True))
    if both:
        prod2 = gram_product_split(Lt)
        ss2 = sum((t - p) ** 2 for t, p in zip(target, prod2.tolist(), strict=True))
        assert ss2 == ss, "the two exact products disagree"
    bound = tau - Fraction(isqrt(ss) + 1, 1 << (2 * s))
    return bound if bound >= lam else None


def check(cert: dict, B3: np.ndarray, RL: dict, both: bool = False) -> tuple[bool, str]:
    scale = cert["scale"]
    QI = expand_gram(B3, cert["values"])
    r = residual(B3, QI, cert["rows"], cert["c"], RL)
    r1 = Fraction(int(np.abs(r).sum()), 2**scale)
    c = Fraction(cert["c"], 2**scale)
    margin = cert["lam"] * len(B3) + c
    if c < 0 or margin <= r1:
        return (
            False,
            f"c = {float(c):.4e}, ||r||_1 = {float(r1):.4e} >= lambda |B_3| + c = {float(margin):.1f}",
        )
    lam = psd_bound(QI, scale, cert["lam"], both)
    if lam is None:
        return False, f"lambda_min(Q) >= {cert['lam']} not certified"
    ok = True
    return ok, (
        f"lambda_min(Q) >= {float(lam):.6f}, c = {float(c):.4e}, ||r||_1 = {float(r1):.4e}, "
        f"lambda |B_3| + c = {float(margin):.1f}"
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("cert", nargs="?", type=Path, default=DEFAULT_CERT)
    ap.add_argument("--no-mutations", action="store_true")
    args = ap.parse_args()
    cert = read_certificate(args.cert)
    B3 = b3()
    RL = autocorrelations()
    log(
        f"certificate: {len(cert['values'])} Gram orbit values, {len(cert['rows'])} multiplier orbits, |B_3| = {len(B3)}"
    )
    ok, msg = check(cert, B3, RL, both=True)
    log(("PASS: " if ok else "FAIL: ") + msg)
    if not ok:
        return 1
    if not args.no_mutations:
        rejected = 0
        bad = dict(cert, values=cert["values"].copy())
        bad["values"][0] += 10**5 << cert["scale"]
        rejected += not check(bad, B3, RL)[0]
        bad = dict(cert, rows=cert["rows"].copy())
        bad["rows"][0, 10] += 1000 << cert["scale"]
        rejected += not check(bad, B3, RL)[0]
        bad = dict(cert, lam=Fraction(12))
        rejected += not check(bad, B3, RL)[0]
        log(f"mutations rejected: {rejected}/3")
        if rejected != 3:
            return 1
    print(
        "RESULT: PASS (no common zero of R_1..R_13 on T^9, so H_27 has no flat point and B_27 < sqrt 27)"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
