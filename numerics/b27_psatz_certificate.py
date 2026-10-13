# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2", "scipy>=1.13"]
# ///
r"""No flat point on the orbit closure H_27, hence B_27 < sqrt 27 (paper, Theorem bn-saturation (4)).  Finds a
degree-3 Positivstellensatz certificate on the 9-torus and writes it for the exact checks
numerics/verify_b27_psatz.py and exact/pari/b27_psatz_check.gp.

PARAMETRIZATION.  mu_r = 2 cos(2 pi r/27) has integer coordinates c(r) in the basis 1, mu_1, ..., mu_8 of
Q(cos 2pi/27): c(0) = 2 e_0, c(r) = e_r (1 <= r <= 8), c(9) = -e_0, c(9+i) = -(e_i + e_{9-i}) (1 <= i <= 4)
and c(-r) = c(r).  With X(t) = (e^{it}, e^{it mu_1}, ..., e^{it mu_8}), which is dense in T^9,
x_r(t) = e^{-it lambda_r} = X(t)^{c(r) - c(0)}, so H_27 = {(X^{c(r) - c(0)})_r : X in T^9}.  A point is flat iff
its periodic autocorrelations R_l(X) = sum_k X^{c(k+l) - c(k)}, l = 1, ..., 13, all vanish.

CERTIFICATE.  Laurent polynomials g_l, a real symmetric matrix Q >= lambda I indexed by the set B_3 of
monomials of word length <= 3 in the generators +-(c(r) - c(0)), and c >= 0 with
    sigma + sum_l g_l R_l + c = r,      sigma = m^* Q m,  m = (X^alpha)_{alpha in B_3},
where r is a small rounding residual.  On T^9 every |X^alpha| = 1, so sigma >= lambda |B_3|.  At a common zero
of the R_l this gives lambda |B_3| + c <= |r(X)| <= ||r||_1, and the certificate has ||r||_1 far below that.

METHOD.  G = <multiplier 2> x <inversion> (order 18) acts on all the data.  The SDP is reduced to character
blocks (10 blocks of size 168..190, 19993 orbit equations, 14075 multiplier orbits) and solved by
Douglas-Rachford splitting between the affine constraint (one dense Cholesky factorization of A A^*) and the
shifted cone {Q_b >= I}, until every block has lambda_min >= 1/2 at a point of the affine set.  That point is
unfolded to the full Gram matrix on B_3 (3139 monomials), averaged over G, rounded to multiples of 2^-12, and
written as one value per G-orbit of index pairs, together with the multiplier orbits and c.

    uv run --script numerics/b27_psatz_certificate.py [--out FILE]                       (CPU, scipy)
    uv run --python 3.13 --with cupy-cuda13x --script numerics/b27_psatz_certificate.py --gpu     (CUDA, CuPy)

Both paths run the same iteration (they agree to the printed digits at iteration 50); about 1400 iterations
are needed.  The output is checked, not trusted: a rerun can differ in the last digits (floating-point
reduction order), and the exact checks decide its validity.
"""

from __future__ import annotations

import argparse
import gzip
import sys
import time
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp

N = 27
DEGREE = 3
SCALE = 12
DEFAULT_OUT = Path(__file__).resolve().parent.parent / "exact" / "pari" / "b27_psatz_cert.gp.gz"
T0 = time.time()


def log(msg: str) -> None:
    print(f"{msg}  ({time.time() - T0:.1f}s)", flush=True)


def coords() -> np.ndarray:
    c = np.zeros((N, 9), dtype=np.int64)
    c[0, 0] = 2
    for r in range(1, 9):
        c[r, r] = 1
    c[9, 0] = -1
    for i in range(1, 5):
        c[9 + i, i] = c[9 + i, 9 - i] = -1
    for r in range(14, N):
        c[r] = c[N - r]
    return c


C = coords()
GENS = [tuple(int(v) for v in C[r] - C[0]) for r in range(1, 14)]
GENS += [tuple(-v for v in g) for g in GENS]


def ball(radius: int) -> dict[tuple[int, ...], int]:
    """Word length <= radius in GENS, by breadth-first search."""
    zero = (0,) * 9
    length = {zero: 0}
    layer = [zero]
    for d in range(1, radius + 1):
        nxt = []
        for e in layer:
            for g in GENS:
                f = tuple(a + b for a, b in zip(e, g, strict=True))
                if f not in length:
                    length[f] = d
                    nxt.append(f)
        layer = nxt
    return length


def autocorrelation(ell: int) -> dict[tuple[int, ...], int]:
    out: dict[tuple[int, ...], int] = {}
    for k in range(N):
        e = tuple(int(v) for v in C[(k + ell) % N] - C[k])
        out[e] = out.get(e, 0) + 1
    return out


RL = {ell: autocorrelation(ell) for ell in range(1, 14)}


def fold(ell: int) -> int:
    ell %= N
    return min(ell, N - ell)


def multiplier_matrix(u: int) -> np.ndarray:
    """M with M c(r) = c(u r) for every r."""
    m = np.zeros((9, 9), dtype=np.int64)
    m[0, 0] = 1
    for b in range(1, 9):
        m[:, b] = C[(u * b) % N]
    return m


M2 = multiplier_matrix(2)
assert all(np.array_equal(M2 @ C[r], C[(2 * r) % N]) for r in range(N))
GROUP = []  # (a, b, M): X^e -> X^{M e} for the multiplier 2^a followed by inversion if b = 1
_m = np.eye(9, dtype=np.int64)
for _a in range(9):
    GROUP.append((_a, 0, _m.copy()))
    GROUP.append((_a, 1, -_m))
    _m = M2 @ _m
OFFSET = 40


def encode(a: np.ndarray) -> np.ndarray:
    a = a + OFFSET
    assert a.min() >= 0 and a.max() < 2 * OFFSET
    out = np.zeros(a.shape[:-1], dtype=np.int64)
    for i in range(9):
        out = out * (2 * OFFSET) + a[..., i]
    return out


def canon(a: np.ndarray) -> np.ndarray:
    return np.min(np.stack([encode(a @ m.T) for _, _, m in GROUP]), axis=0)


def chi(j: int, s: int, k: int) -> complex:
    a, b, _ = GROUP[k]
    return np.exp(2j * np.pi * j * a / 9) * (-1) ** (s * b)


@dataclass(frozen=True, slots=True)
class Block:
    """Gram block of the character (j, s) on the orbits `valid`; Bt maps vec(Q_b) to the orbit equations."""

    j: int
    s: int
    valid: list[int]
    Bt: sp.csr_matrix

    @property
    def m(self) -> int:
        return len(self.valid)

    @property
    def cplx(self) -> bool:
        return self.j != 0


class Problem:
    """The G-symmetric Positivstellensatz of order DEGREE: Gram blocks per character, Macaulay rows per orbit."""

    def __init__(self) -> None:
        d = DEGREE
        length = ball(2 * d)
        self.pts = np.array(list(length), dtype=np.int64)
        nrm = np.array(list(length.values()), dtype=np.int64)
        code_all = canon(self.pts)
        orb_codes, orb_cnt = np.unique(code_all, return_counts=True)
        orbit_of = dict(zip(orb_codes.tolist(), range(len(orb_codes)), strict=True))
        stab_size = 18 // orb_cnt
        self.n_orb = len(orb_codes)
        self.zero = orbit_of[int(canon(np.zeros((1, 9), dtype=np.int64))[0])]
        self.reps = []
        for code in np.unique(code_all[nrm <= d]).tolist():
            a = self.pts[int(np.nonzero(code_all == code)[0][0])]
            cand = np.stack([a @ m.T for _, _, m in GROUP])
            a = cand[int(np.argmin(encode(cand)))]
            stab = [k for k, (_, _, m) in enumerate(GROUP) if np.array_equal(m @ a, a)]
            self.reps.append((a, stab))
        A_O = np.stack([r[0] for r in self.reps])
        nO = len(self.reps)
        pair_o = np.zeros((18, nO, nO), dtype=np.int64)
        for k, (_, _, m) in enumerate(GROUP):
            E = (A_O @ m.T)[:, None, :] - A_O[None, :, :]
            pair_o[k] = np.array([orbit_of[c] for c in canon(E.reshape(-1, 9)).tolist()]).reshape(
                nO, nO
            )
        # block of character (j, s): Gram entry (p, q) is the coefficient vector of psi_p conj(psi_q),
        # psi_p = sum_k chi(k) X^{k a_p}; j = 1..4 stands for the conjugate pair j, 9 - j (weight 2)
        self.blocks = []
        for j in range(5):
            for s in (0, 1):
                valid = [
                    o
                    for o, (_, st) in enumerate(self.reps)
                    if all(abs(chi(j, s, k) - 1) < 1e-9 for k in st)
                ]
                m = len(valid)
                rows, cols, vals = [], [], []
                for p, o1 in enumerate(valid):
                    for q, o2 in enumerate(valid):
                        acc: dict[int, complex] = {}
                        for k in range(18):
                            o = int(pair_o[k, o1, o2])
                            acc[o] = acc.get(o, 0) + chi(j, s, k) * stab_size[o]
                        for o, v in acc.items():
                            if abs(v) > 1e-12:
                                rows.append(o)
                                cols.append(p * m + q)
                                vals.append(v)
                v = np.array(vals) * (1.0 if j == 0 else 2.0)
                if j == 0:
                    v = v.real
                A = sp.csr_matrix((v, (rows, cols)), shape=(self.n_orb, m * m))
                pq = np.arange(m * m)
                Bt = (
                    (A + A.conj()[:, (pq % m) * m + pq // m]) / 2
                ).tocsr()  # Hermitian part of each row
                Bt.eliminate_zeros()
                self.blocks.append(Block(j=j, s=s, valid=valid, Bt=Bt))
        log(f"blocks {[b.m for b in self.blocks]}, orbit equations {self.n_orb}")
        # Macaulay rows: sum_{k in G} X^{k e} R_{2^a j}, e over Stab_j-orbit representatives of B_{2d-2}
        hpts = self.pts[nrm <= 2 * d - 2]
        mr, mc, mv, self.rows = [], [], [], []
        for j in (1, 3, 9):
            stab = [k for k, (a, _, _) in enumerate(GROUP) if fold(2**a * j) == j]
            sc = np.min(np.stack([encode(hpts @ GROUP[k][2].T) for k in stab]), axis=0)
            _, first = np.unique(sc, return_index=True)
            reps_j = hpts[np.sort(first)]
            base = len(self.rows)
            self.rows += [(j, e) for e in reps_j]
            for a, _, m in GROUP:
                R = RL[fold(2**a * j)]
                T = np.array(list(R), dtype=np.int64)
                cf = np.array(list(R.values()), dtype=np.float64)
                S = (reps_j @ m.T)[:, None, :] + T[None, :, :]
                oid = np.array([orbit_of[c] for c in canon(S.reshape(-1, 9)).tolist()])
                mr.append(np.repeat(np.arange(len(reps_j)) + base, len(T)))
                mc.append(oid)
                mv.append(np.tile(cf, len(reps_j)) / orb_cnt[oid])
        self.Mac = sp.csr_matrix(
            (np.concatenate(mv), (np.concatenate(mr), np.concatenate(mc))),
            shape=(len(self.rows), self.n_orb),
        )
        log(f"multiplier orbits {len(self.rows)}")

    def normal_terms(self) -> Iterator[sp.coo_matrix]:
        """Sparse summands of A A^* (real parts) for A: (Q, y, c) -> sum_b Re(Bt_b vec Q_b) + Mac^T y + c e_0."""
        for b in self.blocks:
            yield (b.Bt @ b.Bt.conj().T).tocoo()
        MacT = self.Mac.T.tocsr()
        for i0 in range(0, self.Mac.shape[0], 500):
            yield (MacT[:, i0 : i0 + 500] @ self.Mac[i0 : i0 + 500]).tocoo()
        yield sp.coo_matrix(([1.0], ([self.zero], [self.zero])), shape=(self.n_orb, self.n_orb))


def douglas_rachford(
    prob: Problem, gpu: bool, max_iter: int
) -> tuple[list[np.ndarray], np.ndarray, float]:
    """A point of {A(Q, y, c) = 0} with every Q_b >= I/2, by DR splitting with the cone {Q_b >= I, c >= 0}."""
    n = prob.n_orb
    if gpu:
        import cupy as cp  # ty: ignore[unresolved-import]
        import cupyx  # ty: ignore[unresolved-import]
        import cupyx.scipy.sparse as xsp  # ty: ignore[unresolved-import]
        from cupy.cuda import device  # ty: ignore[unresolved-import]
        from cupy_backends.cuda.libs import cusolver  # ty: ignore[unresolved-import]

        # assembled on the device: the dense matrix takes 3.2 GB, and a host copy would double that
        S_dev = cp.zeros((n, n))
        flat = S_dev.ravel()
        for P in prob.normal_terms():
            index = cp.asarray(P.row.astype(np.int64) * n + P.col)
            cupyx.scatter_add(flat, index, cp.asarray(np.ascontiguousarray(P.data.real)))
        handle = device.get_cusolver_handle()
        lwork = cusolver.dpotrf_bufferSize(handle, 0, n, S_dev.data.ptr, n)
        work = cp.empty(lwork, dtype=cp.float64)
        info = cp.zeros(1, dtype=cp.int32)
        cusolver.dpotrf(handle, 0, n, S_dev.data.ptr, n, work.data.ptr, lwork, info.data.ptr)
        assert int(info[0]) == 0
        del work

        def solve(r):
            x = cp.array(r, dtype=cp.float64)
            cusolver.dpotrs(handle, 0, n, 1, S_dev.data.ptr, n, x.data.ptr, n, info.data.ptr)
            return x

        xp = cp
        dev = xsp.csr_matrix
        asnumpy = cp.asnumpy
    else:
        xp = np
        S = np.zeros((n, n), order="F")  # Fortran order, so that the Cholesky factor overwrites S
        for P in prob.normal_terms():
            S[P.row, P.col] += P.data.real
        L = sla.cholesky(S, lower=True, overwrite_a=True, check_finite=False)

        def solve(r):
            return sla.cho_solve((L, True), r, check_finite=False)

        def dev(a):
            return a

        def asnumpy(a):
            return a

    log("normal matrix factored")
    blocks = prob.blocks
    Bt = [dev(b.Bt) for b in blocks]
    BtT = [dev(b.Bt.T.tocsr()) for b in blocks]
    Mac, MacT = dev(prob.Mac), dev(prob.Mac.T.tocsr())
    zero = prob.zero

    def fwd(Q, y, c):
        r = MacT @ y
        for bt, q in zip(Bt, Q, strict=True):
            r += (bt @ q.ravel()).real
        r[zero] += c
        return r

    def proj_affine(Q, y, c):
        u = solve(fwd(Q, y, c))
        uc = u.astype(xp.complex128)
        out = []
        for b, btt, q in zip(blocks, BtT, Q, strict=True):
            if b.cplx:
                out.append(q - (btt @ uc).reshape(b.m, b.m).conj())
            else:
                out.append(q - (btt @ u).reshape(b.m, b.m))
        return out, y - Mac @ u, c - u[zero]

    def proj_cone(Q, y, c):
        out = []
        for q in Q:
            w, V = xp.linalg.eigh(q)
            out.append((V * xp.maximum(w, 1.0)) @ V.conj().T)
        return out, y, xp.maximum(c, 0.0)

    zQ = [xp.eye(b.m, dtype=xp.complex128) if b.cplx else xp.eye(b.m) for b in blocks]
    zy = xp.zeros(len(prob.rows))
    zc = xp.asarray(0.0)
    for it in range(1, max_iter + 1):
        xQ, xy, xc = proj_cone(zQ, zy, zc)
        wQ, wy, wc = proj_affine(
            [2 * a - b for a, b in zip(xQ, zQ, strict=True)], 2 * xy - zy, 2 * xc - zc
        )
        zQ = [z + w - x for z, w, x in zip(zQ, wQ, xQ, strict=True)]
        zy = zy + wy - xy
        zc = zc + wc - xc
        if it % 50 == 0:
            lam = min(float(xp.linalg.eigvalsh(q)[0]) for q in wQ)
            log(f"DR iteration {it}: lambda_min over the blocks {lam:.4f}, c {float(wc):.3e}")
            if lam >= 0.5:
                return [asnumpy(q) for q in wQ], asnumpy(wy), float(wc)
    raise SystemExit(f"no certificate after {max_iter} iterations")


def b3_order() -> np.ndarray:
    """B_3 in lexicographic order of the exponent vectors."""
    return np.array(sorted(e for e, d in ball(DEGREE).items()), dtype=np.int64)


def pair_orbits(B3: np.ndarray) -> np.ndarray:
    """Canonical key i*n + j (i <= j, the least pair of its orbit) of the G-orbit of every index pair (i, j)."""
    n = len(B3)
    index = {tuple(v): i for i, v in enumerate(B3.tolist())}
    perms = [np.array([index[tuple(v)] for v in (B3 @ m.T).tolist()]) for _, _, m in GROUP]
    key = np.full((n, n), n * n, dtype=np.int64)
    ii, jj = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
    for p in perms:
        a, b = p[ii], p[jj]
        key = np.minimum(key, np.minimum(a, b) * n + np.maximum(a, b))
    return key


def unfold(prob: Problem, Q_blocks: list[np.ndarray], B3: np.ndarray) -> np.ndarray:
    """Full Gram matrix on B_3: sum over blocks of Phi^* Q_b Phi, Phi_p = conj(psi_p) in the monomial basis."""
    index = {tuple(v): i for i, v in enumerate(B3.tolist())}
    Qf = np.zeros((len(B3), len(B3)))
    for b, Qb in zip(prob.blocks, Q_blocks, strict=True):
        Phi = np.zeros((b.m, len(B3)), dtype=complex)
        for p, orb in enumerate(b.valid):
            a = prob.reps[orb][0]
            for k, (_, _, m) in enumerate(GROUP):
                Phi[p, index[tuple((-(m @ a)).tolist())]] += np.conj(chi(b.j, b.s, k))
        H = Phi.conj().T @ Qb @ Phi
        Qf += 2 * H.real if b.cplx else H.real
    return (Qf + Qf.T) / 2


def write_certificate(
    path: Path, values: np.ndarray, rows: list, y: np.ndarray, c: int, lam: str
) -> None:
    lines = [
        "\\\\ B_27 < sqrt 27: degree-3 Positivstellensatz certificate on H_27, written by "
        "numerics/b27_psatz_certificate.py and checked by numerics/verify_b27_psatz.py and "
        "exact/pari/b27_psatz_check.gp.",
        f"\\\\ Header [N, degree, scale, c, claimed lambda, #pair orbits, #rows]; all values are integers / 2^{SCALE}.",
        "\\\\ [1, v...]: Gram values, one per G-orbit of index pairs (i <= j) of B_3 in lexicographic order,",
        "\\\\ in the order of the least pair of each orbit.  [2, j, e_0..e_8, y]: y * sum_{g in G} X^{g e} R_{g j}.",
        f"[{N}, {DEGREE}, {SCALE}, {c}, {lam}, {len(values)}, {len(rows)}]",
    ]
    for i0 in range(0, len(values), 1000):
        lines.append("[1, " + ", ".join(str(int(v)) for v in values[i0 : i0 + 1000]) + "]")
    for (j, e), v in zip(rows, y, strict=True):
        lines.append(f"[2, {j}, " + ", ".join(str(int(t)) for t in e) + f", {int(v)}]")
    # mtime 0 and no file name in the header, so that equal certificates give equal files
    path.write_bytes(gzip.compress(("\n".join(lines) + "\n").encode("ascii"), mtime=0))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--gpu", action="store_true", help="run the splitting on a CUDA device through CuPy"
    )
    ap.add_argument("--max-iter", type=int, default=20000)
    ap.add_argument(
        "--lambda-claim", default="5", help="PSD lower bound written into the certificate"
    )
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args()
    prob = Problem()
    Q_blocks, y, c = douglas_rachford(prob, args.gpu, args.max_iter)
    B3 = b3_order()
    Qf = unfold(prob, Q_blocks, B3)
    key = pair_orbits(B3)
    up = np.triu_indices(len(B3))
    keys, inv = np.unique(key[up], return_inverse=True)
    mean = np.bincount(inv, weights=Qf[up]) / np.bincount(inv)
    values = np.round(mean * 2**SCALE).astype(np.int64)
    QI = values[np.searchsorted(keys, key)]
    lam = np.linalg.eigvalsh(QI / 2**SCALE)[0]
    log(f"|B_3| = {len(B3)}, pair orbits {len(values)}, lambda_min(Q) = {lam:.4f} (floating point)")
    if c < 0:
        raise SystemExit(f"c = {c} < 0")
    yI = np.round(y * 2**SCALE).astype(np.int64)
    cI = round(c * 2**SCALE)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    write_certificate(args.out, values, prob.rows, yI, cI, args.lambda_claim)
    log(f"wrote {args.out} ({args.out.stat().st_size} bytes)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
