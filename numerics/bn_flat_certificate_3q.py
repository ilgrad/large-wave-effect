# /// script
# requires-python = ">=3.11"
# dependencies = ["python-flint>=0.8"]
# ///
r"""No flat point on the orbit closure H_N for N = 3q (q = 5, 7), hence B_15 < sqrt 15 and B_21 < sqrt 21
(paper, Theorem bn-saturation (4)).  Finds an exact Nullstellensatz certificate, checks it in exact rational
arithmetic, and writes it for the independent PARI/GP check exact/pari/c3q_flat_check.gp.

PARAMETRIZATION.  Write Z/3q = Z/3 x Z/q (CRT) and x(b, a) for x_n with n = b mod 3, n = a mod q.  By
Theorem bn-saturation (1), H_N = {e^{-if} : f real, even, f(0) = 0, f^ supported on {0} and the units}.  For
N = 3q the support condition says that, after subtracting the mean, the 3 x q array f(b, a) has zero row and
column sums; with evenness f(-b, -a) = f(b, a) this leaves h(a) = f(1, a) - mean, sum_a h(a) = 0.  So, up to
a unimodular factor, every point of H_N is
    x(1, a) = y_a,   x(2, a) = y_{-a},   x(0, a) = 1/(y_a y_{-a}),      y in T^q,  prod_a y_a = 1,
and conversely.  A point is flat iff its periodic autocorrelations vanish at every nonzero shift (b, s); the
array is even, so they are real, and with z_a = conj(y_a) = 1/y_a they are the Laurent polynomials
    I_s  = sum_a y_{a+s} y_a y_{-a} + z_{a+s} z_a z_{-a} + y_{-a-s} z_a        (shift (1, s)),  s in Z/q,
    II_s = sum_a y_a y_{-a} z_{a+s} z_{-a-s} + y_{a+s} z_a + y_a z_{a+s}    (shift (0, s)),  1 <= s <= (q-1)/2.

CERTIFICATE.  They live in the group ring Q[L], L = Z^q / Z(1, ..., 1) (prod y_a = 1).  We find
    1 = g_0 I_0 + sum_{u in (Z/q)^*} u.(g_1 I_1) + sum_{u = 1}^{(q-1)/2} u.(h_1 II_1)
with u acting by a -> u a (it maps I_s to I_{us} and II_s to II_{+-us}), g_0 invariant under (Z/q)^* and
under inversion y -> 1/y, g_1 inversion-invariant, h_1 invariant under a -> -a and inversion, all supported
on the shifts that keep every product in the ball  ||e|| = sum_a |e_a - median(e)| <= D  (D = 6 for q = 5,
D = 8 for q = 7).  Such an equivariant certificate exists whenever any certificate in the ball does (average
over the group).  Every I_s, II_s vanishes at a flat point, so no flat point exists.

METHOD.  The equivariant linear system (one equation per orbit of monomials in the ball, one unknown per
orbit of admissible shifts) is row-reduced modulo P to choose a nonsingular square block, which is solved over
Q by FLINT; all equations, and then the expanded identity in Q[L], are checked exactly.

    uv run --script numerics/bn_flat_certificate_3q.py [q] [--out FILE]      (q = 7 by default; ~3 min)
"""

from __future__ import annotations

import argparse
import itertools
import sys
import time

import flint

P = 1048573
DEGREE = {5: 6, 7: 8}


def canon(e: tuple[int, ...]) -> tuple[int, ...]:
    m = min(e)
    return tuple(x - m for x in e)


def norm(e: tuple[int, ...]) -> int:
    med = sorted(e)[len(e) // 2]
    return sum(abs(x - med) for x in e)


def generators(q: int) -> list[dict[tuple[int, ...], int]]:
    """I_0, ..., I_{q-1}, II_1, ..., II_{(q-1)/2} as dicts canonical exponent -> coefficient."""

    def mono(pos=(), neg=()):
        e = [0] * q
        for a in pos:
            e[a % q] += 1
        for a in neg:
            e[a % q] -= 1
        return canon(tuple(e))

    def poly(terms):
        out: dict[tuple[int, ...], int] = {}
        for t in terms:
            out[t] = out.get(t, 0) + 1
        return out

    gens = []
    for s in range(q):
        gens.append(
            poly(
                t
                for a in range(q)
                for t in (
                    mono(pos=(a + s, a, -a)),
                    mono(neg=(a + s, a, -a)),
                    mono(pos=(-a - s,), neg=(a,)),
                )
            )
        )
    for s in range(1, (q - 1) // 2 + 1):
        gens.append(
            poly(
                t
                for a in range(q)
                for t in (
                    mono(pos=(a, -a), neg=(a + s, -a - s)),
                    mono(pos=(a + s,), neg=(a,)),
                    mono(pos=(a,), neg=(a + s,)),
                )
            )
        )
    return gens


class Group:
    def __init__(self, q: int):
        self.q = q
        self.units = list(range(1, q))

    def act(self, u: int, e: tuple[int, ...]) -> tuple[int, ...]:
        out = [0] * self.q
        for a in range(self.q):
            out[(u * a) % self.q] = e[a]
        return tuple(out)

    @staticmethod
    def inv(e: tuple[int, ...]) -> tuple[int, ...]:
        return canon(tuple(-x for x in e))

    def orbit(self, e: tuple[int, ...]) -> set[tuple[int, ...]]:
        out = set()
        for u in self.units:
            f = self.act(u, e)
            out |= {f, self.inv(f)}
        return out


def times(Pd: dict, Qd: dict) -> dict:
    out: dict[tuple[int, ...], int] = {}
    for e, c in Pd.items():
        for f, d in Qd.items():
            k = canon(tuple(x + y for x, y in zip(e, f, strict=True)))
            out[k] = out.get(k, 0) + c * d
    return out


def build(q: int, D: int):
    """Equivariant system: orbit representatives of the ball, and per unknown (kind, shift, column poly)."""
    G = Group(q)
    gens = generators(q)
    I0, I1, II1 = gens[0], gens[1], gens[q]
    B = [e for e in itertools.product(range(D + 1), repeat=q) if min(e) == 0 and norm(e) <= D]
    Bset = set(B)
    rep: dict[tuple[int, ...], int] = {}
    reps = []
    for e in B:
        if e not in rep:
            orb = G.orbit(e)
            for f in orb:
                rep[f] = len(reps)
            reps.append(min(orb))

    def admissible(R):
        return [
            m
            for m in B
            if all(canon(tuple(x + y for x, y in zip(m, e, strict=True))) in Bset for e in R)
        ]

    unknowns = []  # (kind, representative shift, set of shifts in the multiplier, generator, u-range)
    seen: set = set()
    for m in admissible(I0):
        if m not in seen:
            orb = G.orbit(m)
            seen |= orb
            unknowns.append((0, min(orb), orb, I0, [1]))
    seen = set()
    for m in admissible(I1):
        if m not in seen:
            orb = {m, G.inv(m)}
            seen |= orb
            unknowns.append((1, min(orb), orb, I1, G.units))
    seen = set()
    for m in admissible(II1):
        if m not in seen:
            orb = {m, G.inv(m), G.act(q - 1, m), G.inv(G.act(q - 1, m))}
            seen |= orb
            unknowns.append((2, min(orb), orb, II1, list(range(1, (q - 1) // 2 + 1))))
    columns = []
    for _kind, _m, orb, R, urange in unknowns:
        base = times({f: 1 for f in orb}, R)
        tot: dict[tuple[int, ...], int] = {}
        for u in urange:
            for e, c in base.items():
                f = G.act(u, e)
                tot[f] = tot.get(f, 0) + c
        columns.append(tot)
    return reps, rep, unknowns, columns


def pivots(Mat, nrows: int, ncols: int) -> list[int]:
    R, r = Mat.rref()
    out, row = [], 0
    for j in range(ncols):
        if row < r and R[row, j] != 0:
            out.append(j)
            row += 1
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("q", nargs="?", type=int, default=7, choices=sorted(DEGREE))
    ap.add_argument(
        "--out", default=None, help="write the certificate for exact/pari/c3q_flat_check.gp"
    )
    args = ap.parse_args()
    q, D = args.q, DEGREE[args.q]
    t0 = time.time()
    reps, rep, unknowns, columns = build(q, D)
    zero = rep[(0,) * q]
    nr, nc = len(reps), len(unknowns)
    entries = [{rep[e]: c for e, c in col.items() if c and e == reps[rep[e]]} for col in columns]
    print(f"N = {3 * q}: ball radius {D}, {nr} orbit equations, {nc} orbit unknowns")
    A = flint.nmod_mat(nr, nc, P)
    for j, col in enumerate(entries):
        for i, c in col.items():
            A[i, j] = c % P
    piv_cols = pivots(A, nr, nc)
    piv_rows = pivots(A.transpose(), nc, nr)
    r = len(piv_cols)
    assert len(piv_rows) == r and zero in piv_rows
    print(f"  rank {r} mod {P}  [{time.time() - t0:.0f}s]")
    pos = {i: k for k, i in enumerate(piv_rows)}
    M = flint.fmpz_mat(r, r)
    for jj, j in enumerate(piv_cols):
        for i, c in entries[j].items():
            if i in pos:
                M[pos[i], jj] = c
    rhs = flint.fmpz_mat(r, 1)
    rhs[pos[zero], 0] = 1
    X = M.solve(rhs)
    sol = [flint.fmpq(0)] * nc
    for jj, j in enumerate(piv_cols):
        sol[j] = X[jj, 0]
    print(f"  solved over Q  [{time.time() - t0:.0f}s]")
    acc = [flint.fmpq(0)] * nr
    for j, col in enumerate(entries):
        if sol[j] != 0:
            for i, c in col.items():
                acc[i] += c * sol[j]
    eq_ok = all(acc[i] == (1 if i == zero else 0) for i in range(nr))
    tot: dict[tuple[int, ...], flint.fmpq] = {}
    for j, col in enumerate(columns):
        if sol[j] != 0:
            for e, c in col.items():
                tot[e] = tot.get(e, flint.fmpq(0)) + c * sol[j]
    nonzero = {e: v for e, v in tot.items() if v != 0}
    id_ok = list(nonzero) == [(0,) * q] and nonzero[(0,) * q] == 1
    den = flint.fmpz(1)
    for v in sol:
        if v != 0:
            den = den * v.q // den.gcd(v.q)
    nnz = [(j, v) for j, v in enumerate(sol) if v != 0]
    print(
        f"  {len(nnz)} nonzero orbit coefficients, common denominator of {len(str(den))} digits; "
        f"all {nr} orbit equations exact: {eq_ok}; expanded identity sum = 1 in Q[L]: {id_ok}"
    )
    if args.out:
        with open(args.out, "w", encoding="utf-8") as fh:
            fh.write(f"[{q}, {D}, {den}]\n")
            for j, v in nnz:
                kind, m = unknowns[j][0], unknowns[j][1]
                num = v.p * (den // v.q)
                fh.write(f"[{kind}, [{', '.join(map(str, m))}], {num}]\n")
        print(f"  certificate written to {args.out}")
    ok = eq_ok and id_ok
    print(
        "RESULT:",
        f"PASS -- 1 lies in the C_{3 * q} flatness ideal, so H_{3 * q} has no flat point"
        if ok
        else "FAIL",
    )
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
