# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2", "scipy>=1.11"]
# ///
"""Spectral-transfer probe for the omega(m)=2 excess lemma: does a Bedert/Chowla-style
least-eigenvalue (or relation-lattice Gram) quantity track the obstruction `dependent_net`?

This is a TRANSFER TEST, not a proof attempt.  Bedert (arXiv:2509.05260) and
Jin-Milojevic-Tomon-Zhang (arXiv:2509.03490) both bound the NEGATIVITY
    K(n) = - min_x f_A(x),     f_A(x) = sum_{a in A} cos(a x),
of a cosine sum from below (>= n^{1/7-o(1)} resp. n^{1/10-o(1)}).  The JMTZ route is spectral:
for the Cayley graph G = Cay(Z/n, A u -A) the eigenvalues are  hat 1_{AuA-}(xi) = 2 f_A(2 pi xi/n),
so lambda_min(G) = 2 min_x f_A = -2 K(n).  A small |lambda_min| forces a large clique (an AP / Sidon
obstruction), which forces K large.

The excess lemma needs an UPPER bound on  A_N = max_psi sum_r b_r sin(phi_r)  over the relation
subtorus, localized as  A_N <= L_pre + O(1).  At the optimizer this is the SIGNED quantity
    dependent_net = sum_{r>M} b_r sin(phi_r).

We test three concrete dictionaries A <-> dependent band:

  (S1) Cayley/Chowla on the DEPENDENT frequency set.  The dependent modes r=M+1..floor(N/2) are the
       relation-LOCKED ones.  Form the cosine sum over their *indices* and report its min (the direct
       Chowla object), and the least eigenvalue of Cay(Z/N, D u -D), D = {dependent indices}.

  (S2) Relation-lattice GRAM spectrum.  G = B B^T for the LLL-reduced relation basis B (rows = relations
       sum_r k_r omega_r = 0).  Bedert/JMTZ negativity is controlled by smallest eigenvalues of a Gram /
       adjacency matrix; we report lambda_min(Gram), the spectral gap, and the weighted form
       B diag(b) B^T (b = ring weights) whose quadratic form is the natural energy of the band.

  (S3) Weighted dependent-block Gram on the coordinate matrix C (phi_r = sum_s C[r,s] psi_s):
       the dependent rows C[M:] carry the locked phases; M_dep = C[M:] diag over s, and the relevant
       spectral object is the smallest singular value / least eigenvalue of C[M:] C[M:]^T scaled by b.

For each N in {15,21,33,35,45,55,...} we print dependent_net (the target) next to every spectral
candidate, then report the best |Pearson rho| and rank correlation.  A candidate "tracks" only if it
is monotone in dependent_net with |rho| > ~0.9 across the set; anything weaker is reported as NOT a
controlling quantity (honest negative).
"""

from __future__ import annotations

import sys
from fractions import Fraction
from math import gcd, lcm, pi
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from large_wave_effect.cyclotomic import sin_coords

RNG_SEED = 0


def totient(n: int) -> int:
    return sum(1 for k in range(1, n + 1) if gcd(k, n) == 1)


def odd_part(n: int) -> int:
    while n % 2 == 0:
        n //= 2
    return n


def prime_factors(n: int) -> list[int]:
    n = odd_part(n)
    fs, d = [], 3
    while d * d <= n:
        if n % d == 0:
            fs.append(d)
            while n % d == 0:
                n //= d
        d += 2
    if n > 1:
        fs.append(n)
    return fs


def weights(n: int) -> tuple[float, np.ndarray, np.ndarray]:
    m = n // 2
    r = np.arange(1, m + 1)
    omega = 2.0 * np.sin(pi * r / n)
    b = 2.0 / (n * omega)
    if n % 2 == 0:
        b[-1] = 1.0 / (n * omega[-1])
    return float(b.sum()), b, omega


def integer_coord_matrix(n: int) -> np.ndarray:
    """C (m x M): phi_r = sum_s C[r,s] psi_s, exact via Fraction elimination (mirrors verify_excess_omega2)."""
    m, big_m = n // 2, totient(2 * n) // 2
    co = [[Fraction(x) for x in row] for row in sin_coords(2 * n, list(range(1, m + 1)))]
    dim = len(co[0])
    aug = [[co[s][j] for s in range(big_m)] + [co[r][j] for r in range(m)] for j in range(dim)]
    pr = 0
    for col in range(big_m):
        sel = next((i for i in range(pr, dim) if aug[i][col] != 0), None)
        if sel is None:
            continue
        aug[pr], aug[sel] = aug[sel], aug[pr]
        p = aug[pr][col]
        aug[pr] = [v / p for v in aug[pr]]
        for i in range(dim):
            if i != pr and aug[i][col] != 0:
                f = aug[i][col]
                aug[i] = [a - f * b for a, b in zip(aug[i], aug[pr], strict=True)]
        pr += 1
        if pr == big_m:
            break
    c = np.zeros((m, big_m))
    for s in range(big_m):
        for r in range(m):
            val = aug[s][big_m + r]
            assert val.denominator == 1, f"non-integer coordinate at N={n}"
            c[r, s] = float(val)
    return c.astype(int)


def relation_basis(n: int) -> np.ndarray:
    """Integer basis of the relation lattice { k in Z^m : sum_r k_r omega_r = 0 }.

    omega_r is proportional to the cyclotomic coord row coord[r]; a relation is a left-null vector
    of the m x dim integer coordinate matrix, i.e. a nullspace vector of coord^T.  Exact via Fractions.
    """
    m = n // 2
    coord = [[int(x) for x in r] for r in sin_coords(2 * n, list(range(1, m + 1)))]
    at = [[Fraction(coord[r][c]) for r in range(len(coord))] for c in range(len(coord[0]))]
    nr, nc = len(at), len(at[0])
    pivcols, pr = [], 0
    for col in range(nc):
        sel = next((r for r in range(pr, nr) if at[r][col] != 0), None)
        if sel is None:
            continue
        at[pr], at[sel] = at[sel], at[pr]
        piv = at[pr][col]
        at[pr] = [v / piv for v in at[pr]]
        for r in range(nr):
            if r != pr and at[r][col] != 0:
                f = at[r][col]
                at[r] = [a - f * b for a, b in zip(at[r], at[pr], strict=True)]
        pivcols.append(col)
        pr += 1
        if pr == nr:
            break
    free = [c for c in range(nc) if c not in pivcols]
    basis = []
    for fc in free:
        v = [Fraction(0)] * nc
        v[fc] = Fraction(1)
        for i, pc in enumerate(pivcols):
            v[pc] = -at[i][fc]
        den = 1
        for x in v:
            den = lcm(den, x.denominator)
        basis.append([int(x * den) for x in v])
    return np.array(basis, dtype=float) if basis else np.zeros((0, m))


def optimize_AN(n: int, *, seed: int, starts: int) -> tuple[float, np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    big_m = totient(2 * n) // 2
    c = integer_coord_matrix(n)
    _u, b, _omega = weights(n)
    cf = c.astype(float)

    def negf(psi: np.ndarray) -> tuple[float, np.ndarray]:
        ph = cf @ psi
        return -float(b @ np.sin(ph)), -(cf.T @ (b * np.cos(ph)))

    best, bpsi = -np.inf, None
    for _ in range(starts):
        res = minimize(negf, rng.uniform(0, 2 * pi, big_m), jac=True, method="L-BFGS-B")
        if -res.fun > best:
            best, bpsi = -res.fun, res.x
    psi = bpsi.copy()
    for _ in range(80):
        ph = cf @ psi
        g = cf.T @ (b * np.cos(ph))
        h = -(cf.T * (b * np.sin(ph))) @ cf
        try:
            step = np.linalg.solve(h, -g)
        except np.linalg.LinAlgError:
            break
        psi = psi + step
        if np.linalg.norm(step) < 1e-15:
            break
    return float(b @ np.sin(cf @ psi)), np.sin(cf @ psi), b


def chowla_min_over_set(n: int, idx: list[int], *, grid: int = 200000) -> float:
    """min_x sum_{r in idx} cos(2 pi r x / n) on a fine grid (the literal Chowla object on a frequency set)."""
    idxa = np.array(idx, dtype=float)
    xs = np.linspace(0, 2 * pi, grid, endpoint=False)
    vals = np.cos(np.outer(xs, idxa)).sum(axis=1)
    return float(vals.min())


def cayley_least_eig(n: int, idx: list[int]) -> float:
    """Least eigenvalue of Cay(Z/n, idx u -idx): eigenvalues = sum_{a in idx} 2 cos(2 pi a xi / n), xi=0..n-1."""
    xis = np.arange(n)
    idxa = np.array(idx, dtype=float)
    eigs = 2.0 * np.cos(2 * pi * np.outer(xis, idxa) / n).sum(axis=1)
    return float(eigs.min())


def pearson(x: np.ndarray, y: np.ndarray) -> float:
    x = x - x.mean()
    y = y - y.mean()
    denom = np.sqrt((x * x).sum() * (y * y).sum())
    return float((x * y).sum() / denom) if denom > 0 else 0.0


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    rx = np.argsort(np.argsort(x)).astype(float)
    ry = np.argsort(np.argsort(y)).astype(float)
    return pearson(rx, ry)


OMEGA2_NS = [15, 21, 33, 35, 55, 45, 30, 60, 42, 66, 70, 90]


def main() -> int:
    print("=" * 110)
    print("SPECTRAL TRANSFER PROBE  (Bedert/JMTZ Chowla negativity  vs  excess-lemma dependent_net)")
    print("=" * 110)
    print("\nDictionary under test:  Bedert lower-bounds K=-min_x f_A; JMTZ via lambda_min(Cay)= 2 min_x f_A.")
    print("Target = dependent_net = sum_{r>M} b_r sin(phi_r) at the A_N optimizer (the obstruction term).\n")

    hdr = (f"  {'N':>4} {'fac':>8} {'M':>3} {'#dep':>4} "
           f"{'dep_net':>9} {'minChowla':>10} {'lamMin_dep':>11} {'lamMinGram':>11} "
           f"{'lamMinWGram':>12} {'sigMin_Cdep':>12}")
    print(hdr)

    target = []
    cand = {"minChowla": [], "lamMin_dep": [], "lamMinGram": [], "lamMinWGram": [], "sigMin_Cdep": []}

    for n in OMEGA2_NS:
        m = n // 2
        big_m = totient(2 * n) // 2
        _a, sinvals, b = optimize_AN(n, seed=RNG_SEED, starts=200)
        dep_net = float((b[big_m:] * sinvals[big_m:]).sum())

        dep_idx = list(range(big_m + 1, m + 1))  # 1-based ring frequency indices of the dependent band

        # (S1) literal Chowla object + Cayley least eigenvalue on the dependent index set
        min_chowla = chowla_min_over_set(n, dep_idx)
        lam_dep = cayley_least_eig(n, dep_idx)

        # (S2) relation-lattice Gram spectra
        bvec = relation_basis(n)  # (relrank x m)
        assert bvec.shape[0] == m - big_m, f"relation rank mismatch at N={n}: {bvec.shape[0]} != {m - big_m}"
        if bvec.shape[0] > 0:
            gram = bvec @ bvec.T
            lam_gram = float(np.linalg.eigvalsh(gram).min())
            _u, bw, _om = weights(n)
            wgram = (bvec * bw[None, :]) @ bvec.T  # B diag(b) B^T
            lam_wgram = float(np.linalg.eigvalsh(wgram).min())
        else:
            lam_gram = lam_wgram = float("nan")

        # (S3) dependent block of the coordinate matrix C: smallest singular value of C[M:] (weighted)
        c = integer_coord_matrix(n).astype(float)
        _u, bw, _om = weights(n)
        cdep = c[big_m:] * np.sqrt(bw[big_m:])[:, None]
        sv = np.linalg.svd(cdep, compute_uv=False)
        sig_min = float(sv.min()) if sv.size else float("nan")

        target.append(dep_net)
        cand["minChowla"].append(min_chowla)
        cand["lamMin_dep"].append(lam_dep)
        cand["lamMinGram"].append(lam_gram)
        cand["lamMinWGram"].append(lam_wgram)
        cand["sigMin_Cdep"].append(sig_min)

        print(f"  {n:>4} {prime_factors(n)!s:>8} {big_m:>3} {len(dep_idx):>4} "
              f"{dep_net:>9.4f} {min_chowla:>10.4f} {lam_dep:>11.4f} {lam_gram:>11.4f} "
              f"{lam_wgram:>12.5f} {sig_min:>12.5f}")

    t = np.array(target)
    print("\n" + "-" * 110)
    print("CORRELATION WITH dependent_net  (|rho|>0.9 monotone => 'tracks'; else honest NEGATIVE)")
    print(f"  {'candidate':>14} {'Pearson rho':>12} {'Spearman':>10}   verdict")
    any_track = False
    for name, vals in cand.items():
        v = np.array(vals)
        if np.any(~np.isfinite(v)):
            print(f"  {name:>14} {'n/a':>12} {'n/a':>10}   (non-finite entries)")
            continue
        rho = pearson(v, t)
        sp = spearman(v, t)
        tracks = abs(rho) > 0.9 and abs(sp) > 0.9
        any_track |= tracks
        print(f"  {name:>14} {rho:>+12.4f} {sp:>+10.4f}   {'TRACKS' if tracks else 'does NOT track'}")

    print("\n" + "=" * 110)
    if any_track:
        print("RESULT: at least one spectral candidate tracks dependent_net (|rho|>0.9). See table; this is")
        print("NUMERICAL evidence of a controlling quantity, NOT a proof of the excess lemma.")
    else:
        print("RESULT: NO tested spectral quantity tracks dependent_net across N=15..90.  The Bedert/JMTZ")
        print("least-eigenvalue object does not linearly control the excess-lemma obstruction in this range.")
        print("This is the honest negative the transfer analysis predicts (see report): Bedert bounds the")
        print("NEGATIVITY (min) of an UNCONSTRAINED cosine sum from BELOW; the excess lemma needs an UPPER")
        print("bound on a relation-CONSTRAINED weighted SUP.  Different functional, different direction.")
    print("=" * 110)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
