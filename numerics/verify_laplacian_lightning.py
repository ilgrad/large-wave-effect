# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2", "scipy>=1.11"]
# ///
"""Lightning as energy concentration on the discrete Laplacian: dielectric-breakdown (Laplacian) growth.

Energy concentration is not only a wave phenomenon. A lightning leader is the electrostatic version:
charge separation raises the field until it exceeds the breakdown threshold, and the discharge advances
where the field is strongest. The standard model (Niemeyer-Pietronero-Wiesmann dielectric breakdown) is
LAPLACIAN GROWTH on a lattice -- exactly the discrete Laplacian of this paper, now as the field solver:
the potential satisfies the discrete Laplace equation (L^{(2)} phi)_{ij}=0 off the leader (phi=0 on the
leader, phi=1 on the far electrode), and the leader grows into an adjacent empty site with probability
proportional to the local field |phi|^eta. The result is a branched, fractal discharge (like the
branched flow of refractive caustics), with the field -- and the energy -- concentrated at the tips.

This is a DIFFERENT mechanism from the wave large wave (electrostatic, dissipative breakdown, an elliptic
Laplace problem -- not the conservative dispersive wave dynamics), but it shares the discrete operator,
the branched morphology, and the concentration of energy. We verify:
 (A) the discrete-Laplace solver is correct: phi is harmonic off the leader (max |L^{(2)} phi| ~ 0);
 (B) the discharge is fractal/branched: mass M ~ R_g^D with 1.4 < D < 1.9 (compact would be D=2);
 (C) the field concentrates at a few tips: max candidate field >> mean (the lightning tip).
"""

from __future__ import annotations

import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve

RNG = np.random.default_rng(0)


def solve_potential(leader: np.ndarray) -> np.ndarray:
    """Discrete Laplace: phi=0 on leader, phi=1 on bottom row (ground), periodic in x; harmonic else."""
    ny, nx = leader.shape
    phi = np.zeros((ny, nx))
    phi[-1, :] = 1.0  # ground electrode
    free = ~leader.copy()
    free[-1, :] = False  # ground row fixed
    free[0, :] = False   # top row fixed at 0 (cloud plane)
    idx = -np.ones((ny, nx), dtype=int)
    coords = np.argwhere(free)
    for m, (i, j) in enumerate(coords):
        idx[i, j] = m
    nfree = len(coords)
    a = lil_matrix((nfree, nfree))
    b = np.zeros(nfree)
    for m, (i, j) in enumerate(coords):
        a[m, m] = -4.0
        for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            ii, jj = i + di, (j + dj) % nx
            if 0 <= ii < ny:
                if free[ii, jj]:
                    a[m, idx[ii, jj]] += 1.0
                else:
                    b[m] -= phi[ii, jj] if not leader[ii, jj] else 0.0
    sol = spsolve(a.tocsr(), b)
    for m, (i, j) in enumerate(coords):
        phi[i, j] = sol[m]
    return phi


def grow(ny: int, nx: int, steps: int, eta: float = 1.0) -> dict:
    leader = np.zeros((ny, nx), dtype=bool)
    leader[0, nx // 2] = True  # seed at the cloud (top centre)
    masses, radii, tip_ratio = [], [], []
    laplace_err = 0.0
    reached = False
    for s in range(steps):
        phi = solve_potential(leader)
        if s == steps // 2:  # check harmonicity once, mid-run
            lap = (np.roll(phi, 1, 0) + np.roll(phi, -1, 0) + np.roll(phi, 1, 1)
                   + np.roll(phi, -1, 1) - 4 * phi)
            interior = ~leader
            interior[0, :] = interior[-1, :] = False
            laplace_err = float(np.max(np.abs(lap[interior])))
        # candidate empty sites adjacent to the leader
        cand = []
        for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            nb = np.roll(leader, (di, dj), (0, 1))
            cand.append(nb)
        adj = np.any(cand, axis=0) & ~leader
        adj[0, :] = adj[-1, :] = False
        ci, cj = np.where(adj)
        if len(ci) == 0:
            break
        weight = np.maximum(phi[ci, cj], 0.0) ** eta
        if weight.sum() <= 0:
            weight = np.ones_like(weight)
        tip_ratio.append(float(weight.max() / weight.mean()))
        pick = RNG.choice(len(ci), p=weight / weight.sum())
        leader[ci[pick], cj[pick]] = True
        pts = np.argwhere(leader)
        com = pts.mean(0)
        masses.append(len(pts))
        radii.append(float(np.sqrt(((pts - com) ** 2).sum(1).mean())) + 1e-9)
        if ci[pick] == ny - 2:
            reached = True
            break
    return {"leader": leader, "mass": np.array(masses), "rg": np.array(radii),
            "laplace_err": laplace_err, "tip_ratio": np.array(tip_ratio), "reached": reached}


def main() -> int:
    print("=" * 72)
    print("Lightning as Laplacian growth on the discrete Laplacian (dielectric breakdown)")
    print("=" * 72)
    ok = True
    r = grow(70, 81, steps=900, eta=1.0)

    print(f"\n(A) Discrete-Laplace solver: max |L^(2) phi| off the leader = {r['laplace_err']:.2e}")
    a_ok = r["laplace_err"] < 1e-9
    ok &= a_ok
    print(f"    phi is harmonic (the discrete Laplacian is the field solver): {'OK' if a_ok else 'FAIL'}")

    m, rg = r["mass"], r["rg"]
    win = m > 25  # fit the scaling regime
    d_frac = float(np.polyfit(np.log(rg[win]), np.log(m[win]), 1)[0])
    b_ok = 1.4 < d_frac < 1.9
    ok &= b_ok
    print(f"\n(B) Fractal dimension from M ~ R_g^D: D = {d_frac:.2f} (branched; compact would be 2.0): "
          f"{'OK' if b_ok else 'FAIL'} (final mass {m[-1]}, {'reached ground' if r['reached'] else 'grew'})")

    tr = float(np.median(r["tip_ratio"]))
    c_ok = tr > 3.0
    ok &= c_ok
    print(f"\n(C) Field concentration at tips: median max/mean candidate field = {tr:.1f} "
          f"(>> 1 => energy at the tips): {'OK' if c_ok else 'FAIL'}")

    print("\n" + "=" * 72)
    print("RESULT:", "LAPLACIAN-GROWTH LIGHTNING VERIFIED" if ok else "CHECK FAILED")
    print("The discrete Laplacian also concentrates energy electrostatically: it solves the field of a")
    print("dielectric breakdown, growing a branched fractal lightning leader with the field at its tips.")
    print("A different (elliptic, dissipative) mechanism sharing the operator and the branched morphology.")
    print("=" * 72)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
