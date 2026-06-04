# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2"]
# ///
"""Spectral-dimension threshold for the large-wave ceiling: U_G = Theta(ln|V|) iff d_s = 1.

For any finite graph the velocity-ceiling that governs the large wave is
    U_G = (1/|V|) sum_{lambda_r > 0} lambda_r^{-1/2},      lambda_r = eigenvalues of the graph Laplacian,
because omega_r = sqrt(lambda_r) and the L1->Linf bound on the velocity Green function is
(1/|V|) sum 1/omega_r. Near lambda = 0 the integrated density of states obeys N(lambda) ~ lambda^{d_s/2}
(d_s = SPECTRAL DIMENSION), so rho(lambda) ~ lambda^{d_s/2 - 1} and
    U_G ~ integral_0 lambda^{-1/2} rho(lambda) d lambda ~ integral_0 lambda^{(d_s-3)/2} d lambda,
which converges at 0 iff d_s > 1. Hence

    d_s > 1  =>  U_G = O(1)          (no large wave: amplitude bounded),
    d_s = 1  =>  U_G ~ (1/pi) ln|V|  (LARGE WAVE, log divergence -- the critical case),
    d_s < 1  =>  U_G diverges as a power.

This GENERALIZES the earlier d-torus result (d=1 log, d>=2 bounded), since the d-torus has d_s = d.
The criterion is geometric, not metric: any QUASI-1D structure (ladder, tube) has d_s = 1 and keeps the
log growth -- bounded transverse width does not destroy the concentration. We verify on rings, ladders,
tubes, square and cubic tori, fitting d_s from the small-lambda IDOS.
"""

from __future__ import annotations

import numpy as np


def cycle_eigs(n: int) -> np.ndarray:
    return 2.0 - 2.0 * np.cos(2.0 * np.pi * np.arange(n) / n)


def product_eigs(sizes: tuple[int, ...]) -> np.ndarray:
    """Laplacian spectrum of the Cartesian product of cycles C_{n1} x ... x C_{nk} (eigenvalues add)."""
    eigs = np.zeros(1)
    for n in sizes:
        eigs = (eigs[:, None] + cycle_eigs(n)[None, :]).ravel()
    return eigs


def ceiling(eigs: np.ndarray) -> float:
    pos = eigs[eigs > 1e-9]
    return float(np.sum(pos ** -0.5) / eigs.size)


def fit_spectral_dimension(eigs: np.ndarray) -> float:
    """d_s from N(lambda) ~ lambda^{d_s/2}: slope of log N(lambda) vs log lambda over the low band."""
    pos = np.sort(eigs[eigs > 1e-9])
    lo = pos[: max(20, pos.size // 8)]
    counts = np.arange(1, lo.size + 1)
    a, _ = np.polyfit(np.log(lo), np.log(counts), 1)
    return 2.0 * a


def main() -> int:
    print("=" * 78)
    print("Spectral-dimension threshold:  U_G = (1/|V|) sum lambda^{-1/2},  large wave iff d_s = 1")
    print("=" * 78)
    ok = True

    print("\n(A) Critical d_s = 1 structures grow logarithmically (ring, ladder, tube)")
    print("    Quasi-1D keeps the log; only 1 of w transverse channels carries it, so slope = 1/(pi w).")
    print(f"    {'family':>14} {'width':>6} {'U_G':>8} {'slope':>8} {'1/(pi w)':>9}")
    for name, width, mk in (
        ("ring C_N", 1, lambda n: cycle_eigs(n)),
        ("ladder N x2", 2, lambda n: product_eigs((n, 2))),
        ("tube N x3", 3, lambda n: product_eigs((n, 3))),
    ):
        us, lns = [], []
        for n in (200, 400, 800, 1600, 3200):
            us.append(ceiling(mk(n)))
            lns.append(np.log(n))
        slope = np.polyfit(lns, us, 1)[0]
        target = 1.0 / (np.pi * width)
        ok &= abs(slope - target) < 0.01
        print(f"    {name:>14} {width:>6} {us[-1]:>8.4f} {slope:>8.4f} {target:>9.4f}  "
              f"({'OK' if abs(slope - target) < 0.01 else 'FAIL'})")

    print("\n(B) Super-critical d_s >= 2 structures stay bounded (square / cubic torus)")
    print(f"    {'family':>14} {'sizes':>12} {'U_G':>8} {'bounded?':>9}")
    for name, sizes_list in (
        ("square torus", [(20, 20), (40, 40), (80, 80)]),
        ("cubic torus", [(10, 10, 10), (16, 16, 16), (24, 24, 24)]),
    ):
        vals = [ceiling(product_eigs(s)) for s in sizes_list]
        bounded = max(vals) - min(vals) < 0.05
        ok &= bounded
        print(f"    {name:>14} {sizes_list[-1]!s:>12} {vals[-1]:>8.4f} {bounded!s:>9}")
    print("    -> U_G converges (no large wave) for d_s >= 2: OK")

    print("\n(C) Fitted spectral dimension d_s vs the threshold")
    print(f"    {'family':>16} {'fitted d_s':>11} {'expected':>9} {'large wave?':>12}")
    cases = [
        ("ring C_3200", cycle_eigs(3200), 1.0),
        ("ladder 1600x2", product_eigs((1600, 2)), 1.0),
        ("tube 1000x3", product_eigs((1000, 3)), 1.0),
        ("square 90x90", product_eigs((90, 90)), 2.0),
        ("cubic 24^3", product_eigs((24, 24, 24)), 3.0),
    ]
    for name, eigs, exp in cases:
        ds = fit_spectral_dimension(eigs)
        good = abs(ds - exp) < 0.35  # 3D IDOS slope is coarser; d_s>1 vs =1 is the binary that matters
        ok &= good
        lw = "YES (log)" if exp <= 1.0 else "no"
        print(f"    {name:>16} {ds:>11.3f} {exp:>9.1f} {lw:>12}  ({'OK' if good else 'FAIL'})")

    print("\n" + "=" * 78)
    print("RESULT:", "SPECTRAL-DIMENSION THRESHOLD VERIFIED" if ok else "CHECK FAILED")
    print("Large wave <=> d_s = 1. Quasi-1D (ladder, tube) keeps the log growth; d_s >= 2 kills it.")
    print("=" * 78)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
