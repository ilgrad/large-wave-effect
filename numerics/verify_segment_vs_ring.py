# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2"]
# ///
"""Segment (Dirichlet) vs ring (periodic) contrast for the impulse large wave.

Paper #1's thesis: the ring is the HARD open case because its spectrum is twofold
degenerate, omega_r = omega_{N-r} (only ~N/2 distinct frequencies), whereas the
Dirichlet segment tridiag(-1,2,-1) has N DISTINCT frequencies
omega_k = 2 sin(pi k / (2(N+1))).

We make the degeneracy contrast RIGOROUS via a quantity that converges (a time
average, not a sup): the mean square of the impulse Green's function at a node.
With independent phases the time-mean would be (1/2) sum coeff_k^2; degeneracy
doubles the contribution of every locked pair. Hence:

    ring:    mean_t |G(node,t)|^2  =  2 * (naive 1/2 sum coeff^2)   (factor ~2)
    segment: mean_t |G(node,t)|^2  =      (naive 1/2 sum coeff^2)   (factor ~1)

We also run an equal-budget sup-scan (honest lower bounds) for context; the true
sup is not scan-reachable for either geometry (almost-periodic, exp-large times).
"""

from __future__ import annotations

import numpy as np

NS = [32, 64, 128, 256]
RMS_TMAX_MULT = 256.0
RMS_DT = 0.05
SUP_TMAX_MULT = 64.0
SUP_DT = 0.1


# ---- ring (periodic) -------------------------------------------------------
def ring_omega(n: int) -> np.ndarray:
    return 2.0 * np.abs(np.sin(np.pi * np.arange(n) / n))


def ring_g_at_node0_coeffs(n: int) -> tuple[np.ndarray, np.ndarray]:
    """G_ring(0,t) = sum_r coeff_r sin(omega_r t); return (omega_nz, coeff_nz)."""
    om = ring_omega(n)
    nz = om > 1e-12
    coeff = (1.0 / n) / om[nz]
    return om[nz], coeff


def ring_sup_scan(n: int, t_max: float, dt: float, block: int = 16384) -> float:
    om = ring_omega(n)
    nz = om > 1e-12
    ts = np.arange(0.0, t_max, dt)
    best = 0.0
    for s in range(0, len(ts), block):
        tb = ts[s:s + block]
        c = np.zeros((len(tb), n))
        c[:, nz] = np.sin(np.outer(tb, om[nz])) / om[nz]
        u = np.fft.ifft(c, axis=1).real
        best = max(best, float(np.abs(u).max()))
    return best


# ---- segment (Dirichlet) ---------------------------------------------------
def segment_modes(n: int) -> tuple[np.ndarray, np.ndarray]:
    """Return (omega[k], V[j,k]) for tridiag(-1,2,-1), k,j = 1..n (0-indexed arrays)."""
    k = np.arange(1, n + 1)
    omega = 2.0 * np.sin(np.pi * k / (2 * (n + 1)))
    j = np.arange(1, n + 1)
    v = np.sqrt(2.0 / (n + 1)) * np.sin(np.pi * np.outer(j, k) / (n + 1))
    return omega, v


def segment_g_at_node_coeffs(n: int, j0: int) -> tuple[np.ndarray, np.ndarray]:
    """G_seg(j0,t) = sum_k coeff_k sin(omega_k t), coeff_k = v_k(j0)^2 / omega_k."""
    omega, v = segment_modes(n)
    coeff = v[j0 - 1, :] ** 2 / omega
    return omega, coeff


def segment_sup_scan(n: int, j0: int, t_max: float, dt: float, block: int = 8192) -> float:
    omega, v = segment_modes(n)
    src = v[j0 - 1, :] / omega          # v_k(j0)/omega_k
    ts = np.arange(0.0, t_max, dt)
    best = 0.0
    for s in range(0, len(ts), block):
        tb = ts[s:s + block]
        coeff_t = np.sin(np.outer(tb, omega)) * src      # (B, n)
        u = coeff_t @ v.T                                # (B, n) field over all nodes
        best = max(best, float(np.abs(u).max()))
    return best


def time_mean_square(omega: np.ndarray, coeff: np.ndarray, t_max: float, dt: float) -> float:
    """Numerical mean over t of (sum_k coeff_k sin(omega_k t))^2."""
    ts = np.arange(0.0, t_max, dt)
    acc, cnt = 0.0, 0
    for s in range(0, len(ts), 16384):
        tb = ts[s:s + 16384]
        g = np.sin(np.outer(tb, omega)) @ coeff
        acc += float(np.sum(g * g))
        cnt += len(tb)
    return acc / cnt


def main() -> int:
    print("=" * 76)
    print("Segment (Dirichlet) vs ring (periodic): degeneracy contrast")
    print("=" * 76)
    ok = True

    # (A) distinct-frequency count: ring ~ N/2, segment = N
    print("\n(A) Distinct frequencies in the spectrum")
    print(f"  {'N':>5} {'ring distinct':>14} {'~N/2':>6} {'seg distinct':>14} {'=N':>5}")
    for n in NS:
        rd = len(np.unique(np.round(ring_omega(n)[1:], 10)))
        sd = len(np.unique(np.round(segment_modes(n)[0], 10)))
        exp_ring = n // 2
        print(f"  {n:>5} {rd:>14} {exp_ring:>6} {sd:>14} {n:>5}")
        ok &= (rd == exp_ring) and (sd == n)
    print(f"  -> ring has ~N/2 distinct freqs (twofold degeneracy), segment has N: "
          f"{'PASS' if ok else 'FAIL'}")

    # (B) degeneracy factor in the (convergent) time-mean square
    print("\n(B) Degeneracy factor = mean_t|G(node,t)|^2 / (naive 1/2 sum coeff^2)")
    print(f"  {'N':>5} {'ring factor':>13} {'seg factor':>12}")
    ring_ok = seg_ok = True
    for n in NS:
        om_r, co_r = ring_g_at_node0_coeffs(n)
        naive_r = 0.5 * float(np.sum(co_r ** 2))
        mean_r = time_mean_square(om_r, co_r, RMS_TMAX_MULT * n, RMS_DT)
        fac_r = mean_r / naive_r

        j0 = max(1, n // 3)
        om_s, co_s = segment_g_at_node_coeffs(n, j0)
        naive_s = 0.5 * float(np.sum(co_s ** 2))
        mean_s = time_mean_square(om_s, co_s, RMS_TMAX_MULT * n, RMS_DT)
        fac_s = mean_s / naive_s

        print(f"  {n:>5} {fac_r:>13.3f} {fac_s:>12.3f}")
        ring_ok &= 1.85 <= fac_r <= 2.15
        seg_ok &= 0.9 <= fac_s <= 1.1
    print(f"  -> ring factor ~2 (locked pairs): {'PASS' if ring_ok else 'FAIL'}")
    print(f"  -> segment factor ~1 (no degeneracy): {'PASS' if seg_ok else 'FAIL'}")
    ok &= ring_ok and seg_ok

    # (C) equal-budget sup scan (honest lower bounds; true sup not scan-reachable)
    print(f"\n(C) Sup-scan lower bounds at equal budget T={SUP_TMAX_MULT:.0f}N (NOT the true sup)")
    print(f"  {'N':>5} {'A_ring':>9} {'A_ring/lnN':>11} {'A_seg':>9} {'A_seg/lnN':>11}")
    for n in NS:
        a_ring = ring_sup_scan(n, SUP_TMAX_MULT * n, SUP_DT)
        a_seg = segment_sup_scan(n, max(1, n // 3), SUP_TMAX_MULT * n, SUP_DT)
        print(f"  {n:>5} {a_ring:>9.4f} {a_ring / np.log(n):>11.4f} "
              f"{a_seg:>9.4f} {a_seg / np.log(n):>11.4f}")
    print("  (interpretation only: both are truncated lower bounds; see verify_largewave_growth.py)")

    print("\n" + "=" * 76)
    print("RESULT:", "ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED")
    print("Takeaway: the ring's twofold spectral degeneracy omega_r=omega_{N-r} is real and")
    print("quantified (RMS factor 2 vs the segment's 1). This is exactly why the log N lower")
    print("bound is settled on the Dirichlet segment but remains open on the ring (paper #1).")
    print("=" * 76)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
