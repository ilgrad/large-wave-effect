# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2"]
# ///
"""Disorder vs the large wave: unequal masses break the arithmetic resonances that block composite N,
so the Bohr ceiling becomes reachable for EVERY N -- but the gain is realised only over long horizons.

Unequal masses m_j = 1 + eps xi_j turn the chain into the generalized eigenproblem L v = omega^2 M v,
M = diag(m_j). A unit velocity impulse at node 0 gives the node-0 response
    x_0(t) = sum_r c_r sin(omega_r t),   c_r = v_r(0)^2 / omega_r > 0,
so the incoherent ceiling C = sum_r c_r is reached iff all sines peak together (rational independence of
the omega_r). For the ORDERED ring with composite N the frequencies carry an exact rational relation
(the mod-4 obstruction), so the achievable sup A sits strictly below C. Generic masses destroy that
relation almost surely, hence the disordered chain CAN reach its ceiling for any N.

The reviewer's paradox, made precise in three checks:
  (A) the exact resonance residual jumps from 0 to O(eps): the obstruction is destroyed by any disorder;
  (B) over a long horizon the achievable efficiency A/C rises above the ordered (blocked) value;
  (C) over a SHORT horizon there is no gain -- the now-independent frequencies need Diophantine time to
      align, so finite-time amplification is NOT improved. Infinite-time helps; finite-time does not."""

from __future__ import annotations

import numpy as np


def ring_laplacian(n: int) -> np.ndarray:
    lap = 2.0 * np.eye(n) - np.eye(n, k=1) - np.eye(n, k=-1)
    lap[0, -1] = lap[-1, 0] = -1.0
    return lap


def modal_response(n: int, masses: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (omega_r, c_r) for the node-0 velocity response with mass profile `masses` (springs = 1)."""
    inv_sqrt_m = 1.0 / np.sqrt(masses)
    b = inv_sqrt_m[:, None] * ring_laplacian(n) * inv_sqrt_m[None, :]
    w2, wvec = np.linalg.eigh(b)
    keep = w2 > 1e-9
    omega = np.sqrt(w2[keep])
    v0 = wvec[0, keep] * inv_sqrt_m[0]
    return omega, v0**2 / omega


def efficiency(omega: np.ndarray, c: np.ndarray, t_max: float, dt: float = 0.05) -> float:
    t = np.arange(0.0, t_max, dt)
    return float((np.sin(np.outer(t, omega)) * c).sum(axis=1).max() / c.sum())


def main() -> int:
    rng = np.random.default_rng(0)
    print("=" * 80)
    print("Disorder breaks resonances: composite N reaches its ceiling over long (not short) horizons")
    print("=" * 80)
    ok = True

    print("\n(A) The exact resonance residual jumps 0 -> O(eps): obstruction destroyed by any disorder (N=12)")
    print("    C_12 carries the exact relation 2 omega_2 - omega_6 = 0 (modes lambda=1 and lambda=4);")
    print("    residual = |2 f1 - f2|, f1 = freq nearest 1, f2 = freq nearest 2")

    def residual(om: np.ndarray) -> float:
        f1 = om[np.argmin(np.abs(om - 1.0))]
        f2 = om[np.argmin(np.abs(om - 2.0))]
        return abs(2.0 * f1 - f2)

    print(f"    {'eps':>7} {'mean residual (20 draws)':>26}")
    res0 = residual(modal_response(12, np.ones(12))[0])
    print(f"    {0.0:>7.2f} {res0:>26.2e}")
    prev, mono = res0, True
    for eps in (0.02, 0.05, 0.10, 0.20):
        mr = float(np.mean([residual(modal_response(12, 1.0 + eps * rng.uniform(-1, 1, 12))[0]) for _ in range(20)]))
        mono &= mr > prev
        prev = mr
        print(f"    {eps:>7.2f} {mr:>26.2e}")
    ok &= mono and res0 < 1e-9
    print(f"    -> residual is 0 when ordered and grows ~linearly with eps: {'OK' if mono and res0 < 1e-9 else 'FAIL'}")

    print("\n(B) Long horizon (T=3e4): disorder lifts efficiency A/C above the blocked ordered value")
    print(f"    {'N':>4} {'ordered A/C':>12} {'disordered A/C (eps=0.4, mean 8)':>34}")
    for n in (9, 12, 15):
        eff_ord = efficiency(*modal_response(n, np.ones(n)), 30000.0)
        eff_dis = float(np.mean([efficiency(*modal_response(n, 1.0 + 0.4 * rng.uniform(-1, 1, n)), 30000.0)
                                 for _ in range(8)]))
        gain = eff_dis > eff_ord + 0.02
        ok &= gain
        print(f"    {n:>4} {eff_ord:>12.4f} {eff_dis:>28.4f}      ({'OK' if gain else 'FAIL'})")
    print("    -> generic masses raise the achievable amplitude toward the full incoherent ceiling: OK")

    print("\n(C) Short horizon (T=300): the same disorder gives NO gain -- the cost is recurrence time (N=15)")
    eff_ord_s = efficiency(*modal_response(15, np.ones(15)), 300.0)
    eff_dis_s = float(np.mean([efficiency(*modal_response(15, 1.0 + 0.4 * rng.uniform(-1, 1, 15)), 300.0)
                               for _ in range(8)]))
    eff_ord_l = efficiency(*modal_response(15, np.ones(15)), 30000.0)
    eff_dis_l = float(np.mean([efficiency(*modal_response(15, 1.0 + 0.4 * rng.uniform(-1, 1, 15)), 30000.0)
                               for _ in range(8)]))
    no_short_gain = eff_dis_s <= eff_ord_s + 0.02
    yes_long_gain = eff_dis_l > eff_ord_l + 0.02
    ok &= no_short_gain and yes_long_gain
    print(f"    short T=300:  ordered {eff_ord_s:.4f}  vs disordered {eff_dis_s:.4f}  (gain {eff_dis_s - eff_ord_s:+.4f})")
    print(f"    long  T=3e4:  ordered {eff_ord_l:.4f}  vs disordered {eff_dis_l:.4f}  (gain {eff_dis_l - eff_ord_l:+.4f})")
    print(f"    -> the disorder advantage exists only at long T (recurrence cost): "
          f"{'OK' if no_short_gain and yes_long_gain else 'FAIL'}")

    print("\n(D) Physical reading: an infinitesimal perturbation makes the absolute ceiling attainable in")
    print("    principle, yet only over Diophantine-long time. The clean arithmetic of prime / 2^m N is")
    print("    what makes a large wave reachable QUICKLY -- the engineering-relevant statement.")

    print("\n" + "=" * 80)
    print("RESULT:", "DISORDER EXPERIMENT VERIFIED" if ok else "CHECK FAILED")
    print("=" * 80)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
