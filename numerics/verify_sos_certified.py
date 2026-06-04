# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2"]
# ///
"""Certified bounds on the composite-N subtorus maximum (the gap A_N < U_N), rigorously.

For composite N the large wave A_N = max over the orbit-closure subtorus of a trigonometric
polynomial sits STRICTLY below the Bohr ceiling U_N, because rational relations among the frequencies
forbid all sines from peaking at once. We certify this gap rigorously, with no SDP solver, by a
Lipschitz-grid certificate on the low-dimensional subtorus, and contrast a prime N (no gap).

N = 5 (prime, Q-rank 2): independent phases, g_5 = a_1 sin x + a_2 sin y; max = a_1 + a_2 = U_5 (no gap).
N = 6 (composite, Q-rank 2): omega = {1, sqrt3, 2} with the relation omega_3 = 2 omega_1, so the
   third phase is locked to 2x:  g_6 = a_1 sin x + a_3 sin 2x + a_2 sin y. The sin 2x term cannot reach
   +a_3 while sin x reaches +1, forcing max g_6 < a_1 + a_2 + a_3 = U_6.

A Lipschitz certificate on a grid of centers gives  max <= gridmax + (L_x h_x + L_y h_y)/2  (rigorous
upper bound) and  max >= gridmax (an evaluation), bracketing the true value. This is the elementary
member of the Lasserre moment--SOS hierarchy, which yields tighter certified bounds at higher relaxation
order (an SDP; full rigor needs rational rounding of the dual SOS certificate).
"""

from __future__ import annotations

import numpy as np


def ceiling(n: int) -> float:
    r = np.arange(1, n)
    return float(np.sum(1.0 / np.sin(np.pi * r / n)) / (2 * n))


def grid_certify(g, lx: float, ly: float, m: int = 2400) -> tuple[float, float]:
    """Return (lower, upper): lower = max over a grid of cell centers; upper = lower + Lipschitz slack."""
    h = 2 * np.pi / m
    c = (np.arange(m) + 0.5) * h  # cell centers in [0, 2pi)
    gridmax = -np.inf
    for x in c:  # row-blocked to bound memory
        vals = g(x, c)
        gridmax = max(gridmax, float(vals.max()))
    slack = (lx * h + ly * h) / 2.0  # max deviation from center within a cell (per-axis Lipschitz)
    return gridmax, gridmax + slack


def main() -> int:
    print("=" * 70)
    print("Certified bounds on the subtorus maximum: the composite gap A_N < U_N")
    print("=" * 70)
    ok = True

    # ---- N = 5 (prime): no gap, A_5 = U_5 ----
    o1, o2 = 2 * np.sin(np.pi / 5), 2 * np.sin(2 * np.pi / 5)
    a1, a2 = 2 / (5 * o1), 2 / (5 * o2)  # weight 2/N each (degeneracy r and N-r)

    def g5(x, y):
        return a1 * np.sin(x) + a2 * np.sin(y)

    lo5, hi5 = grid_certify(g5, a1, a2)
    u5 = ceiling(5)
    print("\nN=5 (prime, Q-rank 2): g_5 = a1 sin x + a2 sin y")
    print(f"    certified A_5 in [{lo5:.5f}, {hi5:.5f}];  U_5 = {u5:.5f};  a1+a2 = {a1 + a2:.5f}")
    n5_ok = hi5 >= a1 + a2 - 1e-3 and abs(u5 - (a1 + a2)) < 1e-9 and lo5 > (a1 + a2) - 5e-3
    ok &= n5_ok
    print(f"    -> max = a1+a2 = U_5, NO gap (sharp on primes): {'OK' if n5_ok else 'FAIL'}")

    # ---- N = 6 (composite): certified gap A_6 < U_6 ----
    b1, b2, b3 = 2 / (6 * 1.0), 2 / (6 * np.sqrt(3.0)), 1 / (6 * 2.0)  # omega = 1, sqrt3, 2

    def g6(x, y):
        return b1 * np.sin(x) + b3 * np.sin(2 * x) + b2 * np.sin(y)

    lo6, hi6 = grid_certify(g6, b1 + 2 * b3, b2)  # L_x = b1 + 2 b3, L_y = b2
    u6 = ceiling(6)
    trivial = b1 + b2 + b3  # = U_6 contribution at node 0
    print("\nN=6 (composite, Q-rank 2): g_6 = a1 sin x + a3 sin 2x + a2 sin y  (omega_3 = 2 omega_1)")
    print(f"    certified A_6 in [{lo6:.5f}, {hi6:.5f}];  U_6 = {u6:.5f};  trivial sum = {trivial:.5f}")
    print(f"    certified gap  U_6 - A_6 >= {u6 - hi6:.5f}  (> 0: composite value strictly below ceiling)")
    n6_ok = (u6 - hi6) > 0.02 and (hi6 - lo6) < 5e-3 and abs(lo6 - 0.5594) < 5e-3
    ok &= n6_ok
    print(f"    -> certified strict gap below U_6, bracket width {hi6 - lo6:.1e}: "
          f"{'OK' if n6_ok else 'FAIL'}")

    print("\n" + "=" * 70)
    print("RESULT:", "CERTIFIED BOUNDS VERIFIED" if ok else "CHECK FAILED")
    print("Prime N=5: A_5=U_5 (no gap). Composite N=6: A_6 certified strictly below U_6 (Lipschitz grid,")
    print("rigorous). The Lasserre moment-SOS hierarchy scales this to higher dimension and order.")
    print("=" * 70)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
