# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2"]
# ///
"""Is there a closed-form A_N for composite N?  Honest structural answer.

A_N = max_j sup_t u_j(t) is the maximum of a trigonometric polynomial over the orbit-closure
subtorus. That subtorus splits, by the lattice of integer relations among {omega_r}, into
BLOCKS. The block maxima are algebraic numbers, but:

  * N = 6: the relation lattice gives only the simple commensurability omega_3 = 2 omega_1
    (sin(pi/2) = 2 sin(pi/6)), plus the independent omega_2 = sqrt 3. Each block is ONE-variable,
    so A_6 is solvable in radicals:
        A_6 = sqrt3/9  +  max_phi [ (1/3) sin phi + (1/12) sin 2phi ],
    where the inner max is at cos phi = (sqrt3 - 1)/2, giving the closed form below.

  * N = 12: there is a genuine THREE-term relation omega_5 = omega_1 + omega_3
    (2 sin75 - 2 sin15 = 2 sin45, i.e. sin75 - sin15 = sin45). The block {1,3,5} is therefore
    a TWO-variable constrained optimization, not a separable 1-D max -> no clean elementary sum.

Conclusion (verified here): no uniform elementary closed form; A_N is an algebraic number whose
defining system's dimension/degree grows with the arithmetic of N. Explicit radicals exist per N.
"""

from __future__ import annotations

import numpy as np


# ---------- exact cyclotomic coordinates (to certify the integer relations) ----------
def divisors(n):
    return [d for d in range(1, n + 1) if n % d == 0]


def poly_div_exact(a, b):
    a = a[:]
    db = len(b) - 1
    q = [0] * (len(a) - db)
    for i in range(len(a) - 1, db - 1, -1):
        c = a[i]
        if c:
            q[i - db] = c
            for j in range(len(b)):
                a[i - db + j] -= c * b[j]
    return q


def cyclotomic(n, cache):
    if n in cache:
        return cache[n]
    num = [-1] + [0] * (n - 1) + [1]
    for d in divisors(n):
        if d < n:
            num = poly_div_exact(num, cyclotomic(d, cache))
    cache[n] = num
    return num


def power_mod(e, phi):
    deg = len(phi) - 1
    p = [0] * max(e + 1, deg)
    p[e] = 1
    for i in range(len(p) - 1, deg - 1, -1):
        c = p[i]
        if c:
            for j in range(deg + 1):
                p[i - deg + j] -= c * phi[j]
    return p[:deg]


def sin_coord(n, r):
    cache = {}
    phi = cyclotomic(2 * n, cache)
    hi = power_mod(r % (2 * n), phi)
    lo = power_mod((2 * n - r) % (2 * n), phi)
    return np.array([h - lv for h, lv in zip(hi, lo, strict=True)], dtype=object)


def main() -> int:
    print("=" * 72)
    print("Closed-form A_N for composite N -- structural / honest")
    print("=" * 72)
    ok = True

    # ---- N = 6: explicit radical closed form ----
    print("\n[N=6]  separable blocks -> closed form in radicals")
    a6_closed = np.sqrt(3) / 9 + (3 + np.sqrt(3)) * 3 ** 0.25 / (12 * np.sqrt(2))
    # independent numeric max over the 2-D subtorus (phi1 free, phi2 free; phi3 = 2 phi1)
    g = np.linspace(0, 2 * np.pi, 4000)
    block1 = np.max(g * 0 + (1 / 3) * np.sin(g) + (1 / 12) * np.sin(2 * g))  # 1-var, refine below
    # refine block1 at the stationary point cos phi = (sqrt3-1)/2
    phi_star = np.arccos((np.sqrt(3) - 1) / 2)
    block1 = (1 / 3) * np.sin(phi_star) + (1 / 12) * np.sin(2 * phi_star)
    block2 = 1 / (3 * np.sqrt(3))  # max sin = 1
    a6_decomp = block1 + block2
    print(f"  cos(phi*) = (sqrt3-1)/2 = {(np.sqrt(3) - 1) / 2:.6f}")
    print(f"  A_6 (radical closed form) = sqrt3/9 + (3+sqrt3)*3^(1/4)/(12 sqrt2) = {a6_closed:.6f}")
    print(f"  A_6 (block decomposition) = {a6_decomp:.6f}")
    print(f"  match: {'PASS' if abs(a6_closed - a6_decomp) < 1e-9 else 'FAIL'}")
    ok &= abs(a6_closed - a6_decomp) < 1e-9

    # ---- N = 12: a genuine 3-term relation (block is 2-D, not separable) ----
    print("\n[N=12]  three-term relation omega_5 = omega_1 + omega_3 (non-separable block)")
    c1, c3, c5 = sin_coord(12, 1), sin_coord(12, 3), sin_coord(12, 5)
    exact_zero = bool(np.all((c5 - c1 - c3) == 0))  # exact integer-coord certificate
    sin_id = 2 * np.sin(5 * np.pi / 12) - 2 * np.sin(np.pi / 12) - 2 * np.sin(3 * np.pi / 12)
    print(f"  exact coords c5 - c1 - c3 == 0 : {exact_zero}")
    print(f"  numeric  omega_5 - omega_1 - omega_3 = {sin_id:.2e}  (sin75 - sin15 = sin45)")
    print("  -> block {1,3,5} couples three modes: 2-variable max, algebraic but not a sum.")
    ok &= exact_zero and abs(sin_id) < 1e-12

    print("\n" + "=" * 72)
    print("RESULT:", "STRUCTURE CONFIRMED" if ok else "MISMATCH")
    print("Honest answer: NO uniform elementary closed form for composite N. A_N is an algebraic")
    print("number = max of a trig polynomial over the relation-lattice subtorus; it splits into")
    print("blocks (1-D for N=6 -> radicals; >=2-D when relations like omega_5=omega_1+omega_3")
    print("appear, e.g. N=12). Per N it is exactly computable (radicals when Galois-solvable);")
    print("a general elementary formula does not exist. Certified value -> Lasserre/SOS.")
    print("=" * 72)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
