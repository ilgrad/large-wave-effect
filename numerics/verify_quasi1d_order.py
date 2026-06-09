# /// script
# requires-python = ">=3.11"
# dependencies = ["mpmath>=1.3"]
# ///
"""Order criterion (Proposition order-criterion) on a concrete quasi-1D lattice: the NNN ring.

The order theorem's mechanism is graph-agnostic: A_N = Theta(ln N) for any d_s=1 family whose ceiling is
~ c ln N (A) and whose frequencies have a Q-independent subset carrying > 1/2 the ceiling (B).  We
illustrate on the next-nearest-neighbour ring with Laplacian symbol
    lambda_r = 2(1-cos t) + 2g(1-cos 2t),  t = 2 pi r/N,  g = 1/2,
so omega_r = 2 sin(pi r/N) sqrt(2 + cos(2 pi r/N)).

(A) ceiling U_N = sum_{r=1}^{floor N/2} 2/(N omega_r) ~ c ln N with c = 1/(pi sqrt(1+4g)) = 1/(pi sqrt3).
(B) is per-family: at N=15 the omega_r are Q-independent (PSLQ finds no relation) -> the NNN ring SATURATES
    (A_N = U_N), whereas at N=12 there is the degeneracy omega_3 = omega_6 = 2.  So (B) must be decided
    case by case -- exactly what the criterion isolates.
"""

from __future__ import annotations

import mpmath as mp

mp.mp.dps = 80
G = mp.mpf(1) / 2  # next-nearest coupling


def omega(n: int) -> list:
    return [2 * mp.sin(mp.pi * r / n) * mp.sqrt(2 + mp.cos(2 * mp.pi * r / n)) for r in range(1, n // 2 + 1)]


def ceiling(n: int) -> mp.mpf:
    om = omega(n)
    b = [2 / (n * w) for w in om]
    if n % 2 == 0:
        b[-1] /= 2
    return mp.fsum(b)


def main() -> int:
    c_pred = 1 / (mp.pi * mp.sqrt(1 + 4 * G))
    print(f"(A) NNN-ring ceiling U_N ~ c ln N,  c = 1/(pi sqrt3) = {mp.nstr(c_pred, 6)}")
    # U_N/lnN = c + C/lnN converges slowly; test the incremental slope (U_2N - U_N)/ln2 -> c instead.
    ns = [128, 256, 512, 1024]
    us = [ceiling(n) for n in ns]
    for n, u in zip(ns, us, strict=True):
        print(f"    N={n:>5}  U_N={mp.nstr(u, 8)}  U_N/lnN={mp.nstr(u / mp.log(n), 6)}")
    slope = float((us[-1] - us[-2]) / mp.log(2))
    print(f"    incremental slope (U_1024-U_512)/ln2 = {slope:.5f}  (c = {float(c_pred):.5f})")
    a_ok = abs(slope - float(c_pred)) < 0.005

    # genuine relations have small coefficients (the pure ring's is [1,0,0,1,0,-1,0]); cap maxcoeff to
    # reject PSLQ's spurious large-coefficient near-relations.
    print("(B) Q-independence of {omega_r} is per-family (PSLQ, 80 digits, maxcoeff 1000):")
    rel15 = mp.pslq(omega(15), maxcoeff=1000, maxsteps=10**6)
    rel12 = mp.pslq(omega(12), maxcoeff=1000, maxsteps=10**6)
    print(f"    N=15: relation = {rel15}   (None => full rank => SATURATES, A_N=U_N)")
    print(f"    N=12: relation = {rel12}   (= omega_3 - omega_6, a degeneracy)")
    b_ok = rel15 is None and rel12 is not None

    ok = a_ok and b_ok
    print("=" * 64)
    print("RESULT:", "PASS -- d_s=1 ceiling c=1/(pi sqrt3); (B) per-family (N=15 vs N=12)" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
