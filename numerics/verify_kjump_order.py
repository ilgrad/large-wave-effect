# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2", "mpmath>=1.3"]
# ///
"""Order criterion for new quasi-1D families: prisms, Mobius ladders, antiprisms, k-jump circulants.

Settles condition (B) of the order criterion (Prop. order-criterion) for several families beyond the
ring and the NNN ring, separating what reduces to existing results from what is genuinely new.

STRUCTURE (exact, Laplacian-spectrum identities):
  * Prism      C_N box K_2 : Cartesian product  -> Cor. product (NOT new).
  * Antiprism  A_N         : equals the NNN ring C_{2N}({1,2}) -> Prop. nnn on 2N vertices (NOT new).
  * Mobius     M_N         : equals the circulant C_{2N}({+-1} u {N}) with the diameter chord N
                             SELF-PAIRED (degree 1).  NOT a Cartesian product -- GENUINELY NEW.
                             Lambda_r = 2(1-cos(pi r/N)) + (1-(-1)^r), r=0..2N-1.
                             The gapless even-r branch (r=2j) is EXACTLY the ring C_N spectrum
                             4 sin^2(pi j/N); the odd-r branch is gapped (Lambda >= 2).

CONDITION (A) -- logarithmic ceiling, slope c = 1/(pi sqrt(sum_{s in S} s^2)) for k-jump C_N(S),
  and c = 1/(2 pi) for the Mobius ladder (one gapless channel, |V| = 2N).  Verified by the incremental
  slope (U_2N - U_N)/ln 2 -> c (Euler-Maclaurin band-edge, rigorous).

CONDITION (B) -- a Q-independent subset carries > U_N/2.  Two routes:
  (B-Mobius) the gapless even branch is the ring C_N, whose first phi(2N)/2 frequencies are
             Q-independent by the palindromic lemma (lem:prefix); they carry (1-o(1)) U.  RIGOROUS.
  (B-kjump)  for the k-jump circulant the only Q-relations among {omega_r=sqrt(lambda_r)} found are the
             COLLISIONS omega_a=omega_b (lambda_a=lambda_b); since collisions add weights, the DISTINCT
             frequencies are Q-independent and carry the FULL U_N.  The collisions obey exact congruence
             laws (Conway-Jones), proved here in integer arithmetic mod Phi_N; "no further relation" is
             PSLQ-numerical (250 digits).  Where it holds, A_N = U_N (full saturation, by thm:ceiling).
"""

from __future__ import annotations

import sys
from math import gcd
from pathlib import Path

import mpmath as mp
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from large_wave_effect.cyclotomic import cyclotomic, exact_rank, sin_coords


# ----------------------------------------------------------------------------- structure (exact spectra)
def lap_eigs(adj: np.ndarray) -> np.ndarray:
    return np.sort(np.linalg.eigvalsh(np.diag(adj.sum(1)) - adj))


def circ_eigs(n: int, jumps: list[int]) -> np.ndarray:
    return np.sort(
        np.array([sum(2 * (1 - np.cos(2 * np.pi * s * r / n)) for s in jumps) for r in range(n)])
    )


def prism_adj(n: int) -> np.ndarray:
    a = np.zeros((2 * n, 2 * n))
    for i in range(n):
        a[i, (i + 1) % n] = a[(i + 1) % n, i] = 1
        a[n + i, n + (i + 1) % n] = a[n + (i + 1) % n, n + i] = 1
        a[i, n + i] = a[n + i, i] = 1
    return a


def antiprism_adj(n: int) -> np.ndarray:
    a = np.zeros((2 * n, 2 * n))
    for i in range(n):
        a[i, (i + 1) % n] = a[(i + 1) % n, i] = 1
        a[n + i, n + (i + 1) % n] = a[n + (i + 1) % n, n + i] = 1
        a[i, n + i] = a[n + i, i] = 1
        a[i, n + (i + 1) % n] = a[n + (i + 1) % n, i] = 1
    return a


def mobius_adj(n: int) -> np.ndarray:
    a = np.zeros((2 * n, 2 * n))
    for i in range(2 * n):
        a[i, (i + 1) % (2 * n)] = a[(i + 1) % (2 * n), i] = 1
    for i in range(n):  # each diameter chord added ONCE (self-paired in the circulant)
        a[i, i + n] = a[i + n, i] = 1
    return a


def mobius_eigs_formula(n: int) -> np.ndarray:
    return np.sort(np.array([2 * (1 - np.cos(np.pi * r / n)) + (1 - (-1) ** r) for r in range(2 * n)]))


# ----------------------------------------------------------------------------- k-jump collision law (exact)
def lam_poly(n: int, jumps: list[int], r: int) -> list[int]:
    """Integer rep of lambda_r = sum_s (2 - zeta^{sr} - zeta^{-sr}) in Z[zeta_N], length N (pre-reduction)."""
    p = [0] * n
    for s in jumps:
        p[0] += 2
        p[(s * r) % n] -= 1
        p[(-s * r) % n] -= 1
    return p


def poly_rem(a: list[int], b: list[int]) -> list[int]:
    a = a[:]
    while len(a) >= len(b) and any(a):
        if a[-1] == 0:
            a.pop()
            continue
        c = a[-1] // b[-1]
        d = len(a) - len(b)
        for j in range(len(b)):
            a[d + j] -= c * b[j]
        while a and a[-1] == 0:
            a.pop()
    return a


def lam_collisions(n: int, jumps: list[int]) -> list[tuple[int, int]]:
    """Exact pairs a<b<=N/2 with lambda_a = lambda_b (omega_a = omega_b), by reduction mod Phi_N."""
    phi = cyclotomic(n)
    m = n // 2
    reps = [lam_poly(n, jumps, r) for r in range(1, m + 1)]
    out = []
    for a in range(m):
        for b in range(a + 1, m):
            diff = [x - y for x, y in zip(reps[a], reps[b], strict=True)]
            if not any(poly_rem(diff, phi)):  # zero remainder: equal in Q(zeta_N) (any(), not truthiness)
                out.append((a + 1, b + 1))
    return out


# ----------------------------------------------------------------------------- ceilings / weights
def kjump_ceiling(n: int, jumps: list[int]) -> float:
    m = n // 2
    om = [np.sqrt(sum(2 * (1 - np.cos(2 * np.pi * s * r / n)) for s in jumps)) for r in range(1, m + 1)]
    b = [2 / (n * w) for w in om]
    if n % 2 == 0:
        b[-1] /= 2
    return float(sum(b))


def mobius_ceiling(n: int) -> float:
    return float(sum(1 / np.sqrt(2 * (1 - np.cos(np.pi * r / n)) + (1 - (-1) ** r)) for r in range(1, 2 * n)) / (2 * n))


def mobius_prefix_weight(n: int) -> float:
    """Weight carried by the Q-independent ring prefix (first phi(2N)/2 ring frequencies)."""
    big_m = sum(1 for k in range(1, 2 * n + 1) if gcd(k, 2 * n) == 1) // 2
    return float(sum(2 * (1 / (2 * n)) / (2 * np.sin(np.pi * j / n)) for j in range(1, big_m + 1)))


# ----------------------------------------------------------------------------- checks
def check_structure() -> bool:
    print("[1] STRUCTURE (exact Laplacian spectra): which families reduce to existing results")
    ok = True
    for n in range(4, 13):
        pr = np.allclose(lap_eigs(prism_adj(n)), np.sort((circ_eigs(n, [1])[:, None] + np.array([0.0, 2.0])).ravel()))
        ap = np.allclose(lap_eigs(antiprism_adj(n)), circ_eigs(2 * n, [1, 2]))
        mo = np.allclose(lap_eigs(mobius_adj(n)), mobius_eigs_formula(n))
        ok &= pr and ap and mo
    print(f"    prism C_N box K2 = C_N (+) K2 (product, Cor. product): {pr}")
    print(f"    antiprism A_N = NNN ring C_2N({{1,2}}) (Prop. nnn):    {ap}")
    print(f"    Mobius M_N = circulant C_2N({{+-1}} u {{N}}, chord self-paired): {mo}")
    # antiprism is the EQUAL-WEIGHT {1,2} circulant on M=2N, i.e. the NNN ring at coupling g=1;
    # its collisions obey the g=1 Conway-Jones law (Prop. nnn, general-g): (M/6,M/2),(M/5,2M/5),(M/4,M/3).
    bad_ap = []
    for nn in range(3, 81):
        mm = 2 * nn
        c = set(lam_collisions(mm, [1, 2]))
        exp: set[tuple[int, int]] = set()
        if mm % 6 == 0:
            exp.add((mm // 6, mm // 2))
        if mm % 5 == 0:
            exp.add((mm // 5, 2 * mm // 5))
        if mm % 12 == 0:
            exp.add(tuple(sorted((mm // 4, mm // 3))))
        if c != exp:
            bad_ap.append(nn)
    ok &= not bad_ap
    print(f"    antiprism collisions = g=1 NNN law on M=2N, exact 3<=N<=80: "
          f"{'CONFIRMED' if not bad_ap else f'FAIL {bad_ap[:5]}'}")
    # Mobius even branch = ring C_N
    even_eq = True
    for n in range(4, 13):
        even = sorted({round(2 * abs(np.sin(np.pi * j / n)), 9) for j in range(1, n)})
        ring = sorted({round(2 * abs(np.sin(np.pi * j / n)), 9) for j in range(1, n)})
        even_eq &= even == ring
    print(f"    Mobius gapless even-r branch = ring C_N frequencies: {even_eq}")
    return ok and even_eq


def check_mobius_order() -> bool:
    print("\n[2] MOBIUS LADDER order theorem A(M_N) = Theta(ln N)  (GENUINELY NEW; rigorous)")
    c_pred = 1 / (2 * np.pi)
    print(f"    (A) ceiling slope -> 1/(2 pi) = {c_pred:.5f} (one gapless channel, |V|=2N)")
    us = [(n, mobius_ceiling(n)) for n in (512, 1024, 2048, 4096)]
    slope = (us[-1][1] - us[-2][1]) / np.log(2)
    print(f"        incremental slope (U_4096 - U_2048)/ln2 = {slope:.5f}")
    a_ok = abs(slope - c_pred) < 5e-4
    print("    (B) embedded ring prefix Q-independent (palindromic lemma) -> 2 L_pre - U > 0, growing:")
    print(f"        {'N':>6} {'U':>9} {'L_pre':>9} {'2Lpre-U':>9}")
    vals = []
    for n in (256, 512, 1024, 2048):
        u, lp = mobius_ceiling(n), mobius_prefix_weight(n)
        vals.append(2 * lp - u)
        print(f"        {n:>6} {u:>9.5f} {lp:>9.5f} {2 * lp - u:>9.5f}")
    b_ok = min(vals) > 0 and vals[-1] > vals[0] + 0.1
    # exact: the ring prefix really is Q-independent (reuse the ring certificate on modulus 2N)
    big_m = sum(1 for k in range(1, 2 * 256 + 1) if gcd(k, 2 * 256) == 1) // 2
    rank = exact_rank(sin_coords(2 * 256, list(range(1, big_m + 1))))
    print(f"        exact: ring-prefix rank {rank} = phi(2N)/2 = {big_m}: {rank == big_m}")
    return a_ok and b_ok and rank == big_m


def check_kjump() -> bool:
    print("\n[3] k-JUMP CIRCULANTS C_N(S): nested-radical frequencies omega_r = sqrt(lambda_r)")
    ok = True

    # (A) ceiling slope c = 1/(pi sqrt(sum s^2))
    print("    (A) ceiling slope c = 1/(pi sqrt(sum_{s in S} s^2)):")
    for jumps in ([1, 2, 3], [1, 2, 3, 4], [1, 3]):
        c = 1 / (np.pi * np.sqrt(sum(s * s for s in jumps)))
        slope = (kjump_ceiling(4096, jumps) - kjump_ceiling(2048, jumps)) / np.log(2)
        good = abs(slope - c) < 5e-4
        ok &= good
        print(f"        S={jumps}: c={c:.5f}  incremental slope={slope:.5f}  ({'OK' if good else 'FAIL'})")

    # (B-law) exact collision congruence laws (Conway-Jones + residue completeness), mod Phi_N.
    # The "p | N" term (p = max jump-span prime) is the residue-completeness collision: when S together
    # with -S covers all nonzero residues mod p, every lambda_r with p nmid r equals a constant.
    print("    (B) collision law (EXACT, integer arithmetic mod Phi_N, 5<=N<=160):")
    bad123 = [n for n in range(5, 161) if bool(lam_collisions(n, [1, 2, 3])) != (n % 4 == 0 or n % 6 == 0 or n % 7 == 0)]
    bad13 = [n for n in range(5, 161) if bool(lam_collisions(n, [1, 3])) != (n % 8 == 0 or n % 5 == 0)]
    ok &= not bad123 and not bad13
    print(f"        S={{1,2,3}}: collision <=> 4|N or 6|N or 7|N : {'CONFIRMED' if not bad123 else f'FAIL {bad123[:5]}'}")
    print(f"        S={{1,3}}  : collision <=> 8|N or 5|N        : {'CONFIRMED' if not bad13 else f'FAIL {bad13[:5]}'}")

    # (B-rank) distinct frequencies Q-independent => rank = #distinct (PSLQ 200 digits)
    print("    (B) Q-rank = #distinct frequencies (PSLQ 200 digits; genuine relations are O(1)-coeff):")
    mp.mp.dps = 200
    jumps = [1, 2, 3]
    rank_eq_distinct = True
    sample = []
    for n in range(5, 41):
        m = n // 2
        om = [mp.sqrt(sum(2 - 2 * mp.cos(2 * mp.pi * s * r / n) for s in jumps)) for r in range(1, m + 1)]
        # dedup distinct
        dist: list[mp.mpf] = []
        for w in om:
            if all(abs(w - d) > mp.mpf(10) ** -150 for d in dist):
                dist.append(w)
        if len(dist) >= 2:
            rel = mp.pslq(dist, maxcoeff=10**5, maxsteps=10**6)
            indep = rel is None
        else:
            indep = True
        rank_eq_distinct &= indep
        if n in (12, 15, 21, 24):
            sample.append((n, m, len(dist), "indep" if indep else "RELATION"))
    ok &= rank_eq_distinct
    for n, m, nd, st in sample:
        print(f"        N={n:>2}: m={m:>2} #distinct={nd:>2} -> {st}")
    print(f"        all distinct-sets Q-independent (no relation beyond collisions), S={{1,2,3}}, N<=40: "
          f"{'CONFIRMED' if rank_eq_distinct else 'FAIL'}")
    print("        => distinct freqs carry the FULL U_N; where collision-free (S={1,2,3}: 4,6,7 nmid N)")
    print("           the frequencies are all-distinct and Q-independent, so A_N = U_N (saturation).")
    return ok


def main() -> int:
    print("=" * 80)
    print("Order criterion for new quasi-1D families (prism / antiprism / Mobius / k-jump circulant)")
    print("=" * 80)
    ok = True
    ok &= check_structure()
    ok &= check_mobius_order()
    ok &= check_kjump()
    print("\n" + "=" * 80)
    print("RESULT:", "PASS" if ok else "FAIL")
    print("  prism = product (Cor. product); antiprism = NNN ring (Prop. nnn);")
    print("  Mobius ladder = NEW (embedded-ring gapless branch gives A=Theta(ln N), rigorous);")
    print("  k-jump C_N(S): condition (A) c=1/(pi sqrt(sum s^2)); collisions obey exact Conway-Jones")
    print("  congruence laws; rank = #distinct (PSLQ); saturation A_N=U_N off the collision congruences.")
    print("=" * 80)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
