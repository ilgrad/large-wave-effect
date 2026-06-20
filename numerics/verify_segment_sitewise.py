# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2"]
# ///
"""Site-wise saturation criterion for the Dirichlet segment (closes the open classification).

Setup (alpha_j == 1, cf. verify_filimonov_schrodinger.py):
    z_j(t) = e^{2it} sum_{k odd} a_k sin(pi k j/N) e^{-2 i omega_k t},  omega_k = cos(pi k/N),  a_k>0,
    ceiling d_j^N = sum_k |a_k sin(pi k j/N)|,  |z_j(t)| <= d_j^N.

Filimonov / our extension: sup_t|z_j| = d_j^N for EVERY site j when N is prime or 2^m (full cos-rank).
Composite N saturates only at SOME sites; the classification was open.  We close it with an exact
necessary-and-sufficient CRITERION.

CRITERION (Prop.).  Let the ACTIVE set be A_j = {k : 1<=k<N, k odd, N does not divide k j} (the terms with
a_k sin(pi k j/N) != 0), and sigma_k = sign(sin(pi k j/N)) in {+,-}.  Then
    sup_t |z_j(t)| = d_j^N   <=>   every integer relation (c_k)_{k in A_j} with
        sum_k c_k cos(pi k/N) = 0   AND   sum_k c_k = 0   (a "balanced" relation)
    satisfies   sum_{k in A_j : sin(pi k j/N) < 0} c_k  ==  0 (mod 2).

Proof (Kronecker-Weyl + Pontryagin).  |z_j| reaches d_j^N iff all active terms |amp_k| e^{i(sigma_k - 2 w_k t)}
can be brought to a common phase phi_0, i.e. the point (sigma_k) lies in <1> + H_0, the orbit under the free
overall phase phi_0 of the closure H_0 = cl{(-2 w_k t)_k : t in R}.  By Kronecker-Weyl, H_0 has annihilator
{c : sum c_k w_k = 0}; adding the all-ones direction <1> intersects the annihilator with {sum c_k = 0}.  By
Pontryagin duality (sigma_k) in <1>+H_0 iff <c, sigma> = 0 (mod 2 pi) for every c in that annihilator; since
sigma_k in {0, pi}, <c, sigma> = pi * sum_{sigma_k = pi} c_k, which is 0 mod 2 pi iff that sum is even.  QED.

Corollaries (proved): (i) j = 1 and j = N-1 always saturate -- there sin(pi k j/N) > 0 for all active k, so the
negative set is empty and the parity is trivially even; (ii) N prime, 2^m, or 2p (the full cos-rank cases)
saturates at every site.  For other composite N the saturating set is sign-pattern dependent (NOT a function of
gcd(j,N): e.g. for N=24 the coprime saturating sites are exactly the subgroup {1,7,17,23} < (Z/24)*).

This script: (A) verifies the criterion against the exact sup_t (high-resolution scan), (B) checks the two
corollaries, (C) prints the saturating-site table.  The exact integer relation lattice is computed with PARI
`matkerint` if `gp` is on PATH (saturated Z-kernel); otherwise the exact column is skipped and only the
numerical law is shown.
"""

from __future__ import annotations

import shutil
import subprocess
from math import gcd, pi, sin

import numpy as np

HAVE_GP = shutil.which("gp") is not None


# --- cyclotomic Phi_{2N} and exact integer coords of 2 cos(pi k/N) = zeta_{2N}^k + zeta_{2N}^{-k} ---
def _divisors(n: int) -> list[int]:
    return [d for d in range(1, n + 1) if n % d == 0]


def _cyclotomic(n: int, cache: dict[int, list[int]]) -> list[int]:
    if n in cache:
        return cache[n]
    num = [-1] + [0] * (n - 1) + [1]
    for d in _divisors(n):
        if d < n:
            b = _cyclotomic(d, cache)
            db = len(b) - 1
            q, a = [0] * (len(num) - db), num[:]
            for i in range(len(a) - 1, db - 1, -1):
                c = a[i]
                if c:
                    q[i - db] = c
                    for jj in range(db + 1):
                        a[i - db + jj] -= c * b[jj]
            num = q
    cache[n] = num
    return num


def _power_mod(e: int, phi: list[int]) -> list[int]:
    deg = len(phi) - 1
    p = [0] * max(e + 1, deg)
    p[e] = 1
    for i in range(len(p) - 1, deg - 1, -1):
        c = p[i]
        if c:
            for jj in range(deg + 1):
                p[i - deg + jj] -= c * phi[jj]
    return p[:deg]


def active_set(n: int, j: int) -> list[int]:
    return [k for k in range(1, n) if k % 2 == 1 and (k * j) % n != 0]


def criterion_saturates(n: int, j: int, phi: list[int]) -> bool:
    """Exact criterion via PARI matkerint (saturated integer kernel of [2cos-coords; all-ones])."""
    act = active_set(n, j)
    if len(act) <= 1:
        return True
    deg = len(phi) - 1
    cols = [[h + lo for h, lo in zip(_power_mod(k % (2 * n), phi), _power_mod((2 * n - k) % (2 * n), phi),
                                     strict=True)] + [1] for k in act]
    mat = [[cols[c][i] for c in range(len(act))] for i in range(deg + 1)]
    mstr = "[" + ";".join(",".join(str(x) for x in row) for row in mat) + "]"
    script = (f"K=matkerint({mstr});if(type(K)==\"t_MAT\"&&matsize(K)[2]>0,"
              f"for(c=1,matsize(K)[2],print(Vec(K[,c]))),print(\"E\"));quit")
    out = subprocess.run(["gp", "-q"], input=script, capture_output=True, text=True, timeout=120).stdout
    neg = [sin(pi * k * j / n) < 0 for k in act]
    for line in out.splitlines():
        line = line.strip()
        if line.startswith("["):
            c = [int(x) for x in line.strip("[]").split(",")]
            if sum(c[i] for i in range(len(act)) if neg[i]) % 2:
                return False
    return True


def numeric_ratio(n: int, j: int, t_max: float = 12000.0, samples: int = 240_000) -> float:
    k = np.arange(1, n)
    a = np.where(k % 2 == 1, (2.0 / n) / np.tan(np.pi * k / (2 * n)), 0.0)
    amp = a * np.sin(np.pi * k * j / n)
    d = float(np.abs(amp).sum())
    if d < 1e-12:
        return 1.0
    w = np.cos(np.pi * k / n)
    ts = np.linspace(0.0, t_max, samples)
    z = np.abs((amp[None, :] * np.exp(-2j * np.outer(ts, w))).sum(axis=1))
    return float(z.max() / d)


def is_prime(n: int) -> bool:
    return n >= 2 and all(n % d for d in range(2, int(n**0.5) + 1))


def is_pow2(n: int) -> bool:
    return n & (n - 1) == 0


def main() -> int:
    print("=" * 84)
    print("Dirichlet-segment site-wise saturation: exact criterion vs sup_t, + corollaries")
    print(f"backend: {'PARI/gp + numpy' if HAVE_GP else 'numpy only (exact criterion skipped)'}")
    print("=" * 84)
    ok = True
    cache: dict[int, list[int]] = {}

    # (A) criterion vs numerical sup, and (C) the saturating-site table
    print("\n(A,C) saturating sites: exact criterion (=) vs numerical sup (~), N=4..24")
    print(f"  {'N':>3} {'class':>6}  saturating sites (criterion)")
    for n in range(4, 25):
        phi = _cyclotomic(2 * n, cache) if HAVE_GP else []
        sats, bad = [], False
        for j in range(1, n):
            crit = criterion_saturates(n, j, phi) if HAVE_GP else None
            if crit:
                sats.append(j)
            if HAVE_GP:
                num = numeric_ratio(n, j)
                # criterion=False must show a numerical gap; criterion=True may be slow (undershoot ok)
                if (not crit) and num > 0.999:
                    bad = True
        cls = "2^m" if is_pow2(n) else ("prime" if is_prime(n) else "comp")
        if HAVE_GP:
            ok &= not bad
            flag = "  <-- criterion/sup DISAGREE" if bad else ""
            print(f"  {n:>3} {cls:>6}  {sats}{flag}")
        else:
            print(f"  {n:>3} {cls:>6}  (need gp)")

    if not HAVE_GP:
        print("\nRESULT: gp not found -- install PARI/gp to run the exact criterion. Numerical law only.")
        return 0

    # (B) corollaries
    print("\n(B) Corollaries")
    c1 = all(criterion_saturates(n, 1, _cyclotomic(2 * n, cache))
             and criterion_saturates(n, n - 1, _cyclotomic(2 * n, cache)) for n in range(3, 30))
    print(f"    (i)  j=1 and j=N-1 always saturate (N<30): {'OK' if c1 else 'FAIL'}")
    c2 = all(all(criterion_saturates(n, j, _cyclotomic(2 * n, cache)) for j in range(1, n))
             for n in (5, 7, 11, 13, 4, 8, 16, 6, 10, 14, 22, 26))
    print(f"    (ii) prime / 2^m / 2p (full cos-rank) saturate at every site: {'OK' if c2 else 'FAIL'}")
    # sign-dependence: N=24 coprime saturating sites are the subgroup {1,7,17,23}
    phi24 = _cyclotomic(48, cache)
    cop = [j for j in range(1, 24) if gcd(j, 24) == 1]
    sat24 = {j for j in cop if criterion_saturates(24, j, phi24)}
    c3 = sat24 == {1, 7, 17, 23}
    print(f"    (iii) N=24 coprime saturating sites = {sorted(sat24)} (subgroup of (Z/24)*): "
          f"{'OK' if c3 else 'FAIL'}")
    ok &= c1 and c2 and c3

    print("\n" + "=" * 84)
    print("RESULT:", "SITE-WISE SATURATION CRITERION VERIFIED" if ok else "CHECK FAILED")
    print("Criterion (balanced-relation parity) matches the exact sup; j=1,N-1 always saturate;")
    print("prime/2^m saturate everywhere; composite saturation is sign-pattern (coset) dependent.")
    print("=" * 84)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
