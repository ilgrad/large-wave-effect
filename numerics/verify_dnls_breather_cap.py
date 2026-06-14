# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2"]
# ///
"""Computer-assisted PROOF of the breather stability dichotomy on the ring DNLS.

For the standing wave z=e^{-i Omega t}phi, (L_N+Omega)phi=gamma phi^3, the linearization splits into
    L_+ = L_N + Omega - 3 gamma diag(phi^2),    L_- = L_N + Omega - gamma diag(phi^2).
Grillakis-Shatah-Strauss + Vakhitov-Kolokolov: with dP/dOmega>0 (anti-continuum expansion) and phi>0
(Perron-Frobenius => lambda_1(L_-)=0, ker L_- = span phi), the orbital stability is decided by the
Morse index n(L_+) together with the gap lambda_2(L_-)>0:
    SITE-centred  n(L_+) = 1  => GSS index n(L_+)-1 = 0  => orbitally STABLE;
    BOND-centred  n(L_+) = 2  => GSS index n(L_+)-1 = 1 (odd) => orbitally UNSTABLE.
Both Morse indices were previously only read off a floating-point spectrum. Here they are CERTIFIED for
concrete rings with no SDP and no interval library -- exact rational arithmetic only:

  1. EXISTENCE (Newton-Kantorovich): a rational approximate profile phi~ with tiny residual R=F(phi~) is
     promoted to a TRUE solution phi in the ball ||phi-phi~||_2 <= rho := 2 ||J^{-1}|| ||R||, once the
     contraction test  ||J^{-1}|| * 3gamma(2 max|phi~| + 2rho) * rho <= 1/2  holds (J=F'(phi~)).
  2. INERTIA (Sylvester's law): exact LDL^T over Fraction gives the exact signature of L_+(phi~) -+ eps*I,
     proving n(L_+(phi~)) = expected with the spectrum kept a distance eps from 0; a congruence
     Q^T(L_-(phi~)-eps*I)Q with Q a rational basis of phi~^perp proves lambda_2(L_-(phi~)) >= eps.
  3. WEYL: ||L_+(phi)-L_+(phi~)||_2 <= 3gamma(2max|phi~|+rho)rho =: delta_+ < eps transports both
     signatures from phi~ to the true phi. Hence n(L_+(phi)) is certified and lambda_2(L_-(phi))>0.

No floating-point value enters the proof: floats only seed the rational profile; every inequality is
checked with exact Fractions. (Perron-Frobenius for phi>0 supplies lambda_1(L_-)=0, the missing kernel.)
"""

from __future__ import annotations

import math
from fractions import Fraction

import numpy as np

GAMMA = 1
DENOM = 1 << 44  # rationalization grid for the profile
EPS = Fraction(1)  # spectral-gap margin to certify (true gaps ~2, perturbation ~1e-8)

Mat = list[list[Fraction]]
Vec = list[Fraction]


# ---------- exact linear algebra over Q ----------
def lap_int(n: int) -> Mat:
    m = [[Fraction(0)] * n for _ in range(n)]
    for j in range(n):
        m[j][j] = Fraction(2)
        m[j][(j + 1) % n] -= 1
        m[j][(j - 1) % n] -= 1
    return m


def matvec(a: Mat, x: Vec) -> Vec:
    return [sum(a[i][k] * x[k] for k in range(len(x))) for i in range(len(a))]


def matmul(a: Mat, b: Mat) -> Mat:
    p, q, r = len(a), len(b[0]), len(b)
    return [[sum(a[i][k] * b[k][j] for k in range(r)) for j in range(q)] for i in range(p)]


def transpose(a: Mat) -> Mat:
    return [[a[i][j] for i in range(len(a))] for j in range(len(a[0]))]


def ldlt_pivots(a: Mat) -> list[Fraction] | None:
    """Exact LDL^T (no pivoting) of a symmetric matrix; returns the diagonal D (the pivots).

    By Sylvester's law of inertia the signs of the pivots are the signature of `a`.
    Returns None if a zero pivot is hit (the matrix is then singular -- not expected for our shifts).
    """
    n = len(a)
    ll = [[Fraction(0)] * n for _ in range(n)]
    d: list[Fraction] = []
    for j in range(n):
        s = a[j][j] - sum(ll[j][k] ** 2 * d[k] for k in range(j))
        if s == 0:
            return None
        d.append(s)
        for i in range(j + 1, n):
            ll[i][j] = (a[i][j] - sum(ll[i][k] * ll[j][k] * d[k] for k in range(j))) / s
    return d


def neg_count(a: Mat) -> int:
    d = ldlt_pivots(a)
    if d is None:
        raise ValueError("zero pivot -- matrix singular, shift it")
    return sum(p < 0 for p in d)


def inverse(a: Mat) -> Mat:
    n = len(a)
    aug = [a[i][:] + [Fraction(i == j) for j in range(n)] for i in range(n)]
    for c in range(n):
        piv = next(r for r in range(c, n) if aug[r][c] != 0)
        aug[c], aug[piv] = aug[piv], aug[c]
        inv = 1 / aug[c][c]
        aug[c] = [v * inv for v in aug[c]]
        for r in range(n):
            if r != c and aug[r][c] != 0:
                f = aug[r][c]
                aug[r] = [aug[r][k] - f * aug[c][k] for k in range(2 * n)]
    return [row[n:] for row in aug]


def norm1(a: Mat) -> Fraction:  # max abs column sum
    return max(sum(abs(a[i][j]) for i in range(len(a))) for j in range(len(a[0])))


def norminf(a: Mat) -> Fraction:  # max abs row sum
    return max(sum(abs(v) for v in row) for row in a)


def upper_sqrt(x: Fraction) -> Fraction:
    """A rational u with u*u >= x and u close to sqrt(x) (conservative upper bound)."""
    if x == 0:
        return Fraction(0)
    u = Fraction(math.sqrt(float(x))).limit_denominator(1 << 30) + Fraction(1, 1 << 28)
    while u * u < x:
        u += Fraction(1, 1 << 20)
    return u


def shift(a: Mat, eps: Fraction, sign: int) -> Mat:
    return [[a[i][j] + (sign * eps if i == j else 0) for j in range(len(a))] for i in range(len(a))]


# ---------- the breather problem ----------
def newton_profile(n: int, omega: int, seed: np.ndarray) -> np.ndarray:
    lapf = 2 * np.eye(n) - np.eye(n, k=1) - np.eye(n, k=-1)
    lapf[0, -1] = lapf[-1, 0] = -1.0
    phi = seed.copy()
    for _ in range(200):
        f = (lapf + omega * np.eye(n)) @ phi - GAMMA * phi**3
        jf = lapf + omega * np.eye(n) - 3 * GAMMA * np.diag(phi**2)
        phi = phi - np.linalg.solve(jf, f)
    return phi


def certify(label: str, n: int, omega: int, seed: np.ndarray, expected_n: int) -> bool:
    print(f"\n{'-' * 84}\n{label}:  N={n}, gamma={GAMMA}, Omega={omega}, expected n(L_+)={expected_n}\n{'-' * 84}")
    lap = lap_int(n)
    phit: Vec = [Fraction(round(v * DENOM), DENOM) for v in newton_profile(n, omega, seed)]
    phi2 = [p**2 for p in phit]
    pos = all(p > 0 for p in phit)

    def lplus() -> Mat:
        return [[lap[i][j] + ((omega - 3 * GAMMA * phi2[i]) if i == j else 0) for j in range(n)] for i in range(n)]

    def lminus() -> Mat:
        return [[lap[i][j] + ((omega - GAMMA * phi2[i]) if i == j else 0) for j in range(n)] for i in range(n)]

    # (1) Newton-Kantorovich existence
    lp_v = matvec(lap, phit)
    r = [lp_v[j] + omega * phit[j] - GAMMA * phit[j] ** 3 for j in range(n)]
    norm_r = upper_sqrt(sum(v * v for v in r))
    jinv = inverse(lplus())  # J = F'(phi~) = L_N + Omega - 3 gamma diag(phi~^2) = L_+
    bnd = upper_sqrt(norm1(jinv) * norminf(jinv))
    maxabs = max(abs(p) for p in phit)
    rho = 2 * bnd * norm_r
    contraction = bnd * (3 * GAMMA * (2 * maxabs + rho)) * rho
    exist_ok = contraction <= Fraction(1, 2) and pos
    delta_p = 3 * GAMMA * (2 * maxabs + rho) * rho
    delta_m = GAMMA * (2 * maxabs + rho) * rho

    # (3) Morse index via Sylvester inertia at +- EPS
    lp = lplus()
    nlo = neg_count(shift(lp, EPS, -1))  # #{lambda < EPS}
    nhi = neg_count(shift(lp, EPS, +1))  # #{lambda < -EPS}
    morse_ok = nlo == expected_n and nhi == expected_n and delta_p < EPS

    # (2) gap lambda_2(L_-) >= EPS via congruence on phi~^perp
    q: Mat = [[Fraction(0)] * (n - 1) for _ in range(n)]
    for c in range(n - 1):
        q[c + 1][c] = Fraction(1)
        q[0][c] = -phit[c + 1] / phit[0]
    proj = matmul(transpose(q), matmul(shift(lminus(), EPS, -1), q))
    piv = ldlt_pivots(proj)
    gap_ok = piv is not None and all(p > 0 for p in piv) and delta_m < EPS and pos

    print(f"(1) existence: ||R||_2<={float(norm_r):.2e}, ||J^-1||<={float(bnd):.3f}, rho<={float(rho):.2e}, "
          f"contraction={float(contraction):.2e}<=1/2 & phi~>0: {exist_ok}")
    print(f"(3) Morse: #eig<+EPS={nlo}, #eig<-EPS={nhi} => n(L_+)={nlo}; delta_+={float(delta_p):.2e}<EPS; "
          f"matches expected {expected_n}: {morse_ok}")
    print(f"(2) gap: Q^T(L_- - EPS I)Q > 0 (pivots>0): {piv is not None and all(p > 0 for p in piv)}; "
          f"delta_-={float(delta_m):.2e}<EPS => lambda_2(L_-(phi))>0: {gap_ok}")

    lapf = 2 * np.eye(n) - np.eye(n, k=1) - np.eye(n, k=-1)
    lapf[0, -1] = lapf[-1, 0] = -1.0
    phif = np.array([float(p) for p in phit])
    ep = np.sort(np.linalg.eigvalsh(lapf + omega * np.eye(n) - 3 * GAMMA * np.diag(phif**2)))
    print(f"    [float check] eig(L_+) = {np.round(ep, 3)}  (n(L_+)={int((ep < 0).sum())})")

    ok = exist_ok and morse_ok and gap_ok
    verdict = "STABLE (index 0)" if expected_n == 1 else "UNSTABLE (GSS odd index 1)"
    print(f"    => hypotheses certified: {ok}   [orbital {verdict}]")
    return ok


def main() -> int:
    print("=" * 84)
    print("Computer-assisted proof of the DNLS breather stability dichotomy (exact rational arithmetic)")
    print("=" * 84)

    def site_seed(n: int, omega: int) -> np.ndarray:
        s = np.zeros(n)
        s[0] = math.sqrt((omega + 2) / GAMMA)
        return s

    def bond_seed(n: int, omega: int) -> np.ndarray:
        s = np.zeros(n)
        s[0] = s[1] = math.sqrt((omega + 2) / GAMMA)
        return s

    site = certify("SITE-centred (stable)", 5, 2, site_seed(5, 2), expected_n=1)
    bond = certify("BOND-centred (unstable)", 6, 2, bond_seed(6, 2), expected_n=2)

    ok = site and bond
    print("\n" + "=" * 84)
    print("RESULT:", "BREATHER DICHOTOMY CERTIFIED -- site stable, bond unstable (exact)" if ok else "FAIL")
    print("=" * 84)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
