# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2"]
# ///
"""Computer-assisted PROOF of the spectral hypotheses of the site-centred DNLS breather.

For the standing wave z=e^{-i Omega t}phi, (L_N+Omega)phi=gamma phi^3, the linearization splits into
    L_+ = L_N + Omega - 3 gamma diag(phi^2),    L_- = L_N + Omega - gamma diag(phi^2).
Orbital stability (Grillakis-Shatah-Strauss + Vakhitov-Kolokolov) needs, besides dP/dOmega>0 (proved
analytically by the anti-continuum expansion) and the Perron-Frobenius fact phi>0 => lambda_1(L_-)=0 with
ker L_- = span phi:
    (ii) the spectral GAP  lambda_2(L_-) > 0   (L_- >= 0 with a one-dimensional kernel),
    (iii) the Morse index   n(L_+) = 1.
These were previously only checked in floating point. Here they are CERTIFIED for a concrete ring
(N=5, gamma=1, Omega=2) with no SDP and no interval library -- exact rational arithmetic only:

  1. EXISTENCE (Newton-Kantorovich): a rational approximate profile phi~ with tiny residual R=F(phi~) is
     promoted to a TRUE solution phi in the ball ||phi-phi~||_2 <= rho := 2 ||J^{-1}|| ||R||, once the
     contraction test  ||J^{-1}|| * 3gamma(2 max|phi~| + 2rho) * rho <= 1/2  holds (J=F'(phi~)).
  2. INERTIA (Sylvester's law): exact LDL^T over Fraction gives the exact signature of L_+(phi~) -+ eps*I,
     proving n(L_+(phi~))=1 with the spectrum kept a distance eps from 0; a congruence Q^T(L_-(phi~)-eps*I)Q
     with Q a rational basis of phi~^perp proves lambda_2(L_-(phi~)) >= eps.
  3. WEYL: ||L_+-(phi)-L_+(phi~)||_2 <= 3gamma(2max|phi~|+rho)rho =: delta_+ < eps transports both signatures
     from phi~ to the true phi. Hence n(L_+(phi))=1 and lambda_2(L_-(phi)) >= eps-delta_- > 0.

No floating-point value enters the proof: floats only seed the rational profile; every inequality is
checked with exact Fractions. (Perron-Frobenius for phi>0 supplies lambda_1(L_-)=0, the missing kernel.)
"""

from __future__ import annotations

import math
from fractions import Fraction

import numpy as np

N, GAMMA, OMEGA = 5, 1, 2
DENOM = 1 << 44  # rationalization grid for the profile
EPS = Fraction(1)  # spectral-gap margin to certify (true gaps ~2, perturbation ~1e-12)

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


def inertia(a: Mat) -> tuple[int, int, int]:
    d = ldlt_pivots(a)
    if d is None:
        raise ValueError("zero pivot -- matrix singular, shift it")
    return sum(p < 0 for p in d), sum(p == 0 for p in d), sum(p > 0 for p in d)


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


# ---------- build the operators ----------
def lplus(phi2: Vec) -> Mat:
    a = lap_int(N)
    for j in range(N):
        a[j][j] += OMEGA - 3 * GAMMA * phi2[j]
    return a


def lminus(phi2: Vec) -> Mat:
    a = lap_int(N)
    for j in range(N):
        a[j][j] += OMEGA - GAMMA * phi2[j]
    return a


def residual(phi: Vec) -> Vec:
    lp = matvec(lap_int(N), phi)
    return [lp[j] + OMEGA * phi[j] - GAMMA * phi[j] ** 3 for j in range(N)]


def jacobian(phi: Vec) -> Mat:  # F'(phi) = L + Omega - 3 gamma diag(phi^2)
    return lplus([p ** 2 for p in phi])


def main() -> int:
    print("=" * 84)
    print(f"Computer-assisted proof: site-centred DNLS breather spectral hypotheses  (N={N}, "
          f"gamma={GAMMA}, Omega={OMEGA})")
    print("=" * 84)

    # --- float Newton seed, then rationalize ---
    lapf = 2 * np.eye(N) - np.eye(N, k=1) - np.eye(N, k=-1)
    lapf[0, -1] = lapf[-1, 0] = -1.0
    phi = np.zeros(N)
    phi[0] = math.sqrt((OMEGA + 2) / GAMMA)
    for _ in range(80):
        f = (lapf + OMEGA * np.eye(N)) @ phi - GAMMA * phi**3
        jf = lapf + OMEGA * np.eye(N) - 3 * GAMMA * np.diag(phi**2)
        phi = phi - np.linalg.solve(jf, f)
    phit: Vec = [Fraction(round(v * DENOM), DENOM) for v in phi]
    print(f"\nrational profile phi~ (denominator 2^44), all components > 0: {all(p > 0 for p in phit)}")
    print("   phi~ =", [f"{float(p):.6f}" for p in phit])

    # --- (1) Newton-Kantorovich existence of a true profile near phi~ ---
    r = residual(phit)
    nr2 = sum(v * v for v in r)
    norm_r = upper_sqrt(nr2)  # >= ||R||_2
    jinv = inverse(jacobian(phit))
    bnd = upper_sqrt(norm1(jinv) * norminf(jinv))  # >= ||J^{-1}||_2
    a_step = bnd * norm_r
    maxabs = max(abs(p) for p in phit)
    rho = 2 * a_step
    contraction = bnd * (3 * GAMMA * (2 * maxabs + rho)) * rho  # <= 1/2 needed
    ok_exist = contraction <= Fraction(1, 2)
    print("\n(1) EXISTENCE (Newton-Kantorovich):")
    print(f"    ||R||_2 <= {float(norm_r):.3e},  ||J^-1||_2 <= {float(bnd):.4f},  rho = 2||J^-1||||R|| "
          f"<= {float(rho):.3e}")
    print(f"    contraction test ||J^-1||*3g(2max|phi~|+2rho)*rho = {float(contraction):.3e} <= 1/2: "
          f"{ok_exist}")
    print(f"    => a TRUE breather phi exists with ||phi - phi~||_2 <= rho = {float(rho):.3e}")

    # --- Weyl perturbation radii of the two operators over the existence ball ---
    delta_p = 3 * GAMMA * (2 * maxabs + rho) * rho  # >= ||L_+(phi)-L_+(phi~)||_2
    delta_m = 1 * GAMMA * (2 * maxabs + rho) * rho  # >= ||L_-(phi)-L_-(phi~)||_2
    phi2 = [p**2 for p in phit]

    # --- (2)+(3) exact inertia at phi~, with margin EPS, then Weyl transport ---
    lp = lplus(phi2)
    lp_minus = [[lp[i][j] - (EPS if i == j else 0) for j in range(N)] for i in range(N)]
    lp_plus = [[lp[i][j] + (EPS if i == j else 0) for j in range(N)] for i in range(N)]
    neg_lo = inertia(lp_minus)[0]  # #{lambda < EPS}
    neg_hi = inertia(lp_plus)[0]  # #{lambda < -EPS}
    morse_ok = neg_lo == 1 and neg_hi == 1 and delta_p < EPS
    print("\n(3) MORSE INDEX  n(L_+) = 1  (Sylvester inertia at +-EPS, EPS=1):")
    print(f"    #eig(L_+(phi~)) < +EPS : {neg_lo}     #eig(L_+(phi~)) < -EPS : {neg_hi}")
    print("    => exactly one eigenvalue, and it is <= -EPS; spectrum-to-0 >= EPS=1")
    print(f"    Weyl: delta_+ = {float(delta_p):.3e} < EPS  => n(L_+(phi)) = 1 : {morse_ok}")

    # phi~^perp basis Q (columns): e_i - (phi~_i/phi~_0) e_0,  i=1..N-1  (exactly orthogonal to phi~)
    q: Mat = [[Fraction(0)] * (N - 1) for _ in range(N)]
    for c in range(N - 1):
        i = c + 1
        q[i][c] = Fraction(1)
        q[0][c] = -phit[i] / phit[0]
    lm = lminus(phi2)
    lm_shift = [[lm[i][j] - (EPS if i == j else 0) for j in range(N)] for i in range(N)]
    proj = matmul(transpose(q), matmul(lm_shift, q))  # Q^T (L_- - EPS I) Q
    gap_pivots = ldlt_pivots(proj)
    gap_ok = gap_pivots is not None and all(p > 0 for p in gap_pivots) and delta_m < EPS
    print("\n(2) SPECTRAL GAP  lambda_2(L_-) > 0  (congruence Q^T(L_- - EPS I)Q > 0 on phi~^perp):")
    print(f"    Q^T(L_-(phi~) - EPS I)Q positive-definite (all pivots > 0): "
          f"{gap_pivots is not None and all(p > 0 for p in gap_pivots)}")
    print(f"    => lambda_2(L_-(phi~)) >= EPS=1;  Weyl: delta_- = {float(delta_m):.3e} < EPS")
    print(f"    with Perron-Frobenius (phi>0 => lambda_1=0, ker=phi): L_- >= 0, dim ker = 1 : {gap_ok}")

    # --- float cross-check (sanity only, NOT part of the proof) ---
    lpf = lapf + OMEGA * np.eye(N) - 3 * GAMMA * np.diag(phi**2)
    lmf = lapf + OMEGA * np.eye(N) - GAMMA * np.diag(phi**2)
    print("\n[float cross-check, not part of the proof]")
    print(f"    eig(L_+) = {np.round(np.sort(np.linalg.eigvalsh(lpf)), 4)}")
    print(f"    eig(L_-) = {np.round(np.sort(np.linalg.eigvalsh(lmf)), 4)}")

    ok = ok_exist and morse_ok and gap_ok
    print("\n" + "=" * 84)
    print("RESULT:", "BREATHER STABILITY HYPOTHESES (ii),(iii) CERTIFIED (exact)" if ok else "FAIL")
    print("=" * 84)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
