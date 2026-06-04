# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2"]
# ///
"""Operator-norm hierarchy of the wave propagator: one effect, three growth laws.

The velocity propagator S_N(t) = L_N^{-1/2} sin(t sqrt(L_N)) is circulant with kernel G_N(.,t).
Its supremum over time in different operator norms grows at different rates -- this is the
variational/operator-theoretic face of the large wave. With omega_r = 2 sin(pi r/N), zeta(1)=(N^2-1)/12:

    ell^2 -> ell^2     sup_t = max_r 1/omega_r = 1/(2 sin(pi/N))         ~ N/(2pi)   (lowest-mode resonance)
    ell^2 -> ell^inf   sup_t = sqrt( zeta(1)/N ) = sqrt((N^2-1)/(12N))   ~ sqrt(N/12) (energy / variational M_N/sqrt2)
    ell^1 -> ell^inf   sup_t = max_k sup_t |G_N(k,t)| = U_N               ~ (1/pi) ln N (THE large wave A_N)

So the large wave is the *weakest-input / strongest-output* corner (ell^1 -> ell^inf) of a hierarchy;
the energy norm gives sqrt(N) (Theorem 1) and the spectral norm gives N. Each sup is a Bohr supremum
(closed form); a uniform t-scan only confirms the upper bound (it undershoots for the ell^inf outputs).
"""

from __future__ import annotations

import numpy as np


def green(n: int, t: float) -> np.ndarray:
    om = 2.0 * np.sin(np.pi * np.arange(n) / n)
    nz = om > 1e-12
    c = np.zeros(n)
    c[nz] = np.sin(om[nz] * t) / om[nz]
    return np.real(np.fft.ifft(c))


def exact_norms(n: int) -> tuple[float, float, float]:
    """(ell2->ell2, ell2->ellinf, ell1->ellinf) suprema over t, closed forms."""
    n22 = 1.0 / (2.0 * np.sin(np.pi / n))            # max_r 1/omega_r
    n2inf = np.sqrt((n * n - 1) / (12 * n))          # sqrt(zeta(1)/N)
    n1inf = float(np.sum(1.0 / np.sin(np.pi * np.arange(1, n) / n)) / (2 * n))  # U_N
    return n22, n2inf, n1inf


def loglog(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    a, b = np.polyfit(x, y, 1)
    pred = a * x + b
    r2 = 1 - float(np.sum((y - pred) ** 2)) / float(np.sum((y - y.mean()) ** 2))
    return float(a), r2


def main() -> int:
    print("=" * 78)
    print("Operator-norm hierarchy of S_N(t)=L_N^{-1/2} sin(t sqrt L_N)")
    print("=" * 78)
    ok = True

    # (A) scan confirms the closed-form suprema are upper bounds (reachable / approached)
    print("\n(A) sup over a t-grid <= closed form (sanity; ell^inf outputs undershoot via Bohr)")
    print(f"  {'N':>5} {'l2->l2 (=N/2pi)':>16} {'l2->linf(~sqrtN)':>17} {'l1->linf(=U_N)':>15}")
    for n in (16, 64, 256):
        e22, e2i, e1i = exact_norms(n)
        ts = np.linspace(0.01, 40.0 * n, 4000)
        s22 = s2i = s1i = 0.0
        for t in ts:
            g = green(n, t)
            om = 2.0 * np.sin(np.pi * np.arange(1, n) / n)
            s22 = max(s22, float(np.max(np.abs(np.sin(om * t) / om))))
            s2i = max(s2i, float(np.sqrt(np.sum(g * g))))
            s1i = max(s1i, float(np.max(np.abs(g))))
        good = s22 <= e22 + 1e-9 and s2i <= e2i + 1e-9 and s1i <= e1i + 1e-9
        ok &= good
        print(f"  {n:>5} {f'{s22:.3f}<={e22:.3f}':>16} {f'{s2i:.3f}<={e2i:.3f}':>17} "
              f"{f'{s1i:.3f}<={e1i:.3f}':>15}  {'OK' if good else 'FAIL'}")

    # (B) growth rates: N, sqrt(N), ln(N)
    print("\n(B) growth of the three suprema (closed forms)")
    ns = np.array([16, 32, 64, 128, 256, 512, 1024])
    e = np.array([exact_norms(int(n)) for n in ns])
    a22, r22 = loglog(ns.astype(float), e[:, 0])              # linear in N
    a2i, r2i = loglog(np.sqrt(ns), e[:, 1])                   # linear in sqrt N
    a1i, r1i = loglog(np.log(ns), e[:, 2])                    # linear in ln N
    print(f"  ell2->ell2 ~ a*N      : a={a22:.4f} (1/2pi={1 / (2 * np.pi):.4f})  R^2={r22:.4f}")
    print(f"  ell2->ellinf ~ a*sqrtN: a={a2i:.4f} (1/sqrt12={1 / np.sqrt(12):.4f})  R^2={r2i:.4f}")
    print(f"  ell1->ellinf ~ a*lnN  : a={a1i:.4f} (1/pi={1 / np.pi:.4f})  R^2={r1i:.4f}")
    ok &= r22 > 0.999 and r2i > 0.999 and r1i > 0.999
    ok &= abs(a22 - 1 / (2 * np.pi)) < 0.02 and abs(a1i - 1 / np.pi) < 0.03

    print("\n" + "=" * 78)
    print("RESULT:", "HIERARCHY VERIFIED" if ok else "CHECK FAILED")
    print("The large wave is the ell^1->ell^inf corner (ln N); the energy norm gives sqrt(N)")
    print("(variational Theorem 1) and the spectral norm gives N. One propagator, three growth laws.")
    print("=" * 78)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
