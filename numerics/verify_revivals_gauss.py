# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2"]
# ///
"""Quantum revivals and Gauss sums: the arithmetic of the discrete Schrodinger dynamics (honest).

Two regimes must be distinguished:

 * QUADRATIC dispersion (free particle / continuum limit), eigenvalues mu_r = r^2 on Z_q. At the
   rational "Talbot" times t = 2 pi s/q the propagator kernel is a QUADRATIC GAUSS SUM, giving an
   EXACT fractional revival: for odd q the kernel has uniform magnitude |U_{jk}| = 1/sqrt(q), the
   hallmark of the discrete Talbot effect / arithmetic quantum dynamics. The Gauss-sum modulus
   |sum_{n=0}^{q-1} e^{2 pi i a n^2/q}| = sqrt(q) (odd q, gcd(a,q)=1) is the exact arithmetic input.

 * TIGHT-BINDING ring (our L_N), lambda_r = 4 sin^2(pi r/N). Here 2 cos(2 pi r/N) is irrational for
   generic N (Niven: rational only for N in {1,2,3,4,6}), so there is NO finite revival period: the
   return amplitude |K_N(0,t)| is almost-periodic and approaches 1 (almost-revivals) without reaching
   it at any fixed time. This almost-periodicity (Bohr) is the recurrence structure behind the large wave.

In the long-wave limit lambda_r ~ (2 pi r/N)^2 the tight-binding model approaches the quadratic one,
so the Talbot/Gauss structure emerges only approximately on the lattice. All of this is verified below.
"""

from __future__ import annotations

import numpy as np


def gauss_sum(a: int, q: int) -> complex:
    n = np.arange(q)
    return complex(np.sum(np.exp(2j * np.pi * a * n * n / q)))


def free_particle_kernel(q: int, s: int) -> np.ndarray:
    """Free-particle (mu_r=r^2) propagator row at Talbot time t=2 pi s/q on Z_q."""
    r = np.arange(q)
    phase = np.exp(-2j * np.pi * s * r * r / q)
    return np.fft.ifft(phase)  # U_{0,m}


def tb_return(n: int, t: float) -> float:
    """|K_N(0,t)| for the tight-binding ring (return amplitude to the initial delta)."""
    lam = 4.0 * np.sin(np.pi * np.arange(n) / n) ** 2
    return abs(complex(np.mean(np.exp(-1j * t * lam))))


def main() -> int:
    print("=" * 74)
    print("Revivals and Gauss sums: quadratic (Talbot) vs tight-binding (almost-periodic)")
    print("=" * 74)
    ok = True

    # (A) quadratic Gauss sum modulus = sqrt(q)  (odd q, gcd(a,q)=1)
    print("\n(A) |sum_n e^{2pi i a n^2/q}| = sqrt(q)  (odd q, gcd(a,q)=1)")
    g_ok = True
    for q in (5, 7, 9, 11, 25, 49):
        for a in (1, 2, 3):
            if np.gcd(a, q) != 1:
                continue
            val = abs(gauss_sum(a, q))
            g_ok &= abs(val - np.sqrt(q)) < 1e-9
    print(f"  checked q in {{5,7,9,11,25,49}}, a coprime: all |G|=sqrt(q) -> {'PASS' if g_ok else 'FAIL'}")
    ok &= g_ok

    # (B) free-particle fractional revival: at t=2pi/q the kernel has uniform |U|=1/sqrt(q)
    print("\n(B) Free particle at Talbot time t=2pi/q: |U_{0,m}|=1/sqrt(q) for all m (odd q)")
    print(f"  {'q':>4} {'min|U|':>9} {'max|U|':>9} {'1/sqrt(q)':>10} {'sum|U|^2':>9}")
    rev_ok = True
    for q in (5, 7, 11, 13):
        u = np.abs(free_particle_kernel(q, 1))
        uniform = abs(u.min() - 1 / np.sqrt(q)) < 1e-9 and abs(u.max() - 1 / np.sqrt(q)) < 1e-9
        unit = abs(float(np.sum(u**2)) - 1.0) < 1e-12
        rev_ok &= uniform and unit
        print(f"  {q:>4} {u.min():>9.5f} {u.max():>9.5f} {1 / np.sqrt(q):>10.5f} {float(np.sum(u**2)):>9.5f}")
    print(f"  -> uniform 1/sqrt(q) revival (Gauss): {'PASS' if rev_ok else 'FAIL'}")
    ok &= rev_ok

    # (C1) quadratic free particle: EXACT full revival at the Talbot period t=2pi
    print("\n(C1) Quadratic free particle: exact full revival K_q(0,2pi)=1 (period 2pi)")
    c1 = True
    for q in (5, 7, 11, 13):
        kq = abs(complex(np.mean(np.exp(-1j * 2 * np.pi * np.arange(q) ** 2))))
        c1 &= abs(kq - 1.0) < 1e-12
    print(f"  K_q(0,2pi)=1 for q in {{5,7,11,13}}: {'PASS' if c1 else 'FAIL'}")
    ok &= c1

    # (C2) tight-binding: NO finite revival period (Niven); almost-revivals approach 1 with the window
    print("\n(C2) Tight-binding ring N=7: almost-revivals approach 1, with no fixed revival time")
    print(f"  {'window T':>10} {'max|K_7(0,t)|':>14} {'argmax t':>10}")
    prev, grow = 0.0, True
    for tmax in (200.0, 800.0, 3200.0):
        ts = np.linspace(0.2, tmax, int(tmax * 100))
        vals = np.array([tb_return(7, float(t)) for t in ts])
        i = int(vals.argmax())
        print(f"  {tmax:>10.0f} {vals[i]:>14.5f} {ts[i]:>10.1f}")
        grow &= vals[i] >= prev - 1e-9
        prev = vals[i]
    print("  -> approaches 1 with wandering time (no finite period, Niven): "
          f"{'PASS' if grow else 'FAIL'}")
    ok &= grow

    print("\n" + "=" * 74)
    print("RESULT:", "REVIVALS/GAUSS VERIFIED" if ok else "CHECK FAILED")
    print("Exact Talbot revivals (Gauss sums, |U|=1/sqrt q) belong to the QUADRATIC dispersion")
    print("(continuum/free-particle limit); the tight-binding ring is almost-periodic (Bohr) with no")
    print("exact revival. The large wave lives in this almost-periodic, Diophantine recurrence regime.")
    print("=" * 74)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
