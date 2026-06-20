# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2"]
# ///
"""Reproduce Filimonov's discrete-Schrodinger large wave on the DIRICHLET segment (J. Math. Sci. 269, 2023).

Model:  i z'_j = z_{j+1} - 2 z_j + z_{j-1},   z_0 = z_N = 0,   z_j(0) = alpha_j = 1.
Explicit solution (Dirichlet sine basis):
    z_j(t) = e^{2it} sum_{k=1}^{N-1} a_k sin(pi k j/N) e^{-2 i omega_k t},
    omega_k = cos(pi k/N),    a_k = (2/N) sum_{r=1}^{N-1} alpha_r sin(pi k r/N).
For alpha_j = 1:  a_k = (2/N) cot(pi k/2N) for k odd, 0 for k even.
Site-amplitude ceiling:  d_j^N = sum_k |a_k sin(pi k j/N)|,  with |z_j(t)| <= d_j^N always.

Filimonov, Theorem 1:  alpha_j=1, N=2^m  =>  sup_t |z_j(t)| = d_j^N.
Filimonov, Theorem 2:  lim_{N->inf} d_j^N = C ln j + O(1).

This script reproduces all of it AND adds two findings:
  * the constant is identified explicitly:  C = 4/pi^2 ~ 0.4053  (R^2 = 1, plus the limiting integral
    d_j^inf = (2/pi) int_0^{pi/2} cot(theta) |sin(2 j theta)| d theta ~ (2/pi)(2/pi) ln j);
  * cos-rank has the CLOSED FORM rank(N) = phi(2N)/2 - [N=2^m] (PROVED, paper Lemma "cosine rank, closed
    form"; corroborated to N=120). Hence full rank (= (N-1)//2) holds EXACTLY for N prime, 2^m, or 2*odd-prime
    -- BOTH directions proven, unlike the ring (sine) case (prime/2^m only). Proof: the symmetrized powers
    zeta_2N^k + zeta_2N^-k span Q(zeta_2N)^+, so rank = phi(2N)/2 - [1 not in span]; 1 IS in the span iff N
    has an odd prime factor p (since sum_{a=1}^{(p-1)/2} 2cos(2 pi a/p) = -1, terms = cosines k=2aN/p); the
    2^m exception via the Galois involution zeta_2N->zeta_2N^{N+1} + Chebyshev cos(kx)=cos x*V_k(cos^2 x). Full
    rank SUFFICES for sup_t|z_j| = d_j^N at every site, so Filimonov's N=2^m law extends to all primes AND 2p.
    Full rank is NOT necessary, though: composite N can still saturate at individual Dirichlet sites
    (e.g. j=1). The complete site-wise classification is now settled by an exact parity criterion on the
    active relation lattice -- see verify_segment_sitewise.py (Prop. site-wise segment saturation). This is
    separate from the settled ring saturation classification.
"""

from __future__ import annotations

import numpy as np

_P = 2_147_483_647


def a_coeffs(n: int) -> np.ndarray:
    k = np.arange(1, n)
    r = np.arange(1, n)
    return (2.0 / n) * np.sin(np.outer(k, r) * np.pi / n).sum(axis=1)  # alpha_r = 1


def d_jN(n: int) -> np.ndarray:
    k = np.arange(1, n)
    a = a_coeffs(n)
    j = np.arange(1, n)
    return np.abs(a[None, :] * np.sin(np.outer(j, k) * np.pi / n)).sum(axis=1)


# --- cyclotomic rank machinery for the cos(pi k/N) independence (Theorem 1) ---
def divisors(n: int) -> list[int]:
    return [d for d in range(1, n + 1) if n % d == 0]


def cyclotomic(n: int, cache: dict[int, list[int]]) -> list[int]:
    if n in cache:
        return cache[n]
    num = [-1] + [0] * (n - 1) + [1]
    for d in divisors(n):
        if d < n:
            b = cyclotomic(d, cache)
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


def power_mod(e: int, phi: list[int]) -> list[int]:
    deg = len(phi) - 1
    p = [0] * max(e + 1, deg)
    p[e] = 1
    for i in range(len(p) - 1, deg - 1, -1):
        c = p[i]
        if c:
            for jj in range(deg + 1):
                p[i - deg + jj] -= c * phi[jj]
    return p[:deg]


def cos_rank(n: int) -> int:
    """rank_Q{2 cos(pi k/N) = zeta_{2N}^k + zeta_{2N}^{-k} : k=1..N-1} via mod-p Gaussian elimination."""
    m = 2 * n
    cache: dict[int, list[int]] = {}
    phi = cyclotomic(m, cache)
    rows = [[h + lo for h, lo in zip(power_mod(k % m, phi), power_mod((m - k) % m, phi), strict=True)]
            for k in range(1, n)]
    a = np.array(rows, dtype=np.int64) % _P
    nr, nc = a.shape
    rk = 0
    for col in range(nc):
        nz = np.nonzero(a[rk:, col])[0]
        if nz.size == 0:
            continue
        piv = rk + int(nz[0])
        if piv != rk:
            a[[rk, piv]] = a[[piv, rk]]
        inv = pow(int(a[rk, col]), _P - 2, _P)
        a[rk] = (a[rk] * inv) % _P
        cv = a[:, col].copy()
        cv[rk] = 0
        nzr = np.nonzero(cv)[0]
        if nzr.size:
            a[nzr] = (a[nzr] - np.outer(cv[nzr], a[rk])) % _P
        rk += 1
        if rk == nr:
            break
    return rk


def is_prime(n: int) -> bool:
    return n >= 2 and all(n % d for d in range(2, int(n**0.5) + 1))


def main() -> int:
    print("=" * 78)
    print("Filimonov's discrete-Schrodinger large wave (Dirichlet segment, alpha_j=1)")
    print("=" * 78)
    ok = True

    print("\n(A) Closed form  a_k = (2/N) cot(pi k/2N) (k odd), 0 (k even)")
    n = 96
    a = a_coeffs(n)
    k = np.arange(1, n)
    cf = np.where(k % 2 == 1, (2.0 / n) / np.tan(np.pi * k / (2 * n)), 0.0)
    err = float(np.max(np.abs(a - cf)))
    ok &= err < 1e-10
    print(f"    N={n}: max|a_k - closed form| = {err:.2e}; a_k=0 for k even: {np.allclose(a[1::2], 0)}")

    print("\n(B) Theorem 2:  d_j^N = C ln j + O(1) at large N=2^12; identify C")
    n = 4096
    d = d_jN(n)
    j = np.arange(2, 300)
    slope, intc = np.polyfit(np.log(j), d[j - 1], 1)
    r2 = 1 - np.sum((d[j - 1] - (slope * np.log(j) + intc)) ** 2) / np.sum((d[j - 1] - d[j - 1].mean()) ** 2)
    c_pred = 4.0 / np.pi**2
    good = abs(slope - c_pred) < 5e-3 and r2 > 0.999
    ok &= good
    print(f"    fit d_j^N = {slope:.4f} ln j + {intc:.3f},  R^2 = {r2:.4f}")
    print(f"    -> C = {slope:.4f} matches 4/pi^2 = {c_pred:.4f}: {'OK' if good else 'FAIL'}")
    # limiting integral C = (2/pi)*(2/pi): check the integral representation grows like (4/pi^2) ln j
    th = np.linspace(1e-5, np.pi / 2, 200000)
    dint = np.array([(2 / np.pi) * np.trapezoid(1.0 / np.tan(th) * np.abs(np.sin(2 * jj * th)), th)
                     for jj in (10, 100, 1000)])
    cint = (dint[2] - dint[0]) / (np.log(1000) - np.log(10))
    ok &= abs(cint - c_pred) < 0.02
    print(f"    limiting integral (2/pi)int cot|sin 2j th|: slope vs ln j = {cint:.4f} (~4/pi^2): "
          f"{'OK' if abs(cint - c_pred) < 0.02 else 'FAIL'}")

    print("\n(C) Full cos-rank (N-1)//2 holds exactly for prime / 2^m / 2*odd-prime: SUFFICIENT for all sites")
    print(f"    {'N':>4} {'class':>10} {'rank':>5} {'(N-1)//2':>9} {'indep?':>7}")

    def _is2p(n: int) -> bool:  # twice an odd prime
        return n % 2 == 0 and is_prime(n // 2) and (n // 2) % 2 == 1

    for nn in (8, 16, 32, 7, 11, 13, 6, 10, 14, 22, 26, 9, 12, 15, 21):
        rk = cos_rank(nn)
        thr = (nn - 1) // 2
        is2m = (nn & (nn - 1)) == 0
        cls = "2^m" if is2m else ("prime" if is_prime(nn) else ("2*prime" if _is2p(nn) else "composite"))
        indep = rk == thr
        # full cos-rank <=> prime / 2^m / 2*odd-prime (the segment has the extra 2p family vs the ring)
        ok &= indep == (is_prime(nn) or is2m or _is2p(nn))
        print(f"    {nn:>4} {cls:>10} {rk:>5} {thr:>9} {indep!s:>7}")
    print("    -> cos(pi k/N) full rank EXACTLY for prime / 2^m / 2*odd-prime -- NOT just prime/2^m: the")
    print("       segment admits the extra 2p family (ring-p values for even k + primitive cos for odd k,")
    print("       spanning Q(zeta_4p)^+); sufficient for sup_t|z_j|=d_j^N at every site (extends Filimonov).")

    print("\n(C2) Closed form rank(N) = phi(2N)/2 - [N=2^m], and the proof's two key identities")

    def _phi(n: int) -> int:
        r, mm, d = n, n, 2
        while d * d <= mm:
            if mm % d == 0:
                while mm % d == 0:
                    mm //= d
                r -= r // d
            d += 1
        if mm > 1:
            r -= r // mm
        return r

    cf_ok = all(cos_rank(N) == _phi(2 * N) // 2 - (1 if (N & (N - 1)) == 0 else 0) for N in range(2, 49))
    print(f"    rank(N) = phi(2N)/2 - [N=2^m] for N=2..48: {'OK' if cf_ok else 'FAIL'}")
    # odd prime p|N => sum_{a=1}^{(p-1)/2} 2cos(2 pi a/p) = -1, each term a cosine of the family (k=2aN/p)
    mem_ok = True
    for nn, p in [(15, 3), (15, 5), (35, 7), (21, 3), (20, 5), (12, 3)]:
        s = sum(2 * np.cos(2 * np.pi * a / p) for a in range(1, (p - 1) // 2 + 1))
        km = all(2 * a * nn % p == 0
                 and abs(2 * np.cos(2 * np.pi * a / p) - 2 * np.cos(np.pi * (2 * a * nn // p) / nn)) < 1e-10
                 for a in range(1, (p - 1) // 2 + 1))
        mem_ok &= abs(s + 1) < 1e-9 and km
    print(f"    odd factor: sum 2cos(2 pi a/p) = -1, terms = cosines (k=2aN/p) -> 1 in span: "
          f"{'OK' if mem_ok else 'FAIL'}")
    xx = np.linspace(0.1, 0.5, 300)
    cc = np.cos(xx) ** 2
    ch_ok = all(np.max(np.abs(np.polyval(np.polyfit(cc, np.cos(k * xx) / np.cos(xx), (k - 1) // 2), cc)
                              - np.cos(k * xx) / np.cos(xx))) < 1e-9 for k in (1, 3, 5, 7, 9))
    print(f"    2^m: cos(kx)/cos(x) is degree (k-1)/2 in cos^2 x (Chebyshev V_k): {'OK' if ch_ok else 'FAIL'}")

    def _sin_rank(n: int) -> int:  # 2 sin(pi r/N) = 2 cos(pi(N-2r)/(2N)) = zeta_{4N}^|N-2r| + conj
        m = 4 * n
        phip = cyclotomic(m, {})
        rows = [[h + lo for h, lo in zip(power_mod(abs(n - 2 * r) % m, phip),
                                         power_mod((m - abs(n - 2 * r)) % m, phip), strict=True)]
                for r in range(1, n)]
        a = np.array(rows, dtype=np.int64) % _P
        nr, rk = a.shape[0], 0
        for col in range(a.shape[1]):
            nz = np.nonzero(a[rk:, col])[0]
            if nz.size == 0:
                continue
            piv = rk + int(nz[0])
            a[[rk, piv]] = a[[piv, rk]]
            a[rk] = (a[rk] * pow(int(a[rk, col]), _P - 2, _P)) % _P
            cv = a[:, col].copy()
            cv[rk] = 0
            nzr = np.nonzero(cv)[0]
            if nzr.size:
                a[nzr] = (a[nzr] - np.outer(cv[nzr], a[rk])) % _P
            rk += 1
            if rk == nr:
                break
        return rk

    sr_ok = all(_sin_rank(N) == _phi(2 * N) // 2 for N in range(3, 33))
    print(f"    sine companion (ring): rank{{2sin(pi r/N)}} = phi(2N)/2, NO exception, N=3..32: "
          f"{'OK' if sr_ok else 'FAIL'}")
    ok &= cf_ok and mem_ok and ch_ok and sr_ok

    print("\n(D) Ceiling holds: |z_j(t)| <= d_j^N for all sampled t (triangle inequality)")
    n = 256
    d = d_jN(n)
    kk = np.arange(1, n)
    a = a_coeffs(n)
    w = np.cos(np.pi * kk / n)
    ts = np.linspace(0, 3000.0, 30000)
    hold = True
    for jq in (16, 64, 128):
        amp = a * np.sin(np.pi * kk * jq / n)
        z = np.abs((amp[None, :] * np.exp(-2j * np.outer(ts, w))).sum(axis=1))
        hold &= z.max() <= d[jq - 1] + 1e-9
    ok &= hold
    print(f"    N=256: max_t|z_j| <= d_j^N at j=16,64,128: {hold} (finite scan undershoots the sup)")

    print("\n(E) Full rank is NOT necessary: composite N still saturates sup_t|z_j|=d_j^N at some site")
    e_ok = True
    for nn in (6, 9, 10, 12):
        d = d_jN(nn)
        kk = np.arange(1, nn)
        a = a_coeffs(nn)
        w = np.cos(np.pi * kk / nn)
        ts = np.linspace(0.0, 4000.0, 120000)
        ratio = 0.0
        for jq in range(1, nn):
            amp = a * np.sin(np.pi * kk * jq / nn)
            z = np.abs((amp[None, :] * np.exp(-2j * np.outer(ts, w))).sum(axis=1))
            ratio = max(ratio, z.max() / max(d[jq - 1], 1e-15))
        e_ok &= ratio > 0.999
        print(f"    N={nn:>3} (composite): max_j sup_t|z_j|/d_j^N = {ratio:.4f}  "
              f"{'saturates' if ratio > 0.999 else 'GAP'}")
    ok &= e_ok
    print("    -> composite N reaches the ceiling at individual sites; full cos-rank is a SUFFICIENT,")
    print("       not necessary, condition. The exact site-wise classification: verify_segment_sitewise.py.")

    print("\n" + "=" * 78)
    print("RESULT:", "FILIMONOV SCHRODINGER LARGE WAVE REPRODUCED" if ok else "CHECK FAILED")
    print("Thm 1 (sup=d_j^N) + Thm 2 (d_j^N ~ C ln j); identified C = 4/pi^2. Full cos-rank (prime/2^m)")
    print("is SUFFICIENT for saturation at all sites, not necessary -- composite N can saturate too.")
    print("=" * 78)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
