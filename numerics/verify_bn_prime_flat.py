# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2"]
# ///
"""B_p = sqrt(p) for every odd prime p, hence limsup_N B_N/sqrt(N) = 1 (Theorem bn-saturation (2)).

B_N = sup_t ||e^{-it L_N}||_{inf->inf} = sup_t sum_k |K_N(k,t)| <= sqrt(N), the sup over ALL t >= 0.
For an odd prime p the bound is reached in the limit:

  (1) lambda_r = 2 - zeta^r - zeta^-r (zeta = e^{2 pi i/p}), r = 1..(p-1)/2, are Q-independent:
      substitute 1 = -sum_{j=1}^{p-1} zeta^j; zeta^1..zeta^{p-1} is a basis of Q(zeta).
  (2) Kronecker: the flow t -> (t lambda_r)_r is dense in the torus, so e^{-it lambda_r} comes
      arbitrarily close to the chirp zeta^{r^2} (consistent with lambda_r = lambda_{p-r}, lambda_0 = 0).
  (3) Gauss sum: at the chirp |K_p(k)| = |sum_r zeta^{r^2 + r k}|/p = 1/sqrt(p) for every k,
      so sum_k |K_p(k)| = sqrt(p).

This script checks (1) and (3) in exact integer arithmetic, shows (2) numerically by explicit times with
||K_p(.,t)||_1/sqrt(p) within 1e-5 of 1 (p = 5, 7), and records that for small composite N, even and
odd, the long-time sup also exceeds the sup over the ballistic window 0.2N <= t <= 0.6N: the sharp
constants of the Debye analysis are windowed constants (t <= SN/2), not all-time ones.
The exact part is repeated in PARI/GP: exact/pari/bn_prime_flat.gp.
"""

from __future__ import annotations

from fractions import Fraction

import numpy as np

PRIMES = [p for p in range(3, 54) if all(p % d for d in range(2, int(p**0.5) + 1))]


def lambda_coords(p: int) -> list[list[Fraction]]:
    """Rows: lambda_r in the basis zeta^1..zeta^{p-1} of Q(zeta_p) (2 = -2 sum_j zeta^j)."""
    rows = []
    for r in range(1, (p - 1) // 2 + 1):
        row = [Fraction(-2)] * (p - 1)
        row[r - 1] -= 1
        row[p - r - 1] -= 1
        rows.append(row)
    return rows


def rank(rows: list[list[Fraction]]) -> int:
    m = [row[:] for row in rows]
    rk, col = 0, 0
    ncol = len(m[0]) if m else 0
    while rk < len(m) and col < ncol:
        piv = next((i for i in range(rk, len(m)) if m[i][col] != 0), None)
        if piv is None:
            col += 1
            continue
        m[rk], m[piv] = m[piv], m[rk]
        for i in range(len(m)):
            if i != rk and m[i][col] != 0:
                f = m[i][col] / m[rk][col]
                m[i] = [a - f * b for a, b in zip(m[i], m[rk], strict=True)]
        rk += 1
        col += 1
    return rk


def gauss_flat_exact(p: int) -> bool:
    """|G_k|^2 = sum_{r,s} zeta^{e(r)-e(s)} = p, e(r) = r^2 + r k, exactly.

    Writing |G_k|^2 = sum_e n_e zeta^e, the representation is unique up to adding a constant to every
    n_e (sum_e zeta^e = 0), so |G_k|^2 = p iff n_e is the same for all e != 0 and n_0 - n_1 = p.
    """
    for k in range(p):
        e = [(r * r + r * k) % p for r in range(p)]
        n = [0] * p
        for a in e:
            for b in e:
                n[(a - b) % p] += 1
        if len(set(n[1:])) != 1 or n[0] - n[1] != p:
            return False
    return True


def l1_ratio(n: int, t: float) -> float:
    lam = 2.0 - 2.0 * np.cos(2.0 * np.pi * np.arange(n) / n)
    return float(np.abs(np.fft.ifft(np.exp(-1j * t * lam))).sum() / np.sqrt(n))


def kronecker_witness(p: int, nmax: int) -> tuple[float, float]:
    """Pin the r = 1 phase on the chirp, t = 2 pi (n - 1/p)/lambda_1, pick the best n, refine t."""
    m = (p - 1) // 2
    lam = 2.0 - 2.0 * np.cos(2.0 * np.pi * np.arange(1, m + 1) / p)
    n = np.arange(1, nmax + 1, dtype=np.float64)
    t = 2.0 * np.pi * (n - 1.0 / p) / lam[0]
    err = np.zeros_like(t)
    for r in range(2, m + 1):
        err += 1.0 - np.cos(t * lam[r - 1] + 2.0 * np.pi * r * r / p)
    bt = float(t[int(np.argmin(err))])
    best, h = l1_ratio(p, bt), 1e-3
    for _ in range(60):
        up, dn = l1_ratio(p, bt + h), l1_ratio(p, bt - h)
        if up > best:
            bt, best = bt + h, up
        elif dn > best:
            bt, best = bt - h, dn
        else:
            h /= 2
    return bt, best


def long_time_sup(n: int, tmax: float, dt: float) -> float:
    lam = 2.0 - 2.0 * np.cos(2.0 * np.pi * np.arange(n) / n)
    best = 0.0
    for t0 in np.arange(0.0, tmax, 50.0):
        ts = np.arange(t0, min(t0 + 50.0, tmax), dt)
        k = np.abs(np.fft.ifft(np.exp(-1j * np.outer(ts, lam)), axis=1)).sum(axis=1)
        best = max(best, float(k.max()))
    return best / np.sqrt(n)


def main() -> int:
    ok = True
    print("=" * 72)
    print("B_p = sqrt(p) on prime rings; limsup_N B_N/sqrt(N) = 1")
    print("=" * 72)

    r_ok = all(rank(lambda_coords(p)) == (p - 1) // 2 for p in PRIMES)
    print(
        f"[1] exact: rank_Q {{lambda_r}} = (p-1)/2 for {len(PRIMES)} odd primes p <= 53 : "
        f"{'PASS' if r_ok else 'FAIL'}"
    )
    ok &= r_ok

    g_ok = all(gauss_flat_exact(p) for p in PRIMES)
    print(
        f"[2] exact: |sum_r zeta^(r^2+rk)|^2 = p for every k, p <= 53         : "
        f"{'PASS' if g_ok else 'FAIL'}"
    )
    ok &= g_ok

    dev = 0.0
    for p in PRIMES:
        r = np.arange(p)
        chirp = np.exp(2j * np.pi * (r * r % p) / p)
        dev = max(dev, abs(float(np.abs(np.fft.ifft(chirp)).sum()) / np.sqrt(p) - 1.0))
    c_ok = dev < 1e-12
    print(
        f"    float: chirp ||K||_1/sqrt(p) = 1 to {dev:.1e} (FFT), p <= 53       : "
        f"{'PASS' if c_ok else 'FAIL'}"
    )
    ok &= c_ok

    print("[3] Kronecker witnesses (numerical):")
    v3 = l1_ratio(3, 2.0 * np.pi / 9.0)
    print(f"    p = 3: t = 2 pi/9           ratio = {v3:.15f}  (attained)")
    w_ok = abs(v3 - 1.0) < 1e-12
    for p, nmax in ((5, 20000), (7, 400000)):
        t, v = kronecker_witness(p, nmax)
        print(f"    p = {p}: t = {t:>16.6f}  ratio = {v:.12f}")
        w_ok &= 1.0 - 1e-5 < v <= 1.0 + 1e-12
    print(
        f"    B_5/sqrt5 and B_7/sqrt7 within 1e-5 of 1                         : "
        f"{'PASS' if w_ok else 'FAIL'}"
    )
    ok &= w_ok

    print("[4] long times beat the window 0.2N <= t <= 0.6N (t <= 1500, step 0.01):")
    e_ok = True
    for n in (6, 8, 10, 12, 14, 9, 15):
        win = max(l1_ratio(n, t) for t in np.linspace(0.2 * n, 0.6 * n, 2000))
        lt = long_time_sup(n, 1500.0, 0.01)
        print(
            f"    N = {n:>2}: windowed {win:.4f}   sup_(t<=1500) >= {lt:.4f}   gain {lt - win:+.4f}"
        )
        e_ok &= lt >= win - 1e-4
        if n >= 9:
            e_ok &= lt > win + 0.01
    print(
        f"    long-time sup exceeds the windowed sup by > 0.01 for N = 9..15      : "
        f"{'PASS' if e_ok else 'FAIL'}"
    )
    ok &= e_ok

    print("=" * 72)
    print("RESULT:", "PASS -- B_p = sqrt(p), limsup B_N/sqrt(N) = 1" if ok else "FAIL")
    print("=" * 72)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
