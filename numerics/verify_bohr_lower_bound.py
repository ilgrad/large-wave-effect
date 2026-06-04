# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2"]
# ///
"""Bohr closed-form sup -> the log N lower bound, the right tool (no t-scan).

For an almost-periodic f(t) = sum_k b_k sin(omega_k t) with {omega_k} linearly
independent over Q, Bohr's theorem gives the EXACT supremum without scanning time:

    sup_t f(t) = sum_k |b_k|   (Kronecker/Weyl makes all sin(omega_k t) -> sign b_k).

This is exactly the device Myshkis-Filimonov used. We (1) reproduce their segment
(Dirichlet, Zhukovsky forcing) result beta_j(N) ~ (4/pi^2) ln j, and (2) carry it to
the ring (free oscillations, impulse IC), where the degenerate pairs (r, N-r) recombine
and Bohr's sum gives the ring large wave A_N = U_N ~ (1/pi) ln N, attained at the
impulse node. This is the value the uniform t-scan only crept up to from far below.

Caveat (the genuinely open part): the equalities require linear independence over Q of
the frequency set; for N = p (prime) and N = 2^m a cyclotomic-field argument supplies it
(Filimonov, part 2 lecture 4), while composite N can introduce extra rational relations.
"""

from __future__ import annotations

import numpy as np


def ring_U(n: int) -> float:
    """Ring impulse large wave A_N = U_N = (1/2N) sum_{r=1}^{N-1} csc(pi r/N) (Bohr sup at j=0)."""
    r = np.arange(1, n)
    return float(np.sum(1.0 / np.sin(np.pi * r / n)) / (2 * n))


def ring_bohr_sup_per_node(n: int) -> np.ndarray:
    """Per-node Bohr sup_t |u_j(t)| = (1/N) sum_{r=1}^{N-1} |cos(2 pi r j/N)| / omega_r."""
    r = np.arange(1, n)
    omega = 2.0 * np.sin(np.pi * r / n)
    j = np.arange(n)
    return (np.abs(np.cos(2 * np.pi * np.outer(j, r) / n)) / omega).sum(axis=1) / n


def segment_beta(n: int, j: int) -> float:
    """Myshkis-Filimonov segment sup_t S_j(t) = (2/N) sum_k |sin(pi k j/N)| cot(pi k/2N), F=1."""
    k = np.arange(1, n)
    return float(2.0 / n * np.sum(np.abs(np.sin(np.pi * k * j / n)) / np.tan(np.pi * k / (2 * n))))


def loglog_slope(xs: np.ndarray, ys: np.ndarray) -> tuple[float, float]:
    """Fit ys ~ a*ln(xs)+b; return (a, R^2)."""
    lx = np.log(xs)
    a, b = np.polyfit(lx, ys, 1)
    pred = a * lx + b
    ss_res = float(np.sum((ys - pred) ** 2))
    ss_tot = float(np.sum((ys - ys.mean()) ** 2))
    return float(a), 1.0 - ss_res / ss_tot


def main() -> int:
    print("=" * 76)
    print("Bohr closed-form sup: the log N lower bound (no t-scan)")
    print("=" * 76)
    ok = True
    inv_pi = 1.0 / np.pi
    # sup_t S_j = (2F/N) sum|sin|cot -> (2/pi) int_0^pi |sin(jx)| cot(x/2) dx ~ (8/pi^2) ln j.
    # That is 2x Myshkis-Filimonov's quoted (4/pi^2), which is the half-amplitude / (F/N) form.
    seg_const = 8.0 / np.pi**2

    # (A) Ring: A_N = U_N ~ (1/pi) ln N, the EXACT sup (vs the scan's far-below lower bound)
    print("\n(A) Ring free impulse: A_N = U_N (Bohr) ~ (1/pi) ln N")
    print(f"  {'N':>6} {'A_N=U_N':>10} {'A_N/lnN':>9}")
    ns = np.array([17, 31, 61, 127, 251, 509, 1021, 2039, 4093])  # primes (expected lin-indep)
    a_ring = np.array([ring_U(int(n)) for n in ns])
    for n, a in zip(ns, a_ring, strict=True):
        print(f"  {n:>6} {a:>10.4f} {a / np.log(n):>9.4f}")
    slope, r2 = loglog_slope(ns.astype(float), a_ring)
    print(f"  fit A_N ~ a ln N: a={slope:.4f} (1/pi={inv_pi:.4f}), R^2={r2:.4f}")
    print(f"  -> slope matches 1/pi: {'PASS' if abs(slope - inv_pi) < 0.03 else 'FAIL'}")
    ok &= abs(slope - inv_pi) < 0.03

    # (B) Ring: the max over nodes of the Bohr sup is attained at the impulse node j=0
    print("\n(B) Ring: max_j sup_t|u_j(t)| attained at j=0 (so A_N = U_N)")
    argmax_ok = True
    for n in (31, 64, 127, 256):
        per = ring_bohr_sup_per_node(n)
        at0 = int(np.argmax(per)) == 0
        close = abs(per[0] - ring_U(n)) / ring_U(n) < 1e-12
        print(f"  N={n:>4}: argmax_j = {int(np.argmax(per)):>3}  per[0]={per[0]:.4f}  U_N={ring_U(n):.4f}")
        argmax_ok &= at0 and close
    print(f"  -> peak at j=0 and equals U_N: {'PASS' if argmax_ok else 'FAIL'}")
    ok &= argmax_ok

    # (C) Segment (Myshkis-Filimonov): sup_t S_j ~ (8/pi^2) ln j  -- reproduce their mechanism
    print("\n(C) Segment (Dirichlet, Zhukovsky): sup_t S_j(N) ~ (8/pi^2) ln j")
    n_big = 8192
    js = np.array([4, 8, 16, 32, 64, 128, 256, 512])
    betas = np.array([segment_beta(n_big, int(j)) for j in js])
    print(f"  N={n_big}, F=1")
    print(f"  {'j':>5} {'beta_j':>10} {'beta/lnj':>9}")
    for j, bta in zip(js, betas, strict=True):
        print(f"  {j:>5} {bta:>10.4f} {bta / np.log(j):>9.4f}")
    slope_s, r2_s = loglog_slope(js.astype(float), betas)
    print(f"  fit sup_t S_j ~ a ln j: a={slope_s:.4f} (8/pi^2={seg_const:.4f}), R^2={r2_s:.4f}")
    print("  (8/pi^2 = 2x MF half-amplitude 4/pi^2; fit approaches it slowly from below, O(1/ln j))")
    print(f"  -> slope matches 8/pi^2: {'PASS' if abs(slope_s - seg_const) < 0.04 else 'FAIL'}")
    ok &= abs(slope_s - seg_const) < 0.04

    print("\n" + "=" * 76)
    print("RESULT:", "ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED")
    print("Takeaway: with Bohr's sup = sum|coeff| (NOT a t-scan), the ring large wave is")
    print("A_N = U_N ~ (1/pi) ln N -- a genuine logarithmic effect on the ring, conditional")
    print("on linear independence of {2 sin(pi r/N)} over Q (proven for N=p, 2^m; open for")
    print("composite N). This reduces paper #1's lower bound to a clean number-theory lemma.")
    print("=" * 76)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
