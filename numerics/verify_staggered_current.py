# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2", "scipy"]
# ///
"""Staggered-current lemma (paper: Lemma lem:stagger-leading): the leading O(gamma) staggered-momentum
drift of the DNLS self-trapping window vanishes identically.

Along the LINEAR flow  i a' = -L a  (L_n = 2 a_n - a_{n+1} - a_{n-1}) from the single-site seed a = delta_0,
    a_n(t) = exp(2 i t) (-i)^n J_n(2t),   J_n real,
so a_{n+2} = (-J_{n+2}/J_n) a_n is a REAL multiple of a_n, whence Im(conj(a_n) a_{n+2}) = 0 and the
distance-2 current -- the RHS of the staggered-momentum identity (paper, eq:stagger) -- is zero pointwise.
(Closed form solves the ODE: verify_staggered_current.mac, Maxima; real-multiple => zero current:
formal/rocq/StaggeredCurrent.v, Rocq.)  This confirms it on the actual discrete ring flow (equalities hold
up to aliasing ~ machine eps for t <~ N/5), and contrasts with the O(1) distance-1 current.
"""
from __future__ import annotations

import numpy as np
from scipy.special import jv


def linear_flow(lam: np.ndarray, s: float) -> np.ndarray:
    """Exact linear-DNLS propagator e^{iLs} acting on the single-site seed delta_0 (b_hat_k(0) = 1)."""
    return np.fft.ifft(np.exp(1j * lam * s))


def staggered_F(b: np.ndarray) -> float:
    dens = np.abs(b) ** 2
    cur = np.imag(np.conj(b) * np.roll(b, -2))
    return float(-np.sum((dens - np.roll(dens, -2)) * cur))


def main() -> int:
    n_ring = 512
    lam = 2.0 - 2.0 * np.cos(2 * np.pi * np.arange(n_ring) / n_ring)
    idx = ((np.arange(n_ring) + n_ring // 2) % n_ring) - n_ring // 2  # signed index for (-i)^n J_n
    print("=" * 78)
    print("Staggered-current lemma:  Im(conj(a_n) a_{n+2}) = 0 along the linear DNLS flow")
    print("=" * 78)
    print(f"N={n_ring}; line approximation valid for t <~ N/5 = {n_ring // 5}")
    print(f"{'t':>5} {'|a-closed|':>12} {'max|Im(a_n a_{n+2})|':>20} {'max|Im(a_n a_{n+1})|':>20} {'|F[a_lin]|':>11}")
    ok = True
    for s in (1.0, 5.0, 20.0, 50.0, 100.0):
        a = linear_flow(lam, s)
        closed = np.exp(2j * s) * (-1j) ** idx * jv(idx, 2 * s)
        err = float(np.max(np.abs(a - closed)))
        cur2 = float(np.max(np.abs(np.imag(np.conj(a) * np.roll(a, -2)))))
        cur1 = float(np.max(np.abs(np.imag(np.conj(a) * np.roll(a, -1)))))
        f = abs(staggered_F(a))
        print(f"{s:>5.0f} {err:>12.2e} {cur2:>20.2e} {cur1:>20.2e} {f:>11.2e}")
        ok &= (err < 1e-9) and (cur2 < 1e-9) and (f < 1e-9) and (cur1 > 1e-3)
    print("=" * 78)
    print("  closed form matches the discrete flow; the distance-2 current and F vanish (~0);")
    print("  the distance-1 current is O(1) (the chain genuinely spreads).")
    print("RESULT:", "STAGGERED-CURRENT LEMMA CONFIRMED" if ok else "CHECK FAILED")
    print("=" * 78)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
