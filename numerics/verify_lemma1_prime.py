# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2"]
# ///
"""Exact certificate of Lemma 1 for prime N=p: {sin(pi r/p): r=1..(p-1)/2} indep over Q.

Proof encoded (see docs/lower_bound_sketch_claude.md, Thm 1'):
  eta = zeta_{2p} = e^{i pi/p} satisfies eta^p = -1, and Q(eta) has Q-basis
  {eta^0, ..., eta^{p-2}} (degree phi(2p) = p-1). Since eta^{-r} = eta^{2p-r} = -eta^{p-r},

      2i * sin(pi r/p) = eta^r - eta^{-r} = eta^r + eta^{p-r}.

  For r = 1..(p-1)/2 the exponents {r} and {p-r} together cover {1,...,p-1} bijectively,
  and {eta^1,...,eta^{p-1}} is a Q-basis. So the coordinate vectors of the 2i*sin(pi r/p)
  in the basis {eta^0,...,eta^{p-2}} are Q-linearly independent.

This script builds those coordinate vectors as EXACT integer rows (via eta^p=-1 and the
cyclotomic relation Phi_{2p}(eta)=0) and computes their EXACT rank over Q. Rank == (p-1)/2
certifies the lemma rigorously (no floating point in the certificate). A float check then
confirms the integer coordinates really reconstruct eta^r - eta^{-r}.
"""

from __future__ import annotations

from fractions import Fraction

import numpy as np


def primes_up_to(limit: int) -> list[int]:
    out: list[int] = []
    for n in range(3, limit + 1):
        if all(n % d for d in range(2, int(n**0.5) + 1)):
            out.append(n)
    return out


def reduce_power(e: int, p: int) -> np.ndarray:
    """Coordinates of eta^e (eta=zeta_{2p}) in the Q-basis {eta^0,...,eta^{p-2}}, integer vector.

    Uses eta^p = -1 to fold e into [0, p), then for the residue p-1 the cyclotomic relation
    Phi_{2p}(eta) = sum_{j=0}^{p-1} (-1)^j eta^j = 0, i.e. eta^{p-1} = -sum_{j=0}^{p-2}(-1)^j eta^j.
    """
    q, s = divmod(e, p)
    sign = -1 if q % 2 else 1
    coords = np.zeros(p - 1, dtype=object)  # object -> exact Python ints
    if s <= p - 2:
        coords[s] = sign
    else:  # s == p-1
        for j in range(p - 1):
            coords[j] = sign * ((-1) ** (j + 1))
    return coords


def sin_rows(p: int) -> np.ndarray:
    """Integer matrix ((p-1)/2 x (p-1)): row r = coords of 2i*sin(pi r/p) = eta^r - eta^{-r}."""
    half = (p - 1) // 2
    rows = np.zeros((half, p - 1), dtype=object)
    for idx, r in enumerate(range(1, half + 1)):
        rows[idx] = reduce_power(r, p) - reduce_power((2 * p - r) % (2 * p), p)
    return rows


def exact_rank(matrix: np.ndarray) -> int:
    """Exact rank over Q of an integer matrix via fraction-free... here Fraction elimination."""
    m = [[Fraction(int(x)) for x in row] for row in matrix]
    rows, cols = len(m), len(m[0])
    rank, pivot_row = 0, 0
    for col in range(cols):
        sel = next((r for r in range(pivot_row, rows) if m[r][col] != 0), None)
        if sel is None:
            continue
        m[pivot_row], m[sel] = m[sel], m[pivot_row]
        pivot = m[pivot_row][col]
        for r in range(rows):
            if r != pivot_row and m[r][col] != 0:
                factor = m[r][col] / pivot
                m[r] = [a - factor * b for a, b in zip(m[r], m[pivot_row], strict=True)]
        pivot_row += 1
        rank += 1
        if pivot_row == rows:
            break
    return rank


def reconstruct_ok(p: int, rows: np.ndarray) -> float:
    """Max error between integer-coordinate reconstruction and the true 2i*sin(pi r/p)."""
    eta = np.exp(1j * np.pi / p)
    powers = eta ** np.arange(p - 1)               # eta^0 .. eta^{p-2}
    worst = 0.0
    for idx, r in enumerate(range(1, (p - 1) // 2 + 1)):
        recon = complex(np.dot(rows[idx].astype(np.complex128), powers))
        true = eta**r - eta ** (-r)                # = 2i sin(pi r/p)
        worst = max(worst, abs(recon - true))
    return worst


def main() -> int:
    print("=" * 70)
    print("Lemma 1 (prime N=p): {sin(pi r/p)} linearly independent over Q -- exact certificate")
    print("=" * 70)
    print(f"  {'p':>4} {'(p-1)/2':>8} {'exact rank':>11} {'recon err':>11} {'status':>7}")
    ok = True
    for p in primes_up_to(47):
        rows = sin_rows(p)
        rank = exact_rank(rows)
        err = reconstruct_ok(p, rows)
        expected = (p - 1) // 2
        good = (rank == expected) and (err < 1e-9)
        ok &= good
        print(f"  {p:>4} {expected:>8} {rank:>11} {err:>11.1e} {'PASS' if good else 'FAIL':>7}")
    print("=" * 70)
    print("RESULT:", "ALL PRIMES CERTIFIED" if ok else "FAILURE")
    print("Consequence: for odd prime p the ring frequencies {2 sin(pi r/p)} are Q-independent,")
    print("so Bohr/Kronecker give UNCONDITIONALLY A_p = U_p ~ (1/pi) ln p (large wave on the ring).")
    print("=" * 70)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
