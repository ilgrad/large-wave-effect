# 2. The mod-4 ceiling criterion as the characterization of saturation

Status: accepted

## Context

The large wave saturates its ceiling, `A_N = U_N`, on some `N` and falls short on others. The first
characterization was *rational independence* of the frequencies `{2 sin(pi r/N)}` (true exactly for `N`
prime or `2^m`). But independence is only **sufficient**: a dependent `N` could in principle still
saturate. We needed the sharp, necessary-and-sufficient condition.

## Decision

State and prove the **ceiling criterion** (Theorem `thm:ceiling`): `A_N = U_N` iff every integer relation
`k` in the frequency relation lattice `Lambda_N` has coefficient-sum `sum_r k_r ≡ 0 (mod 4)`.

The proof is a Kronecker argument: `A_N <= U_N` always, with the maximum at node 0 where all modal
coefficients are positive, so `U_N` is attained iff the target phase `phi_* = (pi/2, ..., pi/2)` lies in
the orbit-closure subtorus, i.e. `<k, phi_*> = (pi/2) sum_r k_r ≡ 0 (mod 2 pi)` for every relation `k`,
which is exactly `sum_r k_r ≡ 0 (mod 4)`.

`Lambda_N` is computed **exactly** by an integer (Hermite) left-kernel of the cyclotomic-coordinate
matrix (`2 i sin(pi r/N) = zeta_{2N}^r - zeta_{2N}^{-r}`), in pure-Python big-integer arithmetic — no
floating point, no external CAS.

## Consequences

- The criterion is strictly weaker than full independence (independence ⇒ `Lambda_N = {0}` ⇒ vacuous),
  so it *could* admit composite saturators. That possibility is now ruled out by the rotated-root
  obstruction: every composite non-`2^m` has an exact relation with odd coefficient sum, violating both
  mod-4 branches. The finite scans are retained as exact-arithmetic cross-checks, not as the proof.
- The per-`N` integer computation remains useful for auditing certificates, but the global
  prime/`2^m` classification no longer rests on a finite range scan.
- A density-zero corollary follows: the saturating set is exactly `{prime} ∪ {2^m}`, both of density zero.

## Alternatives considered

- **Keep "independence" as the characterization**: rejected — it is only sufficient and hides the genuine
  question of whether a dependent `N` can still saturate.
- **Use GAP/Sage for the lattice / Smith normal form**: rejected — a pure-Python Hermite kernel is exact
  and keeps the one-command `uv run` reproducibility (see ADR 0001); a CAS dependency would not.
