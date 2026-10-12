# Exact-arithmetic layer

Independent, floating-point-free reimplementations of the linear-theory arithmetic, in two
computer-algebra systems, plus a symbolic cross-check of the analytic identities. These corroborate
the Python core (`src/large_wave_effect`) and the paper's Section on exact arithmetic.

| Tool | Script | Checks |
|------|--------|--------|
| **GAP** | `gap/relation_lattice.g` | library: exact `omega_r` coordinates in `Q(zeta_2N)`, `Q`-rank, relation lattice `Lambda_N`, the mod-4 criterion |
| | `gap/scan_classification.g` | finite exact cross-check of `rank = 1/2 phi(2N)` and `A_N = U_N <=> N` prime or `2^m`, for `N <= 800` |
| | `gap/prime_square.g` | the minimal mod-4 obstruction for `N = p^2` has support `p` on `{1} U {kp +- 1}` |
| | `gap/export_rocq.g` | emits the Rocq certificate `formal/rocq/Mod4Criterion.v` |
| **PARI/GP** | `pari/crosscheck.gp` | a second exact `Q`-rank implementation + the cyclotomic-degree identities behind the prime / `2^m` / `2p` proofs |
| | `pari/bn_prime_square.gp` | `B_N = sqrt(N)` for `N = 4, 9, 25, 49`: relation lattices of the eigenvalues (ranks `1, 1, 2, 3`), one point with rational phases that satisfies every relation, and `|K(k)|^2 = 1/N` exactly for every `k`; controls: none for `N = 6`, and the relation test rejects the 12 off-closure flat points of `N = 9` |
| | `pari/defect_sos.gp` | exact SOS certificates `U_N - A_N >= 0.06717197, 0.14754945, 0.11208183` for `N = 9, 15, 21`: subtorus relations in `Q(zeta_2N)`, Sturm enclosures of the weights, fraction-free `LDL*` of the Gram matrix, exact residual; seven mutated certificates must be rejected. Data: `pari/defect_sos_certs.gp`, written by `../numerics/build_defect_sos_certificates.py` |
| | `pari/kjump_relations.gp` | exact `Q`-relation lattices of the `k`-jump circulant frequencies `sqrt(alpha_r)`, `5 <= N <= 126`: square classes in `Q(zeta_N)` (Legendre rows, then `nfroots`), integer kernels inside each class, the mod-4 test. `C_N({1,3})` has relations exactly at `N = 28, 56, 84, 112`, all with coefficient sum `0`; none for `{1,2}`, `{1,2,3}`; control: the ring `{1}` is independent exactly at primes and `2^m` (`N <= 48`) |
| | `pari/kjump_nnn_primes.gp` | `C_p({1,2})` saturates for every prime `p >= 7`: the identities `alpha_r = -(3 + c_r)` mod squares and `prod (3 + c_r) = F_p` in `Z[x]/Phi_p`, `F_p` not a square, `p` not dividing `F_p`, for `7 <= p <= 211`; Legendre-character `F_2`-rank `k` certifies the conclusion; a mutated `alpha_r` must fail |
| | `pari/kjump_kummer.gp` | the same conclusion by the `2^k` subset-product search (`nfeltispower`), primes `7 <= p <= 31` (`2^15 - 1` products at `p = 31`; not in the driver) |
| | `pari/bn_prime_flat.gp` | `B_p = sqrt(p)`: exact `Q`-independence of the Schrödinger eigenvalues `2 - 2cos(2 pi r/p)` and exact Gauss-sum flatness, `p <= 53`, plus numerical Kronecker witnesses for `p = 5, 7` |
| **FriCAS** | `fricas/asymptotics.input` | `csc` Laurent series, the Appendix-B antiderivative, the constants `(2/pi)log(2/pi)` and `-pi/72` |

## Run

```bash
# Fedora: dnf install gap pari-gp fricas rocq-prover
gap  -q --nointeract exact/gap/scan_classification.g
gap  -q --nointeract exact/gap/prime_square.g
gp   -q exact/pari/crosscheck.gp
gp   -q -f exact/pari/bn_prime_flat.gp </dev/null
gp   -q -f exact/pari/bn_prime_square.gp </dev/null
gp   -q -f exact/pari/defect_sos.gp </dev/null
gp   -q -f exact/pari/kjump_relations.gp </dev/null
gp   -q -f exact/pari/kjump_nnn_primes.gp </dev/null
fricas -nosman < exact/fricas/asymptotics.input

# one-shot driver + cross-check against the Python core (skips missing tools):
uv run --script numerics/verify_exact_layer.py
```

## Scope (honest)

GAP and PARI verify the *arithmetic* (ranks, relation lattices, the mod-4 evaluation) in exact
cyclotomic arithmetic. The finite scans are not the proof of the general saturation classification; the
uniform proof is the rotated-root obstruction checked in `../numerics/verify_classification_proof.py`
and formalized in part in `../formal/rocq/Classification.v`. The `prime_square.g` finding is now best
read as a structural cross-check, not as the status of the theorem. See `../formal/rocq/README.md` for
what Rocq does and does not certify.
