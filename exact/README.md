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
| **FriCAS** | `fricas/asymptotics.input` | `csc` Laurent series, the Appendix-B antiderivative, the constants `(2/pi)log(2/pi)` and `-pi/72` |

## Run

```bash
# Fedora: dnf install gap pari-gp fricas rocq-prover
gap  -q --nointeract exact/gap/scan_classification.g
gap  -q --nointeract exact/gap/prime_square.g
gp   -q exact/pari/crosscheck.gp
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
