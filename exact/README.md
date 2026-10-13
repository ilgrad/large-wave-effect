# Exact-arithmetic layer

Independent, floating-point-free reimplementations of the linear-theory arithmetic, in two
computer-algebra systems, exact certificates for `N = 15` (Julia) and `N = 21` (python-flint), both re-checked in
PARI/GP, plus a symbolic
cross-check of the analytic identities. These corroborate
the Python core (`src/large_wave_effect`) and the paper's Section on exact arithmetic.

| Tool | Script | Checks |
|------|--------|--------|
| **GAP** | `gap/relation_lattice.g` | library: exact `omega_r` coordinates in `Q(zeta_2N)`, `Q`-rank, relation lattice `Lambda_N`, the mod-4 criterion |
| | `gap/scan_classification.g` | finite exact cross-check of `rank = 1/2 phi(2N)` and `A_N = U_N <=> N` prime or `2^m`, for `N <= 800` |
| | `gap/prime_square.g` | the minimal mod-4 obstruction for `N = p^2` has support `p` on `{1} U {kp +- 1}` |
| | `gap/export_rocq.g` | emits the Rocq certificate `formal/rocq/Mod4Criterion.v` |
| | `gap/segment_mod4.g` | Dirichlet segment: the mod-4 ceiling criterion on the exact relation lattice holds iff `N+1` is prime or `2^m`, `N <= 40` (Cor. `segceiling`) |
| | `gap/excess_omega2.g` | excess lemma, `omega(m) = 2`: rank and unimodular Smith form of `Lambda_N`, chained minimal relations (Remark `defect`) |
| | `gap/schrodinger_saturation.g` | Schrödinger orbit closure: `rank_Q{lambda_r} = 1/2 phi(N)` (`N <= 150`) and the prime-power relation lattice = coset relations (second implementation of `pari/bn_saturation.gp`) |
| **PARI/GP** | `pari/crosscheck.gp` | a second exact `Q`-rank implementation + the cyclotomic-degree identities behind the prime / `2^m` / `2p` proofs |
| | `pari/bn_prime_square.gp` | `B_N = sqrt(N)` for `N = 4, 9, 25, 49`: relation lattices of the eigenvalues (ranks `1, 1, 2, 3`), one point with rational phases that satisfies every relation, and `|K(k)|^2 = 1/N` exactly for every `k`; controls: none for `N = 6`, and the relation test rejects the 12 off-closure flat points of `N = 9` |
| | `pari/defect_sos.gp` | exact SOS certificates `U_N - A_N >= 0.06717197, 0.14754945, 0.11208183` for `N = 9, 15, 21`: subtorus relations in `Q(zeta_2N)`, Sturm enclosures of the weights, fraction-free `LDL*` of the Gram matrix, exact residual; seven mutated certificates must be rejected. Data: `pari/defect_sos_certs.gp`, written by `../numerics/build_defect_sos_certificates.py` |
| | `pari/defect_lower_bound.gp` | Proposition `deficit-lower`: for every non-saturating `N <= 240` and odd prime `p | N`, the class relations (one per orbit of `x mod 2N/p`) checked in `Q(zeta_2N)`, odd coefficient sums, disjoint supports, and the certified rational bound `sum 1/(2 sum k_r^2/b_r)` (rational enclosures of `pi` and `sin`) above the closed form `(s-1) sin(pi/2p)/(8ps) >= 1/(16p^2)`; the triplet bound for `3 | N` (`>= 0.0431`, limit `ln3/(4 pi)`); the `N = 2p` relation (`S = 4p csc(pi/2p)`); five mutations rejected |
| | `pari/kjump_relations.gp` | exact `Q`-relation lattices of the `k`-jump circulant frequencies `sqrt(alpha_r)`, `5 <= N <= 126`: square classes in `Q(zeta_N)` (Legendre rows, then `nfroots`), integer kernels inside each class, the mod-4 test. `C_N({1,3})` has relations exactly at `N = 28, 56, 84, 112`, all with coefficient sum `0`; none for `{1,2}`, `{1,2,3}`; control: the ring `{1}` is independent exactly at primes and `2^m` (`N <= 48`) |
| | `pari/kjump_nnn_primes.gp` | `C_p({1,2})` saturates for every prime `p >= 7`: the identities `alpha_r = -(3 + c_r)` mod squares and `prod (3 + c_r) = F_p` in `Z[x]/Phi_p`, `F_p` not a square, `p` not dividing `F_p`, for `7 <= p <= 211`; Legendre-character `F_2`-rank `k` certifies the conclusion; a mutated `alpha_r` must fail |
| | `pari/bn_prime_flat.gp` | `B_p = sqrt(p)`: exact `Q`-independence of the Schrödinger eigenvalues `2 - 2cos(2 pi r/p)` and exact Gauss-sum flatness, `p <= 53`, plus numerical Kronecker witnesses for `p = 5, 7` |
| | `pari/bn_saturation.gp` | Schrödinger saturation: `rank_Q{lambda_r} = 1/2 phi(N)` (`N <= 200`), the Gauss chirp flat in `H_p` (`p <= 101`), the prime-power relation lattice = coset relations, the even Frank point flat in `H_{p^2}` (`p <= 13`), no symmetric periodic Golay pair for `3 <= M <= 18` |
| | `pari/kjump_kummer.gp` | saturation `A_p = U_p` of the circulant `C_p({1,2})` for every prime `7 <= p <= 97`, and `Q`-independence of the NNN-ring frequencies at `N = 15`: no subset product of the radicands is a square in the cyclotomic field (quadratic-character certificate, cross-checked exhaustively for `p <= 17`; negative control `N = 12`) |
| | `pari/support3_relations.gp` | Lemma `cleancomb`: all support-3 relations `omega_a + omega_b = omega_c` for `3 | N <= 330` (the standard triplets plus `(N/10, N/6, 3N/10)` when `30 | N`), each decided modulo `Phi_2N` |
| | `pari/c15_flat_check.gp` | rebuilds the `N = 15` flatness polynomials `f_1..f_7` and re-checks the Julia certificate `D y^4 = sum N_j m_j f_j` by exact polynomial arithmetic (input: `C15_CERT=<file>`) |
| | `pari/c3q_flat_check.gp` | Theorem `bn-saturation` (4) for `N = 3q`: rebuilds the `Z_3 x Z_q` parametrization of `H_N` (exact DFT in `Q(zeta_N)`, dimension of `V`) and the autocorrelations from their definitions, then re-checks the certificate `1 = sum (orbit multiplier) * A(shift)` exactly in `Z[L]`, with a negative control (input: `C3Q_CERT=<file>` written by `../numerics/bn_flat_certificate_3q.py --out`; `N = 21` in about 7 s) |
| | `pari/b27_psatz_check.gp` | Theorem `bn-saturation` (4) for `N = 27`: checks the `mu`-coordinates `c(r)` exactly in `Q(zeta_27)` (and the `Q`-independence of `1, mu_1, ..., mu_8`), rebuilds the autocorrelations `R_1..R_13` on the 9-torus from them, expands the Gram matrix from its pair-orbit values and recomputes the residual `rho = m*Qm + sum g_l R_l + kappa` exactly; `lambda |B_3| + kappa > ||rho||_1` for the claimed `lambda`, two mutations rejected. The bound `lambda_min(Q) >= lambda` itself is checked by `../numerics/verify_b27_psatz.py` (exact Cholesky residual, FLINT). Data: `pari/b27_psatz_cert.gp.gz`, written by `../numerics/b27_psatz_certificate.py` |
| **Julia** | `julia/c15_flat_certificate.jl` | Theorem `bn-saturation` (4): no flat point on `H_15`, so `B_15 < sqrt 15`. Finds the Nullstellensatz certificate `D y^4 = sum N_j m_j f_j` (2159 terms, `D` of 253 digits) by `P`-adic lifting of the degree-18 Macaulay system and checks it over `Z`; Base Julia only, about 30 s |
| **FriCAS** | `fricas/asymptotics.input` | `csc` Laurent series, the Appendix-B antiderivative, the constants `(2/pi)log(2/pi)` and `-pi/72` |

## Run

```bash
# Fedora: dnf install gap pari-gp fricas rocq-prover; Julia from juliaup
gap  -q --nointeract exact/gap/scan_classification.g
gap  -q --nointeract exact/gap/prime_square.g
gap  -q -A --quitonbreak exact/gap/schrodinger_saturation.g
gp   -q exact/pari/crosscheck.gp
gp   -q -f exact/pari/bn_prime_flat.gp </dev/null
gp   -q -f exact/pari/bn_prime_square.gp </dev/null
gp   -q -f exact/pari/defect_sos.gp </dev/null
gp   -q -f exact/pari/defect_lower_bound.gp </dev/null
gp   -q -f exact/pari/kjump_relations.gp </dev/null
gp   -q -f exact/pari/kjump_nnn_primes.gp </dev/null
gp   -q -f exact/pari/bn_saturation.gp </dev/null   # prints RESULT: PASS; also check for ***
gp   -q -f exact/pari/support3_relations.gp </dev/null
gp   -q -f exact/pari/kjump_kummer.gp </dev/null
gap  -q -A --quitonbreak exact/gap/segment_mod4.g
gap  -q -A --quitonbreak exact/gap/excess_omega2.g
julia --startup-file=no exact/julia/c15_flat_certificate.jl /tmp/c15_certificate.gp
C15_CERT=/tmp/c15_certificate.gp gp -q -f exact/pari/c15_flat_check.gp </dev/null
uv run --script numerics/bn_flat_certificate_3q.py 7 --out /tmp/c21_certificate.txt   # about 3 min; 5 for N = 15
C3Q_CERT=/tmp/c21_certificate.txt gp -q -f exact/pari/c3q_flat_check.gp </dev/null
gp   -q -f exact/pari/b27_psatz_check.gp </dev/null   # identity part of the N = 27 certificate
uv run --script numerics/verify_b27_psatz.py           # PSD bound and identity, exact (FLINT)
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
