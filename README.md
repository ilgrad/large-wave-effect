# large-wave-effect

[![DOI](https://zenodo.org/badge/1258683042.svg)](https://doi.org/10.5281/zenodo.20795861)

Rigorous study of the **large-wave effect** on periodic discrete chains: a localized impulse on a
homogeneous one-dimensional *lattice* amplifies into a peak that grows with the number of nodes — a
genuinely discrete resonance with no continuum counterpart.

One paper (`paper/`), in two parts built on a shared spectral core (the ring Laplacian `L_N`,
Toeplitz/circulant methods):

- **Part I — Linear theory** (rigorous core). Amplitude `A_N = sup_t max_j |u_j(t)|` of the velocity
  Green's function on the ring; cyclotomic independence, the mod-4 ceiling criterion, the discrete
  Schrödinger contrast, spectral zeta, operator norms, statistics, and the dimensional threshold.
- **Part II — Nonlinear extension**. The discrete NLS (Peregrine soliton, modulational instability,
  rogue-wave statistics), discrete breathers and self-trapping, FPUT recurrence, and the route to
  chaos (Lyapunov exponent).

Every quantitative claim is reproduced by a standalone, version-controlled script.

## Main results

Proved results are kept strictly separate from numerical evidence and the (sharply reduced) open problems.

**Linear theory (Part I).**

- **Order law (theorem, all `N`).** `A_N ~ (1/π) ln N` for *every* `N` — the central result, formerly
  a conjecture. The lower bound comes from a *palindromic argument*: the first `φ(2N)/2`
  frequencies `2 sin(πr/N)` are `ℚ`-independent, so a `(1−o(1))` fraction of the ceiling can
  always be Kronecker-aligned. Precisely `A_N = (1/π) ln N + O(ln ln ln N)`.
- **Exact amplitude / saturation (theorem).**
  `A_N = U_N = (1/π) ln N + (1/π)(γ + ln(2/π)) − π/(72 N²) + ⋯` **iff** `N` is
  prime or a power of two; equivalently every integer frequency relation has coefficient-sum
  `≡ 0 (mod 4)` (antipodal variant for even `N`). Via an explicit root-of-unity obstruction,
  `{N : A_N = U_N} = {primes} ∪ {2^m}` (density zero). The orbit-closure subtorus has
  dimension exactly `φ(2N)/2`. The theorem is cross-checked by exact Python/GAP scans
  (`verify_classification_proof.py`, `exact/gap/scan_classification.g`).
- **Composite deficit.** The order theorem gives the general upper envelope
  `0 ≤ U_N − A_N = O(ln ln ln N)`. The strict gap `A_N < U_N` for every composite non-power-of-two follows from
  the saturation classification, and it is quantitative (Proposition `deficit-lower`): one relation per residue
  class gives `U_N − A_N ≥ 1/(16p²)` for `p` the least odd prime factor, the vertex-disjoint triplets give
  `U_N − A_N > 0.043` whenever `3 | N` (limit `ln 3/(4π)`), and for `N = 2p` the deficit is
  `(π³/16) N⁻² (1 + O(1/N))` (`verify_defect_lower_bound.py`, exact check `exact/pari/defect_lower_bound.gp`).
  Sharper small-`N` gaps are certified by a Lipschitz grid
  (`N = 6, 9, 12, 15`) and by exact rational sum-of-squares certificates (`exact/pari/defect_sos.gp`):
  `U_N − A_N ≥ 0.06717197, 0.14754945, 0.11208183` for `N = 9, 15, 21`, within `10⁻⁷` of the optimizer. The sharper excess law
  `U_N − A_N = (1/π) ln(m/φ(m)) + O(1)` for odd part `m` remains open beyond prime-power odd part. The
  per-family limits `defect_∞(m) = lim_a (U − A)(m · 2^a)` are computed to `10⁻⁸`
  (`0.10368…`, `0.03824…`) as numerical constants with no closed form.
- **Discrete Schrödinger.** A localized state is never amplified (unitarity), yet
  `B_N = sup_t ‖e^{−it L_N}‖_{ℓ^∞ → ℓ^∞} = Θ(√N)`, with `liminf B_N/√N ≥ c_0/√2 = 0.861` proved
  (`c_0 = (2/π)^{3/2} · B(1/2, 3/4)`). **Infinite-time saturation (theorem).** `B_N = √N` iff the
  orbit closure — a subtorus of dimension exactly `φ(N)/2` — contains a flat phase vector
  (ε-uniform mixing of the quantum walk on `C_N`). This holds for `N = p` (the Gauss chirp;
  Godsil–Mullin–Roy) and for `N = p²` (an even Frank sequence; Cao–Zhu, arXiv:2607.04207), so
  `limsup B_N/√N = 1`; it fails for every even `N ≥ 6` (a symmetric periodic Golay-pair obstruction: no
  ε-uniform mixing on even cycles; for `N` not an even square this also follows from Monterde's bipartite
  conditions, arXiv:2603.20977, so the even squares `16, 36, 64, …` are the new cases). It fails for `N = 15` as well (computer-assisted): an exact
  Nullstellensatz certificate `D y⁴ = Σ N_j m_j f_j` (2159 terms) puts `y⁴` in the flatness ideal, `y` a coordinate
  of the orbit-closure torus (`exact/julia/c15_flat_certificate.jl`, re-checked in `exact/pari/c15_flat_check.gp`);
  the Gröbner basis in Cao–Zhu's Thm 3.2 implies it too, as it lives on the whole closure. It fails for
  `N = 21` (computer-assisted, new): in coordinates adapted to `Z_21 = Z_3 × Z_7` the orbit closure is the
  6-torus `x(1,a) = y_a, x(2,a) = y_{-a}, x(0,a) = 1/(y_a y_{-a})`, `Π y_a = 1`, and an exact invariant
  certificate `1 = g_0 I_0 + Σ_u u·(g_1 I_1) + Σ_u u·(h_1 II_1)` over the autocorrelations (4005 orbit
  coefficients, 614-digit denominator) shows they have no common zero (`numerics/bn_flat_certificate_3q.py`,
  re-checked from the definitions in `exact/pari/c3q_flat_check.gp`; the same code reproves `N = 15`). So `C_21`
  admits no uniform mixing, exact or ε-approximate — the case Gray–Patel–Young (arXiv:2607.28345) name as the
  next open one. That `B_N = √N` iff `N` is a prime or a prime square is a conjecture: the other odd composites
  that are not prime squares (from `N = 27` on) are open; numerically `B_N/√N < 1` there, closest at `N = 27`
  (`0.99992`)
  (`verify_bn_saturation.py`, `exact/pari/bn_saturation.gp`; a second implementation of the prime and
  prime-square cases in `verify_bn_prime_flat.py`, `exact/pari/bn_prime_flat.gp`,
  `exact/pari/bn_prime_square.gp`). **Windowed constant.** On the
  dispersive window `t = O(N)` the constant splits by parity — numerically even `N → c_0/√2`, odd
  `N → β_odd = 0.928…`, an elliptic-integral constant for which no elementary closed
  form was found. For odd `N` and each fixed `s = 2t/N ≤ 1` the limiting `ℓ¹`
  profile is the *explicit* elliptic-integral function `F(s)` (one-copy term `A(s) = ₂F₁` in closed form),
  and its maximizer `max F = F(1) = β_odd` is **proved by validated interval arithmetic** on all of
  `[0, 1]` (`verify_bn_profile_rigorous.py`: a chord minorant `β_odd ≥ 0.92735` from AGM-enclosed elliptic
  integrals, analytic bounds on `[0, 0.575]`, an enclosed tangent majorant on `[0.575, 0.90]`, and
  `F′(s) > 0` on `[0.90, 1]`). For `t > N/2` the whole `N → ∞` reduction at each fixed `s` is
  **rigorous** as well — joint equidistribution of the multi-copy Debye phases *and* its uniform-in-`u`
  transfer to the `ℓ¹` constant `F(s)` (panels away from the band edges, `O(ε^{3/4})` edge strips; no
  rate) — leaving two steps: convergence uniform in `s` (on `t ≤ N/2` as on `t > N/2`), and a **single**
  validated residual, the Gaussian-excess inequality `E|Z| − E|G| ≤ K⋆ · Σ_2^{−3/2} Σ_ℓ a_ℓ⁴`
  (`K⋆ = (h(1) − √π/2)/2`) for the profile's **balanced** Debye configurations (the live set is an
  interval of consecutive integers, so the even/odd channel counts differ by `≤ 1`; sharp at the two-copy
  pair — imbalanced configs violate the unrestricted bound) whose integrated majorant clears
  `β_odd` — and whose residual is now **proved for every configuration with at most three copies** by
  rigorous Arb certificates (the singleton-channel band, `dagger_rigorous.py`: 111 boxes, worst
  `R ≤ 0.9205 K⋆`; the near-pair basin, `dagger_nearpair.py`: 1231 boxes, convexity around the pair; coverage
  checked in exact rational arithmetic by `verify_dagger_rigorous.py` and `verify_dagger_nearpair.py`), and
  validated beyond three copies (see open problem (1)). So the windowed constant
  `limsup_N B_N(SN/2)/√N = β_odd` (fixed `S`) holds modulo the rest of this inequality and the `s`-uniform
  transfer. On the **Dirichlet segment**
  `d_j^N ~ (4/π²) ln j` for `N+1` prime/`2^m` (the constant `C = 4/π²` is the published
  Myshkis–Filimonov 2003 value). Three large-wave laws on one Laplacian.
- **Dimension, products & circulants (theorem).** The large wave lives exactly at spectral dimension
  `d_s = 1` (`O(1)` for `d_s ≥ 2`). A graph-agnostic order criterion plus the embedded ring spectrum give
  `A = Θ(ln N)` on *every* Cartesian product `C_N □ H` (ladders, prisms, tubes) and, beyond
  products, on the **Möbius ladder** `C_{2N}({1, N})` via its embedded-ring branch. Multi-jump circulants
  `C_N(S)` (`|S| ≥ 2`) instead **saturate** (`A_N = U_N`): for `S = {1, 2}`, `{1, 3}`, `{1, 2, 3}` this is
  decided exactly for every `5 ≤ N ≤ 126` (`exact/pari/kjump_relations.gp`). The radical frequencies are
  not always `ℚ`-independent — `C_N({1, 3})` has two relations at each of `N = 28, 56, 84, 112` — but
  every relation has coefficient sum `0`, so the mod-4 test passes. For the nearest case `S = {1, 2}`
  saturation is **proved for every prime `p ≥ 7`**: up to squares `α_r = −(3 + c_r)`, and
  `∏ (3 + c_r) = F_p` is not a square (Cohn), so by Kummer theory no subset product of the `α_r` is a
  square in `ℚ(ζ_p)` (`exact/pari/kjump_nnn_primes.gp` checks the identities for `p ≤ 211`;
  `exact/pari/kjump_kummer.gp` certifies the conclusion independently for `7 ≤ p ≤ 97` by quadratic characters). It is an
  *infinite-time* ceiling: the finite-time amplitude is Diophantine-limited (recurrence time
  exponential in `N`).

**Nonlinear extension (Part II), focusing DNLS — now a rigorous chain.** The linear large wave **persists**
under weak nonlinearity (Duhamel); the **modulational-instability band**
`σ_Q = √(λ_Q (2γ A² − λ_Q))` is read off the `L_N` spectrum; localized
**breathers exist** (a finite-dimensional variational argument, no MacKay–Aubry continuation); in the
anti-continuum regime the site-centred breather is **orbitally stable** (Vakhitov–Kolokolov / Grillakis–
Shatah–Strauss) and the bond-centred one unstable, with both verdicts certified at `Ω = 2` in exact rational
arithmetic (`verify_dnls_breather_cap.py`, `N = 5` site, `N = 6` bond). Full dispersal of a single-site seed
is **energetically forbidden for `γP > 4`** (it would force `⟨L_N u, u⟩ → 2P − (γ/2) P² < 0`;
`4 = 2⟨L_N δ, δ⟩` is the infinite-lattice bandwidth), and for `γP > 8` the excitation stays on its site
(Shi–Wong prove this above `9`). The constant `4` is energy-sharp, not dynamically sharp: dispersal is proved
for `γP < 0.43`, and numerically the seed disperses up to `γP ≈ 3.4` and keeps a localized remnant from
`≈ 3.5`. The remaining strongly nonlinear regime (saturation, rogue-wave focusing), FPUT recurrence and the
route to chaos remain computational/exploratory.

**Open problems** (each reduced to a clean statement, with the proven part marked). (1) The Schrödinger
constant. Over infinite time: is `B_N < √N` for every odd composite `N` that is not a prime square (known for `N = 15, 21`), and is
`liminf B_N/√N = 1`? On the dispersive window: the windowed constant
`limsup_N sup_{t ≤ SN/2} ‖K_N(·,t)‖_1/√N = β_odd` for each fixed `S`. For each fixed `s ≤ 1` (odd `N`) the
profile maximizer is **proved** (`max F = β_odd`, validated interval arithmetic), and for `t > N/2` the entire
fixed-`s` `N → ∞` reduction is **rigorous** (joint equidistribution of the Debye phases *and* its uniform-in-`u`
transfer to `F(s)`); a supremum over `t` further needs convergence uniform in `s`. What remains beyond that is a
**single** validated inequality: the Gaussian-excess
(expected-modulus) bound `E|Z| − E|G| ≤ K⋆ · Σ_2^{−3/2} Σ_ℓ a_ℓ⁴` for the
profile's **balanced** Debye configs (live set = consecutive integers ⇒ even/odd counts differ by `≤ 1`;
sharp at the two-copy pair — imbalanced configs violate it), whose integrated majorant clears
`β_odd`. Whether polynomial windows `t ≤ N^C` keep `β_odd` is open. The continuum certificate is operational
(`verify_dagger_continuum.py`). The equal-amplitude singleton zone near `a_1 = a_3 ≈ 0.73` that the θ-ball
middle could not reach (its rigorous radius blows up past `t = 12`; the true ratio there is only `≈ 0.88`, so
it was enclosure-limited, not tight) is now covered by a **θ-Taylor + t-collocation** enclosure
(`middle_far_point`, validated against the reference in `verify_middle_far.py`): the previously-failing
boxes certify (e.g. `(0.70, 0.70) → 0.9915`; the hardest boxes sit at the high-λ edge `a_3 ≈ 0.61` and close
under one refinement, worst leaf `0.9999`). That validated 136-box tiling (`run_dagger_campaign.py`) took
the far `t`-collocation remainder and the stencil's third derivatives from finite differences and maximized
the quadratic model on a float grid. It is now **superseded by a rigorous certificate**
(`dagger_rigorous.py`, `run_dagger_campaign.py --rigorous`): in `u = a_1²`, `v = a_3²` the signed excess on
`[0, 24]` is an exact degree-10 Taylor polynomial of a Gauss–Legendre × periodic-midpoint quadrature with
analytic (Bernstein-ellipse, strip) error bounds and a Cauchy remainder, its box maximum minus `K⋆ D` is bounded
exactly, the tail uses `|J_0(x)| ≤ min(1, √(2/(πx)))`, and every bound is an Arb ball. All 111 boxes meeting the
band prove (worst `R ≤ 0.9205 K⋆`, 33 s on 8 cores); `verify_dagger_rigorous.py` re-derives every verdict from
the stored balls, recomputes the worst boxes bit for bit, checks the band coverage in exact rational arithmetic
and pins hashes and versions in `numerics/dagger_rigorous_manifest.json`. The near-pair
basin `λ ≥ 0.82` of the same three-copy family is proved by `dagger_nearpair.py`: the exact 1-D form
`E|Z| = 2/π + (1/π) ∫ ĝ(t) (1 − J_0(a_1 t) J_0(a_3 t)) t^{−2} dt` (slice `b = 1`) with the odd-copy kernel
`ĝ(s) = 2 − 2s I_1(s/2) K_0(s/2)`, a closed-form leading part, Gauss–Legendre with Bernstein-ellipse bounds,
certified convexity of `Φ` on a rectangle around the pair (308 Hessian boxes) and 923 value boxes elsewhere
(1231 boxes, ~75 s on 12 workers); `verify_dagger_nearpair.py` re-decides every box, recomputes the tightest
ones bit for bit, checks the coverage in exact rational arithmetic and pins hashes in
`numerics/dagger_nearpair_manifest.json`. So the inequality is **proved** (computer-assisted) for every
configuration with at most three copies; it stays validated on configurations with more than three copies. The
window part of problem (1) is reduced to the `s`-uniform transfer plus that piece. (2) The **excess lemma**
`A_N ≤ L_pre + O(1)` (equivalent to a rigorous *growing* deficit bound) — proved when the odd
part of `N` is a prime power (`N = 2p` sharp at `1/(2N)`, `N = 4p`), and **unconditional for any fixed `ω(m)`**
via the Kronecker sandwich `|A_N − L_pre| ≤ U_N − L_pre = (1/π) ln(m/φ(m))`. The open
content is an **absolute** constant **uniform as `ω(m) → ∞`**. A proved building block is the
**per-relation bound** (Lemma; `verify_excess_decomposition.py`): a single support-3 relation
`ω_c = ω_a + ω_b` contributes excess `g(b_a, b_b, b_c) ≤ b_c/2`, from
`g(α, β, γ) ≤ (γ²/2)(1/α + 1/β)` and `1/b_a + 1/b_b = 1/b_c`. The obstruction is the
uniform coupling of `Θ(N)` relations (comb generators + shared phases), the sup-side counterpart of
McGehee–Pigno–Smith/Konyagin in the Chowla cosine-extremum circle (where Jin–Milojević–Tomon–Zhang and
Bedert recently proved the first polynomial bounds, `n^{1/10−o(1)}` and `n^{1/5−o(1)}`, in the *opposite*
direction; Chowla's conjectured `√n` is still open). (3) A rigorous dispersal proof on `γP ∈ [0.43, 4)`: the constant `4` is
energy-sharp, not dynamically sharp — numerically the seed disperses only up to `γP ≈ 3.4` and keeps a
localized remnant from `≈ 3.5` (Johansson–Hörnquist–Riklund: `≈ 3.45`), so a proof can reach at most that
crossover, while the rigorous drop is proven only to `0.43`. The **leading-order drift now vanishes provably** (along the
linear flow the distance-2 current `Im(ā_n a_{n+2}) ≡ 0`, so the staggered momentum drifts
only at `O(γ²)` — Maxima + Rocq `StaggeredCurrent.v`, `verify_modified_energy_window.py`); what remains
is the `O(γ²)` remainder and the nonperturbative approach to the crossover. (4) Condition (B) of the order criterion for quasi-1D
lattices *without* an embedded ring (the Möbius ladder and multi-jump circulants are now settled; the
next-nearest-neighbour degeneracy law holds for all couplings `g` via Conway–Jones). **Priority note
(resolved):** Filimonov's 1992 C. R. Acad. Sci. note has been obtained and compared on the ring — the
**rational-independence dichotomy** (peaks for `N` prime/`2^m`, dependence for even composite, by his
cyclotomic argument) and the unbounded ceiling are Filimonov's; the rate `(1/π) ln N`, the explicit
constant, the **mod-4 saturation classification**, and the all-`N` rank `φ(2N)/2` are new here.

## Reproduce

Uses [`uv`](https://docs.astral.sh/uv/). Each `numerics/verify_*.py` is a PEP 723 standalone script
(inline dependencies), so run it with `uv run --script`.

```bash
uv run pytest                                  # unit tests for the core package
uv run --script numerics/verify_theorem1.py    # any verify_*.py re-checks one claim
uv run --script numerics/make_figures.py       # regenerate paper/figures/*.pdf
cd paper && latexmk -pdf large-wave-effect.tex # build the paper (run twice to resolve refs)
```

Headline checks: `verify_order_theorem.py` (the order law + palindromic prefix lemma),
`verify_classification_proof.py` (saturation classification), `verify_bn_saturation.py` (infinite-time
Schrödinger saturation), `verify_bn_parity.py` / `verify_beta_odd.py` / `verify_bn_profile_max.py` (the
windowed Schrödinger constant, its parity split, and the explicit `t ≤ N/2` profile maximizer), `verify_defect_certified.py` (certified `A_N < U_N`),
`verify_excess_smallcases.py` (excess lemma for prime-power odd part), `verify_product_order.py` and
`verify_kjump_order.py` (`C_N □ H`, Möbius ladder, `k`-jump saturation), the DNLS chain
`verify_dnls_persistence.py` / `verify_dnls_mi.py` / `verify_dnls_breather_stability.py`, and
`verify_selftrapping_transition.py` (the energy-sharp `γP = 4` bound and the onset of trapping near `3.5`).

For the Schrödinger continuum cover: `verify_middle_far.py` validates the θ-Taylor enclosure that closes
the `a_1 = a_3 ≈ 0.73` zone, `prove_box_parallel.py` certifies a single box on all cores, and the full band is
tiled by the adaptive, checkpointed campaign:
`uv run --script numerics/run_dagger_campaign.py --jobs 14` (writes each box to `dagger_campaign.jsonl`).

**Exact-arithmetic + formal layer** (Fedora: `dnf install gap pari-gp fricas rocq-prover`; Julia from juliaup). The
linear-theory arithmetic is reproduced floating-point-free in **GAP** and **PARI/GP**, symbolically in
**FriCAS**, the `N = 15` exclusion `B_15 < √15` is certified in **Julia** and re-checked in PARI/GP, the
`N = 21` exclusion `B_21 < √21` is certified with python-flint and re-checked in PARI/GP, and the
mod-4 criterion is machine-checked in **Rocq** — see [`exact/`](exact/) and
[`formal/rocq/`](formal/rocq/). One driver runs the lot and cross-checks against the Python core (missing
tools are skipped, not failed):

```bash
uv run --script numerics/verify_exact_layer.py   # GAP + PARI + Rocq + Python cross-check
gap -q --nointeract exact/gap/scan_classification.g   # finite exact cross-check to N<=800
rocq compile formal/rocq/Mod4Criterion.v              # machine-checked criterion certificate
```

## Layout

```
src/large_wave_effect/   shared spectral core (ring, cyclotomic, Schrödinger)
numerics/                verify_*.py (one per claim) + make_figures.py
exact/                   GAP / PARI / Julia / FriCAS exact + symbolic checks (gap/, pari/, julia/, fricas/)
formal/rocq/             Rocq machine-checked mod-4 criterion certificate
tests/                   pytest suite for the core package
paper/
  large-wave-effect.tex  root: preamble, front matter, \part dividers, bibliography
  sections/linear.tex    Part I  — linear theory
  sections/nonlinear.tex Part II — nonlinear extension
  figures/               generated by make_figures.py
  large-wave-effect.pdf  compiled paper
```

## Provenance

This work continues the author's 2015 undergraduate diploma on the large-wave effect; the original
C simulation that first observed the peaks is not included here.

## License

Apache License 2.0 — code, paper text, and figures. See `LICENSE` and `NOTICE`.

## Citation

Ilia Gradina, *Large-wave amplification on periodic discrete lattices: linear theory and a nonlinear
extension* (2026). DOI: [10.5281/zenodo.20795861](https://doi.org/10.5281/zenodo.20795861) (concept DOI —
always resolves to the latest version). See `CITATION.cff`.
