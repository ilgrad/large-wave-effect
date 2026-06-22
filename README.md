# large-wave-effect

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
  the saturation classification; quantitative small-`N` gaps are certified by Lipschitz-grid and
  Lasserre moment–SOS checks. The sharper excess law
  `U_N − A_N = (1/π) ln(m/φ(m)) + O(1)` for odd part `m` remains open beyond prime-power odd part. The
  per-family limits `defect_∞(m) = lim_a (U − A)(m · 2^a)` are computed to `10⁻⁸`
  (`0.10368…`, `0.03824…`) as numerical constants with no closed form.
- **Discrete Schrödinger.** A localized state is never amplified (unitarity), yet
  `B_N = sup_t ‖e^{−it L_N}‖_{ℓ^∞ → ℓ^∞} = Θ(√N)`. Here
  `liminf B_N/√N ≥ c_0/√2 = 0.861` is proved (`c_0 = (2/π)^{3/2} · B(1/2, 3/4)`; attained
  numerically on even `N`), and the constant splits by parity — numerically even `N → c_0/√2`, odd
  `N → β_odd = 0.928…`, an elliptic-integral constant for which no elementary closed
  form was found — so `limsup ≥ β_odd`. For `t ≤ N/2` the limiting `ℓ¹` profile is the
  *explicit* elliptic-integral function `F(2t/N)` (one-copy term `A(s) = ₂F₁` in closed form), and its
  maximizer `max F = F(1) = β_odd` is now **fully proved** — analytically on `[1/2, 0.937]`
  (Jensen, sharpened by the tangent to a concave `h`) and by **validated interval arithmetic** on
  `[0.937, 1]` (`F′(s) > 0`, worst case `0.051`). So the `t ≤ N/2` contribution to `B_N` is settled. For
  `t > N/2` the whole `N → ∞` reduction is now **rigorous** — joint equidistribution of the multi-copy
  Debye phases *and* its uniform-in-`u` transfer to the `ℓ¹` constant `F(s)` (a mollified Koksma–Hlawka /
  cubic van der Corput estimate; the modulus has divergent Hardy–Krause variation) — leaving a **single**
  validated residual, the Gaussian-excess inequality `E|Z| − E|G| ≤ K⋆ · Σ_2^{−3/2} Σ_ℓ a_ℓ⁴`
  (`K⋆ = (h(1) − √π/2)/2`) for the profile's **balanced** Debye configurations (the live set is an
  interval of consecutive integers, so the even/odd channel counts differ by `≤ 1`; sharp at the two-copy
  pair — imbalanced configs violate the unrestricted bound) whose integrated majorant clears
  `β_odd` — and whose residual is now **validated across the whole band** (the θ-Taylor enclosure `middle_far_point`
  plus a 130-box rigorous-enclosure tiling, 0 failures), so `limsup B_N/√N = β_odd` holds modulo this one
  numerically-validated inequality. On the **Dirichlet segment**
  `d_j^N ~ (4/π²) ln j` for `N+1` prime/`2^m` (the constant `C = 4/π²` is the published
  Myshkis–Filimonov 2003 value). Three large-wave laws on one Laplacian.
- **Dimension, products & circulants (theorem).** The large wave lives exactly at spectral dimension
  `d_s = 1` (`O(1)` for `d_s ≥ 2`). A graph-agnostic order criterion plus the embedded ring spectrum give
  `A = Θ(ln N)` on *every* Cartesian product `C_N □ H` (ladders, prisms, tubes) and, beyond
  products, on the **Möbius ladder** `C_{2N}({1, N})` via its embedded-ring branch. Multi-jump circulants
  `C_N(S)` (`|S| ≥ 2`) instead **saturate** (`A_N = U_N`) off explicit Conway–Jones collision congruences —
  the nested-radical frequencies destroy the ring's mod-4 obstruction; for the nearest case `S = {1, 2}` this
  is **proved exactly** (Besicovitch–Kummer: no subset product of the `α_r` is a square in
  `ℚ(ζ_p)`) for collision-free primes `p ≤ 23` (`exact/pari/kjump_kummer.gp`). It is an
  *infinite-time* ceiling: the finite-time amplitude is Diophantine-limited (recurrence time
  exponential in `N`).

**Nonlinear extension (Part II), focusing DNLS — now a rigorous chain.** The linear large wave **persists**
under weak nonlinearity (Duhamel); the **modulational-instability band**
`σ_Q = √(λ_Q (2γ A² − λ_Q))` is read off the `L_N` spectrum; the on-site
**breather exists** (a finite-dimensional
variational argument, no MacKay–Aubry continuation) and is **stable** (Vakhitov–Kolokolov / Grillakis–
Shatah–Strauss), while the bond-centred breather is unstable. The **self-trapping threshold is proved
energy-sharp** at the exact constant `γP = 4 = ‖L_N‖`: above it, full dispersal of the single-site seed
is energetically forbidden (it would force `⟨L_N u, u⟩ → 2P − (γ/2) P² < 0`). The remaining
strongly nonlinear regime (saturation, rogue-wave focusing), FPUT recurrence and the route to chaos remain
computational/exploratory.

**Open problems** (each reduced to a clean statement, with the proven part marked). (1) The *sharp* `B_N`
upper bound. The `t ≤ N/2` case is **fully proved** (`max F = β_odd`, analytic + validated
interval arithmetic), and for `t > N/2` the entire `N → ∞` reduction is now **rigorous** (joint
equidistribution of the Debye phases *and* its uniform-in-`u` transfer to `F(s)` — a mollified Koksma–Hlawka /
cubic van der Corput estimate). What remains is a **single** validated inequality: the Gaussian-excess
(expected-modulus) bound `E|Z| − E|G| ≤ K⋆ · Σ_2^{−3/2} Σ_ℓ a_ℓ⁴` for the
profile's **balanced** Debye configs (live set = consecutive integers ⇒ even/odd counts differ by `≤ 1`;
sharp at the two-copy pair — imbalanced configs violate it), whose integrated majorant clears
`β_odd`. The continuum certificate is operational
(`verify_dagger_continuum.py`). The equal-amplitude singleton zone near `a_1 = a_3 ≈ 0.73` that the θ-ball
middle could not reach (its rigorous radius blows up past `t = 12`; the true ratio there is only `≈ 0.88`, so
it was enclosure-limited, not tight) is now covered by a **θ-Taylor + t-collocation** enclosure
(`middle_far_point`, validated against the reference in `verify_middle_far.py`): the previously-failing
boxes certify (e.g. `(0.70, 0.70) → 0.9915`; the hardest boxes sit at the high-λ edge `a_3 ≈ 0.61` and close
under one refinement, worst leaf `0.9999`). The full singleton-band
cover is now **validated across the whole band** — the adaptive, node-parallel, checkpointed tiling `run_dagger_campaign.py`
certified all 130 leaf boxes with 0 failures — so `limsup B_N/√N = β_odd` holds at the
dagger certificate's validated-enclosure standard, leaving problem (1) reduced to upgrading that one numerically-validated inequality to a proof. (2) The **excess lemma**
`A_N ≤ L_pre + O(1)` (equivalent to a rigorous *growing* deficit bound) — proved when the odd
part of `N` is a prime power (`N = 2p` sharp at `1/(2N)`, `N = 4p`), and **unconditional for any fixed `ω(m)`**
via the Kronecker sandwich `|A_N − L_pre| ≤ U_N − L_pre = (1/π) ln(m/φ(m))`. The open
content is an **absolute** constant **uniform as `ω(m) → ∞`**. A proved building block is the
**per-relation bound** (Lemma; `verify_excess_decomposition.py`): a single support-3 relation
`ω_c = ω_a + ω_b` contributes excess `g(b_a, b_b, b_c) ≤ b_c/2`, from
`g(α, β, γ) ≤ (γ²/2)(1/α + 1/β)` and `1/b_a + 1/b_b = 1/b_c`. The obstruction is the
uniform coupling of `Θ(N)` relations (comb generators + shared phases), the sup-side counterpart of
McGehee–Pigno–Smith/Konyagin in the Chowla cosine-extremum circle (recently resolved in the *opposite*
direction by Jin–Milojević–Tomon–Zhang / Bedert). (3) A rigorous dispersal proof on the **proof gap** `γP ∈ [0.43, 4)`: the
threshold `4` is now energy-sharp and dispersal is the numerically established behaviour up to `≈ 4`,
but the rigorous drop is proven only to `0.43`. The **leading-order drift now vanishes provably** (along the
linear flow the distance-2 current `Im(ā_n a_{n+2}) ≡ 0`, so the staggered momentum drifts
only at `O(γ²)` — Maxima + Rocq `StaggeredCurrent.v`, `verify_modified_energy_window.py`); what remains
is the `O(γ²)` remainder and the nonperturbative approach to `4`. (4) Condition (B) of the order criterion for quasi-1D
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
`verify_classification_proof.py` (saturation classification), `verify_bn_parity.py` /
`verify_beta_odd.py` / `verify_bn_profile_max.py` (the Schrödinger constant, its parity split, and the
explicit `t ≤ N/2` profile maximizer), `verify_defect_certified.py` (certified `A_N < U_N`),
`verify_excess_smallcases.py` (excess lemma for prime-power odd part), `verify_product_order.py` and
`verify_kjump_order.py` (`C_N □ H`, Möbius ladder, `k`-jump saturation), the DNLS chain
`verify_dnls_persistence.py` / `verify_dnls_mi.py` / `verify_dnls_breather_stability.py`, and
`verify_selftrapping_transition.py` (the energy-sharp `γP = 4` threshold).

For the Schrödinger continuum cover: `verify_middle_far.py` validates the θ-Taylor enclosure that closes
the `a_1 = a_3 ≈ 0.73` zone, `prove_box_parallel.py` certifies a single box on all cores, and the full band is
tiled by the adaptive, checkpointed campaign:
`uv run --script numerics/run_dagger_campaign.py --jobs 14` (writes each box to `dagger_campaign.jsonl`).

**Exact-arithmetic + formal layer** (Fedora: `dnf install gap pari-gp fricas rocq-prover`). The
linear-theory arithmetic is reproduced floating-point-free in **GAP** and **PARI/GP**, symbolically in
**FriCAS**, and the mod-4 criterion is machine-checked in **Rocq** — see [`exact/`](exact/) and
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
exact/                   GAP / PARI / FriCAS exact + symbolic cross-checks (gap/, pari/, fricas/)
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
extension* (2026). See `CITATION.cff`. A Zenodo DOI will be minted on the first release (connect the
repository to Zenodo, then publish a GitHub release).
