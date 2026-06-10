# large-wave-effect

Rigorous study of the **large-wave effect** on periodic discrete chains: a localized impulse on a
homogeneous one-dimensional *lattice* amplifies into a peak that grows with the number of nodes — a
genuinely discrete resonance with no continuum counterpart.

One paper (`paper/`), in two parts built on a shared spectral core (the ring Laplacian $L_N$,
Toeplitz/circulant methods):

- **Part I — Linear theory** (rigorous core). Amplitude $A_N = \sup_t \max_j |u_j(t)|$ of the velocity
  Green's function on the ring; cyclotomic independence, the mod-4 ceiling criterion, the discrete
  Schrödinger contrast, spectral zeta, operator norms, statistics, and the dimensional threshold.
- **Part II — Nonlinear extension**. The discrete NLS (Peregrine soliton, modulational instability,
  rogue-wave statistics), discrete breathers and self-trapping, FPUT recurrence, and the route to
  chaos (Lyapunov exponent).

Every quantitative claim is reproduced by a standalone, version-controlled script.

## Main results

Proved results are kept strictly separate from numerical evidence and the (sharply reduced) open problems.

**Linear theory (Part I).**

- **Order law (theorem, all $N$).** $A_N \sim \tfrac1\pi\ln N$ for *every* $N$ — the central result, formerly
  a conjecture. The lower bound comes from a *palindromic argument*: the first $\tfrac12\varphi(2N)$
  frequencies $2\sin(\pi r/N)$ are $\mathbb{Q}$-independent, so a $(1-o(1))$ fraction of the ceiling can
  always be Kronecker-aligned. Precisely $A_N = \tfrac1\pi\ln N + O(\ln\ln\ln N)$.
- **Exact amplitude / saturation (theorem).**
  $A_N = U_N = \tfrac1\pi\ln N + \tfrac1\pi(\gamma+\ln\tfrac2\pi) - \tfrac{\pi}{72N^2}+\cdots$ **iff** $N$ is
  prime or a power of two; equivalently every integer frequency relation has coefficient-sum $\equiv 0
  \pmod 4$ (antipodal variant for even $N$). Via an explicit root-of-unity obstruction,
  $\{N : A_N = U_N\} = \{\text{primes}\}\cup\{2^m\}$ (density zero). The orbit-closure subtorus has
  dimension exactly $\tfrac12\varphi(2N)$. Certified to $N\le600$ (Python) / $N\le800$ (exact GAP).
- **Composite deficit.** $U_N - A_N = O(\ln\ln\ln N)$, depending asymptotically only on the odd part of $N$;
  the strict gap $A_N < U_N$ is rigorously certified for small $N$ (Lipschitz-grid and Lasserre moment–SOS).
- **Discrete Schrödinger.** A localized state is never amplified (unitarity), yet
  $B_N = \sup_t\|e^{-itL_N}\|_{\ell^\infty\to\ell^\infty} = \Theta(\sqrt N)$. Here
  $\liminf B_N/\sqrt N \ge c_0/\sqrt2 = 0.861$ is proved ($c_0=(2/\pi)^{3/2}B(\tfrac12,\tfrac34)$; attained
  numerically on even $N$), and the constant splits by parity — numerically even $N \to c_0/\sqrt2$, odd
  $N \to \beta_{\mathrm{odd}} = 0.928\ldots$, an elliptic-integral constant for which no elementary closed
  form was found — so $\limsup \ge \beta_{\mathrm{odd}}$. On the **Dirichlet segment**
  $d_j^N \sim \tfrac{4}{\pi^2}\ln j$ for $N+1$ prime/$2^m$ (the constant $C=4/\pi^2$ is the published
  Myshkis–Filimonov 2003 value). Three large-wave laws on one Laplacian.
- **Dimension & products (theorem).** The large wave lives exactly at spectral dimension $d_s = 1$ ($O(1)$
  for $d_s\ge2$). A graph-agnostic order criterion plus the embedded ring spectrum give $A = \Theta(\ln N)$
  on *every* Cartesian product $C_N\,\square\,H$ (ladders, prisms, tubes). It is an *infinite-time* ceiling:
  the finite-time amplitude is Diophantine-limited (recurrence time exponential in $N$).

**Nonlinear extension (Part II), focusing DNLS — now a rigorous chain.** The linear large wave **persists**
under weak nonlinearity (Duhamel); the **modulational-instability band** $\sigma_Q=\sqrt{\lambda_Q(2\gamma
A^2-\lambda_Q)}$ is read off the $L_N$ spectrum; the on-site **breather exists** (a finite-dimensional
variational argument, no MacKay–Aubry continuation) and is **stable** (Vakhitov–Kolokolov / Grillakis–
Shatah–Strauss), while the bond-centred breather is unstable. The strongly nonlinear regime (saturation,
the precise self-trapping threshold $\approx 4$, rogue-wave focusing), FPUT recurrence and the route to
chaos remain computational/exploratory.

**Open problems** (each reduced to a clean statement). (1) The *sharp* $B_N$ upper bound on $t>N/2$
(proved for $t\le N/2$): an anti-flatness statement for the non-quadratic $\sin^2$ chirp / joint
equidistribution of the multi-copy Debye phases, which would give $\limsup B_N/\sqrt N=\beta_{\mathrm{odd}}$.
(2) A rigorous *growing* lower bound on the deficit, equivalent to the **excess lemma** $A_N\le L_{\mathrm{pre}}+O(1)$.
(3) The strongly nonlinear DNLS (above). (4) Condition (B) of the order criterion for quasi-1D lattices
*without* an embedded ring factor. **Priority note:** Filimonov's 1992 C. R. Acad. Sci. note must be obtained
and compared on the ring before claiming priority (its splash table is reproduced in Andrianov–Awrejcewicz–
Danishevskyy 2021).

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
`verify_classification_proof.py` (saturation classification), `verify_bn_parity.py` and
`verify_beta_odd.py` (the Schrödinger constant and its parity split), `verify_defect_certified.py`
(certified $A_N<U_N$), `verify_product_order.py` ($C_N\,\square\,H$), and the DNLS chain
`verify_dnls_persistence.py` / `verify_dnls_mi.py` / `verify_dnls_breather_stability.py`.

**Exact-arithmetic + formal layer** (Fedora: `dnf install gap pari-gp fricas rocq-prover`). The
linear-theory arithmetic is reproduced floating-point-free in **GAP** and **PARI/GP**, symbolically in
**FriCAS**, and the mod-4 criterion is machine-checked in **Rocq** — see [`exact/`](exact/) and
[`formal/rocq/`](formal/rocq/). One driver runs the lot and cross-checks against the Python core (missing
tools are skipped, not failed):

```bash
uv run --script numerics/verify_exact_layer.py   # GAP + PARI + Rocq + Python cross-check
gap -q --nointeract exact/gap/scan_classification.g   # A_N=U_N <=> prime/2^m, verified to N<=800
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
