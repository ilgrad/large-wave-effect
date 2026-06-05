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

Proved results are kept strictly separate from numerical evidence and the one open conjecture.

**Linear theory (Part I).**

- **T1 — exact amplitude (saturating $N$).**
  $A_N = U_N = \tfrac1\pi\ln N + \tfrac1\pi\bigl(\gamma + \ln\tfrac2\pi\bigr) - \tfrac{\pi}{72N^2} + \cdots$,
  proved for $N$ prime and $N = 2^m$, with the same sharp constant $\tfrac1\pi$ on the composite family
  $N = 2p$ (cyclotomic independence + Bohr almost-periodicity).
- **T2 — subtorus dimension (theorem).** The orbit-closure subtorus has dimension exactly
  $\tfrac12\varphi(2N)$ for *every* $N$.
- **T3 — ceiling criterion (theorem).** $A_N = U_N$ **iff** every integer frequency relation has
  coefficient-sum $\equiv 0 \pmod 4$. On $N \le 160$ the saturating set is exactly
  $\{\text{primes}\}\cup\{2^m\}$ (density zero); whether any composite $N$ ever saturates is decided,
  per $N$, by the criterion.
- **T4 — discrete Schrödinger (two settings).** A localized state is never amplified (unitarity), yet the
  $\ell^\infty\!\to\!\ell^\infty$ amplification is $B_N = \Theta(\sqrt N)$ — **both bounds rigorous**
  (upper from unitarity, lower from a Bessel-front estimate at $t = N/8$). On the **Dirichlet segment**
  (uniform initial state) the amplitude at site $j$ grows as $d_j^N \sim \tfrac{4}{\pi^2}\ln j$ for $N$
  prime or $2^m$ (Filimonov 2023, reproduced — we identify the constant $C = 4/\pi^2$). Three large-wave
  laws thus live on one Laplacian: $A_N \sim \tfrac1\pi\ln N$ (ring, system size),
  $B_N \sim \sqrt N$ (ring, $\ell^\infty$ norm), $d_j^N \sim \tfrac{4}{\pi^2}\ln j$ (segment, site).
- **T5 — dimensionality & reachability.** The large wave lives exactly at spectral dimension $d_s = 1$
  (incl. quasi-1D ladders/tubes; $O(1)$ for $d_s \ge 2$). It is an *infinite-time* ceiling: the
  finite-time amplitude is Diophantine-limited (recurrence time grows exponentially in $N$), and mass
  disorder trades the arithmetic obstruction for an exploding recurrence time.
- **Open conjecture.** $A_N = \Theta(\ln N)$ for *all* $N$ (the sharp constant $\tfrac1\pi$ is proved only
  on the classes above). Numerically $A_N/\ln N \in [0.28, 0.34]$ (including the hardest $N = pq$,
  $p,q \sim \sqrt N$), and the alignable independent low-mode prefix carries a budget $\sim 0.26\ln N$ —
  strong evidence, but **not** a proof: the out-of-prefix modes may interfere. The Dirichlet segment is
  full-rank ($A_N = U_N$) iff $N+1$ is prime or $2^m$.

**Nonlinear extension (Part II).** The same $L_N$ is the linear part of the focusing DNLS and the FPUT
lattice. Peregrine factor $3$, Benjamin–Feir modulational instability, heavy-tailed (leptokurtic)
statistics; self-trapping into discrete breathers above the band top ($\gamma \approx 4$); FPUT
recurrence and a Lyapunov transition to chaos; and a third, refractive route — Airy-fold and cusp
(Pearcey) caustics, branched flow in random currents, and lightning as Laplacian growth — all on the
same $L_N$. Part II is
primarily computational/exploratory; its rigorous statements (Hamiltonian conservation laws, the
modulational-instability band) are flagged as classical, while the hard nonlinear results (breather
existence/stability, the precise self-trapping threshold, the route to chaos) are simulated, not proved.

## Reproduce

Uses [`uv`](https://docs.astral.sh/uv/). Each `numerics/verify_*.py` is a PEP 723 standalone script
(inline dependencies), so run it with `uv run --script`.

```bash
uv run pytest                                  # unit tests for the core package
uv run --script numerics/verify_theorem1.py    # any verify_*.py re-checks one claim
uv run --script numerics/make_figures.py       # regenerate paper/figures/*.pdf
cd paper && latexmk -pdf large-wave-effect.tex # build the paper (run twice to resolve refs)
```

Headline checks: `verify_ceiling_criterion.py` (T3), `verify_qrank_formula.py` (T2),
`verify_schrodinger_lower.py` (T4), `verify_conj_order_map.py` (order conjecture),
`verify_relation_lattices.py` (the mod-4 obstruction made explicit).

## Layout

```
src/large_wave_effect/   shared spectral core (ring, cyclotomic, Schrödinger)
numerics/                verify_*.py (one per claim) + make_figures.py
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
