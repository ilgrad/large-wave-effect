# A1 — the sharp-deficit ("excess") problem for the discrete-ring large wave

*A self-contained open-problem brief: statement, the reduction we have proved, the numerical evidence, an
honest log of every attack that failed, the literature placement, and an alternative ergodic-optimization
formulation. The rest of the* large-wave-effect *project is proved; this is the single hard kernel. It is a
sup-side ("primal") analogue of the McGehee–Pigno–Smith / Konyagin theorem for a specific cyclotomic comb,
and (literature search through 2026, §6) no off-the-shelf tool supplies it.*

---

## 1. Setup

Discrete ring $C_N$ (cycle on $N$ vertices), nearest-neighbour. The node-$0$ velocity response to a unit
impulse is
$$ u_0(t)=\sum_{r=1}^{\lfloor (N-1)/2\rfloor} b_r\,\sin(\omega_r t),\qquad
   \omega_r = 2\sin\!\frac{\pi r}{N},\qquad b_r=\frac{1}{N\sin(\pi r/N)}>0 . $$
Node $0$ is the focusing site, so the **large-wave amplitude** is $A_N=\sup_{t\ge 0} u_0(t)$.

- **Ceiling** $U_N=\sum_r b_r=\frac1{2N}\sum_{r=1}^{N-1}\csc\frac{\pi r}{N}=\frac1\pi\ln N+O(1)$ — the value
  if every $\sin(\omega_r t)$ could equal $1$ simultaneously.
- By **Kronecker–Weyl**, the orbit closure of $t\mapsto(\omega_r t\bmod 2\pi)_r$ is the subtorus
  $$ \mathbb T_N=\Big\{\theta:\ \textstyle\sum_r c_r\theta_r\equiv 0\ (2\pi)\ \text{ whenever } \sum_r c_r\omega_r=0,\ c\in\mathbb Z\Big\},$$
  so $A_N=\max_{\theta\in\mathbb T_N}\sum_r b_r\sin\theta_r$.
- The $\mathbb Q$-rank of $\{\omega_r\}$ is $M:=\tfrac12\varphi(2N)$ (proved). The $M$ largest $b_r$ (indices
  $r=1,\dots,M$, the **prefix**) carry $\mathbb Q$-independent frequencies and free phases; the rest (the
  **dependent band**, $\Theta(N)$ modes for composite $N$) are fixed integer combinations of the prefix.
- **Prefix budget** $L_{\mathrm{pre}}=\sum_{r=1}^{M} b_r$ (the free modes can each reach $\sin=1$), and
  $$ U_N-L_{\mathrm{pre}}=\sum_{r>M} b_r=\tfrac1\pi\ln\frac{m}{\varphi(m)}+O(1),\qquad m=\text{odd part of }N. $$
  Along primorials $m/\varphi(m)\sim e^{\gamma}\ln\ln m$, so $U_N-L_{\mathrm{pre}}=\Theta(\ln\ln\ln N)$ —
  unbounded, but slower than any iterate of $\log$.

Saturation $A_N=U_N$ holds **iff** $N$ is prime or $2^m$ (the $\omega_r$ are then $\mathbb Q$-independent);
this and the order law $A_N\sim\frac1\pi\ln N$ are proved.

## 2. The problem (A1)

Define the **excess** $E_N:=A_N-L_{\mathrm{pre}}\ge 0$. Trivially $0\le E_N\le U_N-L_{\mathrm{pre}}$.

> **Question A1.** As $\omega(m)\to\infty$ (number of distinct primes of the odd part $m$), which holds?
> - **(O1)** *bounded*: $E_N\le C$ absolute — equivalently the sharp deficit law
>   $U_N-A_N=\frac1\pi\ln\frac{m}{\varphi(m)}+O(1)$;
> - **(L)** *logarithmic*: $E_N=c\,(U_N-L_{\mathrm{pre}})+o(\cdot)$ with $c\in(0,1)$ — then the sharp deficit
>   constant is $(1-c)/\pi$, **not** $1/\pi$.

Equivalently: **find the sharp constant $\kappa$ in
$U_N-A_N=\frac{\kappa}{\pi}\ln\frac{m}{\varphi(m)}+O(1)$** ($\kappa=1$ under (O1); $\kappa=1-c<1$ under (L)).
Which alternative holds is open (§4, §9): the numerics are inconclusive, and the disorder of the optimizer
(§9) keeps **both** alternatives live — (O1) via a global energy–entropy bound, (L) only as a spin-glass-type
pressure.

It is a **comparison on the relation-constrained subtorus**: by how much can the constrained maximum
$\max_{\theta\in\mathbb T_N}\sum b_r\sin\theta_r$ exceed the free-prefix maximum $L_{\mathrm{pre}}$? It is the
**sup-side / primal analogue of MPS–Konyagin** — an *upper* bound on a *sup*, driven by the relation lattice.

## 3. What is proved (reduction to a comb-cancellation kernel)

- **Bounded $\omega(m)$ is trivial.** Kronecker gives
  $|E_N|\le U_N-L_{\mathrm{pre}}=\frac1\pi\ln\frac m{\varphi(m)}+O(1)$, bounded for fixed $\omega(m)$. Content
  appears only as $\omega(m)\to\infty$.
- **Per-relation bound.** For $\alpha,\beta,\gamma>0$,
  $\max_{u,v}[\gamma\sin(u+v)-\alpha(1-\sin u)-\beta(1-\sin v)]\le\frac{\gamma^2}{2}(\frac1\alpha+\frac1\beta)$.
  For a support-$3$ relation $\omega_c=\omega_a+\omega_b$, the identity $1/b_a+1/b_b=1/b_c$ makes the
  per-relation excess $\le b_c/2$.
- **Clean relations for $3\mid N$ (via Conway–Jones).** The only $\omega_a+\omega_b=\omega_c$
  ($a<b<c\le N/2$) are $a+b=N/3$, $c=N/3+a$; the triplets $(a,N/3-a,N/3+a)$ are vertex-disjoint.
- **Chain-split.** For $3\mid N$, for any $\lambda\in[0,1)$,
  $$ E_N\le \frac1{1-\lambda}\sum_{c\in T}\frac{b_c}2 \;+\; F_{\mathrm{comb}}(\lambda),\qquad
     F_{\mathrm{comb}}(\lambda)=\max_\psi\Big[\sum_{r\in K}b_r\sin\phi_r-\lambda\sum_{s\le M}b_s(1-\sin\psi_s)\Big], $$
  where $\phi_r=(C\psi)_r$ are the dependent phases (integer combinations of the free $\psi$), $T$ = the
  disjoint clean triplets (their sum is an $O(1)$ Riemann sum $\to\frac{\ln3}{4\pi}$ for square-free $N=3p$),
  and $K$ = the **higher-support comb generators**.
- **Square-free $N=3p$ is closed**: the comb block is a single generator (period-$3$ telescoping of
  $1+\zeta_3+\zeta_3^2=0$), giving $E_N\le\frac{\ln3}{4\pi}+o(1)$ unconditionally.

> **Reduced kernel.** For $3\mid N$, (O1) holds **iff** $F_{\mathrm{comb}}(\lambda)\le C'$ for some
> $\lambda$ — a *signed-cancellation upper bound on the comb generators*, with content only once
> $\omega(m)\ge 3$. This is the whole of A1 (the primal MPS–Konyagin analogue).

The combs are explicit: for $5\mid N$ they are support-$5$ "dilated triplets" from $2\sin(\pi/10)=(\sqrt5-1)/2$;
the second prime $p_2$ supplies the longer combs; $\#K=O(p_2)$.

## 4. Numerical evidence (two-edged; the dichotomy is computationally undecidable)

- $E_N\in[0.05,0.12]$ across **every** tested $N$ ($\omega(m)\le4$, up to $N=1155$, prefix dim $240$) → looks
  bounded, favouring **(O1)**.
- But the **ratio** $E_N/(U_N-L_{\mathrm{pre}})$ (against the *finite* $U_N-L_{\mathrm{pre}}$) is roughly
  stable $\approx0.31$ across $\omega(m)=2,3,4$ (it would $\to0$ under (O1)) → favours **(L)**. *Caveat:*
  against the *asymptotic* $\frac1\pi\ln\frac m{\varphi(m)}$ the ratio is **not** stable (it ranges
  $0.33$–$0.88$ and sits above $1/\pi$); at reachable $N$ the excess is $O(1)$-dominated, so the asymptotic
  $\kappa$ is not cleanly extractable.
- The rigorous split upper bound $\min_\lambda[\frac1{1-\lambda}\sum_T\frac{b_c}2+F_{\mathrm{comb}}(\lambda)]$
  stays $\approx0.19$ flat to $N=1155$ → favours (O1). The two signals point opposite ways and neither turns
  over by $\omega(m)=4$.
- **Undecidable at any feasible $N$**: $U_N-L_{\mathrm{pre}}=\frac1\pi\ln\frac m{\varphi(m)}=\Theta(\ln\ln\ln N)$
  grows so slowly ($0.28$ at $\omega=4$, only $0.31$ at $\omega=6$) that $O(1)$ vs $\Theta(\ln\ln\ln)$ cannot be
  separated; $\omega=5$ ($N=15015$) already puts the exact cyclotomic relation matrix out of reach.
- The optimizer is **non-perturbative**: prefix phases sit at $|\theta-\pi/2|\sim0.8$ (far from the aligned
  point), and $E_N$ is a *balance* of a positive dependent net ($\approx0.17$) against a negative prefix
  deficit ($\approx-0.08$) — a large-amplitude coherent rearrangement, not near-total sign cancellation and
  not a small perturbation.

## 5. What was tried and failed (so it is not repeated)

**A. Duality / spectral / geometry-of-numbers.**
1. *Per-relation Cauchy–Schwarz + Bessel/frame SDP packing* ($\sum_m\lambda_m\hat u_m\hat u_m^\top\preceq I$).
   FAILS: $\mathrm{OPT}_{\mathrm{SDP}}/(U_N-L_{\mathrm{pre}})$ **decays** with $\omega(m)$ ($1.0\to0.4\to0.15$);
   per-relation C–S discards phase coherence (modulus loss); heavy overlap (max degree $\approx271$) makes the
   frame constraint bind.
2. *Signed quadratic projection $\to$ CVP.* The quadratic-model deficit equals
   $\min_{n\in\mathbb Z^r}(w+2\pi n)^\top M^{-1}(w+2\pi n)$, $M=K\operatorname{diag}(1/b)K^\top$ the relation
   Gram, $K=[-C\mid I]$, $w$ = mod-$4$ defect vector; $(2/\pi^2)\,w^\top M^{-1}w\le$ deficit. Numerically
   $w^\top M^{-1}w/(U_N-L_{\mathrm{pre}})$ is stable $\sim1.5$–$2.5$ (suggesting (L)), but: (i) only
   qualitative, not the sharp constant; (ii) the true optimizer is non-perturbative so the quadratic model
   overshoots $\sim2\times$; (iii) it reduces to a **structured spectral lower bound on the small eigenvalues
   of the cyclotomic Gram $M$ in the direction $w$** — genuine analytic number theory, not hand-provable.
3. *Transference (Banaszczyk).* Bounds the **covering radius** (an *upper* bound on $\mathrm{dist}(w,\cdot)$) —
   wrong direction; we need a *lower* bound on the distance of the *specific* point $w$.

**B. Chowla / harmonic-analysis circle.**
4. *MPS / Konyagin / Bloom–Green.* These are $L^1$ **lower** bounds (Littlewood) — wrong side; A1 needs an
   $L^\infty$ (sup) **upper** bound.
5. *Chowla cosine problem* (Bedert `arXiv:2509.05260`; Jin–Milojević–Tomon–Zhang `arXiv:2509.03490`, spectral
   $\to$ cut via Cayley-graph eigenvalues). LOWER-bounds the magnitude of the negative dip (small eigenvalue
   $\Rightarrow$ clique) — wrong direction; and $\omega_r=2\sin(\pi r/N)$ are **not characters** of a cyclic
   group, so there is no Cayley transcription.
6. *Cut-norm / Grothendieck (Alon–Naor).* The matching tool for Chowla; same wrong direction. Also the
   dependent band's **own** supremum is its full weight ($\max_\psi\sum_{r>M}b_r\sin\phi_r=\sum_{r>M}b_r$), so
   the bound is forced by the prefix **trade-off**, not by bounding any single sup.

**C. Signed-cancellation / structural guesses (all falsified numerically).**
7. *Möbius/Mertens divisor decomposition.* Conjectured the comb excess inherits a signed divisor cancellation
   $\to(\varphi(m)/m)\sum_{p\mid m}\frac{\log p}{p-1}\to e^{-\gamma}$ (bounded). FALSIFIED: $E_N/M(m)$ creeps
   up ($0.131\to0.184$) while $E_N/(U_N-L_{\mathrm{pre}})$ stays flat — the excess tracks the **growing**
   envelope, not the bounded Mertens constant.
8. *Linear-response normal form.* Expand around $x_0=\frac\pi2\mathbf1$: $E_N\approx C_0+\frac12 g^\top W^{-1}g$,
   $g_i=\sum_{\rho>M}b_\rho\cos\phi_\rho\,\alpha_{\rho,i}$, $W=\operatorname{diag}(b_{\mathrm{pre}})$. Predicted
   $g^\top W^{-1}g/(U-L)\to2/\pi$. FALSIFIED: it **explodes** ($\approx98$ at $N=1155$) via small-$b$ boundary
   amplification; the rigorous along-gradient lower bound captures a *shrinking* fraction ($10\%$ at
   $\omega=4$). Confirms the optimizer is non-perturbative.
9. *$L^1$ back-projection.* $E_N\le\|v_{\mathrm{pre}}\|_1=\sum_r|\sum_{k>M}b_k c_{k,r}|$ (heuristic $L^1$
   duality). The bound is valid but **loose** ($5$–$13\times$) and **grows** $\sim0.3\,\omega(m)=O(\ln\ln N)$ —
   exponentially **worse** than the trivial sandwich $O(\ln\ln\ln N)$. No $O(1)$ collapse.

**D. Per-relation / combinatorial / dynamical.**
10. *Support-$k$ generalization of the triplet bound.* For higher-support combs the per-comb quadratic bound
    $\sum_j k_j^2/b_{a_j}=\Theta(N)$ **diverges**. (The triplet works only via the collapse $2\sin(\pi/6)=1$
    giving $1/b_a+1/b_b=1/b_c$; for $5\mid N$ the analogue $2\sin(\pi/10)$ is golden $\ne1$, no collapse.)
11. *Ablowitz–Ladik conserved functional* (staggered-momentum analogue). Does **not** transfer — the combs
    align in isolation; the cancellation is the joint prefix trade-off, not a conserved invariant.

**E. Computational bypass.**
12. *Direct GPU $\sup_t$* (compute $A_N$ for large $N$ by FFT, avoiding the cyclotomic machinery). DEAD: the
    excess ($\sim0.1$) lies **below the finite-time sampling floor** (the gap $U_N-\sup_{t\le T}$ is $\sim0.57$,
    $6\times$ the excess, already at $\omega=3$); approaching $A_N$ needs recurrence time $\exp(\Theta(m))$.

**F. Explicit lower-bound (feasible-point) constructions — the right direction, but the optimum resists.**
13. *mod-$4$ phase-CSP.* Discretize $\psi_s=\frac\pi2+\frac\pi2 x_s$, $x_s\in\mathbb Z/4$: the dependent reward
    becomes $\sigma(S_r+(Cx)_r)\in\{0,\pm1\}$ ($S_r=\sum_s C_{rs}$, $\sigma(1){=}1,\sigma(3){=}{-}1,\sigma(0){=}\sigma(2){=}0$)
    and the prefix cost $q(x_s)\in\{0,1,2\}$, a clean weighted Max-CSP whose value is a *rigorous lower bound* on
    $E_N$. But rounding the true optimizer to the grid and amplitude-refining ($\psi=\frac\pi2+\alpha\frac\pi2 x$)
    captures only $20$–$50\%$ of $E_N$ ($0.23,0.49,0.20$ at $\omega=2,3,4$); at the raw grid the value is
    *negative*. Reason: the optimum sits at $|\theta-\pi/2|\sim0.8\approx\pi/4$, a half-step *between* mod-$4$
    nodes; the rounded pattern is sparse with no clean residue/$\gcd$ structure. (A finer grid is faithful but
    loses the clean $\{0,\pm1\}$ arithmetic — no clean *analyzable* discretization is faithful.)
14. *Random-phase mean-field certificate.* For independent $\psi_s=\frac\pi2+\eta_s$ one has the rigorous bound
    $E_N\ge\mathbb E_\rho F_N$ with the *exact* Fourier-moment form
    $\mathbb E_\rho F_N=\sum_{r>M}b_r\,\mathrm{Im}\!\big(i^{S_r}\textstyle\prod_s\hat\rho_{\tau(s)}(C_{rs})\big)-\sum_{s\le M}b_s\big(1-\mathrm{Re}\,\hat\rho_{\tau(s)}(1)\big)$.
    Optimizing a three-atom $\eta\in\{-\alpha,0,\alpha\}$ (uniform, or typed by $b_s$-quantile) captures
    $\le40\%$ at $\omega=2$ and only $0$–$19\%$ at $\omega=3$, and the optimal law *degenerates to a point mass*
    (a deterministic shift). Reason (structural): the reward needs *coherent* alignment, and the product
    $\prod_s\hat\rho(C_{rs})$ of *independent* phases decays in the comb support — independence destroys the
    coherence. Only a *correlated* (block / cavity) ensemble could capture it, which is non-rigorous.

These two pin the obstacle precisely. The maximizing configuration is simultaneously **non-perturbative**
($|\theta-\pi/2|\sim0.8$, killing linear/quadratic/normal-form models — those *explode*), **coherent /
correlated** (killing product / mean-field certificates), and **disordered / generic-angle** (killing mod-$4$
and any clean residue pattern). A proof needs an object that is correlated *and* continuous *and*
non-structured at once — which is precisely why no off-the-shelf tool fits.

## 6. Literature placement (bottom line)

An exhaustive 2026 search found **no published theorem that upper-bounds this sup-side constrained maximum**.
The MPS–Konyagin / Bloom–Green circle is $L^1$/lower-bound; the Chowla resolution (Bedert; JMTZ) is min-side;
Conway–Jones / Lam–Leung classify the relation lattice but carry **no metric/alignment content**;
CVP-in-cyclotomic-lattices and transference give covering-radius (wrong-direction) bounds. The one analytic
tool that fits the weighted-cosecant *object* (Blagouchine–Moreau, `arXiv:2312.16657`) handles a free single
phase, not a max over the cyclotomic phase lattice. So A1 is, as far as the literature goes, **open and
untooled**.

## 7. The clean question to pose

> Fix $3\mid N$. Let $K$ be the support-$>3$ "comb" generators of the integer relation lattice of
> $\{2\sin(\pi r/N)\}_{r=1}^{\lfloor N/2\rfloor}$ (explicit: convolutions of period-$p$ telescoping kernels),
> $b_r=1/(N\sin(\pi r/N))$, prefix size $M=\tfrac12\varphi(2N)$, and $C$ the integer matrix giving the
> dependent phases $\phi=C\psi$. Is
> $$ \inf_{\lambda\in[0,1)}\ \max_\psi\Big[\sum_{r\in K}b_r\sin\big((C\psi)_r\big)
>    -\lambda\sum_{s\le M}b_s(1-\sin\psi_s)\Big] $$
> bounded **uniformly in $N$** (alternative (O1)), or does it grow like $c\cdot\frac1\pi\ln\frac m{\varphi(m)}$
> (alternative (L))? Equivalently, find the sharp $\kappa\in(0,1]$ in
> $U_N-A_N=\frac\kappa\pi\ln\frac m{\varphi(m)}+O(1)$.

A proof in **either** direction is a new theorem. A clean general method would be a genuine sup-side analogue
of McGehee–Pigno–Smith / Konyagin (an $L^\infty$ upper bound on a relation-constrained nonnegative-weight
sine sum).

## 8. Reformulation as ergodic optimization (an alternative framing, for dynamicists)

A1 is naturally an **ergodic-optimization / max-plus** problem — a framing that may reach a different
community (weak-KAM, Mañé–Mather, max-plus spectral theory) than harmonic analysts. The comb objective is a
reward–cost functional
$$ R(\psi)=\sum_{r\in K}b_r\sin\big((C\psi)_r\big)-\lambda\sum_{s\le M}b_s(1-\sin\psi_s), $$
maximised over the prefix phases $\psi$, with the relation matrix $C$ acting as deterministic dynamics that
ties the dependent (reward) phases to the free (cost) phases. The two alternatives of §2 are **exactly** the
ergodic-optimization dichotomy:

> **(O1) $\iff$ bounded coboundary.** A potential ("subaction" / Bellman function) $B$ exists with
> $R(\psi)\le B(T\psi)-B(\psi)+O(1)$ along the comb chain, so $F_{\mathrm{comb}}$ telescopes to $O(1)$.
>
> **(L) $\iff$ positive max-plus cycle.** A calibrated periodic phase pattern of positive mean reward
> $\Gamma>0$ exists and repeats along the chain, forcing $F_{\mathrm{comb}}\gtrsim\Gamma\cdot(\text{length})$.

This is the right *primal* language (an upper bound is a Bellman potential; a lower bound is a calibrated
cycle), and it explains every failed route in §5 at once: the mass-based methods ($L^2$ / SDP / CVP / $L^1$)
cannot see the *orbit/phase* structure that creates the maximum. The closed case **$N=3p$ is literally a
coboundary** — the period-$3$ telescoping $1+\zeta_3+\zeta_3^2=0$ gives $R=B_{\mathrm{next}}-B_{\mathrm{prev}}+o(1)$.

**Two caveats, both tested numerically, both failing as stated:**
1. *The single-prime-cell decomposition does not work.* A single odd prime gives **bounded, vanishing**
   excess: for $\omega(m)=1$ (e.g. $N=p^2$), $E_N\le\frac1\pi\ln\frac{p}{p-1}\sim\frac1{\pi p}\to0$ (the
   Kronecker sandwich is nearly tight; measured $0.045,0.034,0.030$ at $p=7,11,13$). A single prime is **not
   growth-capable**, so the excess growth is irreducibly **multi-prime** and a sum $\sum_{p\mid m}\Gamma_p$ of
   single-prime pressures cannot produce it — the same prime-decoupling error that sank the Möbius route (§5C).
2. *The sharp constant is not $1/\pi$.* The ratio $E_N\big/\big(\frac1\pi\ln\frac m{\varphi(m)}\big)$ is **not**
   stable at $1/\pi=0.318$: it ranges $0.33$–$0.88$, depends strongly on $N$ (not merely $\omega$), and sits
   systematically **above** $1/\pi$. At reachable $N$ the excess is $O(1)$-dominated and the asymptotic $\kappa$
   is not cleanly extractable.

**Net.** Ergodic optimization **correctly repackages** A1 (and may bring weak-KAM / max-plus machinery to
bear), but it does **not reduce** it: the genuine task is a Bellman potential / calibrated cycle on the **full
multi-prime comb** — the single-prime shortcut cannot replace it, and the value $\kappa=1-1/\pi$ is not
supported by the data. The non-perturbative optimiser admits a *positive-cycle* (L) reading, but the
diagnostic of §9 shows the optimum is disordered, so (L) — if true — would have to be a *spin-glass-type*
pressure, not a clean cycle, while (O1) via a global energy–entropy bound is equally open; the dichotomy is
genuinely unresolved.

## 9. Diagnostic verdict: the extremum is disordered (why no clean certificate exists)

Rather than guess an ansatz from above, we extracted structure *from* the true optimizer $\psi^*$ via its
stationarity (KKT) equation $b_s\cos\psi^*_s=-h_s$, $h_s=\sum_{r>M}b_r C_{rs}\cos\phi_r$ (the local field). Two
diagnostics on $\psi^*$ for $N=105,315,1155$ (optimizer reaching the known $E_N=0.094,0.094,0.116$):

- **Compressibility.** Out-of-sample (train-half / predict-half) regression of the displacement
  $\psi^*_s-\frac\pi2$ on arithmetic features ($s/M$, $\gcd(s,N)$, $s\bmod 3,5,7$, $b_s$, comb-degree), then
  reconstructing $\hat\psi$ and measuring the captured excess $G(\hat\psi)/E_N$: a random forest gives
  $-0.17,\,0.23,\,-0.12$ ($\approx0$ or *negative*), and linear regression $-2.3,\,-1.4,\,-0.04$. The
  arithmetic features show **no low-complexity structure under these diagnostics** — predicting the phases
  from them is no better than noise.
- **Local correlations.** $\mathrm{Corr}(\psi^*_s,\psi^*_t)$ for pairs $s,t$ sharing a comb is
  $0.03,\,0.06,\,-0.02$ — *indistinguishable* from random pairs ($-0.03,-0.01,0.04$). No local correlation
  structure; the configuration looks fully frustrated. (The displacements are non-perturbative: mean
  $|\cos\psi^*_s|\approx0.35$, individual modes ranging up to $\sim0.8$ from $\pi/2$.)

So under these diagnostics the maximizing configuration **appears disordered and frustrated** — empirically
incompressible, with no detected local correlations — and **none of the tested low-complexity mechanisms**
(perturbative, residue-class, product / mean-field, typed mean-field, or local-block) captures it. *This is a
strong numerical/structural diagnosis, not a theorem that no structure exists*: a hidden cyclotomic structure
invisible to these tests is not excluded. **A universality test settles part of this:** replacing the true
$C$ by randomized relation matrices (shuffled support / random signs / random sparse, matched in degree
profile and weights) and recomputing the excess, at the optimizable sizes $N=105,315$ the *randomized* $C$
yields a $1.4$–$2\times$ **larger** excess — so the cyclotomic arithmetic is **not** irrelevant; it
*suppresses* the excess below the generic frustrated level (a collective structure invisible to the per-mode
tests, and a mild lean toward \textup{(O1)}). At $N=1155$ ($M=240$) the random landscape is too rough for the
optimizer to solve (spurious negative values even with $150$ restarts), so no clean comparison there — but
that roughness, against the readily-optimized cyclotomic optimum, is itself evidence the cyclotomic case is
special, not generic spin-glass. So A1 is best regarded not as a standard harmonic-analysis problem but as a
**frustrated variational problem with special (suppressing) cyclotomic structure** over the relation lattice,
with the dichotomy genuinely open and localized to: prove either a *bounded energy–entropy upper bound*
(giving \textup{(O1)}, now mildly favoured by the suppression) or a *positive spin-glass-type pressure*
(giving \textup{(L)}, no longer expected from any clean arithmetic cycle). This is the sharpest available
statement of *why* A1 is hard.

**Multiplicity / overlap probe — tilts the dichotomy toward \textup{(O1)}.** Collecting near-optimal local
optima (multi-start) and measuring overlaps $q_{ab}=\frac1M\sum_s\cos(\psi^{(a)}_s-\psi^{(b)}_s)$: at $N=105,315$
they concentrate at $q\approx1$ ($N=105$: $13$ near-optima, all $q=1.000$; $N=315$: $q\in[0.84,0.97]$), so the
maximiser is **essentially unique** — *not* a spin-glass with many competing states (no broad overlap
distribution / replica-symmetry-breaking; at $N=105$ even randomized $C$ has a unique optimum, only larger).
This **refutes the \textup{(L)}-via-spin-glass-pressure mechanism**: with no competing states there is no source
for an extensive positive pressure. Combined with the suppression, the structural evidence now **leans
\textup{(O1)}** (a unique, suppressed optimum); the residual \textup{(L)} signal is only the
$O(1)$-contaminated finite-$N$ ratio. The natural forward route is thus an **energy–entropy / inverse-theorem
upper bound** (a unique optimum means low entropy, so the entropy term of a free-energy bound is small) — but
the free energy itself resists naive Monte-Carlo (the partition function at the relevant inverse temperature is
dominated by *rare* high-$F$ samples — the same rarity wall as the direct $\sup_t$), so it needs the analytic
cluster expansion or a genuine inverse theorem. The dichotomy remains formally open, now leaning \textup{(O1)}.

**Curvature / unique-peak Laplace probe — the energy–entropy route hits an $M$-scaling wall.** Given the unique
optimum, the natural \textup{(O1)} attack is a Laplace bound around the single peak,
$E_N\le\frac1\beta\log Z(\beta)+\frac1{2\beta}\log\det A+\frac{M}{2\beta}\log\beta+O(1/\beta)$, where
$A=-\nabla^2F(\psi^*)=C^\top\mathrm{diag}(b_r\sin\phi^*_r)C\succeq0$. Computing the Hessian spectrum at the optimum
(true vs randomized $C$; reliable at $N=105,315$ — at $N=1155$, $M=240$, the optimizer cannot resolve the peak and
the spectrum is an artifact): the peak is genuinely **non-degenerate** (no soft modes), confirming the
isolated-unique-peak premise — but $\frac1M\log\det A$ is a **stable per-mode constant** ($\approx-2.9,-4.0$
cyclotomic vs $\approx-2.1,-3.0$ random), so $\log\det A\propto M$, and **both** the $\frac1{2\beta}\log\det A$ and
$\frac{M}{2\beta}\log\beta$ terms scale with $M$. They are $O(1)$ only at $\beta\sim M$, where
$\frac1\beta\log Z(\beta)\to E_N$ and the bound collapses into the very quantity it was meant to bound. The unique,
low-entropy peak does **not** remove this — the Gaussian fluctuation determinant itself scales with $M$. The same
probe refutes a tempting suppression mechanism: the cyclotomic curvature is *softer*, not stiffer, than random
($\log\det A/M$ more negative; typical eigenvalue $0.055$ vs $0.124$ at $N=105$), so the suppression is a *low,
flat* peak, not a rigid one. **Consequence:** the free-energy / Laplace / naive-cluster family of \textup{(O1)}
attacks all hit the same wall — to win, one must show $\log Z(\beta)$ lies *below* the naive estimate by more than
the $M$-scale, which is exactly the missing sup-side (Gowers/inverse-type) input. This pins down precisely what a
proof of \textup{(O1)} must supply.

---

*All claims are backed by code in the repository (`numerics/verify_excess_*.py`, `verify_finite_time_*.py`,
the `exact/` PARI + `formal/rocq/` layers) and by the paper's Remark "the defect" + the per-relation and
clean-comb lemmas + the "excess dichotomy" conjecture in `paper/sections/linear.tex`.*
