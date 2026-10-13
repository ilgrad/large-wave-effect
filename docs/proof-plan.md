# Proof plan

Per-task plan: what to prove, the method, the tools, and the **current status in this repo**
(proved / tried-and-failed-this-session / open + what is actually needed). Cross-references the
external review against the session's own findings.

**Unifying finding (updated 2026-06-15).** These two were thought to be the **same** obstruction — a
**sup-side signed-cancellation estimate** (primal McGehee–Pigno–Smith / Konyagin). That holds for the
composite excess lemma (1) = **A1**, whose signed core is reduced to a concrete **weighted-packing bound
on the chained relation graph**, uniform over odd `m` (A1 Mechanism); seven routes are ruled out
(modulus/triangle, multivariate concavity, variance-floor, monotone-tail, Berry–Esseen, Bedert's Chowla
machinery — arXiv:2509.05260, naive Arb ball). But the β_odd continuum cover (5) = **A2** turned out
NOT to be that wall: direct computation showed it is **enclosure-limited** (true pointwise margin `~0.12`,
never near-tight), the binding loss being the capped-`|M|` tail under a θ-BALL middle that blows up past
`t=12`. It is now **resolved** by a θ-Taylor + t-collocation enclosure (`middle_far_point`), reducing A2
to an executable tiling (see A2). So only **A1** (and the self-trapping window, §C) remains the signed wall.

---

## Snapshot (closed vs open)

**Closed this session (rigorous).**
- Dirichlet segment mod-4 ceiling criterion (`cor:segceiling`) — closed the "rank-deficient ⇒ A<U" gap.
- Disorder: generic-independence via eigenvalue-map submersion (`F_k≠0`); finite-T continuity
  (`lem:disorder-finiteT`, gain needs `T≳1/eps`).
- Exact-`A_N` table: `A_N<U_N` certified with explicit gap (Lipschitz grid `N=6,9,12,15`; exact rational
  SOS `N=9,15,21`, within `1e-7` of the optimizer, `exact/pari/defect_sos.gp`).
- Deficit lower bounds for every `N` (`prop:deficit-lower`, 2026-10-13): `U_N-A_N >= 1/(16p^2)` from one
  relation per residue class (`lem:rel-defect`), `> 0.043` for `3|N` (triplets + `lem:triplet`), and
  `(pi^3/16) N^-2 (1+O(1/N))` for `N=2p`. Constant along primorials, so the growth question is untouched.
- Breather **dichotomy** CAP (exact rational): site `n(L_+)=1` STABLE, bond `n(L_+)=2` UNSTABLE.
- Finite-time **upper** bound `T_{0.2}(N)≤t~` certified (interval) for concrete primes `N=5,7,11,13`.
- (Earlier in session: classification, `A_N~(1/π)ln N`, mod-4, β_odd profile/transfer/two-copy — see A2.)

**Resolved this session (2026-06-15).**
- **A2** β_odd continuum cover — the θ-ball middle blow-up past `t=12` is fixed by `middle_far_point`
  (θ-Taylor + t-collocation; `verify_middle_far.py`); the previously-failing boxes now certify
  (e.g. `(0.70,0.70)→0.9915`; hardest boxes at the high-λ edge `a_3≈0.61`, closed by one refinement, worst
  leaf `0.9999`). Full band tiling **DONE**: `run_dagger_campaign.py` certified 136 leaf boxes, 0 fail (~7.2 h; 6 boxes
  added 2026-10-13 after `--coverage` found a gap at `a_1≈1.11–1.13`)
  ⇒ A2 reduced to ONE numerically-validated ≥3-copy inequality. **Corrected 2026-10-12/13:** the target is
  the *windowed* constant (see A2), not the all-time `limsup B_N/√N`, which is `1` (`B_N = √N` for `N = p, p²`);
  and the cover is validated, not rigorous modulo one estimate — see the further non-rigorous steps in A2.

**Open — research wall (needs new mathematics).**
- **A1** excess lemma, uniform in `m` = weighted-packing bound on the chained graph — sup-side
  signed-cancellation estimate (7 routes ruled out).
- Self-trapping window `gP∈[0.43,4)` — same focusing-sign motif (§C).

**Open — execution / resource (idea exists, not a one-shot).**
- β_odd full band tiling — **DONE, RIGOROUS since 2026-10-13** (`dagger_rigorous.py` via `run_dagger_campaign.py --rigorous`: 111 boxes, 0 fail, worst `0.9205 K*`, 33 s on 8 cores; `verify_dagger_rigorous.py` re-decides, recomputes, checks coverage exactly, writes the manifest; the old 136-box validated tiling is superseded — see A2);
  finite-time **lower** bound (+ damping, conditional on it);
  FPUT small-energy (Nekhoroshev, 4-wave resonances) and chaos (CAPD); arbitrary graphs (non-cyclotomic).

**Blocks publication, not the mathematics.**
- Filimonov 1992 (C.R. Acad. Sci.) — **RESOLVED**: note obtained & compared. His Thms 1–4 = ceiling under
  Q-independence + the prime/2^m-vs-even-composite dichotomy (cyclotomic) + unbounded ceiling (no rate);
  ours = rate (1/π)ln N, constant, mod-4 saturation classification, all-N rank ½φ(2N).

(Detail + methods per task in A–E below; this snapshot is the index.)

---

## A. Central problems (deeply attacked this session)

### A1. Excess lemma / composite defect `U_N - A_N`  (review #1)

**Goal.** For `N = 2^a m`, `m` odd: `U_N - A_N = (1/pi) ln(m/phi(m)) + O(1)`; unbounded,
`Theta(ln ln ln N)` along primorials. Reduces to the excess lemma `A_N <= L_pre + O(1)`.

**Method (review + session).** Rewrite as a max of a trig polynomial on the relation subtorus:
`A_N = max_x [ sum_{r<=M} b_r sin x_r + sum_{r>M} b_r sin(l_r . x) ]`, `M = phi(2N)/2`,
`b_r = 1/(N sin(pi r/N))`. Seek a **dual / majorant certificate**
`sum_{r>M} b_r sin(l_r . x) <= C + sum_{r<=M} b_r (1 - sin x_r)`. Per-N via trig-SOS / Lasserre; the
FAMILY bound must be uniform over odd `m` (Scope correction below) — a **weighted-packing estimate on
the chained relation graph**, NOT the `a -> infinity` limit (that direction is trivial).

**Status.**
- PROVED: `omega(m) <= 1` (`N=2p` sharp `1/(2N)`, `N=4p`); rank `phi(2N)/2`; mod-4 criterion;
  classification `A_N=U_N <=> N` prime/`2^m`; order `A_N ~ (1/pi) ln N`.
- PROVED 2026-06-20 (rank/saturation circle, closed forms): cosine rank `rank{2cos(pi k/N)} = phi(2N)/2 - [N=2^m]`
  (both directions; full <=> prime/2^m/2p) and sine rank `rank{2sin(pi r/N)} = phi(2N)/2` (no exception; full <=>
  prime/2^m = classification) -- the cos `-1`-on-`2^m` vs sin no-exception is why the segment has the extra `2p`
  family; Dirichlet-segment SITE-WISE saturation criterion (`prop:segment-site`, Kronecker+Pontryagin parity);
  finite-time amplified-set density `rho(eps) ~ V_m(2 eps U)^{m/2}/sqrt(prod b)/(2pi)^m` (Weyl; measure-side of
  the open worst-case recurrence); `g=1` NNN collision law via Conway-Jones 5-term sums (was verified `N<=120`).
  Commits 678b896,9e49f37,43ea733,1af324f,6982178,d38096d,9cb53df. Walls re-confirmed untooled (2026 lit):
  A1 sup-side MPS (Bedert min-side), cubic DNLS scattering (only modified), Baker for `2sin` vector (only
  qualitative; lead Calegari-Dimitrov-Tang arXiv:2510.04156).
- DONE this session (`omega(m)=2`): relation lattice characterized — unimodular SNF, **chained**
  relations, support set by the smallest prime (`3|N`: `w_a+w_b=w_c`, `a+b=N/3`, from `2 sin(pi/6)=1`;
  `5|N`: golden); per-N certificates `N=15` (elementary Lipschitz grid, `>=0.081`; exact rational SOS,
  `>=0.14754945`), `N=21` (exact rational SOS, `>=0.11208183`); `N=33,35` floating-point moment-SOS
  only. `exact/gap/excess_omega2.g`, `numerics/verify_excess_omega2.py`.
- TRIED → blocked: the **uniform** dual certificate is exactly the signed bound that overshoots.
  The obstruction is a **modest O(1) balance** — verified this session across `omega(m)=2..4`: as
  `dep_wt` grows `0.21 -> 0.34`, BOTH `dependent_net ~ +0.16` and `prefix_deficit ~ -0.05..-0.09` stay
  individually O(1) and flat (NOT a near-cancellation of growing terms; earlier docstring corrected).
  Bedert does NOT transfer (4 mismatches; `verify_excess_omega2_spectral.py`).
- **MECHANISM (this session, numerically pinned).** The dual majorant MUST overshoot: near `x=pi/2` a
  dependent gain is LINEAR (`sin(phi_a+phi_b) ~ u+v`) but the prefix penalty is QUADRATIC (`1-sin ~ u^2/2`),
  so no `C=0` certificate exists. A single support-3 triplet `w_a+w_b=w_c` balances at `x_a=x_b=pi/3` with
  contribution EXACTLY `3*sqrt3/2 - 2 = 0.598*b_c` (verified). For DISJOINT triplets the excess is additive
  (`N=15`: support-3 sum `0.0556` ~ full excess `0.0666`, comb mode `+0.011`). Boundedness of the signed
  net comes from the CHAINED graph forcing dependent modes to SHARE prefix phases. ==> the excess lemma
  REDUCES to a **weighted-packing bound on the chained relation graph**.
- **SCOPE correction.** For FIXED `m` the sandwich `U-L_pre = (1/pi)ln(m/phi(m)) + O(1)` is a CONSTANT, so
  excess (hence defect) converges trivially — the family `2^a pq` (fixed `p,q`) does NOT exhibit the
  unbounded defect (the earlier "NEXT" was off). Unboundedness is along PRIMORIALS (`m` grows,
  `omega(m)->inf`); the excess lemma needs `C` uniform **in `m`**, not in `a`. Realistic target: the
  weighted-packing bound, uniform over odd `m`. (`a -> infinity` is the trivial direction.)

### A2. Windowed Schrödinger constant `limsup_N sup_{t<=SN/2} ||K_N(.,t)||_1/sqrt(N) = beta_odd`  (review #5)

**Correction (2026-10-12/13).** As first stated, A2 asked for `limsup B_N/sqrt(N) = beta_odd` with
`B_N = sup_t ||e^{-itL_N}||_{inf->inf}` over ALL `t`. That is false. `B_N = max` over the orbit closure `H_N`
(dim `phi(N)/2`), and `B_N = sqrt(N)` iff `H_N` holds a flat phase vector (epsilon-uniform mixing on `C_N`):
true for `N = p` (Gauss chirp; Godsil–Mullin–Roy) and `N = p^2` (even Frank point), false for even `N >= 6`
(symmetric periodic Golay-pair obstruction) — `thm:bn-saturation`, `verify_bn_saturation.py`,
`exact/pari/bn_saturation.gp`; a second implementation of the prime and prime-square cases in
`verify_bn_prime_flat.py`, `exact/pari/bn_prime_flat.gp`, `exact/pari/bn_prime_square.gp`. So
`limsup B_N/sqrt(N) = 1`. Also false for `N = 15`: exact Nullstellensatz certificate `D y^4 = sum N_j m_j f_j`
(`exact/julia/c15_flat_certificate.jl`, re-checked in PARI); Cao–Zhu's Thm 3.2 Gröbner basis implies it too.
Also false for `N = 21` (2026-10-13): `Z_21 = Z_3 x Z_7` coordinates make `H_21` the 6-torus
`x(1,a) = y_a, x(2,a) = y_{-a}, x(0,a) = 1/(y_a y_{-a})`, `prod y = 1`; an invariant Laurent Nullstellensatz
certificate in `Q[Z^7/Z(1..1)]` at median-norm radius 8 (6707 orbit equations, 5896 unknowns, rank 5484 mod p;
4005 coefficients, 614-digit denominator; FLINT solve ~3 min) — `numerics/bn_flat_certificate_3q.py`, re-checked
from the definitions in `exact/pari/c3q_flat_check.gp`; `q = 5` (radius 6) reproves `N = 15`. GPY
(arXiv:2607.28345) named `C_21` the next open case for exact uniform mixing. `N = 33` (`q = 11`) is out of reach
this way so far: msolve F4 mod p reached degree 7 at 8.7 GB and was stopped before exhausting RAM, and the
invariant system at radius 8 would be ~63k x 25k (the radius actually needed is unknown).
Open (paper `conj:bn-psq`): the other odd composite non-prime-square `N`, from `N = 27` on (numerics `< 1`,
closest `N = 27`, `0.99992`). `thm:joint-equi` and `thm:uniform-transfer` hold for each FIXED `s = 2t/N`. The statement the
analysis supports is the windowed one (paper `conj:bn-window`): for each fixed `S`, it needs (i) the ≥3-copy
dagger inequality and (ii) the transfer uniform in `s` on `[1, S]` (and on `t <= N/2` as well). Polynomial
windows `t <= N^C`: open.

**Goal.** `beta_odd = 0.92801930480793...` (odd subsequence, windowed); even subsequence `c0/sqrt(2)`
(windowed, numerical).

**Method (= the review's plan, fully executed).** `K_N(k,t)` via periodized Bessel images;
`t=sN/2`; uniform Debye for `J_n(nz)`, Airy at the front, stationary phase for tails; interval
arithmetic on a compact `s`-range; Gaussian majorant for the multi-copy overlap.

**Status.**
- PROVED: `t<=N/2` profile maximizer `max F = F(1) = beta_odd` (analytic + validated interval `F'>0`).
- PROVED (this session): joint Debye equidistribution (`thm:joint-equi`, van der Corput); the
  uniform-in-u transfer to the `l^1` limit (`thm:uniform-transfer`) — required fixing a real gap (the
  modulus has DIVERGENT Hardy–Krause variation; mollified Koksma–Hlawka + cubic van der Corput).
- PROVED: the two-copy core `R(1,rho) <= K_star` (interval, gap-convexity at the tangency).
- COVERED POINTWISE (`>=3` copies): near-pair basin `lambda>=0.82` (Phi convex in amplitude coords);
  singleton-channel band `lambda in [0.65,0.82)` (three-region Arb bound, worst `0.99 K_star`); bulk
  (`>=2`/channel) crude. All `< K_star`. `verify_dagger_{extremal2,decomp,singleton,concavity,tail}.py`.
- CONTINUUM COVER (2026-06-15, mechanism resolved). The equal-amplitude boxes near `lambda≈0.73` that
  failed (`(0.73,0.73,h=0.015)→1.068`, `(0.76,0.76)→1.031`) failed for a now-pinned reason: the binding
  loss is the capped-`|M|` tail (`~0.047 K*`, almost pure — the true SIGNED tail is `≈-0.003`), and the
  signed acb-series `middle_point` could NOT be extended past `T2=12` to recover it because its rigorous
  THETA-BALL radius blows up `~amp*t/nth` (precision-independent; `≈300` at `t=20`; the true pointwise
  margin is a healthy `~0.12` — it was enclosure-limited, never near-tight).
  FIXED by `middle_far_point` (`verify_dagger_continuum.py`): a THETA-TAYLOR + t-COLLOCATION enclosure of
  the signed middle for `t>=8` — theta is a Taylor series variable (no ball inflation), t is sampled at
  POINTS (the 0F1 J0-series only blows the radius under an INTERVAL t), the t-integral is a degree-2
  collocation quadrature with a finite-difference Lagrange remainder. Validated vs the scipy reference
  (containment, radius `~7e-7`; `verify_middle_far.py`). Recipe `mid[1.5,8] + far[8,24] + tail[24,inf)`
  via `prove_box(..., far_slices=...)` certifies the previously-failing boxes (e.g.
  `(0.70,0.70)→0.9915`; the hardest boxes are the high-λ edge `a_3≈0.61`, closed by one refinement, worst
  leaf `0.9999`). Node-parallel (`prove_box_parallel.py`); the full band tiling is the adaptive, checkpointed
  `run_dagger_campaign.py`, which certified all 136 leaf boxes (0 fail, ~7.2 h; `--coverage` checks the band is
  covered) — the cover **validates** the residual ≥3-copy inequality across the band (NOT a proof). It is not
  "rigorous modulo one estimate": besides (a) the far `t`-collocation remainder from finite differences
  (`verify_dagger_continuum.py` `middle_far_point`), further steps are non-rigorous (2026-10-12 audit):
  (b) `build_stencil` takes the amplitude third derivatives from finite differences, not a derivative bound,
  and ignores the error propagated through the least-squares fit; (c) `StencilModel.sup` maximizes the
  quadratic on a 41×41 float grid instead of exactly over the box; (d) `run_dagger_campaign.py` reads the
  per-box `lam_lo/lam_hi` and `bmin` from 3×3 samples — since 2026-10-13 band coverage is checked separately
  with interval bounds (`--coverage`), and `bmin` is attained at the sampled corner because
  `b = sqrt(2 − a1² − a3²)` decreases in both amplitudes (up to float rounding); (f) no manifest — library
  versions, precision, hashes.
- RIGOROUS COVER (2026-10-13; `dagger_rigorous.py`, `verify_dagger_rigorous.py`): (a)–(f) are all removed. In
  `u = a1², v = a3²` the signed Kluyver integral on `[0, 24]` is a Gauss–Legendre (t) × periodic midpoint (θ)
  sum with analytic error bounds (Trefethen ATAP Thm 19.3 on Bernstein ellipses; strip aliasing bound), carried
  as an EXACT degree-10 Taylor polynomial in `(u, v)` (Arb power series through every `J0 = 0F1`) with a
  bivariate Cauchy remainder; `sup_box (P − K_lo D)` is bounded exactly (quadratic part over the rectangle,
  higher terms by absolute values, 4^d recentring); the tail uses Watson 13.74's `|J0(x)| ≤ sqrt(2/(πx))` in
  an upper Riemann sum; `K*` is enclosed by `acb.integral` (cross-checked in Maxima: `0.03593223661504624`).
  111 boxes (h = 0.015, a1 ≥ a3) prove, worst `0.92050` at `(0.85, 0.61)`; the checker re-decides every row
  from its stored balls, recomputes the 3 worst bit for bit, checks band coverage in exact rationals, and pins
  source/data SHA-256, parameters and versions in `numerics/dagger_rigorous_manifest.json`. Assumed, not
  re-proved: (K), Watson 13.74, ATAP 19.3, the a1 ↔ a3 symmetry, Arb / python-flint 0.9.0.
- STILL VALIDATED (not touched by the rigorous cover): the near-pair basin `λ ≥ 0.82` (grid-Lipschitz Hessian
  floor per copy-count chart, `verify_dagger_decomp.py`: no interval J0, no bound uniform in the copy count)
  and the configurations with more than three copies (crude modulus bound where both channels carry
  comparable copies; elsewhere a high-precision search).
- FIXED: the `beta_odd` literal (was wrong past digit 9).

**A1 remains the signed-cancellation wall; A2's singleton-band (2-vs-1) cover is now RIGOROUS (111 boxes, 0 fail,
exact coverage; above), reducing the windowed A2 to the `s`-uniform transfer plus the near-pair basin and the
>3-copy configurations, which stay validated.**

---

## B. Closable now (review "closable now")

| Task | Method | Status |
|---|---|---|
| DNLS weak persistence `‖z-z_lin‖_inf <= |g|T P^{3/2}` | Duhamel | ALREADY proved (`prop:dnls-persist`) |
| `N=2p, 4p` excess `O(1/N)` | sharp two-sided | ALREADY proved |
| No dispersal for `gP>4` (energy-sharp, not dynamically sharp), on-site trapping for `gP>8`; quantitative dispersal for `gP<0.43` | energy / Duhamel + Landau | ALREADY proved (`prop:threshold`, `prop:dispersal`) |
| Dirichlet segment mod-4 reachability | port the ring criterion | **DONE** (`cor:segceiling`, `exact/gap/segment_mod4.g`; criterion <=> `N+1` prime/`2^m`, `N<=40`) |
| Exact `A_N` table, small N | orbit-closure subtorus max + certified upper | **DONE** — Lipschitz-grid certifies `A_N<U_N` w/ explicit gap `N=6,9,12,15` (`verify_defect_certified.py`); exact rational SOS pins `N=9,15,21` to `1e-7` (`exact/pari/defect_sos.gp`); floating-point moment-SOS `N=33,35`; orbit-closure values `verify_exact_AN.py`. |
| Generic disorder independence | eigenvalue map submersion; `F_k != 0` | **DONE** (Jacobian `d lambda_r/dm_j=-lambda_r v_r(j)^2` full row rank `N-1`, `verify_disorder.py` (A'); linear.tex Disorder) |
| Fixed-T continuity under disorder | matrix-function Lipschitz `<= C_N(1+T)eps` | **DONE** (`lem:disorder-finiteT`; gain needs `T >~ 1/eps`; `verify_disorder.py` (A'')) |
| **Damping criterion** (`T_eps >> 1/eta` => ceiling unreachable) | modal `e^{-eta t}` | OPEN — **conditional on `T_eps`** (the open finite-time lower bound); uniform envelope `e^{-gamma t}` kills the build-up only if `T_eps(N) >~ 1/gamma`. Not unconditionally closable. |

---

## C. Harder / open

- **Finite-time `A_N(T)`, `T_eps(N)`** (review #3): **PARTLY DONE this session** — for concrete primes
  `N=5,7,11,13` the **upper** bound `T_{0.2}(N) <= t~` is **certified in interval arithmetic** (a witness
  `t~` where `u(t~) >= 0.8 U_N`, mpmath.iv; `verify_finite_time_certified.py`), upgrading the numerical
  recurrence times to a proof. **Quantified this session** (`verify_finite_time_lower_bound.py`, GPU/CuPy,
  primes `N<=71`): `T_eps(N)` fits the simultaneous-approximation law `T_eps ~ eps^(-a m)`, `a~0.32`,
  `m=floor(N/2)` (R^2=0.88), i.e. `T_0.2 ~ exp(0.245 N)` over `N in [11,47]`; the fixed-horizon efficiency
  `A_N(T)/U_N` decays monotonically in N (at `T=1e3`: 1.0 at N=5 -> 0.58 at N=71). [Corrected the paper's
  `~0.63 per unit N`: that was the pre-asymptotic slope over `N<=19`, inflated by the O(1) times at N=5,7.]
  **STILL OPEN (rigorous):** the matching family LOWER bound `T_eps >~ eps^(-c m)` does NOT follow from the
  homogeneous Dirichlet pigeonhole (the target `pi/2` is inhomogeneous); it needs a Baker-type Diophantine
  lower bound on linear forms in the cyclotomic frequencies `omega_r=2 sin(pi r/N)` — a genuine
  number-theoretic wall, like A1. (Upper Dirichlet `~delta^{-m}` bound + per-N certified witnesses: done.)
- **Self-trapping window `gP in [0.43,4)`** (review #8): reduced to a bound on the signed
  staggered-momentum drift (`eq:stagger`); bounding it by its modulus fails, since the `L^4` spacetime norm
  diverges like `(ln T)^2` even on the free flow (`verify_stagger_morawetz.py`). Numerically the seed
  disperses only up to `gP≈3.4` and keeps a remnant from `≈3.5` (`verify_selftrapping_transition.py`; JHR95
  `≈3.45`), so the target is that crossover, not `4`. **Leading order now PROVED**
  (`lem:stagger-leading`): along the linear flow from the single-site seed `a_n=e^{2it}(-i)^n J_n(2t)`,
  the distance-2 current `Im(conj(a_n) a_{n+2})=0` identically (phase `(-i)^2=-1` is real => `a_{n+2}` is a
  real multiple of `a_n`), so the O(g) drift cancels: `sup_t|ΔS|~(gP)^2` (exp 2.08, N-uniform),
  `inf<v^2>≈0.49>0`. Verified Maxima (ODE residual 0, phase -1) + Rocq (`formal/rocq/StaggeredCurrent.v`,
  real-multiple=>zero-current) + numerics (`verify_staggered_current.py`, `verify_modified_energy_window.py`,
  GPU). OPEN: N-uniform bound on the `O(g^2)` remainder + nonperturbative approach to the crossover `gP≈3.5` (focusing sign
  survives there — still a signed-cancellation motif).
- **Breather stability** (review #4): VK/GSS stated; (i) `dP/dOmega>0` for large `Omega` (IFT in `eps=1/Omega`), and enclosed at `Omega=2`
  (site `N=5`: `[0.9184911, 0.9184913]`; bond `N=6`: `[2.0204989, 2.0204992]`). **DONE this session** —
  the full site/bond **dichotomy** is **CAP-certified in exact rational arithmetic**: Newton--Kantorovich
  existence of the true profile (`rho<=4e-9`) + Sylvester-inertia `LDL^T` + Weyl transport
  (`verify_dnls_breather_cap.py`). Site `N=5,Omega=2`: `n(L_+)=1`, `lambda_2(L_-)>0` => GSS index 0, STABLE.
  Bond `N=6,Omega=2`: `n(L_+)=2` => GSS odd index 1, UNSTABLE. Both verdicts upgraded from numerical to
  computer-assisted. REMAINING: the full `Omega`-branch (continuation) stays numerical
  (`verify_dnls_breather_stability.py`); a family CAP would re-run the certificate on an `Omega`-grid.
- **FPUT/chaos** (review #9): small energy via Birkhoff/Nekhoroshev — **reassessed this session as NOT a
  quick win**: the ring has no THREE-wave resonances (`sin a + sin b = sin c` forces `b=0`) but DOES have
  four-wave ones (the FPUT paradox itself), so a rigorous Nekhoroshev bound is real work, not "closable".
  Chaos via CAPD computer-assisted (separate project).
- **Multidimensional / graphs** (review #10): `d_s=1` critical (done for products, Möbius);
  arbitrary graphs lose the cyclotomic structure — confirmed this session that even circulant
  `Cay(Z_N, S)` with `|S|>1` already fails (`omega_j = sqrt(sum_{s in S} 2-2cos(2*pi*j*s/N))` is not a
  cyclotomic integer, so the relation-lattice machinery does not port; only `|S|=1`, the ring, keeps it).

---

## D. Priority / publication

- **Filimonov 1992** (C.R. Acad. Sci.): **RESOLVED** — note obtained (filimonov/Scan1–5.pdf) and compared.
  Priority: the ring Q-independence dichotomy (Thm 2 prime/2^m independent via cyclotomic argument; Thm 3
  even-composite dependent via φ(2N)=N−2^m), the Bohr ceiling under independence (Thm 1), and the unbounded
  ceiling b_j→∞ (Thm 4, no rate) are Filimonov's. New here: the rate (1/π)ln N, the constant −π/72N², the
  mod-4 saturation classification (composite *lowers* the ceiling), the all-N rank ½φ(2N), B_N=Θ(√N), etc.
  Paper "Relation to prior work" + Remark rem:f92 + bib updated accordingly.
- PDF/TeX: kept in sync (`l^1 -> l^inf`, rebuilt every commit).

---

## E. Recommended order

1. ~~**B-block** (segment mod-4 + disorder continuity/independence + exact-`A_N` table)~~ — **DONE**.
   All rigorously closed this session except the damping criterion (conditional on the open
   finite-time lower bound). This is the self-contained "finite-time / disorder" content.
2. **Excess lemma, uniform in `m`** — a **weighted-packing bound on the chained relation graph** (the
   signed-cancellation core; mechanism pinned this session, see A1). A real new theorem. (NB the
   `2^a pq` fixed-`p,q` family is the trivial direction — it converges for free.)
3. **beta_odd offline tiling** — execution, off the interactive session.
4. **Sup-side MPS** (closes A1 and A2 at once) — the research breakthrough.
