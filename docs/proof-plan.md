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
- Exact-`A_N` table: `A_N<U_N` certified with explicit gap (Lipschitz-grid + moment-SOS).
- Breather **dichotomy** CAP (exact rational): site `n(L_+)=1` STABLE, bond `n(L_+)=2` UNSTABLE.
- Finite-time **upper** bound `T_{0.2}(N)≤t~` certified (interval) for concrete primes `N=5,7,11,13`.
- (Earlier in session: classification, `A_N~(1/π)ln N`, mod-4, β_odd profile/transfer/two-copy — see A2.)

**Resolved this session (2026-06-15).**
- **A2** β_odd continuum cover — the θ-ball middle blow-up past `t=12` is fixed by `middle_far_point`
  (θ-Taylor + t-collocation; `verify_middle_far.py`); the previously-failing boxes now certify
  (`(0.70,0.70)→0.948`, `(0.73,0.73)→0.943`, `(0.76,0.76)→0.922`). Full band tiling **COMPLETE**:
  `run_dagger_campaign.py` certified 130 leaf boxes, 0 fail (~7.2 h) ⇒ A2 closed (`limsup B_N/√N = β_odd`
  at the dagger certificate's validated-enclosure standard).

**Open — research wall (needs new mathematics).**
- **A1** excess lemma, uniform in `m` = weighted-packing bound on the chained graph — sup-side
  signed-cancellation estimate (7 routes ruled out).
- Self-trapping window `gP∈[0.43,4)` — same focusing-sign motif (§C).

**Open — execution / resource (idea exists, not a one-shot).**
- β_odd full band tiling — **DONE** (`run_dagger_campaign.py`: 130 leaf boxes, 0 fail, ~7.2 h; A2 closed);
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
  qualitative; lead Calegari-Huang arXiv:2510.04156).
- DONE this session (`omega(m)=2`): relation lattice characterized — unimodular SNF, **chained**
  relations, support set by the smallest prime (`3|N`: `w_a+w_b=w_c`, `a+b=N/3`, from `2 sin(pi/6)=1`;
  `5|N`: golden); per-N certificates `N=15` (elementary Lipschitz grid, `>=0.082`), `N=21,33,35`
  (moment-SOS). `exact/gap/excess_omega2.g`, `numerics/verify_excess_omega2.py`.
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

### A2. Schrödinger constant `limsup B_N/sqrt(N) = beta_odd`  (review #5)

**Goal.** `beta_odd = 0.92801930480793...` (odd subsequence); even subsequence `c0/sqrt(2)`.

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
  via `prove_box(..., far_slices=...)` certifies the previously-failing boxes:
  `(0.70,0.70,h=0.005)→0.948`, `(0.73,0.73)→0.943`, `(0.76,0.76)→0.922` (the worst true value is at
  `a≈0.70`). Node-parallel (`prove_box_parallel.py`); the full band tiling is the adaptive, checkpointed
  `run_dagger_campaign.py`, which certified all 130 leaf boxes (0 fail, ~7.2 h) — the cover is complete.
- FIXED: the `beta_odd` literal (was wrong past digit 9).

**A1 remains the signed-cancellation wall; A2's continuum cover is complete (130 boxes, 0 fail; above).**

---

## B. Closable now (review "closable now")

| Task | Method | Status |
|---|---|---|
| DNLS weak persistence `‖z-z_lin‖_inf <= |g|T P^{3/2}` | Duhamel | ALREADY proved (`prop:dnls-persist`) |
| `N=2p, 4p` excess `O(1/N)` | sharp two-sided | ALREADY proved |
| Self-trapping threshold `gP=4` energy-sharp; converse to `0.43` | energy / virial | ALREADY proved |
| Dirichlet segment mod-4 reachability | port the ring criterion | **DONE** (`cor:segceiling`, `exact/gap/segment_mod4.g`; criterion <=> `N+1` prime/`2^m`, `N<=40`) |
| Exact `A_N` table, small N | orbit-closure subtorus max + certified upper | **DONE** — Lipschitz-grid certifies `A_N<U_N` w/ explicit gap `N=6,9,12,15` (`verify_defect_certified.py`); moment-SOS `N=21,33,35`; orbit-closure values `verify_exact_AN.py`. Tight two-sided pin still needs SOS/SDP. |
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
- **Self-trapping window `gP in [0.43,4)`** (review #8): reduced to a focusing interaction-Morawetz
  `L^4` bound = a bound on the staggered-momentum drift (`eq:stagger`). **Leading order now PROVED**
  (`lem:stagger-leading`): along the linear flow from the single-site seed `a_n=e^{2it}(-i)^n J_n(2t)`,
  the distance-2 current `Im(conj(a_n) a_{n+2})=0` identically (phase `(-i)^2=-1` is real => `a_{n+2}` is a
  real multiple of `a_n`), so the O(g) drift cancels: `sup_t|ΔS|~(gP)^2` (exp 2.08, N-uniform),
  `inf<v^2>≈0.49>0`. Verified Maxima (ODE residual 0, phase -1) + Rocq (`formal/rocq/StaggeredCurrent.v`,
  real-multiple=>zero-current) + numerics (`verify_staggered_current.py`, `verify_modified_energy_window.py`,
  GPU). OPEN: N-uniform bound on the `O(g^2)` remainder + nonperturbative top `gP->4` (focusing sign
  survives there — still a signed-cancellation motif).
- **Breather stability** (review #4): VK/GSS stated; (i) `dP/dOmega>0` analytic. **DONE this session** —
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
