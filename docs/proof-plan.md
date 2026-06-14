# Proof plan

Per-task plan: what to prove, the method, the tools, and the **current status in this repo**
(proved / tried-and-failed-this-session / open + what is actually needed). Cross-references the
external review against the session's own findings.

**Unifying finding.** The two central open targets — the composite excess lemma (1) and the
Schrödinger constant `limsup B_N/sqrt(N) = beta_odd` (5) — reduce to the **same** obstruction: a
**sup-side signed-cancellation estimate** (the primal analogue of McGehee–Pigno–Smith / Konyagin).
Seven routes have been ruled out for it: modulus/triangle bound, multivariate concavity, the
variance-floor, monotone-tail, Berry–Esseen (sharp Esseen constant 1/2 overshoots), Bedert's Chowla
machinery (min-side, does not transfer — verified from arXiv:2509.05260), and naive Arb ball
arithmetic (destroys the cancellation). A breakthrough on the signed estimate closes both. This session
sharpened the (1) side: its signed core is now reduced to a concrete **weighted-packing bound on the
chained relation graph**, uniform over odd `m` (A1 Mechanism).

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

**Open — research wall (needs new mathematics).**
- **A1** excess lemma, uniform in `m` = weighted-packing bound on the chained graph ‖ **A2** β_odd
  continuum cover — the SAME sup-side signed-cancellation estimate (7 routes ruled out).
- Self-trapping window `gP∈[0.43,4)` — same focusing-sign motif (§C).

**Open — execution / resource (idea exists, not a one-shot).**
- β_odd offline tiling (~20–35 min/box); finite-time **lower** bound (+ damping, conditional on it);
  FPUT small-energy (Nekhoroshev, 4-wave resonances) and chaos (CAPD); arbitrary graphs (non-cyclotomic).

**Blocks publication, not the mathematics.**
- Filimonov 1992 (C.R. Acad. Sci.) — physical access for priority attribution.

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
- OPEN: the **continuum** cover between grid nodes. Method built (multivariate amplitude Taylor model
  kills the 4x ball inflation; `verify_dagger_continuum.py`); remaining = sharp signed tail
  (near-degenerate ring beat) + a theta-Taylor model + an **offline tiling** (~20-35 min/box).
  This is execution, not a missing idea.
- FIXED: the `beta_odd` literal (was wrong past digit 9).

**Both A1 and A2 are the same signed-cancellation wall.**

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
  recurrence times to a proof. OPEN: the **family rate** `~ delta^{-m}` does NOT follow from the
  homogeneous Dirichlet pigeonhole — the target `pi/2` is inhomogeneous, needing a Baker-type Diophantine
  lower bound on linear forms in the `omega_r`. The **lower** bound (wait is exponentially long) is the
  real content, still open (LLL / badly-approximable certificates).
- **Self-trapping window `gP in [0.43,4)`** (review #8): reduced this session to a focusing
  interaction-Morawetz `L^4` bound (staggered-momentum identity). Blocked by the focusing sign —
  another signed-cancellation motif. Try sharp discrete Gagliardo–Nirenberg (SOS per N).
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

- **Filimonov 1992** (C.R. Acad. Sci.): obtain and compare on the ring before any "first proof"
  claim. Blocks publication, not the mathematics.
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
