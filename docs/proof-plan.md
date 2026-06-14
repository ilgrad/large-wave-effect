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
arithmetic (destroys the cancellation). A breakthrough on the signed estimate closes both.

---

## A. Central problems (deeply attacked this session)

### A1. Excess lemma / composite defect `U_N - A_N`  (review #1)

**Goal.** For `N = 2^a m`, `m` odd: `U_N - A_N = (1/pi) ln(m/phi(m)) + O(1)`; unbounded,
`Theta(ln ln ln N)` along primorials. Reduces to the excess lemma `A_N <= L_pre + O(1)`.

**Method (review + session).** Rewrite as a max of a trig polynomial on the relation subtorus:
`A_N = max_x [ sum_{r<=M} b_r sin x_r + sum_{r>M} b_r sin(l_r . x) ]`, `M = phi(2N)/2`,
`b_r = 1/(N sin(pi r/N))`. Seek a **dual / majorant certificate**
`sum_{r>M} b_r sin(l_r . x) <= C + sum_{r<=M} b_r (1 - sin x_r)`. Per-N via trig-SOS / Lasserre;
families via cyclotomic block decomposition (p-block + q-block) and the limit `a -> infinity`.

**Status.**
- PROVED: `omega(m) <= 1` (`N=2p` sharp `1/(2N)`, `N=4p`); rank `phi(2N)/2`; mod-4 criterion;
  classification `A_N=U_N <=> N` prime/`2^m`; order `A_N ~ (1/pi) ln N`.
- DONE this session (`omega(m)=2`): relation lattice characterized — unimodular SNF, **chained**
  relations, support set by the smallest prime (`3|N`: `w_a+w_b=w_c`, `a+b=N/3`, from `2 sin(pi/6)=1`;
  `5|N`: golden); per-N certificates `N=15` (elementary Lipschitz grid, `>=0.082`), `N=21,33,35`
  (moment-SOS). `exact/gap/excess_omega2.g`, `numerics/verify_excess_omega2.py`.
- TRIED → blocked: the **uniform** dual certificate is exactly the signed bound that overshoots.
  The obstruction is a **modest O(1) balance** (dependent_net ~ +0.15 vs prefix_deficit ~ -0.07,
  neither growing), not a near-cancellation. Bedert does NOT transfer (4 mismatches; spectral
  quantities do not track dependent_net; `verify_excess_omega2_spectral.py`).
- NEXT (realistic): `N = 2^a pq` for fixed `p,q` via a **signed** block-cancellation lemma per
  cyclotomic block + SOS + the `a -> infinity` limit. A genuine new result if it closes.

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
- **FPUT/chaos** (review #9): small energy via Birkhoff/Nekhoroshev (closable); chaos via CAPD
  computer-assisted (separate project).
- **Multidimensional / graphs** (review #10): `d_s=1` critical (done for products, Möbius);
  arbitrary graphs lose the cyclotomic structure.

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
2. **Excess lemma `2^a pq`** (signed block-cancellation, fixed `p,q`) — a real new theorem.
3. **beta_odd offline tiling** — execution, off the interactive session.
4. **Sup-side MPS** (closes A1 and A2 at once) — the research breakthrough.
