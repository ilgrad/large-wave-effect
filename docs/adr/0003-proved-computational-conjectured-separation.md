# 3. Strict separation of proved / computational / conjectured

Status: accepted

## Context

The work mixes rigorous theorems (cyclotomic independence, the rank `phi(2N)/2`, the mod-4 criterion,
`B_N = Theta(sqrt N)`), claims verified only numerically over a finite range, open conjectures, and an
exploratory nonlinear Part II. An external review caught real overclaims where numerical results were
phrased as theorems.

Status update: the former order conjecture `A_N = Theta(ln N)` for all `N` has since been upgraded to
the order theorem (`verify_order_theorem.py`), and the saturation classification `A_N = U_N iff N` is
prime or a power of two has since been upgraded from finite-range evidence to a theorem
(`verify_classification_proof.py`, with exact/formal cross-checks). The remaining live open problems are
now narrower: the sharp Schrödinger limsup constant, the uniform composite excess lemma, the self-trapping
proof gap, and the quasi-1D non-embedded-ring condition.

## Decision

Every claim carries its epistemic status explicitly and consistently:

- **Proved** statements are `theorem`/`proposition`/`corollary` with proofs.
- **Computational** statements say "verified for `N <= ...`, conjectured in general" and name the script.
- **Conjectures** are `conjecture` environments, referenced as such everywhere they appear.
- Attribution is split: borrowed results are cited to their authors (e.g. Filimonov's segment theorems to
  `\cite{F23}`); our additions on top (the constant `C = 4/pi^2`, the prime extension) are marked "we
  identify / we show".
- A consolidated **Limitations** paragraph collects all caveats in one place.

## Consequences

- The headline `A_N = Theta(ln N)` is now stated as a theorem, but only where the palindromic prefix
  proof is invoked. Older exploratory scripts are retained as corroboration and must call it the former
  conjecture, not a live open problem.
- Finite-range computations remain explicitly finite-range unless a separate uniform argument is supplied
  (`verify_classification_proof.py` is such an argument; a raw scan is not).
- The paper is honest at preprint standard and survives a referee's "is this proved or computed?" pass.
- Cost: more hedged language; some readers may prefer cleaner (over)statements. Accepted.

## Alternatives considered

- **State numerical regularities as facts before proof** ("`A_N = Theta(ln N)`"): rejected — it was an
  overclaim at the time and the exact failure mode the review flagged. The statement is now allowed only
  because the proof has been added.
- **Drop the unproved material**: rejected — the conjectures and the finite-range evidence are part of the
  contribution (they define the research program) and are valuable when honestly labeled.
