# 3. Strict separation of proved / computational / conjectured

Status: accepted

## Context

The work mixes rigorous theorems (cyclotomic independence, the rank `phi(2N)/2`, the mod-4 criterion,
`B_N = Theta(sqrt N)`), claims verified only numerically over a finite range (no composite saturator
`N <= 160`, the segment rank), open conjectures (`A_N = Theta(ln N)` for all `N`, the exact constant
`lim B_N/sqrt N`), and an exploratory nonlinear Part II. An external review caught real overclaims where
numerical results were phrased as theorems.

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

- The headline `A_N = Theta(ln N)` is never stated as proved for all `N`; only the upper bound and the
  prime/`2^m`/`2p` classes are. The prefix budget `L_pre` is labeled heuristic, not a lower bound.
- The paper is honest at preprint standard and survives a referee's "is this proved or computed?" pass.
- Cost: more hedged language; some readers may prefer cleaner (over)statements. Accepted.

## Alternatives considered

- **State numerical regularities as facts** ("`A_N = Theta(ln N)`"): rejected — it is an overclaim and was
  the exact failure mode the review flagged.
- **Drop the unproved material**: rejected — the conjectures and the finite-range evidence are part of the
  contribution (they define the research program) and are valuable when honestly labeled.
