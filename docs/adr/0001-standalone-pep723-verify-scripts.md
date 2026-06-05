# 1. Standalone PEP 723 verify scripts over a shared library

Status: accepted

## Context

Every quantitative claim in the paper (asymptotic constants, ranks, operator norms, the ceiling
criterion, the Schrödinger bounds, the nonlinear simulations) must be reproducible by a reader who clones
the repository. There are two natural ways to organize the ~50 verification programs in `numerics/`:

- import a shared, tested core package (`src/large_wave_effect/`) — DRY, one source of truth;
- make each `verify_*.py` a self-contained PEP 723 script with inline dependencies, runnable with
  `uv run --script numerics/verify_X.py` and no package install.

The core mathematics is genuinely small (ring spectrum, cyclotomic rank, Schrödinger kernel — ~230 lines
in `src/`), while the verification surface is large and one-claim-per-file.

## Decision

Each `verify_*.py` is a **standalone PEP 723 script** with its own inline dependency block. The scripts do
**not** import `large_wave_effect`; primitives that several scripts need (`cyclotomic`, `power_mod`,
mod-`p` rank, `U_N`, the integer Hermite kernel) are duplicated per script. `src/large_wave_effect/` stays
as the **tested reference core** exercised only by `tests/` (pytest).

## Consequences

- A reviewer reproduces any single claim in isolation — clone, `uv run --script <one file>`, done; no
  build, no install, no cross-file state. This is the property that matters for a preprint companion.
- Cost: the core (`cyclotomic`, `power_mod`, rank, ceiling) is copied across ~13 scripts; a bug fix must
  be applied in each place (this happened during development). Mitigated by the pytest-covered `src/`
  acting as the canonical implementation to diff against.
- `src/` is therefore small (~3.5% of the code) and load-bearing only for the tests, not for the
  paper's numerics.

## Alternatives considered

- **Scripts import `src/`** (run via `uv run`, project env): eliminates duplication but loses the
  one-command-no-install reproducibility, and a PEP 723 `--script` env cannot see the local package
  without packaging it. Rejected: reproducibility-first wins for a paper artifact.
- **One monolithic `verify_all.py`**: rejected — claims would no longer be checkable in isolation, and a
  single failure would obscure which claim broke.
