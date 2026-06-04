# 4. GitHub-first + Zenodo, arXiv deferred

Status: accepted

## Context

The work is ready as a preprint. arXiv submission (math.SP / nlin.SI) currently requires endorsement the
author does not yet have. A citable, timestamped, reproducible artifact is wanted now.

## Decision

Publish **GitHub-first**: the repository (paper source, figures, the `verify_*.py` reproduction suite, the
tested `src/` core) is the primary artifact, with a **Zenodo DOI** minted on the first GitHub release for
citability. arXiv stays open for later — a public GitHub/Zenodo record does not burn the arXiv option.

## Consequences

- A reader can cite the Zenodo DOI and reproduce every number with one command; the artifact is the paper
  *and* its checks, not just a PDF.
- `CITATION.cff` is present; on the first release add `version`, `date-released`, `repository-code`, and
  the Zenodo DOI.
- The paper is structured to stand alone (self-contained bibliography, reproducibility appendix) so it can
  be moved to arXiv unchanged later.

## Alternatives considered

- **Wait for arXiv endorsement first**: rejected — delays a finished, citable artifact for an
  administrative step; GitHub+Zenodo provides the timestamp and DOI now.
- **PDF-only release**: rejected — discards the reproduction suite, which is a core part of the value.
