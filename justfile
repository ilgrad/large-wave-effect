# Verification loop for this repository.
#
# `just check` is the fast gate (format, lint, types, tests): run it before calling a change done.
# The certificate layers take minutes and run separately: `just exact`, `just dagger`, `just paper`.

default:
    @just --list

# Everything fast, cheapest first. Stops at the first failure.
check: fmt lint types types-scripts test

fmt:
    uv run ruff format --check .

lint:
    uv run ruff check .

# The package and its tests, in the project environment.
types:
    uv run ty check

# Every numerics/ script in its own PEP 723 environment; `types` excludes them.
types-scripts:
    #!/usr/bin/env bash
    set -euo pipefail
    one() {
        uv sync --quiet --script "$1"
        py=$(uv python find --script "$1")
        # the interpreter uv resolved, not the declared floor: the resolved numpy's stubs need >= 3.12
        ver=$("$py" -c 'import sys; print("%d.%d" % sys.version_info[:2])')
        out=$(uv run --quiet ty check --output-format concise --python "$py" --python-version "$ver" \
            --extra-search-path src --extra-search-path numerics "$1" 2>&1) && return 0
        printf '%s\n' "$out"
        return 1
    }
    export -f one
    git ls-files -z 'numerics/*.py' | xargs -0 -P 8 -n 1 bash -c 'one "$0"'

test:
    uv run pytest -q

fix:
    uv run ruff check --fix .
    uv run ruff format .

# GAP, PARI/GP, Julia, FLINT and Rocq certificates, cross-checked against the Python core.
exact:
    timeout 3600 uv run --script numerics/verify_exact_layer.py

# Re-decide, recompute and check the band coverage of the rigorous dagger cover.
dagger:
    timeout 1800 uv run --script numerics/verify_dagger_rigorous.py --recompute 3 --crosscheck

# Re-decide, recompute and check the coverage of the rigorous near-pair basin certificate.
nearpair:
    timeout 1800 uv run --script numerics/verify_dagger_nearpair.py --recompute 3 --crosscheck

# Build the paper and fail on unresolved references (latexmk exits 0 with "??" in the PDF).
paper:
    #!/usr/bin/env bash
    set -euo pipefail
    cd paper
    timeout 600 latexmk -pdf -interaction=nonstopmode -synctex=1 large-wave-effect.tex
    ! rg -n 'Undefined control sequence|! LaTeX Error|Reference .* undefined|Citation .* undefined|multiply-defined' large-wave-effect.log
