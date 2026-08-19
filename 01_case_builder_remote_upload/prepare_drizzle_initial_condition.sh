#!/bin/bash
set -euo pipefail

# Run this from the case's run directory before launching simexec.
cd "$(dirname "$0")"

# The available Python module differs between clusters.  Loading these is
# harmless when a module is unavailable, and lets the search below find the
# cluster's NumPy/SciPy installation.
if command -v module >/dev/null 2>&1; then
    module load anaconda3/2023.09 >/dev/null 2>&1 || true
    module load python/3.8.10 >/dev/null 2>&1 || true
fi

find_python() {
    local candidate
    if [[ -n "${RAINY_DRIZZLE_PYTHON:-}" ]]; then
        printf '%s\n' "$RAINY_DRIZZLE_PYTHON"
        return 0
    fi
    for candidate in python3 python \
        /share/apps/anaconda3/bin/python3 \
        /share/apps/anaconda3/bin/python \
        /public/software/apps/anaconda3/2023.09/bin/python; do
        if command -v "$candidate" >/dev/null 2>&1 || [[ -x "$candidate" ]]; then
            if "$candidate" -c "import numpy, scipy" >/dev/null 2>&1; then
                printf '%s\n' "$candidate"
                return 0
            fi
        fi
    done
    return 1
}

DRIZZLE_PYTHON="$(find_python || true)"
if [[ -z "$DRIZZLE_PYTHON" ]]; then
    echo "ERROR: Python with NumPy/SciPy not found." >&2
    echo "Set RAINY_DRIZZLE_PYTHON to a Python executable that can import numpy and scipy." >&2
    exit 21
fi

"$DRIZZLE_PYTHON" -c "import numpy, scipy" >/dev/null

tmpdir="$(mktemp -d .drizzle_generate.XXXXXX)"
trap 'rm -rf "$tmpdir"' EXIT

"$DRIZZLE_PYTHON" ./generate_drizzle_initial_condition.py \
    --bou-in ./bou.in \
    --profile "$tmpdir/drizzle_init.dat" \
    --metadata "$tmpdir/drizzle_init_meta.json" \
    --profile-n 401 \
    --perturb-amp 1e-4 \
    --saturation-width 1e-8 \
    > "$tmpdir/generator_output.json"

"$DRIZZLE_PYTHON" ./check_drizzle_before_submit.py \
    --bou-in ./bou.in \
    --profile "$tmpdir/drizzle_init.dat" \
    --expected-perturb 1e-4

mv "$tmpdir/drizzle_init.dat" ./drizzle_init.dat
mv "$tmpdir/drizzle_init_meta.json" ./drizzle_init_meta.json
mv "$tmpdir/generator_output.json" ./drizzle_init_last_generation.json

echo "Ready: drizzle_init.dat was regenerated from the current bou.in"
