#!/bin/bash
set -euo pipefail

# Run this from the case's run directory before launching simexec.
cd "$(dirname "$0")"

DRIZZLE_PYTHON="${RAINY_DRIZZLE_PYTHON:-/share/apps/anaconda3/bin/python3}"
if [[ ! -x "$DRIZZLE_PYTHON" ]]; then
    echo "ERROR: Python with NumPy/SciPy not found: $DRIZZLE_PYTHON" >&2
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
