#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="$SCRIPT_DIR/platform_config.sh"
[[ -f "$CONFIG" ]] || { echo "Missing $CONFIG"; exit 2; }
source "$CONFIG"

NEW_BINARY="$REFERENCE_BINARY_ROOT/simexec_new_platform"
[[ -s "$NEW_BINARY" ]] || {
  echo "Missing $NEW_BINARY; run compile_on_new_platform.sh first."
  exit 2
}

APPLY=0
[[ "${1:-}" == "--apply" ]] && APPLY=1
count=0
while IFS= read -r -d '' run_dir; do
  count=$((count + 1))
  if [[ "$APPLY" -eq 0 ]]; then
    echo "DRY RUN: install new binary -> $run_dir/simexec"
    continue
  fi
  if [[ -f "$run_dir/simexec" && ! -f "$run_dir/simexec.old_platform_20260808" ]]; then
    cp -p "$run_dir/simexec" "$run_dir/simexec.old_platform_20260808"
  fi
  cp -p "$NEW_BINARY" "$run_dir/simexec"
  chmod +x "$run_dir/simexec" "$run_dir/subjob.sh" 2>/dev/null || true
  mkdir -p "$run_dir"/{fact,flowmov,data,contstr,movie,diagnostics}
  echo "INSTALLED: $run_dir/simexec"
done < <(find "$CASE_ROOT" -mindepth 2 -maxdepth 2 -type d -name run -print0 | sort -z)

[[ "$count" -eq 14 ]] || {
  echo "Expected 14 cases but found $count. No submission was attempted."
  exit 3
}
if [[ "$APPLY" -eq 0 ]]; then
  echo "Dry run complete for $count cases. Re-run with --apply after review."
else
  echo "Installed the new-platform binary into $count cases."
fi

