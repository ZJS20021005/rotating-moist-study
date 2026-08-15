#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="$SCRIPT_DIR/platform_config.sh"
[[ -f "$CONFIG" ]] || { echo "Missing $CONFIG"; exit 2; }
source "$CONFIG"

DO_SUBMIT=0
[[ "${1:-}" == "--submit" ]] && DO_SUBMIT=1
bash "$SCRIPT_DIR/check_cases.sh"

count=0
while IFS= read -r -d '' run_dir; do
  count=$((count + 1))
  if [[ "$DO_SUBMIT" -eq 0 ]]; then
    if [[ "$SUBMIT_STYLE" == "csub" ]]; then
      echo "cd '$run_dir' && $SUBMIT_COMMAND < subjob.sh"
    else
      echo "cd '$run_dir' && $SUBMIT_COMMAND subjob.sh"
    fi
    continue
  fi
  cd "$run_dir"
  if [[ "$SUBMIT_STYLE" == "csub" ]]; then
    "$SUBMIT_COMMAND" < subjob.sh
  else
    "$SUBMIT_COMMAND" subjob.sh
  fi
done < <(find "$CASE_ROOT" -mindepth 2 -maxdepth 2 -type d -name run -print0 | sort -z)

if [[ "$DO_SUBMIT" -eq 0 ]]; then
  echo "Previewed $count submissions. Re-run with --submit only after review."
else
  echo "Submitted $count cases."
fi

