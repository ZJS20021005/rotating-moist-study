#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="$SCRIPT_DIR/platform_config.sh"
if [[ -f "$CONFIG" ]]; then
source "$CONFIG"
else
  BUNDLE_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
  CASE_ROOT="$BUNDLE_ROOT/03_current_cases"
fi

if [[ -n "${INTEL_SETVARS:-}" && -f "$INTEL_SETVARS" ]]; then
  set +eu
  source "$INTEL_SETVARS" >/dev/null 2>&1 || true
  set -eu
fi
if [[ -n "${HDF5_ROOT:-}" ]]; then
  export PATH="$HDF5_ROOT/bin:$PATH"
  export LD_LIBRARY_PATH="$HDF5_ROOT/lib:${LD_LIBRARY_PATH:-}"
fi
if [[ -n "${FFTW_ROOT:-}" ]]; then
  export LD_LIBRARY_PATH="$FFTW_ROOT/lib:${LD_LIBRARY_PATH:-}"
fi

CHECK_HDF5=0
[[ "${1:-}" == "--check-hdf5" ]] && CHECK_HDF5=1
required=(continua_dsal.h5 continua_q1.h5 continua_q2.h5 continua_q3.h5 continua_qvap.h5 continua_grid.dat bou.in simexec subjob.sh)
fail=0
count=0

while IFS= read -r -d '' run_dir; do
  count=$((count + 1))
  case_name="$(basename "$(dirname "$run_dir")")"
  missing=0
  for name in "${required[@]}"; do
    if [[ ! -s "$run_dir/$name" ]]; then
      echo "MISSING: $case_name/$name"
      missing=1
      fail=1
    fi
  done
  nread="$(awk '/Restart flag/{getline; getline; print $1; exit}' "$run_dir/bou.in" 2>/dev/null || true)"
  if [[ "$nread" != "1" ]]; then
    echo "BAD NREAD: $case_name has NREAD=$nread"
    fail=1
  fi
  if [[ "$CHECK_HDF5" -eq 1 && "$missing" -eq 0 ]]; then
    if command -v h5dump >/dev/null 2>&1; then
      for name in continua_dsal.h5 continua_q1.h5 continua_q2.h5 continua_q3.h5 continua_qvap.h5; do
        h5dump -H "$run_dir/$name" >/dev/null || {
          echo "BAD HDF5: $case_name/$name"
          fail=1
        }
      done
    else
      echo "NOTICE: h5dump is unavailable; size checks only for $case_name"
    fi
  fi
  if [[ "$missing" -eq 0 && "$nread" == "1" ]]; then
    echo "OK: $case_name"
  fi
done < <(find "$CASE_ROOT" -mindepth 2 -maxdepth 2 -type d -name run -print0 | sort -z)

if [[ "$count" -ne 14 ]]; then
  echo "Expected 14 cases but found $count"
  fail=1
fi
[[ "$fail" -eq 0 ]] || exit 1
echo "All $count continuation cases passed."
