#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="$SCRIPT_DIR/platform_config.sh"
if [[ ! -f "$CONFIG" ]]; then
  echo "Missing $CONFIG"
  echo "Copy platform_config.sh.example to platform_config.sh and edit it first."
  exit 2
fi
source "$CONFIG"

if [[ -n "${MODULE_INIT:-}" && -f "$MODULE_INIT" ]]; then
  source "$MODULE_INIT"
fi
if command -v module >/dev/null 2>&1; then
  module purge
  [[ -n "${MPI_MODULE:-}" ]] && module load "$MPI_MODULE"
  [[ -n "${FFTW_MODULE:-}" ]] && module load "$FFTW_MODULE"
  [[ -n "${HDF5_MODULE:-}" ]] && module load "$HDF5_MODULE"
fi
if [[ -n "${INTEL_SETVARS:-}" && -f "$INTEL_SETVARS" ]]; then
  source "$INTEL_SETVARS" >/dev/null 2>&1
fi

if [[ -n "${HDF5_ROOT:-}" ]]; then
  export PATH="$HDF5_ROOT/bin:$PATH"
  export LD_LIBRARY_PATH="$HDF5_ROOT/lib:${LD_LIBRARY_PATH:-}"
fi
if [[ -n "${FFTW_ROOT:-}" ]]; then
  export LD_LIBRARY_PATH="$FFTW_ROOT/lib:${LD_LIBRARY_PATH:-}"
fi

cd "$SOURCE_ROOT"
make clean
make FFTW_ROOT="${FFTW_ROOT:-$HOME/lib/fftw-3.3.7}"
test -s simexec
chmod +x simexec

mkdir -p "$REFERENCE_BINARY_ROOT"
cp -p simexec "$REFERENCE_BINARY_ROOT/simexec_new_platform"
sha256sum simexec "$REFERENCE_BINARY_ROOT/simexec_new_platform"
echo "Build complete: $REFERENCE_BINARY_ROOT/simexec_new_platform"
