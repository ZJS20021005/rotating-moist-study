#!/usr/bin/env bash
#PBS -N __JOB_NAME__
#PBS -q xhacnormalc
#PBS -l nodes=1:ppn=32
#PBS -j oe

set -euo pipefail

WORKDIR="${PBS_O_WORKDIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
cd "$WORKDIR"
export OMP_NUM_THREADS=1

# The executable was built with these Intel MPI/HDF5/FFTW installations.
# Intel's setup script is not guaranteed to be nounset-clean.
set +eu
source "$HOME/software/intel_hpckit/setvars.sh" >/dev/null 2>&1
intel_status=$?
set -eu
if [[ "$intel_status" -ne 0 ]]; then
  echo "Failed to load the Intel oneAPI environment." >&2
  exit "$intel_status"
fi
export PATH="$HOME/software/hdf5-1.10.5-new/install/bin:$PATH"
export LD_LIBRARY_PATH="$HOME/software/hdf5-1.10.5-new/install/lib:$HOME/software/fftw-3.3.7/lib:${LD_LIBRARY_PATH:-}"

mkdir -p fact flowmov data contstr movie diagnostics
mpirun -np "${PBS_NP:-32}" ./simexec
