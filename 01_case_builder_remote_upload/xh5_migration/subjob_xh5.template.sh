#!/usr/bin/env bash
#SBATCH -J __JOB_NAME__
#SBATCH -N 1
#SBATCH --ntasks-per-node=64
#SBATCH -p xhacnormalc
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.err

set -euo pipefail

WORKDIR="${SLURM_SUBMIT_DIR:-$PWD}"
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
# xh5's srun route layer fails for a one-node hostlist.  Run Hydra locally;
# all current continuation cases request exactly one node.
unset I_MPI_PMI_LIBRARY || true
unset I_MPI_FABRICS || true
export I_MPI_HYDRA_BOOTSTRAP=fork

mkdir -p fact flowmov data contstr movie diagnostics
mpirun -launcher fork -np "${SLURM_NTASKS:-64}" ./simexec
