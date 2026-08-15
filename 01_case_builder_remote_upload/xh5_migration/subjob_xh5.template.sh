#!/usr/bin/env bash
#SBATCH --job-name=__JOB_NAME__
#SBATCH --partition=xhacnormalc
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --ntasks-per-node=64
#SBATCH --time=30-00:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

set +eu
source "$HOME/software/intel_hpckit/setvars.sh" >/dev/null 2>&1 || true
set -eu
export PATH="$HOME/software/hdf5-1.10.5-new/install/bin:$PATH"
export LD_LIBRARY_PATH="$HOME/software/hdf5-1.10.5-new/install/lib:$HOME/software/fftw-3.3.7/lib:${LD_LIBRARY_PATH:-}"

mkdir -p fact flowmov data contstr movie diagnostics
mpirun -np "${SLURM_NTASKS:-64}" ./simexec
