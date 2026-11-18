#!/usr/bin/env bash
set -euo pipefail

# Compile the staged, multi-day diagnostic executable on the Shuguang system.
# Run this script from the source directory on the remote platform.
module purge
module load compiler/intel/2017.5.239
module load mpi/intelmpi/2017.4.239
module load mathlib/fftw/3.3.8-intel-2017-double

export HDF5_ROOT="${HDF5_ROOT:-$HOME/software/hdf5-1.10.5-impi2017}"
export FFTW_ROOT="${FFTW_ROOT:-/public/software/mathlib/fftw/3.3.8/double/intel}"
export PATH="$HDF5_ROOT/bin:$PATH"
export LD_LIBRARY_PATH="$HDF5_ROOT/lib:$FFTW_ROOT/lib:${LD_LIBRARY_PATH:-}"

# The source Makefile has incomplete parallel dependencies; compile serially.
make -f Makefile.shuguang clean PROGRAM=simexec_diag05
make -f Makefile.shuguang PROGRAM=simexec_diag05

chmod 755 simexec_diag05
echo "built: $(pwd)/simexec_diag05"
