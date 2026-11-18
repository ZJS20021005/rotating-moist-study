#!/usr/bin/env bash
set -euo pipefail

# Canonical production build for the Shuguang platform.
module purge
module load compiler/intel/2017.5.239
module load mpi/intelmpi/2017.4.239
module load mathlib/fftw/3.3.8-intel-2017-double

export HDF5_ROOT="${HDF5_ROOT:-$HOME/software/hdf5-1.10.5-impi2017}"
export FFTW_ROOT="${FFTW_ROOT:-/public/software/mathlib/fftw/3.3.8/double/intel}"
export PATH="$HDF5_ROOT/bin:$PATH"
export LD_LIBRARY_PATH="$HDF5_ROOT/lib:$FFTW_ROOT/lib:${LD_LIBRARY_PATH:-}"

# Compile serially because the original Makefile does not encode every module
# dependency needed by a parallel make.
make -f Makefile.shuguang clean
make -f Makefile.shuguang

chmod 755 simexec
sha256sum simexec
echo "built production executable: $(pwd)/simexec"
