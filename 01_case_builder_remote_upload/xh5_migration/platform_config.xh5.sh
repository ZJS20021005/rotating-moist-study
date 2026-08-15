#!/usr/bin/env bash

# Configuration for the xh5 / hpccube PBS-compatible qsub platform.
BUNDLE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SOURCE_ROOT="$BUNDLE_ROOT/01_dns_program/source_merged"
CASE_ROOT="$BUNDLE_ROOT/03_current_cases"
REFERENCE_BINARY_ROOT="$BUNDLE_ROOT/01_dns_program/reference_binaries"

MODULE_INIT=""
MPI_MODULE=""
FFTW_MODULE=""
HDF5_MODULE=""
INTEL_SETVARS="$HOME/software/intel_hpckit/setvars.sh"
HDF5_ROOT="$HOME/software/hdf5-1.10.5-new/install"
FFTW_ROOT="$HOME/software/fftw-3.3.7"

RUN_COMMAND="srun --nodes=1 --ntasks=\${SLURM_NTASKS:-64} --mpi=pmi2 ./simexec"
SUBMIT_STYLE="pbs"
SUBMIT_COMMAND="qsub"
PBS_QUEUE="xhacnormalc"
PBS_NP="64"
