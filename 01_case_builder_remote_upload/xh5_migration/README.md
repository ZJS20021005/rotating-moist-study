# xh5 migration workflow

This folder records the Slurm migration used for the 14 current `Ra=8e6`
continuation cases on the `xh5` cluster.

Remote destination:

`/work/home/jiasenzhang/rotating_moist_migration_bundle_20260815`

The migration uses:

- Slurm partition `xhacnormalc`;
- one node and 64 MPI ranks per case;
- 500 additional physical time units with `NREAD=1`;
- `$HOME/software/intel_hpckit/setvars.sh`;
- parallel HDF5 under `$HOME/software/hdf5-1.10.5-new/install`;
- FFTW under `$HOME/software/fftw-3.3.7`.

The upload and preparation stage does not call `sbatch`. Preview submissions
with `bash 04_platform_scripts/submit_all.sh`; submit only with the explicit
`--submit` option.
