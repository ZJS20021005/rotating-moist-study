# xh5 migration workflow

This folder records the PBS/qsub migration used for the 14 current `Ra=8e6`
continuation cases on the `xh5` cluster.

Remote destination:

`/work/home/jiasenzhang/rotating_moist_migration_bundle_20260815`

Deployment status, verified on 2026-08-15:

- The source in `01_dns_program/source_merged` was compiled successfully on
  xh5 with Intel MPI 2021.6, parallel HDF5 1.10.5, and FFTW 3.3.7.
- The deployed `simexec_new_platform` SHA-256 is
  `f4504bda4e9f7524664fb615f61a879a0d8eb5fa5c2d24f1abfd4b263ee8b371`.
- The binary, `NREAD=1`, `TMAX=500d0`, and the PBS `subjob.sh` have been
  installed and checked for all fourteen cases: `Ek1e-1`, `Ek5e-2`,
  `Ek3e-2`, `Ek1e-2`, `Ek7e-3`, `Ek5e-3`, `Ek3e-3`, `Ek2e-3`, `Ek1e-3`,
  `Ek7e-4`, `Ek5e-4`, `Ek2e-4`, `Ek1p5e-4`, and `norotating`.
- This stage only uploads and verifies the continuation setup. No `qsub`
  command has been issued; the queue was empty at verification time.

The migration uses:

- the PBS-compatible `qsub` entry point with queue `low`;
- one node and 32 MPI ranks per case;
- 500 additional physical time units with `NREAD=1`;
- `$HOME/software/intel_hpckit/setvars.sh`;
- parallel HDF5 under `$HOME/software/hdf5-1.10.5-new/install`;
- FFTW under `$HOME/software/fftw-3.3.7`.

The upload and preparation stage does not call `qsub`. Preview submissions
with `bash 04_platform_scripts/submit_all.sh`; submit only with the explicit
`--submit` option.
