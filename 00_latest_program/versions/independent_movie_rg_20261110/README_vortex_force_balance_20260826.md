# Vortex/Q and conditional force-balance DNS version

This source tree contains the 2026-08-26 diagnostic extension for the
rotating moist Rainy–Bénard DNS.  The previous source tree is preserved in
`E:\moist RB\rotating_case_inventory\05_source_reference\source_before_vortex_force_20260826`.

## What was added

1. `force_balance_operators.f90` extracts the staggered-grid momentum stencil
   from the full-field diagnostic.  It is shared by the full-field and
   vortex-conditioned calculations, so the definitions of inertia, Coriolis,
   pressure (including `-dp/dz`), viscosity, buoyancy, and the finite-difference
   transient term cannot silently diverge.
2. `dynamical_vortex_diagnostics.f90` implements the notebook's 2-D z-cut
   criterion.  At the grid layer nearest the physical height `vortex_q_z`
   (default `0.7`), it forms cell-centered `u,v`, periodic
   centered derivatives, `Q_2D`, and `omega_z`.  It identifies periodic
   8-neighbour components using both `Q_2D>0` and `Q_2D>Q_rms`, tracks component
   IDs across frames, and writes births, merges, breakups, and deaths.
3. `vortex_force_balance_online.f90` computes conditional plane RMS forces
   inside the current `Q_2D>Q_rms` components.  The horizontal mean of each
   force component on the selected plane is removed before the RMS.  The output
   is a plane-conditioned diagnostic; it is not a 3-D vortex-volume average.

## New outputs

```text
diagnostics/vortex_q/vortex_q_positive_history.dat
diagnostics/vortex_q/vortex_q_core_history.dat
diagnostics/vortex_q/vortex_q_count.dat
diagnostics/vortex_q/tracker_q_positive.restart
diagnostics/vortex_q/tracker_q_core.restart
diagnostics/vortex_q/vortex_force_balance.dat
diagnostics/vortex_q/vortex_force_balance_summary.dat
```

`vortex_q_count.dat` columns are:

```text
time z Q_rms N(Q>0) Rmax(Q>0) N(Q>Q_rms) Rmax(Q>Q_rms)
```

`vortex_force_balance.dat` has one row per current Q-core.  After time, ID,
status, area, and radius, the force groups are ordered as:

```text
F_I F_C F_P F_V F_B F_T
F_I,h F_C,h F_P,h F_V,h F_B,h F_T,h
F_I,z F_C,z F_P,z F_V,z F_B,z F_T,z
R_G R_z R_momentum
```

The status code is `0=alive`, `1=birth`, `2=merge`, `3=breakup`, `4=death`.

## Build on Shuguang/xh5

The binary currently deployed for new Shuguang cases is
`simexec_vortex_force_balance_20260826` (also copied as `simexec`).  Its
SHA256 is:

```text
4df94de4ed9b5fe227889887f65542a91a571925c636711c0965dcd3f56530aa
```

Load the Shuguang modules and use the platform-specific makefile.  On the
current `xh5` environment, the successful build used the parallel-HDF5
`h5pfc` wrapper from the user's local Intel-MPI HDF5 installation, together
with the system FFTW and MKL libraries.  The exact command was:

```bash
module purge
module load compiler/intel/2017.5.239
module load mpi/hpcx/intel-2017.5.239
module load mathlib/hdf5/1.8.20-intel-2017
module load mathlib/fftw/3.3.8-intel-2017-double
export PATH=/work/home/jiasenzhang/software/hdf5-1.10.5-impi2017/bin:$PATH
export LD_LIBRARY_PATH=/work/home/jiasenzhang/software/hdf5-1.10.5-impi2017/lib:$LD_LIBRARY_PATH
make -f Makefile.shuguang clean
make -f Makefile.shuguang -j1
```

The platform-specific makefile uses `h5pfc` and does not submit a job.

The case builder is already configured to copy the deployed remote binary
from:

```text
/work/home/jiasenzhang/rotating_case/latest_program/source/simexec
```

Therefore a newly prepared Shuguang case receives the Q/vortex and
conditional-force diagnostics automatically.  The builder writes the final
`bou.in` movie row as
`mov_zcut_k  tframe_me  stat_me  vortex_lc  vortex_q_z`; `vortex_lc` is
still a manually supplied value and the Q-component filter is
`radius > 2*vortex_lc`.  `vortex_q_z` defaults to `0.7`.

## Interpretation boundary

The Q diagnostic is deliberately called `Q_2D`: a single z cut cannot provide
the vertical derivatives required by a full 3-D Q criterion.  Use the existing
MSE diagnostics for wet/aggregation structures and use `vortex_q` for the
dynamical rotational structures.  Do not compare their component counts as if
they were the same object population.
