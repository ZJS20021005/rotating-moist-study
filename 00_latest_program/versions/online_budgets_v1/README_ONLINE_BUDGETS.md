# DNS online-budgets-v1

Baseline: Hainan `lm_merger_beta1p4_program_20261108_mpi_safe`, matching the
`conti4/run/simexec` SHA-256
`426abd61da0d38b091d9800181bec81199b6314a0fe63681ad8b58c3f4aac821`.

The executable enables the three budgets at `TSTAT_DIAG`. It skips the
MSE and Q vortex trackers and the Q-conditioned vortex force output.
The other output schedules, input format, restart format, and solver
equations remain as in the baseline.

CSV files appear in `diagnostics/budget/`:

- `lmse_budget.csv`: variance and transport terms use the scalar mesh;
  `D_m` is minus the complete discrete diffusion contribution. Condensation
  cancels and unequal diffusivities are included through `Sm`.
- `ke_budget.csv`: buoyancy and viscous work use the native velocity
  locations and control-volume weights.
- `ke2d_budget.csv`: native depth-weighted horizontal velocities, Reynolds
  stress transfer, horizontal dissipation, and top-minus-bottom wall work.
  Additional columns reproduce the single-snapshot, cell-center energy ratio
  used by the offline movie analysis.

The two kinetic RHS columns include signed linear drag work when `alpha`
is nonzero. `closure` is the change per unit time minus the trapezoidal
mean of both endpoint RHS values, evaluated at actual diagnostic times.
The first row of each execution has NaN for the interval quantities.
`production_residual` is `P_H+P_V-P_grad`; a stretched grid need not obey
the continuous scalar product rule exactly.

Build with the baseline compiler and library settings. Use serial `make`
on a clean tree because the inherited Makefile has dependencies between
Fortran module files that are implicit in the serial build order.

## Hainan release

Executable: `~/rotating_case/program_versions/online-budgets-v1/simexec`.
SHA-256: `8bd5d712db47df33c47c29844cd63f88189e6b66008dc6ca1310898ade32c031`.
This release is based on the named continuation executable. Use its explicit
path when preparing a new run from that line of cases.

`budget_smoke_results.json` records a 40-step, 33×33×17 DNS check with a
no-slip bottom and free-slip top. Baseline and new 2-rank checkpoints have
identical values in all five state arrays. In this early-stage smoke case,
peak interval closure divided by peak interval RHS was 0.00318%, 0.0122%,
and 0.106% for LMSE, full KE, and 2D KE, respectively. The 2D RHS was very
small in this short run. The independent movie-based legacy energy ratio
agreed within `3.8e-25` in absolute value. The reported 0.027 s budget time
covers 40 samples on 2 ranks at this small resolution; it is not a
production-grid performance measurement.
