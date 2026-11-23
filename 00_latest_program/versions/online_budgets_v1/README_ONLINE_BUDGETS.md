# DNS online-budgets-v2

Baseline: Hainan `lm_merger_beta1p4_program_20261108_mpi_safe`, matching the
`conti4/run/simexec` SHA-256
`426abd61da0d38b091d9800181bec81199b6314a0fe63681ad8b58c3f4aac821`.

The executable enables the four budgets at `TSTAT_DIAG`. It skips the
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
- `gmse_budget.csv`: horizontal variance of the column-integrated moist
  static energy anomaly, with horizontal advection, horizontal diffusion,
  and the signed upper-minus-lower wall diffusion contribution.

## GMSE discretization

The diagnostic forms `M = sum_k dz(k) * (b + gamma*q - mbar(k))`, using
the scalar cell widths `zcr(k+1)-zcr(k)`, and writes `mean_xy(M**2)`.
The column integral has no height normalization or factor of one half.
It reuses the LMSE face fluxes and horizontal diffusion stencil. The
column-wise boundary flux uses the same one-sided scalar gradients and
wall values as the scalar solver. Vertical advection integrates to zero
at the impermeable walls. The buoyancy and vapor condensation sources
cancel in `b + gamma*q`; the prognostic scalar equations have no separate
background-gradient source.

The scalar diffusivities are `kappa_b = 1/sqrt(Ra*Prs)` and
`kappa_q = Sm*kappa_b`. The code uses the combined diffusion of
`b + gamma*Sm*q` for both the horizontal and wall terms. For `Sm=1`,
the horizontal term `GMSE_DH` equals the diffusivity times the squared
horizontal face gradients of `M`, averaged horizontally. For `Sm!=1`,
it is the corresponding cross-gradient expression and may have either
sign. The wall term `GMSE_SB` is a signed covariance with the top minus
bottom diffusive gradient flux, including when wall values are uniform.

`GMSE_rhs = 2*(GMSE_PH - GMSE_DH + GMSE_SB)`. `GMSE_PH_ibp` independently
pairs each integrated horizontal scalar face flux with the matching
discrete gradient of `M`, using horizontal periodicity. `GMSE_PH_residual`
is `GMSE_PH - GMSE_PH_ibp`. The scalar grid is decomposed solely in `z`.
One packed MPI reduction sums the five integrated horizontal arrays onto
rank zero, which computes the final averages and writes the CSV.

The GMSE columns are `time,GMSE_variance,GMSE_PH,GMSE_DH,GMSE_SB,GMSE_rhs,
d_GMSE_variance_dt,rhs_interval,closure,closure_abs,closure_rel,
GMSE_PH_ibp,GMSE_PH_residual`. `closure_abs` is `abs(closure)` and
`closure_rel` divides it by the larger magnitude of the measured rate and
interval RHS, with a machine-tiny floor. Files use the existing budget
restart and append behavior.

## Interval closure and build

The two kinetic RHS columns include signed linear drag work when `alpha`
is nonzero. `closure` is the change per unit time minus the trapezoidal
mean of both endpoint RHS values, evaluated at actual diagnostic times.
The first row of each execution has NaN for the interval quantities.
`production_residual` is `P_H+P_V-P_grad`; a stretched grid need not obey
the continuous scalar product rule exactly.

Build with the baseline compiler and library settings. Use serial `make`
on a clean tree because the inherited Makefile has dependencies between
Fortran module files that are implicit in the serial build order.

## Release executables

Each cluster has a build for its own MPI and library stack. Use the
corresponding `~/rotating_case/program_versions/online-budgets-v2/simexec`.

| Platform | Executable SHA-256 |
| --- | --- |
| Hainan | `ebcafcc199d2b41b8cacae94ccba022480314bf8aaebe0b80a35de0245f69ddf` |
| xh5 | `48e2c9a7a94190cdb473fb83bde1e9ec09bf6ec1dc754fd1b5c10df5cf4a1b15` |
| Local Shuguang | `c2068a6106ea71a508a341967801174a2da29ba2f644fdf2e2fb1c5c353e2d0c` |

The v2 source remains in `versions/online_budgets_v1` to preserve the
existing source entry path. Check `VERSION` or the Git tag for the release.
Use `verify_gmse_smoke.py --one-rank RUN1 --two-rank RUN2 --baseline RUN_V1
--offline` to reproduce the checks on run directories started from the same
checkpoint and sampled at the same times. The results are in
`gmse_smoke_results.json`.

## Previous release (v1)

Executable: `~/rotating_case/program_versions/online-budgets-v1/simexec`.
SHA-256: `8bd5d712db47df33c47c29844cd63f88189e6b66008dc6ca1310898ade32c031`.
This release is based on the named continuation executable. Use its explicit
path when reproducing the previous release. The baseline source remains
available at Git tag `online-budgets-v1`.

`budget_smoke_results.json` records a 40-step, 33×33×17 DNS check with a
no-slip bottom and free-slip top. Baseline and new 2-rank checkpoints have
identical values in all five state arrays. In this early-stage smoke case,
peak interval closure divided by peak interval RHS was 0.00318%, 0.0122%,
and 0.106% for LMSE, full KE, and 2D KE, respectively. The 2D RHS was very
small in this short run. The independent movie-based legacy energy ratio
agreed within `3.8e-25` in absolute value. The reported 0.027 s budget time
covers 40 samples on 2 ranks at this small resolution; it is not a
production-grid performance measurement.
