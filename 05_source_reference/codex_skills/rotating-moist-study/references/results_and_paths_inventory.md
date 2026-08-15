# 当前结果与路径清单

## 主路径

Remote hosts used recently: `c01n0034`, `c01n0037`; connection node may change.

High-resolution and main root:

`/share/org/SHUTUANL/shu_zhangjs/rainy model/ns/transition_study/beta1`

Low-resolution sweep:

`/share/org/SHUTUANL/shu_zhangjs/rainy model/ns/transition_study/beta1/lowrestest`

Aspect-ratio study:

`/share/org/SHUTUANL/shu_zhangjs/rainy model/ns/transition_study/beta1/lowrestest/aspect_ratio_study`

Remote postprocessor index:

`/share/org/SHUTUANL/shu_zhangjs/rainy model/ns/transition_study/beta1/postprocess_scripts_index`

Local organized result root:

`E:/moist RB/moist/result/spectrum/ns/transition_study/RaT1e8pr07/beta1_highres/comparison_outputs/organized_results_20260705`

## Canonical new-program case scope (2026-08-13)

Unless the user explicitly names another dataset, all phrases such as “current
cases”, “new-program cases”, or “all cases” refer only to these 14 `Ra=8e6`
cases:

`Ek1e-1`, `Ek5e-2`, `Ek3e-2`, `Ek1e-2`, `Ek7e-3`, `Ek5e-3`,
`Ek3e-3`, `Ek2e-3`, `Ek1e-3`, `Ek7e-4`, `Ek5e-4`, `Ek2e-4`,
`Ek1p5e-4`, and `norotating`.

Do not silently mix old-program runs, other Rayleigh numbers, other moisture
boundary conditions, AR10 studies, or legacy beta sweeps into these analyses.
Current primary local storage is `G:/moist convection` for 13 cases and
`H:/rotating_case/Pr0p7/Ra8e6/Ek1e-1` for `Ek1e-1`. The intermittent plume-burst
cases are `Ek2e-4`, `Ek5e-4`, and `Ek7e-4`; omit their steady scalar points when
the user requests non-intermittent comparisons and mark their Ek positions.

## 主要结果集

- `profiles`: standard/profile collections for `m`, `q`, `b`, RH, velocity, flux and correlations.
- `ra1e6_research_logic_profiles`: latest remote Ra1e6 profile reduction and boundary-layer comparison.
- `Ek_scaling`: Ek-scaling figures including the legacy `l_h(z=0.75)` convention.
- `funnel_transport_decomposition`: Ra1e6 and Ra1e8 funnel/environment MSE transport decomposition.
- `funnel_identification_visualization`: direct mask QA at selected heights.
- `height_spectra_gaussian_flux`: Ra1e8 high-resolution height spectra and coarse-grained energy flux.
- `height_shell_to_shell_transfer`: height-resolved shell-to-shell transfer products.
- `vortex_geometry`: Q/vorticity-based vortex radius and spacing diagnostics.
- `kolmogorov_resolution`: dissipation-based resolution checks.
- `E:/moist RB/rotating_case_inventory/04_outputs_and_figures/kolmogorov_resolution_audit_20260811`: 2026-08-11 restart-field Kolmogorov/Batchelor, time-step, and near-wall resolution audit for the 14 current `Ra=8e6` migration cases, plus the disk allocation plan. This is a single-snapshot screen; no case currently passes all strict spatial, temporal, and boundary-layer checks.
- `parameter_audit_20260712`: remote parameter audit.
- `thermal_wind`: selected-case thermal-wind diagnostic.
- `time_series_avgvar`: kinetic-energy and legacy `l0` time series.
- `aggregation_dns_module_20260715`: modified online MSE-aggregation DNS module archive.
- `aspect_ratio_study/AR10_kinetic_energy_timeseries`: current AR10 kinetic-energy histories.
- `aspect_ratio_study/AR10_self_aggregation_Ek3e-3_current`: current `Ra1e8, AR10, Ek3e-3` MSE aggregation result.
- `aspect_ratio_study/AR10_Voronoi_latest10_tge800_20260715`: latest-10-field periodic Voronoi products for mature AR10 cases.
- `aspect_ratio_study/AR10_mrms_scales_timeseries_20260715`: current online `sqrt(A_m)` and field-derived `L_peak/L_integral` time series for the AR10 cases using the new MSE diagnostic.
- `aspect_ratio_study/AR10_REMOTE_MSE方差与聚集尺度更新_AR10_remote_MSE_variance_scales_updated_20260715`: authoritative remote-only update of Ra1e8 AR10 `A_m`, `L_peak`, `L_integral`, and periodic connected-cluster radii, with all three cases extended through `t=800`.
- `aspect_ratio_study/AR10_256不旋转MSE方差与聚集尺度_AR10_res256_norotating_MSE_scales_current_20260715`: current remote-only Ra1e8 AR10 nonrotating 256x256x65 MSE variance and aggregation-scale histories.
- `aspect_ratio_study/AR10_MSE算例动能时间序列_AR10_MSE_cases_kinetic_energy_timeseries_20260716`: paired kinetic-energy histories for every current AR10 case that has an online MSE time series.
- `aspect_ratio_study/AR10中层垂直速度谱长度10帧平均_AR10_midplane_w_spectral_length_last10_20260716`: frame-first last-10-field midplane vertical-velocity spectral length for all remote AR10 datasets with complete fields reaching `t>=800`.

## 当前AR10 Voronoi筛选

Processing date: 2026-07-15. Mature criterion used: latest physical time `t>=800`.

- `Ra=1e6`: `Ek=1e-2`, `9e-4`, using `t=710–800`.
- `Ra=1e8`: no rotation and `Ek=1e-2`, `7e-3`, `5e-3`, `2e-3`, `1e-4`, using `t=710–800`.
- `Ra=1e8, Ek=9e-4`: continuation latest ten fields, `t=1510–1600`.
- `Ra=1e8, Ek=3e-3`: only reached `t=220` at this check, so excluded from the mature set.

Important provenance: the remote movie files for `Ra=1e8, Ek=1e-4, AR10` were reinitialized immediately after the first reduction, but the reduced latest-ten-field NPZ was successfully extracted before the movie directory was cleared.

## Current AR10 MSE-rms and scale time series

Processing snapshot: 2026-07-15.

- `Ra=1e8, Ek=1e-4`: online MSE to `t=296.0`, complete fields to `t=290`.
- `Ra=1e8, Ek=3e-3`: online MSE to `t=553.8`, complete fields to `t=540`.
- `Ra=1e8`, no rotation: online MSE to `t=483.5`, complete fields to `t=470`.

These were the only AR10 cases with `run/data/mse_aggregation.out` at this
check. The current local result folder is

`E:/moist RB/moist/result/spectrum/ns/transition_study/RaT1e8pr07/beta1_highres/comparison_outputs/organized_results_20260705/aspect_ratio_study/AR10_mrms_scales_timeseries_20260715`.

Remote reduced arrays are at

`/share/org/SHUTUANL/shu_zhangjs/rainy model/ns/transition_study/beta1/lowrestest/aspect_ratio_study/remote_profile_exports/ar10_mrms_scales_timeseries_20260715`.

## Updated remote-only AR10 MSE variance and aggregation scales

Update snapshot: 2026-07-15. This is the authoritative update; no local field
or local `mprime` value is used.

Current output:

`E:/moist RB/moist/result/spectrum/ns/transition_study/RaT1e8pr07/beta1_highres/comparison_outputs/organized_results_20260705/aspect_ratio_study/AR10_REMOTE_MSE方差与聚集尺度更新_AR10_remote_MSE_variance_scales_updated_20260715`.

Remote provenance root:

`/share/org/SHUTUANL/shu_zhangjs/rainy model/ns/transition_study/beta1/lowrestest/aspect_ratio_study` on host `c01n0037`.

The online `mse_aggregation.out` histories for `Ra=1e8, AR=10` now contain
8000 samples and reach `t=800` for `Ek=1e-4`, `Ek=3e-3`, and the nonrotating
case. Their current endpoint `A_m` values are `0.0062932`, `0.0101872`, and
`0.0125212`, respectively.

The full-field aggregation diagnostics were recomputed remotely from all 80
stored snapshots (`t=10,20,...,800`) for all three cases. The authoritative
reduced arrays are in
`remote_profile_exports/ar10_mrms_scales_timeseries_to_t800_20260715`; no local
field or local `mprime` value is used. Endpoint values at `t=800` are:

- `Ek=1e-4`: `L_peak/H=0.90909`, `L_integral/H=0.097556`;
- `Ek=3e-3`: `L_peak/H=10`, `L_integral/H=0.914224`;
- nonrotating: `L_peak/H=10`, `L_integral/H=1.185719`.

Late-time means over the last ten snapshots (`t=710--800`) are:

- `Ek=1e-4`: `<L_peak>/H=0.98990`, `<L_integral>/H=0.096289`;
- `Ek=3e-3`: `<L_peak>/H=10`, `<L_integral>/H=1.169994`;
- nonrotating: `<L_peak>/H=9`, `<L_integral>/H=0.940030`.

Remote online/full-field `A_m` agreement at matching times is better than
`5e-11`. The update also retains periodic connected-component mean and maximum
radii for `m'_2D>0` and `m'_2D>sigma_m`; Voronoi is not a primary scale.

## Current AR10 nonrotating 256x256x65 aggregation run

Status snapshot: 2026-07-15 23:49 CST. The remote case
`Ra1e8/AR10/norotating/res25625665` was still running and had reached
approximately `t=253.9`. Twenty-five complete full-field snapshots cover
`t=10,20,...,250`; the reduced scale result deliberately freezes this complete
field set rather than mixing in a newly written field.

Remote reduced output:

`/share/org/SHUTUANL/shu_zhangjs/rainy model/ns/transition_study/beta1/lowrestest/aspect_ratio_study/remote_profile_exports/ar10_res256_norotating_scales_current_20260715`.

Local output:

`E:/moist RB/moist/result/spectrum/ns/transition_study/RaT1e8pr07/beta1_highres/comparison_outputs/organized_results_20260705/aspect_ratio_study/AR10_256不旋转MSE方差与聚集尺度_AR10_res256_norotating_MSE_scales_current_20260715`.

At `t=250`, `L_peak/H=10`, `L_integral/H=0.979535`, the maximum equivalent
radius for `m'_2D>0` is `3.92214 H`, and the maximum equivalent radius for
`m'_2D>sigma_m` is `1.27881 H`. The last-ten-field means are
`<L_peak>/H=5.83333` and `<L_integral>/H=0.781193`. Because `L_integral` is
still rising and `L_peak` only reached the box mode near `t=240`, this run is
not yet a converged aggregation-scale plateau.

## 数据更新协议

1. Re-scan remote physical times; do not rely on this inventory as permanently current.
2. Merge continuation directories by physical time and deduplicate.
3. Update the phase diagram whenever completed cases change.
4. Rebuild all related summary plots from master reduced tables rather than appending points manually.
5. Store only reduced data locally and preserve exact remote source paths.
6. If a definition changes, update `diagnostics_and_data_sources.md` and the relevant method reference in the same turn.

## Ra=1e8 nonrotating betaqs study

Prepared on 2026-07-16; no jobs were submitted. Remote root:

`/share/org/SHUTUANL/shu_zhangjs/rainy model/ns/transition_study/beta_not1/Ra1e8/AR10/norotating`.

Prepared cases are `Beta1p05`, `Beta1p10`, `Beta1p20`, `Beta1p30`,
`Beta1p40`, and `Beta1p50`, corresponding to `betaqs=1.05, 1.10, 1.20,
1.30, 1.40, 1.50`. Every case has `Ra=1e8`, `Pr=0.7`, `invRo=0`,
`129x129x65`, `REXT1=REXT2=10`, `DT=DTMAX=7e-3`, `TPIN=0.1`, `TMAX=800`,
bottom no-slip and top free-slip. `gamma=1.1` is unchanged; only `betaqs`
varies.

Each `run` contains the same executable SHA-256
`62a553db4599b142da91cdbea9ae84c9fdfed7b32f2548e23c68223a3ac35081`
and writes the online z/H=0.75 vertical-velocity and MSE aggregation-scale
histories documented in `diagnostics_and_data_sources.md`. The preparation
script and README are retained locally in
`comparison_outputs/program_updates/beta_not1_online_length_20260716` and at
the remote study root.

### Beta1p20 online-history status and reduced plots

Verified on 2026-07-24. The primary `Beta1p20` run is complete from
`t=0.1` to `t=800.0` with 8,000 online samples on the `129x129x65` grid.
The direct sources are:

- `Beta1p20/run/data/mse_aggregation.out`;
- `Beta1p20/run/data/mse_aggregation_scales.out`;
- `Beta1p20/run/data/w_z075_spectral_length.out`.

Do not merge the separate
`Beta1p20/384384192/run/data` history into this curve. It is a distinct
`385x385x193` calculation and had only reached `t=2.4` at the scan time.

The synchronized reduced tables, plotting program, PNG figures, and PDF
figures are stored locally under:

`comparison_outputs/organized_results_20260705/beta_not1/Beta1p20_mprime_scales_timeseries_20260724`.

The primary figures use the raw online samples without temporal smoothing.
They include `A_m`, `sqrt(A_m)`, `L_peak`, `L_integral`, `2*pi*ell_m`,
`2*pi*ell_w(z=0.75)`, and the two threshold-based cluster-radius histories.

The same folder also contains the Beta1p20 kinetic-energy histories generated
from `avgvar.out`.  The definition is
`K(t)=0.5*avgvar_column_4(t)`, with duplicate physical times removed and
non-finite/uninitialized rows filtered.  Current files include the full
primary `129x129x65` curve through `t=800`, a `t>=100` zoom, and a separate
early-time `385x385x193` high-resolution curve through `t=8.3`.

## Generic rotating-case local builder

Created on 2026-07-25 and upgraded with case-specific drizzle initialization
on 2026-07-28. Current local tool folder:

`E:\moist RB\rotating_case_inventory\01_case_builder_remote_upload`

Main entry points:

- `case_config.json`: default editable parameter file;
- `batch_cases.json`: example batch parameter file with list-valued sweeps;
- `create_rotating_case.py`: prepares one remote case;
- `create_case.ps1`: VS Code / PowerShell wrapper;
- `create_cases_batch.ps1`: VS Code / PowerShell wrapper for batch creation;
- `generate_drizzle_initial_condition.py`: calls the supplied
  `moist_base_state`, reads the current `bou.in`, and rejects
  nonconverged/linear-fallback results;
- `prepare_drizzle_initial_condition.sh`: manual one-command profile
  generation to run before each new `nread=0` calculation;
- `check_drizzle_before_submit.py`: checks profile parameters and endpoints;
- `README.md`: user instructions.

Default remote root:

`/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case`

Generated directory hierarchy:

`Pr/Ra/Ek-or-norotating/AR/Beta/qbot-qtop/grid/run`.

Example:

`rotating_case/Pr0p7/Ra1e8/Ek9e-4/AR4/Beta1p30/qbot1_qtop0p004978/N257x257x129/run`.

Parameter mapping:

- `Ek=norotating` gives `invRo=0`; otherwise `invRo=sqrt(Pr/Ra)/Ek`.
- `alpha` in the local config maps to solver `alphaqs`; solver `alpha` remains
  configurable as `control_alpha` and defaults to zero.
- `beta` maps to `betaqs` and also sets the buoyancy boundary row as
  `dsaltop=beta-1` and `dsalbot=0`.
- `aspect_ratio` sets both `REXT1` and `REXT2`.

The tool copies the configured remote template `run`, installs the latest
`simexec`, places the manual drizzle generator and supplied solver in the
case, writes `case_info.json`, and creates a `submit_after_check.sh` helper.
It neither solves drizzle silently nor runs `csub` automatically. Before
every new `nread=0` calculation, the user runs
`./prepare_drizzle_initial_condition.sh` inside that case's `run` directory;
this rereads the current `bou.in` and overwrites the profile.

At `nread=0`, the current DNS reads the drizzle `b_0(z),q_0(z)`, adds a
`1e-4*sin(pi*z/H)` Gaussian perturbation to buoyancy only, keeps drizzle
moisture unperturbed, and initializes velocity to zero. The submitted linear
stability solver source is retained unchanged under
`linear_stability_reference/stability_solver.py`.

The compiled current executable SHA-256 is
`08527afd2ebb084c84a5e47b1d6f6c63e18f180a439894868b7a35220d7c4d77`.
On 2026-07-28, all 39 existing cases received the manual generator, solver,
submission check, and current executable with zero update failures. No job
was submitted. Audit:

`/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case/manual_drizzle_workflow_update_report_20260728.json`.

Batch mode reads `batch_cases.json`, merges `defaults` with each entry in
`cases`, and expands list-valued entries as Cartesian products. The same
post-write boundary check is applied to every case.

Critical boundary rule: the builder must always write the `bou.in` spatial
modulation row

`A_stopmod k_stopmod A_sbotmod k_sbotmod dsaltop dsalbot qvaptop qvapbot`

with `dsaltop=beta-1`, `dsalbot=0`, `qvaptop=<config qvaptop>`, and
`qvapbot=<config qvapbot>`. After writing, it reads the row back and aborts if
any of these four boundary values differ from the requested values.

Prepared rotating-case sweep on 2026-07-25:

- root:
  `/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case/Pr0p7/Ra8e6`;
- parameters: `Ra=8e6`, `Pr=0.7`, `beta=1.02`, `AR=16`,
  grid `257x257x65`, `qvaptop=0.004978`, `qvapbot=1.0`;
- Ek cases prepared but not submitted:
  `7e-4`, `1e-3`, `1.5e-3`, `2e-3`, `3e-3`, `4e-3`, `5e-3`,
  `7e-3`, and `1e-2`;
- each case passed read-back checks for `invRo=sqrt(Pr/Ra)/Ek`,
  `dsaltop=0.02`, `dsalbot=0`, `qvaptop=0.004978`, and `qvapbot=1.0`.

## Local field diagnostics generator

Created on 2026-07-26. Local tool folder:

`E:/moist RB/field_diagnostics_generator`

Root launcher:

`E:/moist RB/generate_ra8e6_field_diagnostics.ps1`

Default data root:

`E:/moist RB/moist/result/data/ns/transition_study/rotating_self_aggregation/Ra8e6`

The tool scans every `fieldXXXXX.h5` and idempotently generates missing
ParaView-readable pairs:

- `mprimeXXXXX.h5/.xmf`, dataset `MPI_M`;
- `vortzXXXXX.h5/.xmf`, dataset `VORTZ_me`.

Definitions:

- `m = DSAL_me + gamma * QVAP_me`, default `gamma=1.1`;
- `mprime = m - mean_xy(m)` at each height and time;
- `VORTZ_me = d(VY_me)/dx - d(VX_me)/dy`, using horizontal spectral
  derivatives and periodic x/y directions.

The first Ra8e6 local pass generated 34 complete output pairs from 17 source
fields with no errors. Subsequent runs detected all 34 existing pairs and
generated nothing, confirming that newly downloaded fields can be added by
rerunning the same launcher.

## Latest-program Ra8e6 stable MSE-scale selection (2026-08-03 13:16)

Current source:

`E:/moist RB/rotating_case_inventory/04_outputs_and_figures/high_resolution_timeseries_latest_program_20260803/high_resolution_timeseries_long.csv`

Strict stable-only `2*pi*Lm` versus `Ek` figure:

`E:/moist RB/rotating_case_inventory/04_outputs_and_figures/high_resolution_timeseries_latest_program_20260803/lm_stable_scaling/Ra8e6_stable_only_lm_vs_ek_loglog.png`

Included points and averaging windows:

- `Ek=1e-3`: newly stable, latest 500 time units;
- `Ek=3e-3,5e-3,7e-3`: established stable, latest 500 time units;
- `Ek=1e-2`: statistical steady, latest 1000 time units to sample its
  low-frequency fluctuations.

The near-steady `Ek=3e-2` point is excluded because its `2*pi*Lm` still
decreased by about 12% over the latest 500 time units. The nonrotating,
`5e-4`, `7e-4`, `5e-2`, and `1e-1` cases are also excluded because at least
one of kinetic energy, MSE variance, or `2*pi*Lm` remains nonstationary. The
current five-point exploratory fit is

`2*pi*Lm = 65.6887 Ek^0.623046`, with `R^2=0.983442`.

## Latest-program Ra8e6 time-series and scheduler status (2026-08-03 13:44)

Updated figures and reduced histories:

`E:/moist RB/rotating_case_inventory/04_outputs_and_figures/high_resolution_timeseries_latest_program_20260803`

The update contains 13 cases with the latest four diagnostic files. It
merges base/continuation rows by physical time and retains the historical
reduced-data patch for the accidentally deleted `Ek=7e-3`, `t=10--1000`
interval. Current maximum times are: NR `844.7`; `Ek=1.5e-4` `1200`;
`2e-4` `1200`; `5e-4` `800`; `7e-4` `800`; `1e-3` `846.3`; `3e-3`
`2400`; `5e-3` `2400`; `7e-3` `2200`; `1e-2` `2390`; `3e-2` `2400`;
`5e-2` `959.203`; and `1e-1` `1194.102`.

At 13:44 the scheduler retained only one compute job: job `74474575`,
`Ek=5e-2`, which was still actively writing diagnostics. Jobs `74444575`
(`Ek=1e-1`), `74474581` (NR), and `74474650` (`Ek=1e-3`) had stopped
writing for roughly 11--12.5 hours and were then killed by the
user/administrator account at 13:44, exiting by signal 13. They must not be
reported as normally running or normally completed.

Durable status table:

`E:/moist RB/rotating_case_inventory/04_outputs_and_figures/high_resolution_timeseries_latest_program_20260803/latest_program_case_status_20260803.csv`

## Intermittent-plume low-valley MSE scale (2026-08-03)

For the Ra8e6 latest-program cases, `Ek=2e-4`, `5e-4`, and `7e-4`
show a well-separated intermittent high/low `2*pi*Lm` cycle. Low valleys are
selected from a 10-time-unit smoothed history by a two-cluster split after
removing the initial pre-burst low segment. To prevent the falling and
regrowing shoulders from contaminating the plume-scale mean, accepted
samples must additionally satisfy `abs(d(2*pi*Lm)/dt)<=3e-3`, and each
continuous flat interval must last at least 20 time units. Raw online
`2*pi*Lm` samples inside these flat intervals are time averaged. The
`Ek=1.5e-4` case is now excluded from both the plot and all fits.

Selected means:

- `Ek=2e-4`: flat intervals `353.3--432.6`, `622.9--732.1`, and
  `939.4--1029.8`; mean `0.351092`;
- `Ek=5e-4`: flat intervals `390.5--429.0` and `451.8--488.5`; mean
  `0.620867`;
- `Ek=7e-4`: flat intervals `447.9--469.6`, `491.1--523.1`, and
  `580.6--684.6`; mean `0.733414`.

Do not fit one power law across the small-Ek intermittent points and the
later large-vortex points. A free log-log fit using only the three flat
intermittent-valley means gives

`2*pi*Lm = 56.1538 Ek^0.595217`, with `R^2=0.997962`.

For the five large-vortex points at
`Ek=1e-3--1e-2`, keep the previously selected exponent fixed at `0.69`, fit
only the prefactor, and draw a short parallel dashed guide beside rather
than through the markers. The fixed-slope comparison is

`2*pi*Lm = 95.0377 Ek^0.69`, with log-space `R^2=0.972086` for the five
large-vortex points only. The displayed dashed line is shifted upward for
visual separation and must not be interpreted as a connecting curve.

Outputs:

`E:/moist RB/rotating_case_inventory/04_outputs_and_figures/high_resolution_timeseries_latest_program_20260803/lm_low_valley_scaling`

Update at 2026-08-03 22:54: the restarted latest-program histories reached
NR `t=1130.6`, `Ek=1e-3` `t=1096.9`, `Ek=5e-2` `t=1200.997`, and
`Ek=1e-1` `t=1403.003`. The intermittent flat-valley means were unchanged.
Updating the `Ek=1e-3` late-500 mean from `0.851706` to `0.868044` changes
the fixed large-vortex comparison to
`2*pi*Lm = 95.3996 Ek^0.69`, with log-space `R^2=0.969484`.
# Latest-program Nu and force balance (2026-08-05)

- Local output root:
  `E:\moist RB\rotating_case_inventory\04_outputs_and_figures\high_resolution_timeseries_latest_program_20260805`
- All-case normalized Nu time series:
  `num_timeseries\Ra8e6_latest_program_all_cases_Num_timeseries.png`
- Normalized Nu source CSV:
  `num_timeseries\Ra8e6_latest_program_Num_timeseries_normalized.csv`
- Reduced force profiles and bulk summary:
  `force_balance\force_balance_profiles.csv` and
  `force_balance\force_balance_bulk_summary.csv`
- Per-case classifications and explicit pressure caveat:
  `force_balance\force_balance_case_classification.csv`
- Remote reducer:
  `/share/org/SHUTUANL/shu_zhangjs/postprocess_force_balance_20260805`

## Ra8e6 migration and nonrotating provenance correction (2026-08-08)

- The deployable restart bundle is
  `C:\Users\jiasenzhang\Desktop\修改程序\rotating_moist_migration_bundle_20260808`.
- Do **not** use
  `G:\moist convection\norotating\AR16\Beta1p02\qbot0p5_qtop0p004978\N257x257x65\conti1\conti_strict_force_500\run`
  as a nonrotating source. It is a misdirected local download whose
  `bou.in`, AR4 grid, and 33,819,136-byte HDF5 files match `Ek=2e-4`.
- The corrected nonrotating restart was downloaded directly from
  `/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case/Pr0p7/Ra8e6/norotating/AR16/Beta1p02/qbot0p5_qtop0p004978/N257x257x65/conti1/run`.
- Its audited local active copy is
  `C:\Users\jiasenzhang\Desktop\修改程序\rotating_moist_migration_bundle_20260808\03_current_cases\norotating\run`:
  effective grid `385x385x65`, AR16, `invRo=0`, `NREAD=1`, and five
  75,893,248-byte continuation HDF5 files.
- The separate raw recovery copy is
  `C:\Users\jiasenzhang\Desktop\修改程序\remote_recovery_norotating_20260808`.

## xh5 continuation deployment (2026-08-15)

- Remote bundle:
  `/work/home/jiasenzhang/rotating_moist_migration_bundle_20260815`.
- Fourteen cases are present: thirteen rotating Ek cases from `1e-1` through
  `1.5e-4`, plus `norotating`.
- Every case has `NREAD=1`, `TMAX=500d0`, five validated continuation HDF5
  files, and a Slurm `subjob.sh` for partition `xhacnormalc`, one node, and 64 ranks.
- The xh5-compiled binary SHA-256 is
  `f4504bda4e9f7524664fb615f61a879a0d8eb5fa5c2d24f1abfd4b263ee8b371`.
- The source hashes match the current local source under
  `E:\moist RB\rotating_case_inventory\00_latest_program\source`.
- Upload and validation are complete. No continuation job was submitted as
  part of this deployment.
