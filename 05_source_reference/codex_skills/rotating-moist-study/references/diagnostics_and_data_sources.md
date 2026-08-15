# Diagnostics and data sources

This is the living definition sheet for the rotating moist Rayleigh-Benard project. Update it whenever a diagnostic definition, averaging window, data source, or plotting convention changes.

## Global constants and boundary values

- Moist static energy:

\[
m=b+\gamma q,\qquad \gamma=1.1 .
\]

- Current moisture boundary convention:
  - bottom: \(q_{\rm bot}=q_s=1\);
  - top: \(q_{\rm top}=0.1q_s=0.004978\).

- Moist contrast used in reduced parameters and \(Nu_m\):

\[
\Delta_m=\gamma(q_{\rm bot}-q_{\rm top}) .
\]

With the current values, \(\Delta_m\approx1.0945\). If the boundary condition changes, update this number before computing \(Ra_m\), reduced \(Ra\), or \(Nu_m\).

- Default Prandtl number in the current transition data: \(Pr=0.7\).

## Remote data roots

Use remote data directly whenever possible; do not pull full HDF5 field files back to local storage.

- Main remote host: `c01n0034`.
- Current accessible host for the 2026-07-11 resolution check: `c01n0037`.
- High-resolution / original `Ra=1e8` style cases:

`/share/org/SHUTUANL/shu_zhangjs/rainy model/ns/transition_study/beta1`

- Low-resolution sweep cases:

`/share/org/SHUTUANL/shu_zhangjs/rainy model/ns/transition_study/beta1/lowrestest`

- Preferred remote-reduced small-output folder:

`/share/org/SHUTUANL/shu_zhangjs/rainy model/ns/transition_study/beta1/remote_profile_exports`

- Local working output folder:

`E:\moist RB\moist\result\spectrum\ns\transition_study\RaT1e8pr07\beta1_highres\comparison_outputs`

## Remote postprocessing script index

Remote postprocessing/helper scripts under the `beta1` tree are organized on `c01n0020` at:

`/share/org/SHUTUANL/shu_zhangjs/rainy model/ns/transition_study/beta1/postprocess_scripts_index`

Important files:

- `README.md`: organization notes;
- `MANIFEST.tsv`: category, script name, original path, local-folder copy, and central indexed copy;
- `all_scripts_by_source/<category>/`: central copies grouped by diagnostic category.

Additionally, each original script directory has a `postprocess_scripts/` folder containing copies of the scripts found in that directory. Original scripts were not moved or deleted.

Current categories:

- `transport_num`;
- `re_rossby`;
- `profiles_lengths_rh`;
- `ekman_corr_boundary`;
- `vortex_structure`;
- `plotting_visualization`;
- `case_submission_helpers`.

When looking for or reusing a remote postprocessing script, check `postprocess_scripts_index/MANIFEST.tsv` first.

## Averaging rules

- Prefer the latest completed late-time field frames.
- For field-based diagnostics, the common default is the latest 10 movie frames unless the user asks for another window.
- Record `nframes`, `field_first`, `field_last`, `nz`, `ny`, `nx`, and `remote_path` in every reduced CSV.
- For `avgvar.out` time series diagnostics, merge continuation outputs in chronological order and remove duplicated physical times when necessary.

## Phase diagram

The reduced moist Rayleigh number is

\[
Ra_m=Ra\,\Delta_m,\qquad R_m=Ra_m Ek^{4/3}.
\]

Use `R_m` and `Ek` for the reduced phase diagram unless the user explicitly asks for a different phase space.

Data source:

- completed-case CSVs in `comparison_outputs`;
- remote case status should be refreshed from `c01n0034` before adding new points.

Rules:

- The phase diagram is a living summary. Refresh it first when new data finish.
- Add only cases stable enough for late-time averaging.
- Keep regime labels synchronized with `moist_rb_core_refs.md`.

## Moist static energy profiles

\[
m(z)=\langle b+\gamma q\rangle_{x,y,t}.
\]

Field variables:

- `DSAL_me` for \(b\);
- `QVAP_me` for \(q\).

Default data source:

- remote HDF5 movie fields, reduced remotely into small CSV/NPZ files.

Important:

- Do not use old local single-frame exports when the user asks for averaged profiles.
- Profile figures should use one color per case family.

## Moist perturbation field

\[
m'(x,y,z,t)=m(x,y,z,t)-\langle m\rangle_{x,y}(z,t).
\]

Use this for moist-structure visualization and moist-structure length scales. It is not the same diagnostic as the velocity-based \(l_h\).

Generated field products:

- local `mprime` HDF5/XMF files are derived from existing local `field*.h5` files;
- highres and lowres local field folders should both be processed when the user asks for all cases.
- use the fixed model value `gamma=1.1`, read `b=DSAL_me` and `q=QVAP_me`,
  and store the result as `MPI_M` in `mprimeNNNNN.h5`; the matching
  `mprimeNNNNN.xmf` must reference `:/MPI_M`;
- calculate the horizontal mean independently for every height and every
  snapshot.  Do not subtract a volume mean, a time mean, or a single mean
  shared by all heights;
- the maintained local single-frame generator is
  `E:\moist RB\moist\result\spectrum\ns\transition_study\RaT1e8pr07\beta1_highres\generate_mprime_h5_xmf.py`;
  it reads in vertical blocks and writes through temporary files so large-AR
  cases do not exhaust memory or leave false complete outputs after an
  interrupted run;
- the local recursive batch driver and independent audit are
  `batch_generate_all_local_mprime.py` and `audit_all_local_mprime.py` in the
  same directory.  Batch outputs and provenance tables are stored under
  `E:\moist RB\moist\result\spectrum\ns\transition_study\_mprime_batch_reports`.

Local all-case generation audit on 2026-07-15:

- discovered 398 local `field*.h5` sources in 68 directories;
- 397 HDF5/XMF `mprime` pairs are complete after generating 361 missing pairs,
  retaining 35 pre-existing pairs, and recovering one pair from a
  byte-identical complete local copy;
- all 396 directly readable source fields passed the formula check exactly at
  sampled heights, and all available outputs passed the all-height
  zero-horizontal-mean and XMF-reference checks;
- the only unavailable product is for
  `RaT1e8pr07\beta1_highres\3e-4\field00073.h5`, whose local source is
  physically truncated (`119013376` bytes versus the expected `671094784`)
  and has no exact local duplicate.  Do not synthesize this frame from a
  neighboring time; restore the original source before generating it.

Additional local generation on 2026-07-24:

- generated 79 `mprime*.h5`/`mprime*.xmf` pairs for
  `E:\moist RB\moist\result\spectrum\ns\transition_study\lowres\beta1.2\Ra1e8\norotating`;
- source files were the existing local `field*.h5` products in that folder;
  `field00071.h5` was absent locally, so no `mprime00071` product was
  synthesized;
- definition unchanged: `MPI_M = b + 1.1 q - <b + 1.1 q>_{xy}(z,t)`;
- audit samples showed `MPI_M` shape `(64,128,128)` and layer-mean residuals
  at roughly `1e-11`.

## \(Nu_m\), \(Nu_q\), and \(Nu_b\)

## Strict rotating force-balance diagnostics

The maintained pressure-enabled `simexec` writes the pressure movie as
`movie/pressureNNNNN.h5:/PR_me` with matching
`movie/pressureNNNNN.xmf`. The pressure XMF uses the movie/base-grid
dimensions `n1m,n2m,n3m` and the coordinates in
`movie/cordin_info.h5:/x`, `:/y`, and `:/z`.

The strict offline diagnostic is:

`E:\moist RB\rotating_case_inventory\02_inventory_and_plot_scripts\strict_force_balance_from_pressure_movies.py`

and the synchronized remote copy is under:

`/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case/postprocess/strict_force_balance_20260805/`

The continuation staging command is:

`/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case/postprocess/strict_force_balance_20260805/prepare_strict_force_continuations.sh`

It only creates missing `conti_strict_force_500/run` directories and never
calls `csub`. Existing continuation directories are protected from overwrite.

For each snapshot and each height, it evaluates the complete nondimensional
momentum terms using the coefficients read from `bou.in`:

\[
\mathbf F_I=-(\mathbf u\cdot\nabla)\mathbf u,\quad
\mathbf F_C=Ro_C^{-1}(v,-u,0),\quad
\mathbf F_P=-\nabla p,
\]
\[
\mathbf F_V=\sqrt{Pr/Ra}\nabla^2\mathbf u,\quad
\mathbf F_B=(0,0,b-\langle b\rangle_{xy}),\quad
\mathbf F_T=\partial_t\mathbf u.
\]

The vertical pressure gradient `-∂z p` is included. Every component is
demeaned over the horizontal plane separately at each `(z,t)` before
forming total, horizontal, and vertical RMS values. The diagnostic also
outputs the horizontal geostrophic residual `R_G`, the vertical momentum
residual `R_z`, the full momentum closure residual, and the ratio of that
closure residual to `F_T`.

For the current force-versus-`Ek` comparison and the intermittent-case force
time series, the green Coriolis-pressure cancellation curve is **not** `R_G`.
Use the complete three-dimensional vector residual requested by the user,

\[
F_{CP}(z,t)=\left|\mathbf F_C+\mathbf F_P\right|_{\rm rms}
=\sqrt{R_G^2(z,t)+F_{P,z}^2(z,t)},
\]

where all force components have first been horizontally demeaned at each
`(z,t)` and `F_{C,z}=0`. Form this plane RMS before averaging over the chosen
bulk-height interval. Do not reconstruct it from already height-averaged
`R_G` and `F_P_z`; use `force_balance_z.out`, the strict profile table, or
synchronized velocity/pressure fields. Label it
`|\mathbf F_C+\mathbf F_P|`, while reserving "horizontal geostrophic
residual" for `R_G=|\mathbf F_{C,h}+\mathbf F_{P,h}|`.

Output files:

- `strict_force_balance_zt.npz`: all `z-t` arrays and metadata;
- `strict_force_balance_bulk_timeseries.csv`: wide bulk time series;
- `strict_force_balance_timeseries_long.csv`: explicit one-row-per-time,
  one-quantity force/residual time series;
- `strict_force_balance_profiles_long.csv`: long-form `z-t` table.

The default bulk averaging interval is `0.1<z<0.9`; record any changed
interval in the provenance table. The previous
`remote_force_balance_from_movies.py` does not include the complete pressure
gradient or closure residual and must not be used as the strict
geostrophic/CIA判据.

### Online DNS force time series (2026-08-06)

The maintained remote source now includes
`fluid_solver/force_balance_online.f90` and writes
`run/data/force_balance.out` at the same statistics cadence as `avgvar.out`.
It uses the DNS staggered-grid momentum stencils and actual coefficients for
`F_I`, `F_C`, `F_P`, `F_V`, and `F_B`; the configured drag coefficient is
zero and no drag diagnostic is written. `F_P` includes all three pressure
gradients. `F_T` is a finite difference between successive statistics
outputs. The first call initializes the previous velocity and does not write
a row.

Each force component is horizontally demeaned at every height before total,
horizontal, and vertical RMS values are accumulated over `0.1H<=z<=0.9H`.
The file also contains `R_G`, `R_z`, and the complete momentum closure
residual. Exact columns, source hashes, build instructions, originals, and the
compiled local backup are in:

`E:\moist RB\rotating_case_inventory\00_latest_program\force_balance_online_modified_20260806\README_force_balance_online.md`

The compiled remote executable is
`/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case/latest_program/source/simexec`.
It is not automatically copied into existing or running cases. New
continuations must explicitly use this latest executable before submission.

### Current Re and moist Nusselt reduction (2026-08-05)

Remote source:

`/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case/postprocess/re_num_20260805/`

The reducer is:

`E:\moist RB\rotating_case_inventory\02_inventory_and_plot_scripts\reduce_re_num_remote.py`

It merges the base run and all continuation segments for each actual grid,
with later physical-time segments replacing duplicate times. The definitions
are:

\[
Re(t)=
\sqrt{\left\langle u^2+v^2+w^2\right\rangle_V}
\sqrt{\frac{Ra}{Pr}},
\]

\[
Nu_m(t)=
\frac{1}{\Delta_m}
\int_0^1 Nu_m(z,t)\,dz,
\qquad
\Delta_m=(b_{\rm bot}-b_{\rm top})
\gamma(q_{\rm bot}-q_{\rm top}).
\]

For the current \(Ra=8\times10^6\), \(Pr=0.7\), \(\beta=1.02\),
`qbot=0.5`, `qtop=0.004978` cases,
\(\Delta_m=0.5245242\). `Re` is obtained from `avgvar.out` column 4
(one-based) and `Nu_m(z,t)` from column 4 of `data/nu_profiles.out`.

Stability screening uses 10-time-unit block averages over the latest
200 physical-time units. A case is included in an Ek scaling figure only if
both Re and `Nu_m` have absolute linear-fit drift and first-half/second-half
shift no larger than 10 percent.

Local output:

`E:\moist RB\rotating_case_inventory\04_outputs_and_figures\re_num_20260805\`

The current project `Nu_m` is based on moist static energy \(m=b+\gamma q\).

Volume/integral transport form:

\[
Nu_m
=1+\frac{\sqrt{RaPr}}{\Delta_m}
\int_0^1 \langle w b\rangle_{x,y,t}
       +\gamma\langle w q\rangle_{x,y,t}\, dz .
\]

Profile form:

\[
Nu_m(z)=
\frac{
\sqrt{RaPr}\,\langle w(b+\gamma q)\rangle_{x,y,t}
-\partial_z\langle b+\gamma q\rangle_{x,y,t}
}{\Delta_m}.
\]

The bulk value of \(Nu_m(z)\) should be approximately height-independent only if the definition and averaging are consistent.

Related quantities:

\[
Nu_q=1+\frac{\sqrt{RaPr}}{\Delta_q}
\int_0^1 \langle wq\rangle_{x,y,t}\,dz ,
\]

\[
Nu_b \text{ is based on } \langle wb\rangle
\text{ and the corresponding } b \text{ contrast.}
\]

For the current setup, verify the \(b\)-contrast before using \(Nu_b\) as a wall-normalized Nusselt number. Do not rename \(Nu_b\), \(Nu_q\), and \(Nu_m\) interchangeably.

Data sources and current output tables:

- `comparison_outputs/remote_lowrestest_num_master.csv`
- `comparison_outputs/remote_lowrestest_num_summary.csv`
- `comparison_outputs/Ra1e6_num_summary.csv`
- `comparison_outputs/Ra5e6_num_summary.csv`
- `comparison_outputs/Ra5e7_reduced/Ra5e7_num_summary.csv`
- `comparison_outputs/Ra1e8_lowrestest_num_summary.csv`

Field variables:

- `VZ_me`, `DSAL_me`, `QVAP_me`.

Default reduction:

- compute profiles and fluxes remotely from late-time movie fields;
- save small CSVs locally;
- normalize by same-Ra nonrotating \(Nu_{m,0}\) only for plots explicitly labeled \(Nu_m/Nu_{m,0}\).

## Reynolds numbers

Global Reynolds number:

\[
Re_{\rm global}=u_{\rm rms}\sqrt{\frac{Ra}{Pr}} .
\]

The exact `avgvar.out` column used must be checked against the solver output before reuse. In prior scripts, kinetic energy and velocity variance came from `avgvar.out`; do not silently switch columns.

Vertical-velocity Reynolds number:

\[
Re_z(z)=w_{\rm rms}(z)\sqrt{\frac{Ra}{Pr}},
\qquad
w_{\rm rms}(z)=\langle w^2\rangle_{x,y,t}^{1/2}.
\]

Field variable:

- `VZ_me`.

Important:

- The user may ask for full-velocity \(Re\) or vertical-only \(Re_z\). Do not relabel one as the other.
- If normalized by nonrotating data, use same \(Ra\), same diagnostic, and same height.

Current Ra=5e7 and Ra=1e8 all-case global-Re figure:

- definition: full-velocity \(Re=U_{\rm rms}\sqrt{Ra/Pr}\), with
  \(U_{\rm rms}=\langle u^2+v^2+w^2\rangle_{x,y,z,t}^{1/2}\);
- output: `comparison_outputs/organized_results_20260705/Ek_scaling/Ra5e7_Ra1e8_global_Re_vs_Ek_all_cases.png` and `.pdf`;
- audit table: `comparison_outputs/organized_results_20260705/Ek_scaling/Ra5e7_Ra1e8_global_Re_vs_Ek_all_cases.csv`;
- high-resolution source: `comparison_outputs/global_Re_Rez_all_highres_remote/Ra1e8_global_Re_Rez_all_highres_last30_no3e5.csv`, using 30 late fields;
- low-resolution supplements: `comparison_outputs/re_global_remote_summary_allra.csv`, using 10 late fields and only filling Ek values absent from the high-resolution table;
- `Ek=3e-5` supplement: local high-resolution `3e-5/field00038.h5` and `field00080.h5`; it is a two-frame estimate and is less converged than the other points;
- Ra=5e7 source: remote latest 10 full fields for eight cases under `/share/org/SHUTUANL/shu_zhangjs/rainy model/ns/transition_study/beta1`, reduced to `comparison_outputs/organized_results_20260705/Ek_scaling/Ra5e7_global_Re_last10.csv`;
- plot scope: all 17 finite-Ek Ra=1e8 cases plus eight Ra=5e7 cases, no error bars, dashed saturated-blue/red lines with filled circle/square markers.

Current two-Ra normalized Num figures:

- include only `Ra=5e7` and `Ra=1e8`; the former `Ra=5e6` purple series is intentionally excluded;
- Ek figure: `comparison_outputs/num_vs_ek_phase_complete_no_ra1e6.png`, with the legend shown in the upper right;
- convective-Rossby figure: `comparison_outputs/num_vs_roc_phase_complete_no_ra1e6.png`, without a repeated legend;
- active point table: `comparison_outputs/num_phase_complete_points_no_ra1e6.csv`, filtered to these two Ra values.

## Local Rossby number / `Rozl`

The project local Rossby number is a measured diagnostic, not the control convective Rossby number.

\[
Ro_z^\ell(z)
=\frac{Re_z(z)\,Ek}{l_h(z)}.
\]

Use:

- \(Re_z(z)\) from `VZ_me`;
- \(l_h(z)\) from the velocity-based horizontal spectral length at the same height.

Do not call this `Ro_c`.

The control convective Rossby number is

\[
Ro_c=\sqrt{\frac{Ra\,Ek^2}{Pr}} .
\]

## Velocity-based horizontal length \(l_h\)

Default \(l_h\) is velocity-based and uses the vertical velocity field \(w\).

For each frame:

\[
l_h(z,t)=
\frac{\sum_{k_h>0}E_w(k_h,z,t)}
{\sum_{k_h>0}k_hE_w(k_h,z,t)},
\qquad
E_w=|\hat w(k_x,k_y,z,t)|^2 .
\]

Then average lengths in time:

\[
\overline{l_h}(z)=\langle l_h(z,t)\rangle_t.
\]

Default scaling plots use \(z=0.75\):

\[
\overline{l_h}(z=0.75).
\]

This definition was confirmed as the project default on 2026-06-30.

Do not silently use:

- \(1/\langle l_h^{-1}\rangle_t\);
- \(\sum_t\sum E_w/\sum_t\sum k_hE_w\).

Those are useful checks but not the default. The corrected all-case table is:

- `comparison_outputs/lh_length_z075_frame_mean_corrected_all_cases.csv`.

Interpretation:

- This is a spectral velocity scale, suited for comparison with the classical \(Ek^{1/3}\) convective-scale theory.
- It is not the geometric diameter of a funnel or large-scale vortex.

### AR10 midplane vertical-velocity length

When the user asks for the length in the supplied midplane definition, use the
same frame-first spectral moment at \(z=0.5\):

\[
\ell_{w,1/2}(t)=
\frac{\sum_{k_h>0}|\hat w(k_x,k_y,z=0.5,t)|^2}
{\sum_{k_h>0}k_h|\hat w(k_x,k_y,z=0.5,t)|^2}.
\]

For the current even 64-plane output grid, evaluate \(z=0.5\) by linear
interpolation between the two surrounding output planes, rather than silently
choosing the lower or upper plane. Exclude \(k_h=0\), use angular wavenumber,
and do not multiply by an extra \(2\pi\). Compute each frame's length first,
then take the arithmetic mean of the requested frames.

Current AR10 mature-case result (2026-07-16):

- inclusion: latest complete field time \(t\ge800\), at least 10 unique fields;
- averaging: latest 10 physical-time snapshots per case;
- \(Ra=10^8,Ek=9\times10^{-4}\) uses its continuation at \(t=1510--1600\);
- all other included cases use \(t=710--800\);
- remote output:
  `/share/org/SHUTUANL/shu_zhangjs/rainy model/ns/transition_study/beta1/lowrestest/aspect_ratio_study/remote_profile_exports/ar10_midplane_w_length_last10_20260716`;
- local output:
  `comparison_outputs/organized_results_20260705/aspect_ratio_study/AR10涓眰鍨傜洿閫熷害璋遍暱搴?0甯у钩鍧嘷AR10_midplane_w_spectral_length_last10_20260716`.
- current plotting convention: show only the (Ra=10^8) rotating points and
  draw the authoritative current `128x128x65` nonrotating value
  (ell/H=0.0687256) as a horizontal gray dashed reference; do not show the
  (Ra=10^6) points unless the user asks to restore them.

## Moist-structure length \(l_{m'}\)

For visual/moist structure scales use \(m'\), not \(w\):

\[
l_{m'}(z,t)=
\frac{\sum_{k_h>0}E_{m'}(k_h,z,t)}
{\sum_{k_h>0}k_hE_{m'}(k_h,z,t)},
\qquad
E_{m'}=|\widehat{m'}|^2 .
\]

Default time average should match \(l_h\):

\[
\overline{l_{m'}}(z)=\langle l_{m'}(z,t)\rangle_t.
\]

Earlier exploratory `l_mprime` profile files may include an energy-weighted time-spectrum version; regenerate with the frame-first definition before using in final figures.

Use \(l_{m'}\) when comparing to RH or \(m'\) visualizations. Keep the notation distinct from \(l_h\).

Current plotting shorthand:

\[
l_{hm}\equiv l_{m'}
\]

when the figure is explicitly labeled as the moist-structure horizontal length. Use \(l_{hm}\) for \(m'\)-based Ekman-number scaling plots and reserve \(l_h\) for the velocity-based \(w\) scale.

Current Ra=1e8 output:

- `comparison_outputs/lhm_mprime/lhm_mprime_Ra1e8_z075_frame_mean_clean.csv`
- remote source: `c01n0020`, under the `beta1` and `lowrestest/Ra1e8` trees;
- height: `z_target=0.75`;
- averaging: frame-first, latest 10 readable movie fields per case;
- columns include both `lhm_mean` and `lambda_hm_mean=2*pi*lhm_mean`.

Current Ra=1e8 profile comparison for the transitional Ek range:

- local output folder: `comparison_outputs/lhm_mprime/profiles_Ra1e8_Ek1e-4_to_1e-3_last30`;
- CSV: `Ra1e8_Ek1e-4_to_1e-3_lh_lhm_profiles_last30_frame_mean.csv`;
- summary at \(z\approx0.75\): `Ra1e8_Ek1e-4_to_1e-3_z075_summary.csv`;
- remote source: `c01n0020`, under `/share/org/SHUTUANL/shu_zhangjs/rainy model/ns/transition_study/beta1`;
- included Ek: \(10^{-4},2\times10^{-4},3\times10^{-4},4\times10^{-4},5\times10^{-4},6\times10^{-4},7\times10^{-4},8\times10^{-4},9\times10^{-4},10^{-3}\);
- averaging: latest 30 readable movie frames, frame-first length averaging;
- diagnostics: velocity-based \(l_h(z)\), moist-structure \(l_{hm}(z)=l_{m'}(z)\), and \(2\pi l_{hm}(z)\);
- high-resolution cases use the original `Ek*/RaTwithkh` tree; low-resolution intermediate cases use `lowrestest/Ra1e8/Ra1e8Pr07Ek*`.

## Correlation diagnostics

The Ekman-pumping-style correlation is a Pearson correlation computed on horizontal planes:

\[
\mathrm{corr}(a,b)=
\frac{\langle (a-\langle a\rangle)(b-\langle b\rangle)\rangle}
{\sigma_a\sigma_b}.
\]

Commonly used:

\[
\mathrm{corr}(w,\omega_z),
\quad
\mathrm{corr}(w,T'),
\quad
\mathrm{corr}(w,q'),
\quad
\mathrm{corr}(w,b').
\]

Scripts:

- `analyze_ekman_pumping.py`
- `analyze_ekman_pumping_ra5e7.py`

Output tables:

- `comparison_outputs/ekman_pumping_correlations.csv`
- `comparison_outputs/ekman_pumping_correlations_ra5e7.csv`

Interpretation:

- Correlation measures geometric coherence, not transport amplitude.
- A high \(\mathrm{corr}(w,\omega_z)\) does not by itself imply large \(Nu_m\).

## Boundary-layer diagnostics

Moisture or moist-static-energy boundary-layer thicknesses have been exploratory and must be treated cautiously.

Existing files:

- `comparison_outputs/moisture_ekman_boundary_layer_thickness.csv`
- `comparison_outputs/moisture_boundary_layer_flux_crossover.csv`

Preferred future direction:

- define the thermal/moist boundary layer from consistent \(Nu_m(z)\) or \(m(z)\) gradients;
- define Ekman-layer thickness from a theory-checked velocity profile, not an arbitrary RMS maximum, unless explicitly stated.

## Kolmogorov-resolution check

For DNS-resolution screening, compute the velocity-gradient dissipation from the velocity fields:

\[
\epsilon_u(z,t)
=\nu\left\langle
\sum_{i,j}\left(\partial_j u_i\right)^2
\right\rangle_{x,y},
\qquad
\nu=\sqrt{\frac{Pr}{Ra}} .
\]

Then use the local Kolmogorov length:

\[
\eta_K(z,t)=\left(\frac{\nu^3}{\epsilon_u(z,t)}\right)^{1/4}.
\]

The current screening records:

- global and minimum \(\eta_K\);
- \(k_{\max}\eta_K\), using the \(2/3\)-dealiased horizontal cutoff;
- near-wall and all-domain \(\Delta z/\eta_K\);
- the number of grid points inside the nominal Ekman thickness \(\delta_E=\sqrt{Ek}\).

Current interpretation rule:

- `resolved`: \(k_{\max}\eta_{\min}\ge 1\) and near-wall \(\Delta z/\eta \le 2\);
- `marginal`: close to one of those thresholds;
- `under-resolved`: fails the strict small-scale or near-wall spacing check;
- `not_checked`: usually an initial or placeholder field with grid/coordinate mismatch, not a physical resolution verdict.

Current output:

- remote source host: `c01n0037`;
- source tree: `/share/org/SHUTUANL/shu_zhangjs/rainy model/ns/transition_study/beta1`;
- remote output folder: `/share/org/SHUTUANL/shu_zhangjs/rainy model/ns/transition_study/beta1/remote_profile_exports/kolmogorov_resolution_all_cases_20260711`;
- local output folder: `comparison_outputs/organized_results_20260705/kolmogorov_resolution`;
- table: `kolmogorov_resolution_all_cases_latest_field_clean.csv`.

This is a latest-field screening, not a time-averaged dissipation table. For final publication-quality claims, recompute selected cases over several late-time frames.

### 2026-08-11 restart-field resolution audit

The current `Ra=8e6` migration cases were re-screened locally from the restart
velocity fields with the Zhang et al. (2017) strain-tensor definition

\[
\epsilon_u=2\nu S_{ij}S_{ij},\qquad
\eta_K=(\nu^3/\epsilon_u)^{1/4}.
\]

The audit also records the gradient-form dissipation as an incompressibility
and interpolation cross-check, `dt/tau_eta`, the number of cell centers inside
the nominal `delta_E=sqrt(Ek)` layer, and a thermal-scalar Batchelor proxy.
Do not use the Batchelor proxy as the final q/m resolution verdict until the
actual moisture diffusivity or Schmidt number has been applied.

Current local script and output:

- script: `E:/moist RB/rotating_case_inventory/02_inventory_and_plot_scripts/check_kolmogorov_resolution_local.py`;
- output: `E:/moist RB/rotating_case_inventory/04_outputs_and_figures/kolmogorov_resolution_audit_20260811`;
- source fields: the 14 current-case restart folders in the verified 2026-08-08 migration bundle.

Interpretation rule for this audit:

- a low-energy or quiet intermittent snapshot can only show that the snapshot
  is resolved; it cannot establish resolution during a later burst;
- intermittent cases must be checked at several burst peaks and quiet troughs,
  and the worst dissipation scale controls the verdict;
- `kmax*eta >= 1` alone is insufficient when the nominal Ekman layer contains
  fewer than about eight z points or when `dt/tau_eta >= 0.01`;
- the nonrotating boundary-layer verdict remains unconfirmed until an actual
  viscous-layer thickness is measured from the velocity profiles.

## Vertical grid stretching

The current solver reads the grid parameters from `bou.in` at runtime. The relevant line is:

```text
ALX3D     REXT1     REXT2     ISTR3       STR3       Lmax
```

Current production cases usually use:

```text
1.0       2.0d0     2.0d0     1           16.0d0      1
```

In the source `cordin.f90`, `ISTR3=1` activates a clipped Chebyshev-like vertical grid. The parameter `STR3` is used as the number of clipped endpoint levels:

\[
n_{\rm clip}=STR3 \times mref3 .
\]

Smaller `STR3` retains more Chebyshev endpoint clustering and therefore gives finer near-wall spacing; larger `STR3` clips more of the endpoint clustering. Therefore, for new resolution-sufficient sweeps, estimate Ekman-layer resolution from the actual stretched grid, not from a uniform \(\Delta z=1/(N_3-1)\).

With the existing default `STR3=16`, approximate numbers of vertical cell centers inside \(\delta_E=\sqrt{Ek}\) are:

| \(N_3\) | \(Ek=10^{-5}\) | \(3.16\times10^{-5}\) | \(10^{-4}\) | \(3.16\times10^{-4}\) | \(10^{-3}\) | \(3.16\times10^{-3}\) | \(10^{-2}\) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 65 | 0 | 1 | 1 | 2 | 3 | 5 | 9 |
| 129 | 1 | 2 | 3 | 5 | 8 | 13 | 20 |
| 257 | 3 | 5 | 8 | 13 | 20 | 30 | 45 |
| 385 | 6 | 9 | 15 | 23 | 34 | 49 | 70 |
| 513 | 9 | 14 | 22 | 33 | 48 | 68 | 96 |

For small Ek, one can either increase \(N_3\) or reduce `STR3` case-by-case. For example, `N3=385, STR3=8` gives about 8 cell centers inside \(\delta_E\) at \(Ek=10^{-5}\), while `N3=257, STR3=4` gives about 9 centers inside \(\delta_E\) at \(Ek\approx3.16\times10^{-5}\).

## Barotropic / two-dimensionalization diagnostics

Barotropic fraction is intended to measure large-scale depth-independent horizontal flow:

\[
F_{\rm bt}
=
\frac{\frac12\langle |\overline{\boldsymbol{u}}^{\,z}(x,y)|^2\rangle_{x,y}}
{\frac12\langle |\boldsymbol{u}(x,y,z)|^2\rangle_{x,y,z}},
\]

where \(\overline{\boldsymbol{u}}^{\,z}\) is the depth-averaged horizontal velocity.

Existing script/table:

- `compute_barotropic_fraction.py`
- `comparison_outputs/barotropic_fraction.csv`

Use this, or low-wavenumber energy fraction, to discuss large-scale vortex / two-dimensionalization. Do not use \(l_h\) alone as evidence for large-scale vortex size.

## Relative humidity and velocity profiles

RH profile:

\[
RH(z)=\langle RH_{\rm me}\rangle_{x,y,t}.
\]

RH-structure horizontal length uses the relative-humidity field itself, not the fluctuation field. The \(k_h=0\) mode is excluded only because it has zero wavenumber and cannot define a finite horizontal scale.

\[
l_{RH}(z,t)=
\frac{\sum_{k_h>0}E_{RH}(k_h,z,t)}
{\sum_{k_h>0}k_hE_{RH}(k_h,z,t)},
\qquad
E_{RH}=|\widehat{RH}|^2 .
\]

As with \(l_h\) and \(l_{m'}\), the default time average is frame-first:

\[
\overline{l_{RH}}(z)=\langle l_{RH}(z,t)\rangle_t.
\]

Use \(l_{RH}\) when comparing to relative-humidity visualizations; it should not be conflated with the velocity-based \(l_h\).

Velocity profiles:

\[
u_{h,\rm rms}(z)=\langle u^2+v^2\rangle_{x,y,t}^{1/2},
\qquad
w_{\rm rms}(z)=\langle w^2\rangle_{x,y,t}^{1/2}.
\]

When a figure is labeled \(v_{\rm rms}\) in the current Ra=1e6 research-logic profile set, it means full speed rms:

\[
v_{\rm rms}(z)=\left\langle u^2+v^2+w^2\right\rangle_{x,y,t}^{1/2}.
\]

Field variables:

- `RH_me`, `VX_me`, `VY_me`, `VZ_me`.

Recent remote-averaged cell-regime outputs:

- `comparison_outputs/remote_cell_profiles`.

Current Ra=1e8 RH-length profile output:

- local output folder: `comparison_outputs/lhm_mprime/profiles_Ra1e8_Ek1e-4_to_1e-3_last30_rh`;
- direct-RH CSV: `Ra1e8_Ek1e-4_to_1e-3_lhRH_direct_profiles_last30_frame_mean.csv`;
- direct-RH summary at \(z\approx0.75\): `Ra1e8_Ek1e-4_to_1e-3_lhRH_direct_z075_summary.csv`;
- remote source: `c01n0020`, under `/share/org/SHUTUANL/shu_zhangjs/rainy model/ns/transition_study/beta1`;
- included Ek: \(10^{-4},2\times10^{-4},3\times10^{-4},4\times10^{-4},5\times10^{-4},6\times10^{-4},7\times10^{-4},8\times10^{-4},9\times10^{-4},10^{-3}\);
- averaging: latest 30 readable movie frames, frame-first length averaging;
- field variable: `RH_me`;
- figures: `lhRH_direct_profiles_Ra1e8_Ek1e-4_to_1e-3_last30.png` and `lambda_hRH_direct_profiles_Ra1e8_Ek1e-4_to_1e-3_last30.png`.

Current Ra=1e6 research-logic profile output:

- remote source host: `c01n0037`;
- remote source tree: `/share/org/SHUTUANL/shu_zhangjs/rainy model/ns/transition_study/beta1/lowrestest/Ra1e6`;
- remote output folder: `/share/org/SHUTUANL/shu_zhangjs/rainy model/ns/transition_study/beta1/remote_profile_exports/ra1e6_research_logic_profiles_20260712`;
- local output folder: `comparison_outputs/organized_results_20260705/ra1e6_research_logic_profiles`;
- profile table: `Ra1e6_research_logic_profiles_latest10.csv`;
- summary table: `Ra1e6_research_logic_profiles_summary_latest10.csv`;
- averaging: latest 10 readable movie fields per completed case;
- completed/readable cases included: `Ek=2e-4,3e-4,4e-4,5e-4,6e-4,7e-4,8e-4,9e-4,1e-3,2e-3,3e-3,8e-3,1e-2`, plus `norotating`;
- skipped cases with only initial tiny-grid fields, such as current `Ek=1e-4`.

The synchronized correlation and scalar-standard-deviation refresh is stored
remotely in
`/share/org/SHUTUANL/shu_zhangjs/rainy model/ns/transition_study/beta1/remote_profile_exports/ra1e6_corr_m_b_std_20260714`
and locally as `Ra1e6_corr_m_b_standard_deviation_latest10.csv` in the same
Ra=1e6 output folder. It uses the latest 10 complete fields for each case.
For `Ek=1e-3` this is the continuation interval `field00191`--`field00200`;
the other currently included cases use `field00071`--`field00080`.

The computed profiles are:

\[
\langle wm\rangle(z)=\left\langle w(b+\gamma q)\right\rangle_{x,y,t},
\qquad
\langle wq\rangle(z)=\left\langle wq\right\rangle_{x,y,t}.
\]

\[
\mathrm{corr}(w,\omega_z)(z),
\qquad
\omega_z=\partial_x v-\partial_y u,
\]

and

\[
\mathrm{corr}(w,b)(z).
\]

Scalar standard-deviation profiles use the combined horizontal/time variance, for example:

\[
q\ \mathrm{standard\ deviation}(z)
=\left(\left\langle q^2\right\rangle_{x,y,t}
-\left\langle q\right\rangle_{x,y,t}^2\right)^{1/2}.
\]

Use the identical strict definition for (b) and for
(m=b+1.1q):

\[
\sigma_m(z)=
\left(\left\langle m^2\right\rangle_{x,y,t}
-\left\langle m\right\rangle_{x,y,t}^2\right)^{1/2},
\qquad
\sigma_b(z)=
\left(\left\langle b^2\right\rangle_{x,y,t}
-\left\langle b\right\rangle_{x,y,t}^2\right)^{1/2}.
\]

Do not construct (sigma_m) by adding (sigma_b) and
(1.1\sigma_q), because that omits the (b)--(q) covariance. For
correlation profiles, calculate the horizontal-plane Pearson correlation in
each frame first and then average the resulting profile over the 10 frames.

The exploratory moisture/thermal boundary-layer comparison in this Ra=1e6 set uses scalar-standard-deviation peak locations. This is not yet a final thermal boundary-layer definition because many \(b\) or \(T\) standard-deviation profiles peak in the bulk rather than near the lower wall.

Current Ra=1e6 moisture/Ekman boundary-layer comparison for
\(7\times10^{-4}\le Ek\le3\times10^{-3}\):

- local output folder: `comparison_outputs/organized_results_20260705/ra1e6_research_logic_profiles`;
- source profiles: `Ra1e6_research_logic_profiles_latest10.csv`;
- reduced table: `Ra1e6_deltaq_deltaE_Ek7e-4_to_3e-3.csv`;
- averaging: latest 10 readable remote fields, with horizontal statistics followed by time averaging;
- lower moisture thickness \(\delta_q\): first local maximum of the combined \(q\) standard-deviation profile measured from the lower wall;
- Ekman thickness \(\delta_E\): first near-wall local maximum of \(u_{h,\mathrm{rms}}(z)\) measured from the lower wall;
- peak positions are refined by a local three-point quadratic interpolation on the stretched vertical grid;
- the \(Ek^{1/2}\) comparison fixes the exponent at \(1/2\) and determines only one vertical prefactor in log space; it is a reference-slope test, not a free-exponent fit.

For this current six-point interval, \(\delta_q/\delta_E>1\), while a free log-log diagnostic slope of the measured \(\delta_E\) is approximately 0.31. Therefore the data show the expected increase of Ekman thickness with \(Ek\), but do not quantitatively establish an \(Ek^{1/2}\) law over this limited interval.

## Funnel/environment moist-static-energy flux decomposition

The first-pass funnel diagnostic decomposes the instantaneous horizontal moist-static-energy flux after explicitly removing horizontal means:

\[
w'=w-\langle w\rangle_{xy},\qquad
m'=m-\langle m\rangle_{xy},\qquad m=b+\gamma q.
\]

The current default mask is evaluated independently at every height and frame:

\[
I_f=\mathbf{1}_{w'>0}\,\mathbf{1}_{\omega_z>0}.
\]

Then

\[
A_f=\langle I_f\rangle_{xy},\qquad
F_f=A_f\langle w'm'\mid I_f=1\rangle,
\]

\[
F_{env}=(1-A_f)\langle w'm'\mid I_f=0\rangle,qquad
F_{total}=F_f+F_{env}.
\]

The conditional factorization is

\[
F_f=A_f\,C_f\,\sigma_{w,f}\,\sigma_{m,f}.
\]

Current outputs:

- remote host: `c01n0037`;
- remote folder: `remote_profile_exports/funnel_transport_decomposition_20260712`;
- local folder: `comparison_outputs/organized_results_20260705/funnel_transport_decomposition`;
- separate result folders: `Ra1e6` and `Ra1e8`;
- averaging: latest 10 readable fields per case, instantaneous mask/statistics first and time averaging second;
- bulk interval: `0.2 <= z/H <= 0.8`, integrated on the nonuniform z grid;
- representative joint PDFs: use the nearest plane to `z/H=0.5`; remove the
  horizontal mean separately in each frame, compute the All/Funnel/Environment
  density on common bins for each of the latest 10 fields, and then take an
  equal-weight mean of the 10 framewise PDFs. Do not use only the last frame or
  pool unequal conditional sample counts;
- tables: `transport_decomposition_profiles.csv`, `transport_decomposition_timeseries.csv`, `transport_decomposition_summary.csv`, and `transport_decomposition_validation.csv`.

Instantaneous identification QA for the exact sign-only transport mask is
stored locally in
`comparison_outputs/organized_results_20260705/funnel_identification_visualization`
and remotely in
`remote_profile_exports/funnel_identification_visualization_20260713`.
For the four Ra=1e8 high-resolution funnel cases (`Ek1e-4`, `Ek3e-4`,
`Ek7e-4`, and `Ek1e-3`), the reduced planes contain (w'), (omega_z), and
the exact binary mask at the nearest planes to `z/H=0.5`, `0.75`, and `0.9`.
These are instantaneous visual checks at the established representative late
frames, not time-averaged morphology. Keep this sign-mask QA distinct from the
periodic Q-criterion component masks used for vortex radius and spacing.

Important interpretation rule: this sign-only mask is a cyclonic-updraft subset. It should not be treated as a unique morphological funnel in nonrotating, cell, or plume cases. Final morphology claims require sensitivity checks with percentile thresholds and an RH or condensation gate. The implemented module supports `w`, `condensation`, and `combined` modes and configurable absolute/percentile thresholds.

Current mechanism result:

- at Ra=1e6, rotation concentrates transport into the selected subset but environmental transport decreases enough to cancel most of the subset gain;
- at Ra=1e8 in the funnel interval, the subset contribution grows much more strongly than the environment decreases;
- the Ra=1e8 gain is caused by the product of stronger `sigma_w_f`, stronger `sigma_m_f`, and stronger `Corr_f`, while the selected area fraction generally decreases relative to nonrotating.

## Modal kinetic-energy spectra

When the user asks for the previous or standard 2D/3D kinetic-energy spectra, use the modal decomposition from the older cascade-analysis scripts, not a height-by-height 2D horizontal spectrum.

The barotropic 2D mode is computed first by a \(z\)-weighted vertical average:

\[
U(x,y,t)=\langle u\rangle_z,\qquad
V(x,y,t)=\langle v\rangle_z .
\]

The 3D/baroclinic velocity is then

\[
u'(x,y,z,t)=u-U,\qquad
v'(x,y,z,t)=v-V,\qquad
w'(x,y,z,t)=w .
\]

For horizontal shell wavenumber \(k_h\), the spectra are

\[
E_{\rm total}(k_h)=
\frac12\left\langle
\sum_{|\mathbf{k}_h|\in k_h}
\left(|\hat u|^2+|\hat v|^2+|\hat w|^2\right)
\right\rangle_z / \Delta k ,
\]

\[
E_{2D}(k_h)=
\frac12
\sum_{|\mathbf{k}_h|\in k_h}
\left(|\hat U|^2+|\hat V|^2\right) / \Delta k ,
\]

\[
E_{3D}(k_h)=
\frac12\left\langle
\sum_{|\mathbf{k}_h|\in k_h}
\left(|\widehat{u'}|^2+|\widehat{v'}|^2+|\hat w|^2\right)
\right\rangle_z / \Delta k .
\]

Do not replace \(E_{2D}\) with a vertical average of per-height horizontal spectra. That quantity is a different diagnostic and does not isolate the barotropic mode.

### Height-resolved full-field spectrum and Gaussian flux

When the user explicitly asks for the complete field at selected heights
without 2D/3D decomposition, use a distinct diagnostic. At each height,

\[
E(k_h,z)=\frac{1}{2\Delta k}\sum_{|\mathbf{k}_h|\in k_h}
(|\hat u|^2+|\hat v|^2+|\hat w|^2).
\]

For horizontal Gaussian coarse graining,

\[
G_\ell(k_h)=e^{-\ell^2 k_h^2/2},\qquad k_c=1/\ell,
\]

\[
\tau_{ij}^{(\ell)}=\overline{u_i u_j}^{\ell}
-\bar u_i^{\ell}\bar u_j^{\ell},\qquad
\Pi(k_c,z)=-\langle\tau_{ij}^{(\ell)}S_{ij}^{(\ell)}\rangle_{xy}.
\]

Use all three velocity components and the complete 3D strain tensor at the
requested height. Horizontal derivatives are spectral; vertical derivatives
use adjacent planes on the actual stretched z grid. Positive `Pi` is
forward/downscale transfer and negative `Pi` is inverse/upscale transfer. This
quantity must not be labeled 2D or 3D modal flux.

Current Ra=1e8 calculation:

- high-resolution cases: `Ek1e-4`, `Ek3e-4`, `Ek7e-4`, `Ek1e-3`;
- requested heights: `0.70`, `0.90`, `0.95`; actual planes: `0.6957433`,
  `0.9022776`, `0.9503293`;
- averaging: latest 10 fields, frame-first spectra/flux then time average;
- remote host/output: `c01n0037`,
  `remote_profile_exports/height_spectra_gaussian_flux_20260713`;
- local output:
  `comparison_outputs/organized_results_20260705/height_spectra_gaussian_flux`;
- plot only complete isotropic shells with `k_h <= pi/dx` and Gaussian cutoff
  coordinates with `k_c <= 0.9 pi/dx`.

For the matching one-dimensional nonlinear shell-transfer spectrum at a
height, use

\[
\mathbf N=-(\mathbf u\cdot\nabla)\mathbf u,
\]

\[
T(k_p,z)=\frac{1}{\Delta k}\sum_{|\mathbf{k}_h|\in k_p}
\mathrm{Re}[\hat u^*\widehat{N_u}+\hat v^*\widehat{N_v}
+\hat w^*\widehat{N_w}].
\]

Positive `T` means that the shell gains energy from nonlinear advection;
negative `T` means that it loses energy. Compute each frame first and average
the resulting transfer spectra. At a fixed z, the shell integral need not be
zero because vertical advection transports kinetic energy between planes, so
do not equate the sign of `T(k_p,z)` directly with cumulative cascade direction.
Use Gaussian `Pi(k_c,z)` for the direct coarse-grained flux sign.

The legacy files are `height_full_velocity_Tkp_frames.csv`,
`height_full_velocity_Tkp_10frame_summary.csv`, and `height_Tkp_metadata.csv`
in the height spectra/Gaussian flux output folder. Despite the historical
`Tkp` filename, this is `T(k,z)`, not a shell-to-shell `T(K,P)` matrix.

### Height-resolved shell-to-shell kinetic-energy transfer

For a true horizontal shell-to-shell matrix at fixed height, shell `P` is the
donor and shell `K` is the receiver. With all three velocity components,

\[
T_h(K,P;z)=-\frac{1}{\Delta k}
\left\langle u_i^K (u\partial_x+v\partial_y)u_i^P\right\rangle_{xy},
\]

\[
T_v(K,P;z)=-\frac{1}{\Delta k}
\left\langle u_i^K w\partial_z u_i^P\right\rangle_{xy},
\qquad
T(K,P;z)=T_h+T_v.
\]

Positive `T(K,P)` means donor shell `P` supplies kinetic energy to receiver
shell `K`. Horizontal derivatives are spectral. The vertical derivative of
the donor-filtered field uses adjacent planes on the actual stretched grid.
The advecting velocity is the complete unfiltered velocity.

At fixed height the matrix is not generally antisymmetric. In the continuous
equations,

\[
T(K,P)+T(P,K)
=-\frac{1}{\Delta k}\partial_z
\left\langle w\,\boldsymbol{u}^K\cdot\boldsymbol{u}^P\right\rangle_{xy}.
\]

Therefore use

\[
T_{\rm exchange}(K,P)=\frac{T(K,P)-T(P,K)}{2}
\]

for conservative inter-shell exchange, and

\[
T_{\rm vertical}(K,P)=\frac{T(K,P)+T(P,K)}{2}
\]

for the pair-symmetric vertical-transport contribution. The raw `T_v` term is
the contribution from vertical advection `w partial_z`; it is not by itself a
vertical-wavenumber cascade. A vertical spectral cascade would require a full
3D modal expansion compatible with the nonperiodic vertical boundaries.

Current calculation convention:

- `Ra=1e8` high-resolution cases `Ek1e-4`, `Ek3e-4`, `Ek7e-4`, `Ek1e-3`;
- requested heights `z/H=0.70,0.90,0.95`;
- latest 10 readable fields, matrix first in every frame and then arithmetic
  time average;
- horizontal domain `Lx=Ly=2`, `Delta k=pi`, shells `1<=K,P<=75`;
- remote host/output: `c01n0037`,
  `remote_profile_exports/height_shell_to_shell_transfer_20260713`;
- local output:
  `comparison_outputs/organized_results_20260705/height_shell_to_shell_transfer`.

Current corrected local output for the three-case Ra=1e8 comparison:

- local output folder: `comparison_outputs/organized_results_20260705/spectra_flux`;
- spectra table: `old_modal_velocity_spectra_three_cases_local_field00080.csv`;
- summary table: `old_modal_velocity_spectra_three_cases_local_field00080_summary.csv`;
- figures: `Ra1e8_three_cases_old_modal_spectra_field00080_same_panel.png` and `Ra1e8_three_cases_old_modal_spectra_field00080_overlay.png`;
- included cases: \(Ek=3\times10^{-4},7\times10^{-4},10^{-3}\);
- current source: local `field00080.h5` and `cordin_info.h5` in each case folder, because the remote SSH proxy was temporarily unavailable;
- averaging: single-frame field diagnostic for this corrected quick version. When remote access is available, recompute the same formulas over the latest 10 frames before final quantitative interpretation.

## Kinetic energy time series

Current plotting convention:

\[
K(t)=\frac12 A_4(t),
\]

where \(A_4(t)\) is the 4th column of `avgvar.out` using one-based counting, i.e. `avgvar[:,3]` in Python zero-based indexing.

The current production executable emits two horizontal integral wavelengths in
`avgvar.out`:

- column 8, `avgvar[:,7]`: \(l_0\) near \(z=0.5\);
- column 9, `avgvar[:,8]`: \(l_0\) near \(z=0.75\).

They are distinct from the default vertical-velocity \(l_h\). Direct field
validation gives the effective output convention

\[
l_0(z,t)=2\pi
\frac{\sum_{k_h>0}E_{u_h}(k_h,z,t)/k_h}
{\sum_{k_h>0}E_{u_h}(k_h,z,t)},
\qquad
E_{u_h}=|\hat u|^2+|\hat v|^2 .
\]

The archived `RaT/source/fluid_solver/avgvar.f90` snapshot shows one mid-height
moment, while the executable emits two lengths and values matching the
effective `2 pi` wavelength convention. Treat current output plus direct field
validation as authoritative; the source snapshot is not identical to the
production executable.

Use `avgvar.out` for solver time series. At arbitrary heights, recompute the
same horizontal-velocity moment from full fields and label it as
field-reconstructed \(l_0(z,t)\).

Current remote-reduced output:

- remote source host: `c01n0020`;
- source tree: `/share/org/SHUTUANL/shu_zhangjs/rainy model/ns/transition_study/beta1`;
- remote output folder: `/share/org/SHUTUANL/shu_zhangjs/rainy model/ns/transition_study/beta1/remote_profile_exports/avgvar_ke_l0_timeseries_Ra1e6_Ra1e8`;
- local output folder: `comparison_outputs/organized_results_20260705/time_series_avgvar`;
- raw table: `avgvar_ke_l0_timeseries_Ra1e6_Ra1e8.csv`;
- filtered plot table: `avgvar_ke_l0_timeseries_Ra1e6_Ra1e8_filtered_for_plots.csv`.

Current AR=10 aspect-ratio kinetic-energy time series (updated 2026-07-15):

- remote source host: `c01n0037`;
- source tree: `/share/org/SHUTUANL/shu_zhangjs/rainy model/ns/transition_study/beta1/lowrestest/aspect_ratio_study`;
- included cases are exactly the AR10 runs with valid rows in `run/data/avgvar.out`:
  - (Ra=10^6): (Ek=9\times10^{-4}) and (10^{-2}), both through (t\simeq800);
  - (Ra=10^8): (Ek=10^{-4},9\times10^{-4},2\times10^{-3},3\times10^{-3},5\times10^{-3},7\times10^{-3},10^{-2}) plus nonrotating;
  - all listed Ra=10^8 cases currently reach (t\simeq800), except (Ek=9\times10^{-4}), which merges its continuation through (t\simeq1600), and the still-running (Ek=3\times10^{-3}) case, whose 2026-07-15 latest reduced snapshot reaches (t=55.1);
- remote reduced output folder: `/share/org/SHUTUANL/shu_zhangjs/rainy model/ns/transition_study/beta1/remote_profile_exports/ar10_kinetic_energy_timeseries_20260715`;
- local output folder: `comparison_outputs/organized_results_20260705/aspect_ratio_study/AR10_kinetic_energy_timeseries`;
- reduced table: `AR10_kinetic_energy_timeseries.csv`;
- provenance/status table: `AR10_kinetic_energy_timeseries_metadata.csv`;
- definition remains \(K(t)=A_4(t)/2\), where one-based `avgvar.out` column 4 is \(\langle u^2+v^2+w^2\rangle_V\); data are the raw solver time records after removing nonfinite values, duplicate physical times, and uninitialized rows with \(K\ge0.1\). No temporal smoothing or averaging is applied.
- the current `Ek=9e-4` continuation source is
  `Ra1e8/AR10/Ra1e8Pr07Ek9e-4AR10/conti1/run/data/avgvar.out`;
  its physical solver time starts at `800.100`, so concatenate it directly
  with the base segment ending at `800.002` rather than adding a manual time
  offset.

MSE-paired AR10 kinetic-energy update (2026-07-16):

- includes exactly the four cases with `run/data/mse_aggregation.out`:
  `Ek=1e-4`, `Ek=3e-3`, nonrotating `128x128x65`, and nonrotating
  `256x256x65`;
- the first three reach `t=800`; the running 256 case snapshot reaches
  `t=295.6`;
- output folder:
  `comparison_outputs/organized_results_20260705/aspect_ratio_study/AR10_MSE绠椾緥鍔ㄨ兘鏃堕棿搴忓垪_AR10_MSE_cases_kinetic_energy_timeseries_20260716`;
- the definition remains `K=0.5*avgvar[:,3]`; curves are raw solver records
  after finite-value, duplicate-time, and `K<0.1` checks, with no smoothing.

AR=10 Ra=1e8 weak-rotation case preparation (2026-07-14):

- prepared but **not submitted** under
  `/share/org/SHUTUANL/shu_zhangjs/rainy model/ns/transition_study/beta1/lowrestest/aspect_ratio_study/Ra1e8/AR10`;
- complete prepared endpoint/infill set:
  \(Ek=10^{-3},2\times10^{-3},3\times10^{-3},\ldots,9\times10^{-3},10^{-2}\);
- all use \(Ra=10^8\), \(Pr=0.7\), grid `129 129 65`, horizontal aspect ratio `REXT1=REXT2=10`, `DT=DTMAX=7e-3`, `NREAD=0`, and `TMAX=800`;
- `invRo` is generated from
  \(invRo=\sqrt{Pr/Ra}/Ek\) and is recorded case by case in
  `AR10_Ra1e8_weak_rotation_preparation_manifest.csv` at the aspect-ratio-study root;
- the existing `Ra1e8Pr07Ek9e-4AR10` run is continued in
  `Ra1e8Pr07Ek9e-4AR10/conti1/run`, with `NREAD=1`, `TMAX=1600`, and job name
  `Ra1e8Pr07Ek9e-4AR10_conti1600`;
- the continuation preserves the original rounded `invRo=0.092962` exactly and uses byte-identical copies of the six `continua_*` restart files written at `t=800.002`; changing it to the more precise formula value during a continuation would introduce a parameter discontinuity;
- reusable preparation script:
  `comparison_outputs/organized_results_20260705/aspect_ratio_study/remote_case_preparation/prepare_ar10_weak_rotation_cases_and_conti.py`, mirrored remotely in `aspect_ratio_study/scripts/`.

Current Ra=1e8 high-resolution time-height scale audit:

- cases: `Ek1e-4`, `Ek3e-4`, `Ek7e-4`, `Ek1e-3`;
- remote host/output: `c01n0037`,
  `remote_profile_exports/lh_l0_time_height_20260713`;
- local output:
  `comparison_outputs/organized_results_20260705/lh_l0_time_height_Ra1e8_highres`;
- field table: `Ra1e8_highres_lh_l0_time_height_all_cases.csv`;
- 32 stretched-grid heights per frame, frame-first spectral moments;
- `Ek7e-4` `t=10` is a `31x15x15` placeholder and is excluded;
- for visual wavelength comparison save `lambda_h=2*pi*l_h` separately; do
  not change the established no-`2 pi` definition of `l_h` itself;
- save `lambda_peak,w=2*pi/k_peak,w` as a discrete dominant-spacing check,
  not as a vortex-core radius.

Current mechanism result at `z approximately 0.75`: from early `t=20..100`
to mature `t=600..800`, `l0` grows by about `60%` for `Ek3e-4` and `Ek7e-4`
and their dominant `w` wavelength doubles to the box scale, while `l_h` does
not grow. This supports a nonlinearly organized/merged mature funnel scale.
`Ek1e-4` shows no corresponding `l0` or peak-wavelength growth and remains
closer to an instability-constrained convective scale. Keep the finite-box
qualification until the `Gamma=4` and `Gamma=10` runs are available.

Filtering rule for completed-case summary plots:

- merge continuation `avgvar.out` files, remove duplicate times, and keep only long runs with \(t_{\max}\ge300\);
- remove placeholder/uninitialized rows with \(K\ge0.1\). In the current files, these correspond to `avgvar` 4th-column values near 1 and are not physical kinetic-energy data.

For an explicitly requested "all currently running AR10 cases" time-series
update, retain valid short current output as well and mark it
`short_current_output` in the metadata rather than presenting it as
statistically complete. The 2026-07-15 AR10 figure follows this rule for
Ra=1e8, Ek=3e-3.

Before interpreting absolute values across different \(Ra\), verify the exact `avgvar.out` column and nondimensionalization. The nonrotating comparison file is:

- `comparison_outputs/three_nonrotating_kinetic_energy_timeseries.csv`.

## Plot style

Use `jiasen-scientific-plot-style`.

Stable Ra colors:

- \(Ra=10^6\): saturated green;
- \(Ra=5\times10^6\): saturated purple;
- \(Ra=5\times10^7\): saturated red;
- \(Ra=10^8\): saturated blue.

Use filled markers by default. Avoid error bars unless explicitly requested.

## Required update protocol

Whenever a new diagnostic is computed or a definition changes:

1. Add or update the formula in this document.
2. Record the data source path and output CSV/NPZ.
3. Record the averaging window and whether averaging is frame-first, field-first, profile-first, or time-spectrum weighted.
4. If a plot uses a non-default definition, state it in the filename or legend/caption.
5. Update `moist_rb_core_refs.md` only for high-level theory/regime memory; keep detailed diagnostic definitions here.

## Funnel/environment transport case scope

For the Ra=1e8 funnel/environment moist-static-energy transport decomposition,
only the user-defined funnel regime is included as rotating data:

\[
10^{-4}\le Ek\le10^{-3}.
\]

The retained rotating cases are `Ek1e-4`, `Ek2e-4`, `Ek3e-4`, `Ek4e-4`,
`Ek5e-4`, `Ek6e-4`, `Ek7e-4`, `Ek8e-4`, `Ek9e-4`, and `Ek1e-3`.
Keep `norotating` only as the normalization/reference baseline. Do not include
the stronger-rotation cell/unclassified cases below `Ek1e-4` or the plume cases
above `Ek1e-3` in Ra=1e8 funnel-mechanism figures or numerical ranges.

## Cyclonic-vortex radius and spacing

For `Ra=1e8` vortex-geometry scaling in the funnel interval, use the latest 10
readable fields of each case and the nearest planes to `z/H=0.2`, `0.5`,
`0.75`, and `0.9`. Identify geometry in each frame first and then average the
framewise results.

The primary mask is the horizontal Q criterion with cyclonic sign:

\[
Q_{2D}=\frac12(\lVert\Omega_h\rVert^2-\lVert S_h\rVert^2),\qquad
Q_{2D}>Q_{2D,\mathrm{rms}},\quad\omega_z>0.
\]

Use periodic eight-connected components and discard components with equivalent
radius below `0.025H`. For component area `A_i`, define

\[
r_i=\sqrt{A_i/\pi},\qquad
r_v=\frac{\sum_i A_i r_i}{\sum_i A_i},\qquad
\ell_v=\sqrt{\frac{L_xL_y}{N_v}}.
\]

Here `r_v` is the area-weighted equivalent-core radius and `ell_v` is the
number-density spacing. Keep `omega_z>2 omega_z,rms` only as a sensitivity
definition because it includes shear filaments at weak rotation.
Direct periodic nearest-neighbor distance may be plotted as a secondary check,
but do not fit it when many frames contain fewer than two retained vortices;
use `ell_v` for the primary spacing exponent in that situation.

Current source/output:

- remote host: `c01n0037`;
- remote reduced output: `remote_profile_exports/vortex_geometry_20260713`;
- local output: `comparison_outputs/organized_results_20260705/vortex_geometry`;
- frame table: `Ra1e8_funnel_vortex_geometry_frames.csv`;
- 10-frame summary: `Ra1e8_funnel_vortex_geometry_10frame_summary.csv`;
- fit table: `Ra1e8_funnel_vortex_geometry_powerlaw_fits.csv`.

Current primary-Q result: `ell_v ~ Ek^0.36` at `z/H=0.2` (`R^2=0.83`) and
`ell_v ~ Ek^0.40` at `z/H=0.5` (`R^2=0.85`), consistent with an approximate
`Ek^(1/3)` spacing law. The fits at `z/H=0.75` and `0.9` are weak, and the
equivalent-core radius has no robust single power law.

For identification QA, use the latest-frame `z/H鈮?.5` reduced planes in
`vortex_identification_all_cases_z05.npz`. Plot normalized Q or vertical
vorticity as the background and the retained, periodic, radius-filtered mask
as a black contour. Keep these figures explicitly described as instantaneous
visual checks rather than 10-frame statistics.

## Per-frame RMS-standardized visualization fields

For ParaView/XDMF visualization requests described as an RMS field, retain
the full spatial structure by standardizing independently on every horizontal
plane and frame. Do not replace the three-dimensional field by a repeated
one-dimensional RMS profile. Use

\[
\widetilde{\omega}_z=
\frac{\omega_z-\langle\omega_z\rangle_{xy}}
{\sqrt{\langle(\omega_z-\langle\omega_z\rangle_{xy})^2\rangle_{xy}}},
\]

\[
\widetilde b=
\frac{b-\langle b\rangle_{xy}}
{\sqrt{\langle(b-\langle b\rangle_{xy})^2\rangle_{xy}}},\qquad
\widetilde m=
\frac{m-\langle m\rangle_{xy}}
{\sqrt{\langle(m-\langle m\rangle_{xy})^2\rangle_{xy}}},
\quad m=b+1.1q.
\]

For the local Ra=1e6, AR=10, Ek=9e-4 visualization case, the source is
`transition_study/lowres/AR10/Ek9e-4/field*.h5`. The 33 existing source
frames are stored on a `64 x 128 x 128` visualization grid. Derived raw and
standardized HDF5/XMF pairs are in `derived_rms_fields/{omega_z,b,m}`, with
temporal collection files `omega_z_rms_all_frames.xmf`,
`b_rms_all_frames.xmf`, and `m_rms_all_frames.xmf`. Each HDF5 file also stores
the framewise horizontal mean and RMS/standard-deviation profiles used for
normalization. The source sequence has no frame 26 and jumps from frame 40 to
frame 58; process all files that actually exist without inventing missing
frames.

The same standardized-field convention is also used for the local Ra=1e6,
AR=10, Ek=9e-4 visualization case at
`transition_study/lowres/AR10/Ra1e6/Ek9e-4`. Its current local source set
contains frames 5, 8, and 80 (times approximately 50, 80, and 800), and the
derived HDF5/XMF files are stored in that case's `derived_rms_fields` folder.

Two additional local low-resolution visualization cases use the same derived
field convention and output layout:

- `transition_study/lowres/Ra1e6/Ek1e-3`, currently source frame 80 at time 800;
- `transition_study/lowres/Ra1e8/Ek9e-4`, currently source frame 80 at time 800.

For the centered x-z comparison of these two cases, identify one cyclonic
center in each frame from the maximum of a periodic Gaussian-smoothed
\(\widetilde{\omega}_z\) plane at the nearest grid level to \(z/H=0.75\)
(smoothing width: two horizontal grid cells). Use the detected center's fixed
\(y\) coordinate for the x-z cut, then translate periodically in x so that
the center appears at \(x/H\simeq1\). Overlay the projected instantaneous
velocity field \((u,w)\) as streamlines. The current frame-80 centers are
\((x/H,y/H)=(0.1328125,0.9140625)\) for `Ra1e6/Ek1e-3` and
\((1.7421875,0.7109375)\) for `Ra1e8/Ek9e-4`; both were detected at
\(z/H=0.745985\). Use a shared displayed color range of
\(-8\le\widetilde{\omega}_z\le8\) to expose the weaker low-Ra structure;
larger values are only color-clipped, not removed from the data. Current
output: `transition_study/lowres/comparison_outputs/xz_vorticity_streamlines_centered`.

Their `derived_rms_fields` folders contain the raw and standardized
three-dimensional fields, the horizontal mean/RMS profiles, the individual
frame XMF files, and the three temporal-collection XMF entry points.

## Moist convective self-aggregation module (2026-07-15)

Use moist static energy

\[
m(x,y,z,t)=b(x,y,z,t)+\gamma q(x,y,z,t),\qquad \gamma=1.1,
\]

and always remove the instantaneous horizontal mean separately at every
height,

\[
m'(x,y,z,t)=m(x,y,z,t)-\langle m(x,y,z,t)\rangle_{xy}.
\]

The primary self-aggregation amplitude is not a three-dimensional variance
about one volume mean. First compute

\[
{\rm Var}_m(z,t)=\left\langle [m'(x,y,z,t)]^2\right\rangle_{xy},
\]

then vertically average it with the actual stretched-grid cell thicknesses,

\[
A_m(t)=\frac{\int {\rm Var}_m(z,t)\,dz}{\int dz}.
\]

The modified DNS keeps every legacy column of `data/avgvar.out` unchanged and
writes two new files from the same `avgvar` calls:

- `data/mse_aggregation.out`: time, `A_m`, and `sqrt(A_m)`;
- `data/mse_variance_profile.out`: time, z, `mean_xy(m)`, and `Var_xy(m)`.

Modified source and build provenance:

- build/source root: `/share/org/SHUTUANL/shu_zhangjs/rainy model/ns/transition_study/aggregation_mse_dns_20260715`;
- modified files: `source/fluid_solver/avgvar.f90` and `source/fluid_solver/openfi.f90`;
- source template matching the previous production binary: `/share/org/SHUTUANL/shu_zhangjs/rainy model/ns/transition_study/Rotating_Moist_RB3-ns-deltab609`;
- new `simexec` MD5: `ccab677d8b20edffc60f2b20c53e514b`;
- previous `simexec` MD5: `dd6857c04f9ad3e14147e018a326b261`;
- deployed to all 28 existing Ra=1e6 and Ra=1e8 AR10 `run` directories;
- each old binary is retained as `simexec_before_mse_aggregation_20260715`;
- replacing the executable does not modify an already running process. New
  online files begin only after a new start or continuation uses the updated
  executable. Deployment does not submit or restart jobs.

The full-field aggregation postprocessor is intentionally separate from
`avgvar` because spectra, correlations, clusters, and Voronoi cells require
the two-dimensional horizontal fields. Its authoritative remote location is

`/share/org/SHUTUANL/shu_zhangjs/rainy model/ns/transition_study/beta1/lowrestest/aspect_ratio_study/self_aggregation_tools`.

For each physical-time snapshot, it computes diagnostics first and only then
forms temporal summaries. Continuation files are ordered and de-duplicated by
physical time. The default two-dimensional aggregation field is the
stretched-grid-weighted vertical mean of `m'`; a selected-height mode is also
available. Definitions and output rules are:

- peak wavelength: radially sum `|FFT2(m'_2D)|^2`, exclude `k=0`, find the
  shell containing the spectral maximum, and use `L_peak=2*pi/k_peak`;
- integral scale: compute the normalized periodic two-dimensional
  autocorrelation, radially average it, and integrate its positive lobe from
  zero separation to the first zero crossing;
- clusters: use periodic eight-connected components for both `m'>0` and
  `m'>sigma_m`; save area, equivalent radius `sqrt(A/pi)`, and maximum
  periodic extent for every object and time;
- Voronoi: find periodic local maxima above `sigma_m`, construct the tessellation
  with a 3-by-3 periodic tiling, and save central-cell areas and `sqrt(A)`;
- do not identify a single vortex radius with the aggregation scale. MSE
  organization and dynamical-vortex geometry remain separate diagnostics.

Per-case output folder: `aggregation_metrics/`. Required arrays include
`MSE_variance.npy`, `Var_m_z_time.npy`, `L_peak.npy`, `L_integral.npy`,
`cluster_radius_threshold0.npy`, `cluster_radius_threshold_sigma.npy`, and
`voronoi_area.npy`, together with frame provenance, a summary JSON, and the
standard PNG/PDF figures. The cross-case script uses the final ten available
snapshots for the reported final averages.

Desktop archive:

`C:/Users/jiasenzhang/Desktop/淇敼绋嬪簭/Rainy_Benard_self_aggregation_20260715`.

### Current AR10 aggregation products and Voronoi limitation

Current local `Ra=1e8, AR=10, Ek=3e-3` output:

`comparison_outputs/organized_results_20260705/aspect_ratio_study/AR10_self_aggregation_Ek3e-3_current`.

At the 2026-07-15 check, online `A_m` reached `t=174.5` and 17 complete
fields reached `t=170`. Online and full-field `A_m` agreed to a maximum
absolute difference of `4.1e-11` at matching times. The data show coarsening
but do not yet establish a saturated aggregation scale.

Latest-ten-field periodic Voronoi products for all AR10 cases that had reached
physical time `t>=800` are stored locally at

`comparison_outputs/organized_results_20260705/aspect_ratio_study/AR10_Voronoi_latest10_tge800_20260715`

and remotely at

`/share/org/SHUTUANL/shu_zhangjs/rainy model/ns/transition_study/beta1/lowrestest/aspect_ratio_study/remote_profile_exports/voronoi_latest10_tge800_20260715`.

For the requested averaged visualization, construct `m'_2D` separately in
each field, average the latest ten two-dimensional fields, and identify
periodic local maxima of the mean field above its own standard deviation.
Also retain framewise core counts because moving structures can be smeared by
the ten-field mean.

Do not use Voronoi cell area as the primary aggregate-size measure for broad
weak-rotation or nonrotating MSE structures. Voronoi partitions the entire
domain around point maxima and can split one broad aggregate into many cells.
Use it for isolated-core spacing. For aggregate size, use a low-pass MSE mask,
periodic connected components, radius of gyration, spectral-moment length,
and correlation length, with threshold/filter sensitivity checks.

### Online aggregation scales and z/H=0.75 vertical-velocity length

Update date: 2026-07-16. The modified DNS source and executable are stored at

`/share/org/SHUTUANL/shu_zhangjs/rainy model/ns/transition_study/beta_not1/program_online_scales_20260716`.

The compiled executable is
`source/simexec_online_mse_scales_wz075_20260716`, with SHA-256
`62a553db4599b142da91cdbea9ae84c9fdfed7b32f2548e23c68223a3ac35081`.
It was compiled from the existing online-MSE source with Intel Fortran,
MPI, HDF5, FFTW, bounds checking, and traceback enabled. The previous
`avgvar.f90` and `openfi.f90` are retained under
`source/backup_before_online_scales_20260716`.

All new quantities are evaluated at every `avgvar` call, hence at the same
`TPIN` cadence as kinetic energy. The legacy `avgvar.out` columns are not
changed.

`data/w_z075_spectral_length.out` contains `time`, `ell_w`, and
`2*pi*ell_w`, where the vertical velocity is linearly interpolated to
`z/H=0.75` before the horizontal FFT and
`ell_w=sum(|w_hat|^2)/sum(k_h|w_hat|^2)`. Exclude the zero horizontal mode
and restore the Hermitian multiplicity of the real-to-complex transform.
This supersedes the tentative online midplane version; do not use the
midplane label for this executable.

`data/mse_aggregation_scales.out` uses the stretched-grid-weighted field

`m2d(x,y)=integral[m-mean_xy(m)(z)]dz/integral dz`, with `m=b+gamma*q`.

Its columns are: time, `L_peak`, `L_integral`, MSE spectral-moment `ell_m`,
`2*pi*ell_m`, cluster count/mean equivalent radius/maximum equivalent radius
for `m2d>0`, and the same three cluster statistics for `m2d>sigma_m`.
`L_peak=2*pi/k_peak` uses radial shell-summed power. `L_integral` is the
positive-lobe integral of the radially averaged periodic autocorrelation to
its first zero. Clusters use periodic eight-neighbour connectivity. Voronoi
is deliberately not an online primary scale because broad weak-rotation
aggregates can contain several local maxima.

Desktop synchronized source archive:

`C:/Users/jiasenzhang/Desktop/淇敼绋嬪簭/Rainy_Benard_beta_not1_online_scales_wz075_20260716`.

### Complete online scale and moist-structure diagnostics (2026-07-27)

Desktop source:

`C:/Users/jiasenzhang/Desktop/淇敼绋嬪簭/Rotating_Moist_RB3-ns-deltab727`.

The `TPIN` diagnostics are separated from `avgvar.f90`:

- `aggregation_scale_diagnostics.f90` writes `lc_w(z~0.5)`, the
  stretched-grid vertical average of plane-wise `lc_w(z)`, `Lm`, `Lq`,
  `Lzeta`, shell-peak wavelengths for `w`, `m`, and `q`, condensation-cell
  count/fraction, and `mprime^2` at `z/H~0.25,0.5,0.75` plus its
  stretched-grid volume mean.
- `mse_vortex_diagnostics.f90` independently tracks periodic eight-connected
  components of the vertically weighted `mprime` field at thresholds
  `mprime2d>0` and `mprime2d>sigma_mprime2d`. It writes persistent IDs,
  birth/merge/breakup/death status, maximum radii, structure ratios, and
  per-component mean velocity, mean/rms vertical vorticity, and mean
  `mprime`. Binary tracker state files preserve IDs across continuation runs.
- Starting from the `vortex_lc` update, tracked components can be filtered by
  the user-specified linear-instability length scale. The value is read from
  the final `bou.in` movie/diagnostic row as `vortex_lc`; if `vortex_lc>0`,
  both `mprime2d>0` and `mprime2d>sigma_mprime2d` component lists keep only
  components with equivalent radius `radius > 2*vortex_lc`. If
  `vortex_lc=0`, the radius filter is disabled. Do not auto-fill this value
  from `k_peak`; the user currently wants to type it explicitly in `bou.in`.

The integral spectral length is

`L=sum(E)/sum(k_h E)`,

with the zero mode excluded and the real-to-complex Hermitian multiplicity
restored. Peak wavelengths use radially shell-summed spectra and
`Lpeak=2*pi/kpeak`.

Movie output is split into:

- `fieldNNNNN.h5`: `DSAL_me`, `VZ_me`, `RH_me`;
- `horizontal_velocityNNNNN.h5`: `VX_me`, `VY_me`;
- `mprimeNNNNN.h5`: per-height `MPRIME_me`;
- `condensationNNNNN.h5`: `COND_me`.

The Q-criterion movie is disabled. Every HDF5 file has a matching XMF file.

Remote compile and deployment provenance:

- remote source and archive root:
  `/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case/Rotating_Moist_RB3-ns-deltab727`;
- compiled with Intel oneAPI 2023.2, `mpi/latest`,
  `fftw/3.3.10-intel-2023`, and `hdf5/1.10.6-intel-2023`;
- `simexec` SHA-256:
  `194ef95272454916111b9491bee0f142feed989ae206ac883b68389e62ff9c6c`;
- deployed on 2026-07-27 to all 39 discovered `rotating_case` run
  directories containing `bou.in`;
- each prior executable is retained as
  `simexec_before_deltab727_20260727`;
- deployment audit:
  `/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case/deltab727_deployment_report_20260727.json`;
- the portable source ZIP SHA-256 is
  `8b8701e49424fab7823a49af54795ed28177ef084ed0ad77e54432c8f4b391de`.

Replacing `simexec` does not alter an executable already loaded by a running
MPI process. The new diagnostics start when a case is newly launched or
continued with the deployed executable.

### `vortex_lc` executable and case-input refresh (2026-07-28)

The active `rotating_case/latest_program/source` version reads an optional
fourth value from the final `bou.in` `mov_zcut_k` row:

`mov_zcut_k  tframe_me  stat_me  vortex_lc`.

The compiled executable at
`/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case/latest_program/source/simexec`
was deployed to the 39 discovered `rotating_case` run directories on
2026-07-28. Its SHA-256 was
`21b841ae38fb26a1641a05c3c794d7f224c19794581bbc3ae433ac9346777cfc`.
Prior executables were backed up as `simexec_before_latest_20260728_155044`.
Deployment audit:
`/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case/latest_simexec_update_report_20260728_155044.json`.

The same 39 existing `bou.in` files were updated from the three-value
diagnostic row to the four-value row by appending `vortex_lc=0`. Prior inputs
were backed up as `bou.in_before_vortex_lc_20260728_155331`. Audit:
`/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case/bouin_vortex_lc_update_report_20260728_155331.json`.

Local maintenance scripts:

- `E:/moist RB/rotating_case_inventory/update_simexecs.ps1`
- `E:/moist RB/rotating_case_inventory/update_bouin_vortex_lc.ps1`

Use `update_bouin_vortex_lc.ps1 --vortex-lc VALUE --set-existing` only when
the user explicitly wants to assign the same `vortex_lc` to all existing cases.

### Running deltab727 diagnostic snapshot (2026-07-27)

Case:

`/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case/Pr0p7/Ra8e6/Ek5e-3/AR16/Beta1p02/qbot0p5_qtop0p004978/N257x257x65/run`.

The diagnostic tables and tracking histories were reduced at their common
latest time `t=67.2`, giving 672 `TPIN=0.1` samples. The case was still
running, so this snapshot must not be interpreted as a statistically
saturated state.

Remote reduced tables:

`/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case/remote_profile_exports/deltab727_Ra8e6_Ek5e-3_AR16_qbot0p5_current_20260727`.

Desktop plots, reduced CSV files, metadata, and plotting script:

`C:/Users/jiasenzhang/Desktop/淇敼绋嬪簭/Ra8e6_Ek5e-3_AR16_qbot0p5_deltab727_diagnostics`.

At `t=67.2`, representative endpoint values are `lc_w(z~0.5)=0.09488`,
vertically averaged `lc_w=0.08896`, `Lm=0.17963`, `Lq=0.17973`,
`Lzeta=0.06968`, `Rmax_positive=6.52135`, `Rmax_core=0.93759`,
condensation fraction `0.01272`, and volume-mean `mprime^2=0.002127`.

An updated snapshot downloaded while the case was running was aligned to the
common time `t=200.2` (2026-07-27). Its desktop folder is

`C:/Users/jiasenzhang/Desktop/淇敼绋嬪簭/Ra8e6_Ek5e-3_AR16_qbot0p5_deltab727_diagnostics/update_t200p2`.

At that time, `lc_w(z~0.5)=0.11060`, vertically averaged `lc_w=0.10954`,
`Lm=0.27923`, `Lq=0.29525`, `Lzeta=0.09606`, and all three shell-peak
wavelengths were `2.66667`. The full-volume kinetic energy
`K=0.5*<u^2+v^2+w^2>` was `0.0013241`; volume-mean `mprime^2` was
`0.0029859`. The case remained active, so these are running values rather
than final saturated statistics.

### Dry rotating cases (2026-07-28)

Two dry `Ra=8e6`, `Pr=0.7` cases were read to `t=600` from:

- `/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case/Pr0p7/Ra8e6/dry/Ek3p42e-5`;
- `/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case/Pr0p7/Ra8e6/dry/Ek5p66e-5`.

Desktop copies, plots, and fit summaries are in
`C:/Users/jiasenzhang/Desktop/淇敼绋嬪簭/dry_cases_20260728`.
Using `K=0.5*<u^2+v^2+w^2>`, clean exponential-energy fits gave
`140<=t<=195`, `s_E=0.24712` (`R^2=0.99991`) for `Ek=3.42e-5`, and
`6<=t<=14`, `s_E=0.91190` (`R^2=0.99596`) for `Ek=5.66e-5`.

### Ek=5e-3 t=10 high-wavenumber audit (2026-07-28)

The `AR=16`, `256x256x64` moist case frame
`movie/field00001.h5` (`t=10`) was audited at `z/H=0.5094`. Results and
figures are in
`C:/Users/jiasenzhang/Desktop/淇敼绋嬪簭/Ek5e-3_t10_highk_analysis`.

The online value `lc=0.0656365` was reproduced exactly. It corresponds to
`2*pi*lc=0.412406`, only 6.60 horizontal grid intervals, whereas the
shell-peak wavelength was `0.64` (10.24 grid intervals). Modes with radial
index `n_h>32` carried 49.57% of `w` energy and 71.85% of the
`sum(k_h E)` denominator, demonstrating strong high-wavenumber bias.

Historical note for the executable used to create this `t=10` field:
the then-deployed `inqpr.f90` initialized
`dsal=0.1*sin(pi*z/H)*Gaussian_grid_white_noise`, replaces rather than
perturbs a background profile, and then sets `qvap=qs`. For this case the
random buoyancy RMS near mid-height is much larger than the imposed
boundary difference (`0.02`). This is the primary early high-wavenumber
source; saturated initialization and `tau_cond=0.001` can reinforce it.
This is not the current initialization after the drizzle update documented
below.

### Case-specific drizzle initialization (authoritative from 2026-07-28)

For all new `rotating_case` runs with `nread=0`, the active initialization is
now a case-specific one-dimensional conductive/drizzle state generated with
the user-supplied
`transfer_linear_stability_Ra8e6_20260728/stability_solver.py`,
function `moist_base_state`.

The solved base state is

\[
\frac{d^2b_0}{dz^2}+\gamma C_0=0,\qquad
S_m\frac{d^2q_0}{dz^2}-C_0=0,
\]

\[
q_{s0}=\exp[\alpha_{qs}(b_0-\beta_{qs}z)],\qquad
C_0=\frac{H_{\rm smooth}(q_0-q_{s0})(q_0-q_{s0})}{\tau_{\rm cond}},
\]

with the exact `bou.in` scalar boundary values. The nonlinear solve must
converge; the linear-profile fallback built into the supplied helper is
explicitly rejected for DNS initialization.

The DNS initial fields are

\[
b(x,y,z,0)=b_0(z)+10^{-4}\sin(\pi z/H)\,\xi(x,y,z),
\qquad
q(x,y,z,0)=q_0(z),
\qquad
\boldsymbol{u}(x,y,z,0)=0,
\]

where \(\xi\) is zero-mean, unit-variance Gaussian noise. The perturbation is
therefore added to drizzle buoyancy, not to a linear buoyancy profile, and
moisture is not reset to local saturation after perturbing buoyancy.

The authoritative workflow is explicitly two-step. Before every new
`nread=0` launch, enter that case's `run` directory and execute
`./prepare_drizzle_initial_condition.sh`. This script rereads the current
`bou.in`, calls the supplied `moist_base_state`, rejects its linear fallback,
validates the result, and atomically overwrites `drizzle_init.dat` and
`drizzle_init_meta.json`. Only after that should `simexec` be run or
submitted. The case builder does not silently solve drizzle at case-creation
time.

Per-case tools in each `run` directory:

- `prepare_drizzle_initial_condition.sh`: manual one-command generator;
- `generate_drizzle_initial_condition.py`: parses the current `bou.in`;
- `stability_solver.py`: unchanged supplied solver source;
- `check_drizzle_before_submit.py`: compares the profile against `bou.in`;
- `drizzle_init.dat`: generated profile and perturbation amplitude;
- `drizzle_init_meta.json`: nonlinear residuals and provenance.

Both `subjob.sh` and `submit_after_check.sh` validate the existing profile
before launch/submission, but neither recomputes it. `inqpr.f90`
independently verifies all control parameters and scalar endpoints before
reading the profile. A missing or mismatched profile causes a hard stop for
`nread=0`. Existing continuations with `nread=1` still read restart fields.

Active local workflow:

- root: `E:/moist RB/rotating_case_inventory`;
- builder:
  `E:/moist RB/rotating_case_inventory/01_case_builder_remote_upload/create_rotating_case.py`;
- manual generator staged in every case:
  `E:/moist RB/rotating_case_inventory/01_case_builder_remote_upload/prepare_drizzle_initial_condition.sh`;
- DNS initialization source:
  `E:/moist RB/rotating_case_inventory/00_latest_program/source/fluid_solver/inqpr.f90`.

The remote generator and code passed a direct two-step `nread=0` test using
`Ra=8e6`, `Pr=0.7`, `beta=1.02`, `qbot=0.5`, `qtop=0.004978`, and
`17x17x65`. The manual generator solved and checked the current `bou.in`;
the subsequent DNS startup reported 401 profile points, perturbation
amplitude `1e-4`, saturation width `1e-8`, zero initial divergence, and
completed normally. A second `beta=1.2`, `qbot=1` case also converged, with
maximum reported nonlinear residual `2.67e-7`.

Active compiled executable:

- path:
  `/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case/latest_program/source/simexec`;
- SHA-256:
  `08527afd2ebb084c84a5e47b1d6f6c63e18f180a439894868b7a35220d7c4d77`.

All 39 existing `rotating_case` run directories were updated with the manual
generator, supplied solver, checker, and this executable; no job was
submitted. Existing profile files were not treated as authoritative: each
new `nread=0` launch must rerun the manual generator. Remote audit:

`/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case/manual_drizzle_workflow_update_report_20260728.json`.

Cleanup on 2026-07-28 removed the superseded per-case Python wrapper,
renamed reference solver, three generations of `simexec`/`subjob` backups,
the temporary validation cases, and the two superseded update reports.
Current simulation data, active profile files, required four-tool workflow,
active executable, and final audit report were retained.
The final strict cleanup also removed three still older executable copies,
all `bou.in` backups, editor swap files, `fort.97`, and Python bytecode
caches from the 39 cases. It removed 207 additional entries and released
297.61 MiB. Restart fields, diagnostics, scheduler logs, postprocessing
results, and current drizzle provenance were retained.

### Ek=5e-3 vortex-radius snapshot through t=800 (2026-07-28)

The deltab727 moist case vortex histories and online maximum-radius table
were snapshotted at `t=800`. Local data, lifetime tables, filtered tracks,
figures, and the reproducible plotting script are in:

`C:/Users/jiasenzhang/Desktop/淇敼绋嬪簭/Ek5e-3_vortex_radius_t800`.

For individual-ID radius tracks, only positive-area records are used and
the lifetime is defined by
`max(time with area>0)-min(time with area>0)`. IDs with lifetime below 50
are excluded. This leaves 85 `mprime>0` IDs and 107
`mprime>sigma_mprime` IDs. Gaps longer than 1.5 output intervals are drawn
as breaks rather than connected lines. The online `Rmax_positive(t)` and
`Rmax_core(t)` curves remain unfiltered instantaneous maxima over all
structures, consistent with their original definitions.

The matching scale files were also snapshotted through `t=800` and plotted
without smoothing or time averaging. Outputs are under
`C:/Users/jiasenzhang/Desktop/淇敼绋嬪簭/Ek5e-3_vortex_radius_t800/scale_time_series`.
The endpoint values are `lc_w(z~0.5)=0.12358724`, vertically averaged
`lc_w=0.10666333`, `Lm=0.35509426`, `Lq=0.34045042`,
`Lzeta=0.14894177`, and `Lpeak_w=Lpeak_m=Lpeak_q=4.0`.

The matching MSE-anomaly variance and full-volume kinetic-energy series are
under
`C:/Users/jiasenzhang/Desktop/淇敼绋嬪簭/Ek5e-3_vortex_radius_t800/amplitude_time_series`.
At `t=800`,
`<mprime^2>_xy(z~0.25)=0.00491838689`,
`<mprime^2>_xy(z~0.50)=0.00306312916`,
`<mprime^2>_xy(z~0.75)=0.00248005173`, and the vertical average is
`0.00320453823`. With
`K=0.5*<u^2+v^2+w^2>_V`, `K(t=800)=0.00299059579`; the full-series maximum
is `0.00746858150` at `t=10.8`, during the strong initialization transient.

## Current high-resolution time-series collection (2026-08-01)

The current `rotating_case` tree was rescanned on host `c01n0006`. The
selection uses the actual `N1,N2,N3` row in each `bou.in`, not the grid
directory name. A case is included when both horizontal dimensions are at
least 257. If identical physical parameters have several resolutions, retain
the highest actual grid and use the longest physical-time branch only as an
equal-grid tie breaker. Base and continuation segments are stitched by
physical time and duplicate times retain the later segment.

Current local output (refreshed 2026-08-02):

`E:/moist RB/rotating_case_inventory/04_outputs_and_figures/high_resolution_timeseries_latest_program_20260802`

Reproducible script:

`E:/moist RB/rotating_case_inventory/02_inventory_and_plot_scripts/plot_highres_energy_mprime_scales.py`

The 2026-08-02 latest-program-only scan selected 13 moist cases, all at
`Ra=8e6`, `Pr=0.7`, `beta=1.02`, and `qbot=0.5`. They are the nonrotating
case and rotating cases at `Ek=1.5e-4,2e-4,5e-4,7e-4,1e-3,3e-3,5e-3,
7e-3,1e-2,3e-2,5e-2,1e-1`. The two `AR=4` strong-rotation cases remain on
the actual `257x257x65` grid. The nonrotating case and all `AR=16` cases from
`Ek=5e-4` through `1e-1` use the actual `385x385x65` grid recorded in
`bou.in`, even where an older grid-directory label says 257. Base and all
continuation branches are merged by physical time before plotting.

Unified quantities and sources are:

- `K=0.5*avgvar[:,3]` after the established finite and `K<0.1` checks;
- MSE variance `<<mprime^2>xy>z` from the stretched-grid-weighted column of
  `diagnostics/thermo/mprime_square.dat`, or the equivalent `A_m` column of
  `data/mse_aggregation.out` for the earlier executable;
- horizontal-velocity `l0(z~0.5)` and `l0(z~0.75)` from one-based columns 8
  and 9 of `avgvar.out`;
- MSE peak wavelength from `peak_scale.dat` or the equivalent `L_peak`
  column of `mse_aggregation_scales.out`;
- MSE spectral wavelength `2*pi*sum(E_m)/sum(k_h E_m)` from
  `moist_integral_scale.dat` or the equivalent `2*pi*ell_m` column of
  `mse_aggregation_scales.out`;
- the newer `2*pi*lc_w(z~0.5)` and vertically averaged `2*pi*lc_w` are kept
  in separate panels because they are distinct height diagnostics.

The output folder contains PNG/PDF figures, selected-case metadata, excluded
lower-resolution metadata, a long-form reduced CSV, and a JSON manifest. No
full HDF5 fields were transferred locally.

Every refresh also writes two fixed-style single-case figures under
`individual_cases/`: (1) a dual-y-axis plot with red `K(t)` on the left axis
and blue `<<mprime^2>xy>z(t)` on the right axis, both axes remaining black
and the curves identified by a frameless in-panel legend, and (2) a blue
`2*pi*Lm(t)` plot by itself. Single-case figures have no parameter header.
The directory contains an
`individual_case_figure_index.csv` that records parameters, grid, time
coverage, and both PNG paths. These single-case products use the same merged
base/continuation histories and the permanent `Ek=7e-3` repair described
below.

### Permanent Ek=7e-3 history repair

The `Ra=8e6`, `Pr=0.7`, `Ek=7e-3`, `AR=16`, `beta=1.02`, `qbot=0.5`,
actual `385x385x65` case lost its remote diagnostic history between the early
run and the continuation after an accidental file operation. The current
remote histories jump from `t=100` to `t=1000.1`.

For every future refresh by
`E:/moist RB/rotating_case_inventory/02_inventory_and_plot_scripts/plot_highres_energy_mprime_scales.py`,
restore `10<=t<=1000` from the immutable reduced table:

`E:/moist RB/rotating_case_inventory/03_inventory_tables/historical_time_series_patches/Ra8e6_Ek7e-3_t10_1000_latest_program.csv`

The patch contains kinetic energy, MSE variance, both `l0` heights, MSE peak
and spectral scales, and the midheight and vertically averaged convective
scales at `dt=0.1`. Merge the historical rows first and the current remote
rows second, so the current remote value always wins at a duplicate physical
time. Record the patch path, interval, and number of inserted rows in the
output manifest and selected-case metadata. Do not extrapolate or interpolate
this missing interval.

## MSE-variance spectrum and Gaussian scale flux (2026-08-02)

For transfer analyses, use the plane anomaly

\[
m'=m-\langle m\rangle_{xy}(z,t),\qquad m=b+\gamma q,
\]

and the variance density `mprime^2/2`. In the current equal-diffusivity model,
condensation cancels from `m`, giving

\[
D_t m=\kappa\nabla^2m,\qquad \kappa=(RaPr)^{-1/2}.
\]

The approved coarse-grained diagnostic follows the project's previous
Gaussian-flux convention:

\[
G_\ell(k_h)=\exp(-\ell^2k_h^2/2),\qquad k_c=1/\ell,
\]

\[
\Pi_m^\ell=-\langle\tau_i^{(\ell)}\partial_i m'_\ell\rangle_V,
\qquad
\tau_i^{(\ell)}=(u_i m')^\ell-u_i^\ell m'_\ell.
\]

The filter is horizontal, while the contraction includes all three velocity
components and all three spatial derivatives. `Pi_m>0` denotes forward
(downscale) transfer and `Pi_m<0` denotes inverse (upscale) transfer. Keep the
Gaussian cutoff `k_c=1/ell` distinct from shell-center `k_h`; compare trends,
not points one to one.

The first full application used the latest ten snapshots of seven current
`Ra=8e6`, `qbot=0.5` cases: `Ek=1.5e-4, 2e-4, 3e-3, 5e-3, 7e-3, 1e-2,
3e-2`. The AR4 cases used movie grids `256x256x64`; the AR16 cases used
`384x384x64`. Snapshot diagnostics were computed first and then averaged.
No raw HDF5 fields were copied locally.

Remote reduced output:

`/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case/postprocess/mse_variance_gaussian_flux_20260802/final_latest10`

Local data, figures, and interpretation:

`E:/moist RB/rotating_case_inventory/04_outputs_and_figures/Ra8e6_mse_variance_gaussian_flux_20260802`

Reproducible scripts:

- `E:/moist RB/rotating_case_inventory/02_inventory_and_plot_scripts/compute_mse_variance_gaussian_flux_remote.py`;
- `E:/moist RB/rotating_case_inventory/02_inventory_and_plot_scripts/plot_mse_variance_gaussian_flux_latest10.py`.

Current defensible finding: the mature intermediate-rotation cases
`Ek=3e-3--1e-2` show negative low-`k_c` flux followed by positive high-`k_c`
flux. This supports upscale MSE-variance organization at large scales and
downscale mixing at small scales. `Ek=3e-2` is positive at every resolved
cutoff, so its box-scale spectral peak must not be called an inverse cascade
without other evidence. The `Ek=1.5e-4` signal is only order `1e-11`, and
`Ek=2e-4` remains transient; neither is sufficient for a stationary cascade
claim. The Gaussian budget closure residual is about 5.5--8.2 percent for the
mature finite-amplitude cases, and reconstructed `mprime` matches stored
`MPRIME_me` to about `1e-14`.

## Vertically averaged MSE mode and shell-to-shell matrix (2026-08-02)

Do not describe the full-field Gaussian MSE flux above as a vertically
averaged diagnostic. It keeps `mprime(x,y,z)` at every height and volume
averages only after forming the three-dimensional nonlinear transfer.

When a vertically coherent scalar mode is requested, use stretched-grid
cell-thickness weights:

\[
M=\langle m'\rangle_z=\sum_jw_jm'(z_j),\quad
U=\sum_jw_ju(z_j),\quad V=\sum_jw_jv(z_j),
\]

where `sum_j w_j=1`. Never replace these weights by an equal average over the
stored `z` planes. The two-dimensional scalar self-transfer is

\[
T_m^{2D}(K,P)=-\frac{1}{\Delta k}
\langle M_K(U\partial_x+V\partial_y)M_P\rangle_{xy},
\]

with donor `P` on matrix columns and receiver `K` on rows. Positive transfer
above the diagonal (`K>P`) is downscale; positive transfer below the diagonal
(`K<P`) is upscale. Plot the conservative exchange matrix

\[
T_{m,ex}^{2D}=\frac{T_m^{2D}-(T_m^{2D})^T}{2}
\]

and place `sum_P T_m,ex^{2D}(K,P)` in the right-hand marginal panel.

The mode is not closed. The exact vertically averaged nonlinear input is

\[
N_M=-\nabla_h\cdot\langle\boldsymbol u_hm'\rangle_z,
\]

and

\[
F_{3D\rightarrow2D}(K)=N_M(K)-\sum_PT_m^{2D}(K,P)
\]

must be retained. This term is the forcing by vertically varying velocity and
MSE fluctuations; it is not a scalar donor shell in the `T_2D` matrix.

The first calculation used `Delta k=pi/2`, a common physical-wavenumber edge
at `k=59.6903`, and the latest ten fields of the seven current `Ra=8e6`,
`qbot=0.5` cases. Matrices were computed frame first and then averaged. For
the mature `Ek=3e-3--1e-2` cases, the vertically coherent self-transfer is
predominantly downscale, while the maximum `3D-to-2D` forcing is about 19--33
times the maximum self-transfer net input. This does not contradict the
negative low-cutoff full-field Gaussian flux because the two diagnostics act
on different projections.

Remote reduced output:

`/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case/postprocess/mse_shell_to_shell_20260802/final_2dmode_latest10`

Local data, figures, notes, and scripts:

`E:/moist RB/rotating_case_inventory/04_outputs_and_figures/Ra8e6_mse_2dmode_shell_to_shell_20260802`
# Latest-program moist Nusselt and force-balance diagnostics (2026-08-05)

- `data/nu_profiles.out` column 4 is the raw moist-static-energy flux profile,
  not an already normalized Nusselt number. Compute
  `Nu_m(t)=trapz_z(F_m(z,t))/Delta_m`, with
  `Delta_m=(b_bot-b_top)+gamma*(q_bot-q_top)`. For the current Ra8e6,
  beta1.02, qbot0.5 cases, `Delta_m=0.5245242`; the conductive initial value
  is then exactly one.
- RRBC-style force magnitudes are horizontal-plane RMS values after removing
  each force component's horizontal mean. Use
  `F_I=-u dot grad(u)`, `F_C=-invRo zhat cross u`,
  `F_V=sqrt(Pr/Ra) laplacian(u)`, and `F_B=b' zhat`.
- The force local Rossby diagnostic is `Ro_F=F_I/F_C`. `Ro_F<1` establishes
  rotational constraint but does not by itself prove geostrophy. Direct
  geostrophic verification requires pressure to test Coriolis-pressure
  cancellation. Current movie output has no pressure, so use the wording
  "geostrophic-compatible" and preserve this caveat.
- Current reduced force statistics use one latest common 3-D movie snapshot
  and a bulk integral over `0.2 <= z <= 0.8`. Nonstationary cases remain
  provisional. CIA requires geostrophic leading order plus an ageostrophic
  inertia-buoyancy balance with viscosity subdominant; do not infer CIA from
  `Ro_F<1` alone.
