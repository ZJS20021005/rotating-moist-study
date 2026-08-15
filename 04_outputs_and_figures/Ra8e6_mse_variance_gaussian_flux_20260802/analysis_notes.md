# Ra=8e6 moist-static-energy variance transfer

## Scope and provenance

- Cases: `Ek=1.5e-4, 2e-4, 3e-3, 5e-3, 7e-3, 1e-2, 3e-2`.
- Common parameters: `Ra=8e6`, `Pr=0.7`, `beta=1.02`, `gamma=1.1`, `qbot=0.5`.
- Each result is computed snapshot first and then averaged over the latest ten full-field outputs.
- The remote reduction used the full three-dimensional fields; only reduced NPZ/CSV products were copied locally.
- Remote output: `/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case/postprocess/mse_variance_gaussian_flux_20260802/final_latest10`.

Time windows:

| Ek | AR | movie grid | time window |
|---:|---:|---:|---:|
| 1.5e-4 | 4 | 256x256x64 | 1110--1200 |
| 2e-4 | 4 | 256x256x64 | 510--600 |
| 3e-3 | 16 | 384x384x64 | 1520--1610 |
| 5e-3 | 16 | 384x384x64 | 1520--1610 |
| 7e-3 | 16 | 384x384x64 | 1520--1610 |
| 1e-2 | 16 | 384x384x64 | 1520--1610 |
| 3e-2 | 16 | 384x384x64 | 1530--1620 |

## Definitions

The conserved moist scalar is

\[
m=b+\gamma q,
\qquad
m'=m-\langle m\rangle_{xy}(z,t).
\]

For the present equal-diffusivity cases, condensation cancels from the `m`
equation:

\[
D_t m=\kappa\nabla^2m,
\qquad
\kappa=(RaPr)^{-1/2}.
\]

The scalar-variance density is `m'^2/2`. Horizontal Gaussian filtering uses

\[
G_\ell(k_h)=\exp\!\left(-\frac{\ell^2k_h^2}{2}\right),
\qquad k_c=\ell^{-1}.
\]

The direct subfilter MSE-variance flux is

\[
\Pi_m^\ell=-\left\langle\tau_i^{(\ell)}\,\partial_i m'_\ell\right\rangle_V,
\qquad
\tau_i^{(\ell)}=(u_i m')^\ell-u_i^\ell m'_\ell.
\]

Sign convention:

- `Pi_m > 0`: resolved MSE variance is transferred to scales smaller than `ell` (forward/downscale).
- `Pi_m < 0`: MSE variance is transferred to scales larger than `ell` (inverse/upscale).

This is a horizontal Gaussian cutoff but the contraction uses all three
velocity components and all three derivatives. Therefore it is a
three-dimensional MSE-variance transfer diagnosed across horizontal scale.

## Main result

The mature intermediate-rotation cases `Ek=3e-3, 5e-3, 7e-3, 1e-2` have a
continuous negative low-wavenumber interval followed by positive flux at
higher wavenumber. They therefore show simultaneous upscale organization of
large-scale MSE variance and downscale mixing at smaller scales.

The strongest negative value occurs for `Ek=3e-3`,
`min(Pi_m)=-1.64e-5`. The negative interval narrows and weakens toward
`Ek=1e-2`.

For `Ek=3e-2`, `Pi_m` is positive at every resolved cutoff even though the
MSE spectrum peaks at the box fundamental. Its box-scale MSE organization is
therefore not evidence of an inverse MSE-variance cascade; direct low-k
production is the more consistent interpretation.

The `Ek=1.5e-4` flux is only order `1e-11`, so its negative sign is below the
amplitude needed for a physical cascade claim. The `Ek=2e-4` case is still
evolving and has a non-negligible variance tendency relative to its tiny
flux. Neither strong-rotation case is used to establish a steady cascade.

## Numerical checks and caveats

- Reconstructed `m'` agrees with stored `MPRIME_me` to about `1.2e-14` or better.
- Maximum relative Gaussian-budget residual is about 0.9--8.2 percent across the seven cases and about 5.5--8.2 percent for the mature finite-amplitude cases.
- `k_c=1/ell` on the Gaussian-flux figure is not the same coordinate as the shell-center `k_h` on the spectrum. Compare their trends, not individual points one to one.
- Aspect ratio changes the fundamental horizontal wavenumber: AR4 starts at `k_h=pi/2`, whereas AR16 starts at `k_h=pi/8`. This matters when comparing the two strongly rotating cases with the five AR16 cases.
- The normalized-flux figure is supplementary only: normalizing the nearly zero `Ek=1.5e-4` curve amplifies noise.

## Reproducible scripts

- Remote reduction: `E:/moist RB/rotating_case_inventory/02_inventory_and_plot_scripts/compute_mse_variance_gaussian_flux_remote.py`.
- Local publication plotting: `E:/moist RB/rotating_case_inventory/02_inventory_and_plot_scripts/plot_mse_variance_gaussian_flux_latest10.py`.

