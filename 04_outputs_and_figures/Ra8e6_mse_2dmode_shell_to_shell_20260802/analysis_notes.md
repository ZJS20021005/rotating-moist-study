# Vertically averaged MSE mode and shell-to-shell transfer

## What was vertically averaged

The previously calculated Gaussian MSE-variance flux was **not** based on a
vertically averaged scalar. It retained `mprime(x,y,z)` at every height and
then volume averaged the three-dimensional nonlinear transfer.

The present shell-to-shell analysis introduces a separate vertically
coherent scalar mode:

\[
m'(x,y,z,t)=m-\langle m\rangle_{xy}(z,t),
\]

\[
M(x,y,t)=\langle m'\rangle_z=\sum_j w_jm'(x,y,z_j,t).
\]

The weights `w_j` are cell-thickness weights reconstructed from the actual
stretched vertical grid and normalized to sum to one. This is not an equal
average over the stored planes. The corresponding two-dimensional velocity
mode is

\[
U=\sum_jw_ju(z_j),\qquad V=\sum_jw_jv(z_j).
\]

## Shell transfer

Horizontal Fourier modes are grouped into common physical-wavenumber bins
with `Delta k=pi/2` and edges from zero through `k=59.6903`. For scalar donor
shell `P` and receiver shell `K`,

\[
T_m^{2D}(K,P)=-\frac{1}{\Delta k}
\left\langle M_K(U\partial_x+V\partial_y)M_P\right\rangle_{xy}.
\]

Positive values mean that donor `P` supplies MSE variance to receiver `K`.
The heat maps plot the conservative exchange part

\[
T_{m,ex}^{2D}(K,P)=\frac{T_m^{2D}(K,P)-T_m^{2D}(P,K)}{2}.
\]

With donor `P` on the horizontal axis and receiver `K` on the vertical axis:

- positive values above the diagonal (`K>P`) are downscale transfer;
- positive values below the diagonal (`K<P`) are upscale transfer.

The right-hand curve is `sum_P T_m,ex^{2D}(K,P)`, the net conservative
transfer gained by each receiver shell.

## Why this is not a closed two-dimensional scalar

The exact vertically averaged nonlinear scalar flux is

\[
N_M=-\nabla_h\cdot\langle\boldsymbol u_hm'\rangle_z.
\]

Writing `u_h=(U,V)+u_h_doubleprime` and `mprime=M+m_doubleprime` gives

\[
\langle\boldsymbol u_hm'\rangle_z=(U,V)M+
\langle\boldsymbol u_h''m''\rangle_z.
\]

Therefore the `T_2D` matrix describes only self-transfer inside the
vertically coherent modes. The remaining shell forcing

\[
F_{3D\rightarrow2D}(K)=N_M(K)-\sum_PT_m^{2D}(K,P)
\]

is the effect of vertically varying three-dimensional fluctuations on `M`.

## Results

All matrices were calculated independently in each of the latest ten fields
and then averaged. The plotted exchange matrices are almost antisymmetric;
the mean-matrix relative residual is about 0.05--3.7 percent across the seven
cases and below 1.5 percent for `Ek=1.5e-4` through `1e-2`.

For the mature finite-amplitude cases `Ek=3e-3, 5e-3, 7e-3, 1e-2`, the main
positive band lies above the diagonal. Their vertically coherent MSE mode
therefore transfers variance predominantly downscale, not upscale. The
positive-upscale to positive-downscale transfer ratios are respectively
`0.058`, `0.0053`, `0.0136`, and `0.0111`.

More importantly, the maximum `3D-to-2D` forcing is approximately 19--33
times the maximum self-transfer net input in these four cases. Thus `M` is
not dynamically closed like an autonomous barotropic passive scalar. Its
large-scale variance budget is controlled mainly by coupling to vertically
varying fluctuations and by the mean-gradient/diffusive terms not contained
in the self-transfer matrix.

This does not contradict the negative low-cutoff full-field Gaussian flux:
that flux includes every height-dependent MSE mode, whereas this matrix is a
projection onto `M=<mprime>_z`. The projection introduces an explicit
`3D-to-2D` exchange term and can have a different transfer direction.

The `Ek=1.5e-4` and `2e-4` matrix amplitudes are only order `1e-11`; their
normalized shapes must not be treated as finite-amplitude cascade evidence.

## Provenance

- Remote reduced data: `/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case/postprocess/mse_shell_to_shell_20260802/final_2dmode_latest10`.
- Local results: `E:/moist RB/rotating_case_inventory/04_outputs_and_figures/Ra8e6_mse_2dmode_shell_to_shell_20260802`.
- Remote calculation script: `E:/moist RB/rotating_case_inventory/02_inventory_and_plot_scripts/compute_mse_2dmode_shell_to_shell_remote.py`.
- Core transfer routines: `E:/moist RB/rotating_case_inventory/02_inventory_and_plot_scripts/compute_mse_shell_to_shell_remote.py`.
- Plotting script: `E:/moist RB/rotating_case_inventory/02_inventory_and_plot_scripts/plot_mse_2dmode_shell_to_shell_latest10.py`.

