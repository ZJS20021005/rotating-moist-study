# Core references for moist rotating convection

Use these as the default theory/history anchors when interpreting results.

For exact diagnostic definitions, averaging conventions, data-source paths, and current output tables, use and update `diagnostics_and_data_sources.md`.

- Vallis, Parker & Tobias (2019): the original Rainy-Bénard model paper.
  - Key takeaways:
    - introduces the moist Boussinesq model with condensational heating;
    - defines moist static energy `m = b + gamma q`;
    - shows the drizzle solution and its stability;
    - in the non-rotating runs, kinetic energy grows with Rayleigh number in the regime studied there.

- Desktop reference library:
  - `C:\Users\jiasenzhang\Desktop\moist rotating convection\`
  - Treat this folder as the working reference set for rotating-moist convection papers, especially for:
    - rotating Rayleigh-Bénard scaling,
    - Ekman pumping / boundary-layer interpretation,
    - transition stories and figure-style analogs,
    - moist convection literature beyond the original model paper.

## Literature map for the current story

- Chandrasekhar (1961), *Hydrodynamic and Hydromagnetic Stability*:
  linear rotating-convection onset and the classical horizontal onset length
  proportional to `Ek^(1/3)`.
- Julien et al. (2012), *Heat transport in low-Rossby-number
  Rayleigh-Benard convection*, Phys. Rev. E 85, 016313:
  asymptotic/geostrophic framework for rapidly rotating heat transport and
  regime interpretation. Use it as a dry-RRBC mechanism reference, not as a
  direct proof that the same moist boundary-layer criterion is sufficient.
- Guervilly, Hughes & Jones (2014), *Large-scale vortices in rapidly rotating
  Rayleigh-Benard convection*:
  large-scale-vortex formation and inverse-transfer context. Distinguish a
  barotropic dynamical vortex from an MSE aggregation region.
- Aurnou et al. (2015), *Rotating convective turbulence in Earth and planetary
  cores*, PEPI:
  review of cellular, convective-Taylor-column, plume and geostrophic-turbulent
  regimes; useful for presentation background and regime vocabulary.
- Chien et al. (2022), *Hurricane-Like Vortices in Conditionally Unstable Moist
  Convection*:
  close phenomenological precedent for moist rotating vortex structures. The
  gap relevant to this project is the lack of a systematic mechanism/scaling
  account connecting rotation, transport, aggregation, vortex merger and
  regime transitions.

When invoking these references, separate three claims: what the paper
demonstrated, what the present DNS directly measured, and what remains a
mechanistic inference.

When a new result seems to contradict the expected trend, first check:

1. whether the plotted quantity is a raw kinetic-energy proxy or a normalized quantity,
2. whether the same no-rotation baseline was used,
3. whether the averaging window is comparable across cases,
4. whether the case is actually complete or still in a startup / continuation transient.
5. whether profile figures such as `m(z)`, `q(z)`, and `Nu_m(z)` were rebuilt from the same remote reduction pipeline as the `Ra=10^8` reference, rather than mixed with older local profile exports.

## Local Rossby number convention

Use `Ro_z^l` / `Rozl` for the local Rossby number based on a height-dependent vertical velocity scale and horizontal integral length scale:

\[
Ro_z^l(z)=\frac{\langle u_z^2\rangle_{x,y,t}^{1/2}}{2\Omega\,l_h(z)}
       =\frac{Ro_z(z)}{l_h(z)/d}.
\]

In the solver's free-fall nondimensionalization with domain height `d=1`,

\[
Re_z(z)=w_{\rm rms}(z)\sqrt{\frac{Ra}{Pr}},\qquad
Ro_z^l(z)=\frac{Re_z(z)\,Ek}{l_h(z)}.
\]

Here `l_h(z)` is the horizontal integral length at the same height, nondimensionalized by the layer depth. Do not call this quantity the convective Rossby number `Ro_c`; `Ro_c=sqrt(Ra Ek^2/Pr)` is a control parameter, while `Ro_z^l` is a measured local flow diagnostic.

## Integral-length convention

Unless the user explicitly asks for another height or a different diagnostic, the default horizontal length-scale figure should use the legacy `z=0.75` integral-length definition:

- `l_h(z=0.75)` from the inverse-wavenumber spectrum diagnostic based on the vertical velocity field `w` / `VZ_me`;
- do not replace it by a height-averaged `l_h` unless the user explicitly says to change the definition;
- when comparing to the `Ek^{1/3}` guideline, use the same `z=0.75` convention consistently across Ra cases.

The default time averaging must be frame-first, then average the length:

\[
l_h(z,t)=
\frac{\sum_{k_h>0} E_w(k_h,z,t)}
     {\sum_{k_h>0} k_h E_w(k_h,z,t)},\qquad
\overline{l_h}(z)=\langle l_h(z,t)\rangle_t .
\]

Do **not** silently replace this with either of the following unless explicitly requested:

- `1 / mean(l_h^{-1})`, where \(l_h^{-1}(z,t)=\sum k_hE_w/\sum E_w\);
- the energy-weighted time-averaged spectrum length \(\sum_t\sum E_w / \sum_t\sum k_hE_w\).

Those alternatives are useful diagnostics, but they are not the project's default `l_h`. This distinction matters near boundaries and when the spectral energy varies strongly in time.

The detailed provenance and output-table rules for `l_h`, \(l_{m'}\), \(Nu_m\), `Rozl`, correlations, profiles, and time series are maintained in `diagnostics_and_data_sources.md`; update that file whenever a definition or source changes.

If a figure uses vertical bracketing lines to indicate the `Ek^{1/3}`-matching interval, keep the same visual convention across all Ra panels and label the highlighted window as the scaling-consistent range rather than as a new fitted law.

For moist-structure scales, use a separate notation such as \(l_{m'}\), with \(m'=m-\langle m\rangle_{xy}\) and \(m=b+\gamma q\). Do not mix \(l_{m'}\) with the default velocity-based \(l_h\) in scaling plots unless the legend and caption explicitly distinguish them.

## Regime labels for phase diagrams

The user is building a morphology/regime classification for later phase diagrams. Preserve these labels exactly unless the user revises them:

- `cell`: `Ra=1e6, Ek=3e-4`.
- `cell`: `Ra=1e8, Ek=1e-5`.

Planned regime vocabulary includes `cell`, `intermittent plume burst`, `funnel`, and `plume`. When the user assigns more points, append them here and use them in later phase diagrams.

Physical interpretation note for `cell`: near-saturated relative humidity over much of the field does not imply strong motion. In the cell regime, the thermodynamic field can remain close to saturation while rotation/near-onset dynamics strongly constrain vertical velocity and moist-static-energy transport.

## Why Ra=1e6 does not show strong moist-transport enhancement

The ordering \(\delta_E<\delta_q\) is at most a necessary geometric condition, not a sufficient transport criterion. In the current Ra=1e6 analysis, \(\delta_q\) is the first lower-wall peak of q standard deviation. It is a moisture-fluctuation-layer scale, not the conductive mean-gradient thickness of the conserved moist static energy. Water vapour q is not separately conserved because condensation exchanges q and buoyancy; the transport diagnostic is controlled by \(m=b+\gamma q\).

For a horizontal plane with negligible mean vertical velocity,

\[
\langle wm\rangle = \mathrm{corr}(w,m)\,\sigma_w\sigma_m.
\]

Therefore a thick moisture-fluctuation layer or a large geometric correlation cannot enhance \(Nu_m\) when rotation suppresses \(\sigma_w\) and the flux amplitudes \(\langle wq\rangle\) and \(\langle wm\rangle\). At fixed Ek, Ra=1e6 also has a convective Rossby number ten times smaller than Ra=1e8 and is much less supercritical relative to rotating onset. It can form organized cells/columns without producing vigorous plume transport or strong condensation feedback.

For an RRBC-style boundary-layer mechanism test, compare \(\delta_E\) with a thickness derived from the mean m gradient, the diffusive/convective m-flux crossover, or a conserved \(Nu_m(z)\), rather than treating the q-standard-deviation peak as a direct analog of the dry thermal boundary layer.

## Raw kinetic energy across Ra

The solver uses free-fall nondimensionalization. Raw plotted kinetic energy is

\[
K^*=\frac12\langle |\boldsymbol{u}^*|^2\rangle,
\qquad
Re=\sqrt{2K^*}\sqrt{\frac{Ra}{Pr}}.
\]

Thus a larger dimensionless \(K^*\) at lower Ra does not mean larger dimensional velocity or Reynolds number. Because the free-fall velocity used for normalization increases with Ra, \(K^*\) can decrease with Ra even while Re and dimensional kinetic energy increase. Compare Re, or equivalently \(K^*Ra/Pr\), when comparing different Ra.
