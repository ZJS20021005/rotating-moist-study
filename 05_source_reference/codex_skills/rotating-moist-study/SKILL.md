---
name: rotating-moist-study
description: Research and workflow memory for the rotating moist Rainy-Benard project, including physical interpretation, diagnostic definitions, self-aggregation and vortex-scale analysis, phase diagrams, Nusselt/Re/Rossby/profile figures, spectra and energy transfer, remote data provenance, and completed-case bookkeeping. Use for any analysis, plotting, interpretation, remote reduction, case update, or research-story task in this project.
---

# Rotating Moist Study

## Overview

Use this skill for the recurring plotting and bookkeeping loop in the rotating moist RB transition study.

Its main job is to keep the completed-case figures synchronized with new remote output: phase diagram, `Num`/`Nuq`/`Nub`, `Re`, `Rol`, correlation, profiles, `l_h` scaling plots, and moist-convective self-aggregation diagnostics.

Load only the references needed for the current task:

- `references/diagnostics_and_data_sources.md`
- `references/research_story_and_current_findings.md`: physical knowledge, current evidence, open questions, and defensible conclusions.
- `references/analysis_method_selection.md`: choose diagnostics for transport, self-aggregation, vortex geometry, inverse cascade, boundary layers, and numerical resolution.
- `references/results_and_paths_inventory.md`: current local/remote result roots and the major processed result sets.
- `references/remote_connection.md`: verified SSH alias, credential boundary, connection tests, and troubleshooting for the active cluster host.
- `references/moist_rb_core_refs.md`: literature anchors and stable model interpretation.

## Core workflow

1. When new remote results appear, first identify which cases are truly complete.

- Prefer the latest summary CSVs and completed-case point tables.
- Exclude unstable or still-growing cases from the main completed-case phase diagram.
- If a `norotating` baseline exists, use it for normalization when the figure definition requires it.

2. Update the reduced phase diagram whenever completed cases change.

- Keep the diagram as a clean scatter plot unless the user explicitly asks for lines.
- Use the reduced moist Rayleigh number `R_m = Ra_m Ek^{4/3}`, with `Ra_m = Ra * Delta_m` and `Delta_m = gamma (q_bottom - q_top)`.
- Rebuild the figure from the completed-case CSV rather than manually editing points.
- Every time new data are finished, the phase diagram is the first plot to refresh.
- If the definition or source of any plotted quantity changes, update `references/diagnostics_and_data_sources.md` in the same turn.

3. For transport and stability summaries, keep normalization rules explicit.

- If a nonrotating reference exists, normalize `Nu_m`, `Re`, or related quantities by that baseline only when the plot definition calls for it.
- Do not silently change definitions between figures.
- Use the same case-color mapping across all related figures.

4. For profile figures, keep the established conventions.

- Use the project?s scientific plot style.
- Keep one color per case family when plotting profiles.
- Use the user?s requested height selection directly; do not invent a different diagnostic height.
- If the user asks for `m`, plot `m = b + gamma q`.
- If the user asks for `m'`, use `m-<m>_{xy}(z,t)` at every height unless the user explicitly requests a different anomaly.
- For velocity-based `l_h`, use the frame-first length average documented in `references/diagnostics_and_data_sources.md`. Do not silently use inverse-mean or time-spectrum-weighted alternatives.

5. For `Num`-type quantities, respect the definition.

- `Nuq`, `Nub`, and `Num` must be computed from the corresponding flux definitions, not by renaming another quantity.
- If the user asks whether a quantity is conserved, check the governing definition first; do not assume conservation from a similar-looking label.

6. For local Rossby / Reynolds / correlation figures, keep all completed cases in sync.

- `Re` and `Rol` should be built from the full velocity quantity the user requested, not a vertical-only surrogate unless explicitly asked.
- For correlation and Ekman-layer diagnostics, preserve the current analysis convention and annotate the interpretation carefully when it is an inference rather than a directly measured fact.

7. For self-aggregation, use the dedicated MSE-anomaly definitions.

- Use `m'=m-<m>_xy` separately at each height; do not use a volume-mean anomaly.
- The primary amplitude is the vertical average of the per-height horizontal MSE variance `A_m(t)`.
- In figures, do not label the y-axis only as `A_m` unless the user explicitly
  asks for the shorthand. Prefer the direct expression
  `\langle\langle m'^2\rangle_{xy}\rangle_z` for variance and
  `\langle\langle m'^2\rangle_{xy}\rangle_z^{1/2}` for the plotted standard
  deviation.
- Compute spectra, correlation lengths, periodic clusters, and periodic Voronoi cells snapshot-first from the full fields. Keep MSE aggregation scales separate from dynamical vortex radii.
- Treat Voronoi as a point-pattern diagnostic for isolated cores and their spacing. It is not the default size measure for broad weak-rotation aggregates containing many internal maxima. For those, use low-pass MSE, periodic connected components, radius of gyration, spectral length, and correlation length.
- Use the online and offline programs, output names, hashes, and deployment provenance in `references/diagnostics_and_data_sources.md`.

8. Keep claims tied to what the diagnostic can actually measure.

- Distinguish instability scale, vortex radius, funnel spacing, MSE aggregation scale, and box-limited condensate scale.
- A negative coarse-grained kinetic-energy flux supports inverse transfer only with the documented sign convention and over a resolved scale interval; a negative one-dimensional shell-transfer value alone is not sufficient.
- A moisture-standard-deviation peak is a fluctuation-layer thickness, not automatically the conserved MSE conductive boundary layer.
- Record unresolved alternatives and sensitivity tests rather than selecting a favorable threshold or replacement point.

## What to refresh when new data arrive

When a new case finishes, update in this order:

1. Completed-case phase diagram.
2. `Num` / `Re` / `Rol` / correlation summary plots.
3. Profile figures that compare cases across `Ek` or `Ra`.
4. Any scaling figure that depends on the new cases.

## Remote bookkeeping

- Use the remote filesystem to inspect which cases have actually finished before adding them to summary plots.
- Prefer reading the latest CSV summaries over recomputing from raw outputs unless the summary is stale or missing.
- If a plot uses `norotating`, keep it as the reference case consistently across future updates.
- Every reduced CSV should record enough provenance to audit it: remote path, frame range or time window, grid size, and diagnostic definition.

## Plotting rules that matter here

- Use `jiasen-scientific-plot-style` as the authoritative plotting standard. In
  particular, keep the approved frame width `4.5`, axis-title size `24`,
  tick-label size `13`, and fixed axes rectangle identical in with-legend and
  no-legend versions; legends must not resize the plotting box.
- Avoid error bars unless the user explicitly asks for them.
- Avoid overcrowded legends and overlapping labels.
- Keep the existing color mapping stable when adding new cases.
- If a figure is a completed set, update the whole set; do not only append the newest points.
- When the user asks to revise or modify an existing figure and does not explicitly ask to keep old versions, delete or overwrite the superseded/incorrect versions. Do not leave stale variants in the active result folder where they can be confused with the current figure.
- For scalar fluctuation profiles of `q`, `b`, and `m`, label the plotted quantity as standard deviation, not generic rms, unless the user explicitly asks for variance. The strict definition is the square root of the combined horizontal-plane/time variance, e.g. `q standard deviation = sqrt(<(q - <q>_{x,y,t}(z))^2>_{x,y,t})`.

## Current regime labels

Use these user-defined labels when annotating or interpreting the phase diagram and related case-by-case summaries:

- `Ra=1e8`, `1e-4 <= Ek <= 1e-3`: funnel regime.
- `Ra=1e8`, `Ek > 1e-3`: plume regime.
- Previously noted cell cases:
  - `Ra=1e6, Ek=3e-4`: cell regime.
  - `Ra=1e8, Ek=1e-5`: cell regime.

Do not infer labels for other points unless the user explicitly assigns them.

## Formula response preference

When the user asks for a LaTeX formula, provide PowerPoint-compatible LaTeX by default: give copyable equation strings using standard commands such as `\frac{}`, `\sqrt{}`, `\left...\right`, subscripts with `{}`, and `\mathrm{}` for roman text. Avoid code-only explanations unless the user explicitly asks for code. For project diagnostics, prefer the exact working definitions documented in `references/diagnostics_and_data_sources.md`.

## Bundled analysis scripts

Use or adapt these deterministic scripts instead of reconstructing the same
pipeline from memory:

- `scripts/aggregation_analysis.py`: full-field MSE variance, spectra,
  correlation length, periodic clusters, and periodic Voronoi diagnostics.
- `scripts/compute_ar10_voronoi_latest10_remote.py`: scan mature AR10 cases and
  export latest-ten-field averaged periodic Voronoi reduced data remotely.
- `scripts/plot_ar10_voronoi_latest10.py`: produce fixed-style single-case and
  multi-case Voronoi figures from reduced NPZ data.
- `scripts/plot_current_aggregation.py`: plot online `A_m`, height-time MSE
  variance, spectral scales, cluster scales, and Voronoi scale for one case.

Inspect paths and command-line arguments before reuse. Keep the remote
calculation and local publication plotting stages separate so raw HDF5 files
remain remote.

## Typical triggers

Use this skill when the user asks to:

- add newly finished remote cases to an existing plot,
- rebuild the reduced phase diagram,
- update `Num`, `Re`, `Rol`, `corr`, or `l_h` figures,
- normalize by the nonrotating run,
- compute or plot `m = b + gamma q` or its perturbation,
- diagnose self-aggregation from MSE variance, spectra, clusters, or Voronoi statistics,
- keep a consistent plotting story across repeated new data drops.
