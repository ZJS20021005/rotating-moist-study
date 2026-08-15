# Historical time-series patches

This folder stores immutable reduced-data repairs for simulation histories
that are no longer complete on the remote filesystem.

## Ra8e6, Ek7e-3

- File: `Ra8e6_Ek7e-3_t10_1000_latest_program.csv`
- Physical case: `Ra=8e6`, `Pr=0.7`, `Ek=7e-3`, `AR=16`, `beta=1.02`,
  `qbot=0.5`, actual grid `385x385x65`.
- Time range: `10 <= t <= 1000`, sampled every `0.1`.
- Source: the complete reduced table saved on 2026-08-01 at
  `04_outputs_and_figures/high_resolution_timeseries_latest_program_20260801/high_resolution_timeseries_long.csv`.
- Included metrics: kinetic energy, MSE variance, both `l0` heights,
  MSE peak and spectral scales, and both convective scales.
- SHA-256: `3FD26B8E068EBC7FF78F513DE6F3FC53FE4F6455D7755E6138921F789F87A08F`.

The unified plotting script loads this file before current remote rows.
Therefore current remote data take priority at duplicate physical times.
