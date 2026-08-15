# Phase diagram and completed-case rules

- Treat the reduced phase diagram as a living summary of completed runs.
- Add a point only when the case is complete enough for stable averaging.
- Keep the nonrotating case as a reference only when the figure definition needs normalization.
- Use the same color family for the same Ra in related figures.
- When new results appear, refresh the phase diagram before the derivative summary plots.
- Do not silently switch from a full-field measure to a vertical-only surrogate.
# Intermittent low-Ek cases

- For the current Ra8e6 latest-program set, cases with `Ek < 1e-3`
  (`1.5e-4`, `2e-4`, `5e-4`, `7e-4`) are treated as intermittent states.
- Do not describe these cases using a single "stable/unstable" label and do
  not assign a force-balance regime from one snapshot or a full-time mean.
- Diagnose force balance as a time series and distinguish burst/high-energy
  and quiescent/low-energy phases. The current reduced time-series source is
  `force_balance/intermittent_timeseries/force_balance_timeseries.csv` under
  the 2026-08-05 latest-program output root.
