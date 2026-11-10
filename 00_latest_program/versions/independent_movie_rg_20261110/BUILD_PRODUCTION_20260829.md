# Production build, 2026-08-29

This is the canonical Shuguang DNS source used for newly created cases.

- Runtime array-bounds checking is disabled in both `Makefile` and
  `Makefile.shuguang`.
- `-O3`, `-fma`, and `-traceback` remain enabled.
- The diagnostic schedule remains: ordinary statistics every 0.5 time units,
  aggregation scales every 1.0, and full force diagnostics every 5.0 unless
  explicitly overridden by the separate `Diagnostic output intervals` section
  in `bou.in`.
- The movie/vortex row remains five columns. The output intervals are written
  on their own three-value row (`TSTAT_DIAG`, `TSCALE_DIAG`, `TFORCE_DIAG`).
  Legacy five-column and eight-column files remain readable.
- The canonical remote executable is
  `/work/home/jiasenzhang/rotating_case/latest_program/simexec`.

Existing running cases are not overwritten automatically. New cases created by
`create_rotating_case.py` use the canonical executable above.
