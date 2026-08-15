# Strict force-balance post-processing

`strict_force_balance_from_pressure_movies.py` reads synchronized
`fieldNNNNN.h5`, `horizontal_velocityNNNNN.h5`, and `pressureNNNNN.h5`
movie snapshots from one Rainy-Bénard `run` directory.

It uses the actual nondimensional coefficients from `bou.in`:

- `F_I = -(u·∇)u`
- `F_C = invRo (v, -u, 0)`
- `F_P = -∇p`, including `-∂z p`
- `F_V = sqrt(Pr/Ra) ∇²u`
- `F_B = (0, 0, b - <b>xy)`
- `F_T = ∂t u`

At every height and time, each component has its horizontal mean removed
before RMS diagnostics are formed. The default bulk interval is
`0.1 <= z <= 0.9`.

Run locally:

```powershell
python strict_force_balance_from_pressure_movies.py `
  "E:\path\to\run" `
  --output-dir "E:\path\to\run\diagnostics\force_balance_strict"
```

Outputs:

- `strict_force_balance_zt.npz`: all `z-t` profiles and metadata;
- `strict_force_balance_bulk_timeseries.csv`: wide bulk time series;
- `strict_force_balance_timeseries_long.csv`: one row per time/quantity;
- `strict_force_balance_profiles_long.csv`: long-form `z-t` table.

The `F_*` entries include total, horizontal, and vertical RMS values.
`R_G`, `R_z`, `R_momentum`, and `R_momentum_over_FT` are residual diagnostics.

The script requires `h5dump` on `PATH`.

## DNS online force output

The legacy `data/pressure_force.out` from `avgvar.f90` contains
`time F_P,h`, with

`F_P = sqrt(<(dp/dx)^2 + (dp/dy)^2>_V)`.

The updated 2026-08-06 executable additionally writes
`data/force_balance.out` through `force_balance_online.f90`. It contains
strict total, horizontal, and vertical RMS time series for `F_I`, `F_C`,
`F_P`, `F_V`, `F_B`, and `F_T`, followed by `R_G`, `R_z`, the full
momentum closure residual, and separate x/y component RMS time series for
each force. All components are horizontally demeaned at every height and
the per-height RMS values are then averaged over `0.1H<=z<=0.9H`. The drag
coefficient is zero, so no drag term is calculated.

The same quantities at every output height are written to
`data/force_balance_z.out`. This is the online z-t output; it can be used
directly for height-time plots without reconstructing the force terms from
the movie fields.

For already completed pressure-enabled runs, use
`compute_fp3d_from_pressure_movies.py` to recover the full three-dimensional
`F_P` from `PR_me`, including `-dp/dz`. The legacy `pressure_force.out` must
not be substituted for the strict `F_P` in force-balance figures.

## Stage continuation runs

Windows/VScode terminal:

```powershell
.\prepare_strict_force_continuations.ps1 -DryRun
.\prepare_strict_force_continuations.ps1
```

The second command only creates or reports continuation directories. It never
submits a job. On the cluster, the equivalent command is:

```bash
/share/org/SHUTUANL/shu_zhangjs/rainy\ model/rotating_case/postprocess/strict_force_balance_20260805/prepare_strict_force_continuations.sh --dry-run
```
