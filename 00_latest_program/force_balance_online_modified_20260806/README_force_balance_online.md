# Online force-balance output (2026-08-06)

Remote source root:

`/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case/latest_program/source`

The modified executable writes `run/data/force_balance.out` at the same
statistics cadence as `avgvar.out`. The first statistics call initializes the
previous velocity field; output begins at the second statistics call.

## Definitions

The actual DNS coefficients and staggered-grid stencils are used:

- `F_I = -(u . grad)u`
- `F_C = invRo (v, -u, 0)`
- `F_P = -grad(p)`, including the vertical pressure gradient
- `F_V = sqrt(Pr/Ra) laplacian(u)`
- `F_B = (0, 0, b)` at the vertical-velocity location
- `F_T = [u(t)-u(t_previous_statistics)]/(t-t_previous_statistics)`

The configured drag coefficient is zero, so no drag term is calculated or
written.

At every height, each component is horizontally demeaned before RMS. Bulk
values use `0.1H <= z <= 0.9H`.

## Columns

`force_balance.out` has 34 columns:

1. time
2-7. `F_I F_C F_P F_V F_B F_T`
8-13. `F_Ih F_Ch F_Ph F_Vh F_Bh F_Th`
14-19. `F_Iz F_Cz F_Pz F_Vz F_Bz F_Tz`
20. `R_G = rms(F_C,h + F_P,h)`
21. `R_z = rms(F_I,z + F_P,z + F_V,z + F_B,z)`
22. `R_momentum = rms(F_T - F_I - F_C - F_P - F_V - F_B)`
23-28. `F_Ix F_Cx F_Px F_Vx F_Bx F_Tx`
29-34. `F_Iy F_Cy F_Py F_Vy F_By F_Ty`

`run/data/force_balance_z.out` has 35 columns and contains the same
quantities at every output height:

1. time; 2. z
3-8. total RMS `F_I F_C F_P F_V F_B F_T`
9-14. horizontal RMS `F_Ih F_Ch F_Ph F_Vh F_Bh F_Th`
15-20. vertical RMS `F_Iz F_Cz F_Pz F_Vz F_Bz F_Tz`
21-23. `R_G R_z R_momentum`
24-29. x-component RMS
30-35. y-component RMS

## Modified files

- `fluid_solver/force_balance_online.f90`
- `fluid_solver/openfi.f90`
- `gcurv.f90`
- `Makefile`

Original remote versions are retained beside the source with the suffix
`.before_force_balance_20260806`. Original local copies use `.original`.

## Build

Run `compile_remote.sh` in a shell that can access the cluster. The compiled
remote executable is `latest_program/source/simexec`. A synchronized local
copy is `simexec_force_balance_online_20260806` in this directory.

This build does not replace executables in existing case directories and does
not submit or restart jobs. Future continuations must copy the new latest
`simexec` into their `run` directory before submission.
