# Remote rotating moist case builder

This local tool prepares remote Rainy-Benard / rotating moist convection
cases.  It can create one case from `case_config.json`, or create many cases
from `batch_cases.json`.  It copies a remote template run, edits `bou.in`,
updates the job name in `subjob.sh`, installs the selected `simexec`, checks
the written boundary conditions, generates the case-specific drizzle base
state, and stops before submission.

Default remote output root:

```text
/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case
```

The generated layout is:

```text
rotating_case/
  Pr.../
    Ra.../
      Ek... or norotating/
        AR.../
          Beta.../
            qbot..._qtop.../
              N...x...x.../
                run/
```

For example, the default config creates:

```text
/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case/Pr0p7/Ra1e8/norotating/AR10/Beta1p20/qbot1_qtop0p004978/N129x129x65/run
```

## How to use in VS Code: one case

1. Open this folder in VS Code:

```text
E:\moist RB\rotating_case_inventory\01_case_builder_remote_upload
```

2. Edit `case_config.json`.

3. Run:

```powershell
.\create_case.ps1
```

With no arguments, interactive mode asks for `Ra, Pr, Ek, beta, AR, gamma,
alpha, tau, qbot, qtop, drizzle buoyancy perturbation amplitude, vortex_lc,
n1, n2, n3, time/TMAX`. Press Enter to reuse the default shown in brackets.

You can also override a few parameters from the terminal:

```powershell
.\create_case.ps1 --ra 1e8 --pr 0.7 --ek 9e-4 --beta 1.2 --aspect-ratio 10 --n1 129 --n2 129 --n3 65
```

For a nonrotating case:

```powershell
.\create_case.ps1 --ek norotating
```

To check the target path without writing remotely:

```powershell
.\create_case.ps1 --dry-run
```

## How to use in VS Code: batch cases

Edit `batch_cases.json`, then run:

```powershell
.\create_cases_batch.ps1
```

To preview the whole batch without writing remotely:

```powershell
.\create_cases_batch.ps1 --dry-run
```

The example `batch_cases.json` creates one nonrotating case and a small
rotation/beta sweep:

```json
{
  "defaults": {
    "ra": 100000000.0,
    "pr": 0.7,
    "gamma": 1.1,
    "alpha": 3.0,
    "tau": 0.001,
    "qvapbot": 1.0,
    "qvaptop": 0.004978,
    "aspect_ratio": 10.0,
    "n1": 129,
    "n2": 129,
    "n3": 65
  },
  "cases": [
    {
      "ek": "norotating",
      "beta": 1.2
    },
    {
      "ek": [0.001, 0.0009, 0.0007],
      "beta": [1.1, 1.2]
    }
  ]
}
```

If a value in a case is a list, the program creates every combination.  In the
example above, the second entry creates `3 x 2 = 6` rotating cases.

## SSH update rule

Normally only `~/.ssh/config` needs to change when the cluster gives a new
login node or port.  Keep `ssh_alias` as `c01n0006` if your SSH config uses the
same alias.

If you want a different alias, change `ssh_alias` in `case_config.json`.

## Parameter mapping

- `ra` -> `Ra` in `bou.in`
- `pr` -> `Pr`
- `ek` -> `invRo=sqrt(Pr/Ra)/Ek`; use `norotating` for `invRo=0`
- `gamma` -> `gamma`
- `alpha` -> `alphaqs`
- `beta` -> `betaqs`
- `tau` -> `tau_cond`
- `beta` also sets the buoyancy boundary values:
  - top `dsaltop=beta-1`
  - bottom `dsalbot=0`
- `qvaptop`, `qvapbot` -> moisture boundary values:
  - top `qvaptop`
  - bottom `qvapbot`
- `aspect_ratio` -> `REXT1=REXT2`
- `n1`, `n2`, `n3` -> grid size
- velocity boundary conditions are written as `UBCBOT=1`, `UBCTOP=0`
  (bottom no-slip, top free-slip)

The relevant `bou.in` row is:

```text
A_stopmod   k_stopmod  A_sbotmod   k_sbotmod    dsaltop    dsalbot   qvaptop   qvapbot
```

The program always writes:

```text
dsaltop = beta - 1
dsalbot = 0
qvaptop = qvaptop from config
qvapbot = qvapbot from config
```

For example, with `beta=1.2`, `qvaptop=0.004978`, and `qvapbot=1.0`, the
boundary part of that row is:

```text
0.2d0     0     0.004978d0     1d0
```

After writing `bou.in`, the program reads this row back. If `dsaltop`,
`dsalbot`, `qvaptop`, or `qvapbot` does not match the requested values, it
stops with an error.

The program also checks the number of values in the important `bou.in` rows
after writing:

- `N1 N2 N3 NSST`: 4 values
- `ALX3D REXT1 REXT2 ISTR3 STR3 Lmax`: 6 values
- `Ra Pr invRo alpha gamma Sm alphaqs betaqs tau_cond`: 9 values
- `A_stopmod ... qvapbot`: 8 values

Each prepared case also gets:

- `case_info.json`, recording all parameters and the computed `invRo`
- `submit_after_check.sh`, a one-command remote submission helper
- `run/prepare_drizzle_initial_condition.sh`, the manual drizzle generator
- `run/generate_drizzle_initial_condition.py`, which reads the current
  `bou.in`
- `run/stability_solver.py`, the supplied linear-stability/drizzle solver
- `run/check_drizzle_before_submit.py`, a parameter-consistency check

Before every new (`nread=0`) calculation, enter the case's `run` directory
and execute:

```bash
./prepare_drizzle_initial_condition.sh
```

This command rereads the current `bou.in`, solves that case's converged
one-dimensional drizzle state with the supplied
`stability_solver.moist_base_state`, checks it, and atomically overwrites:

- `drizzle_init.dat`
- `drizzle_init_meta.json`
- `drizzle_init_last_generation.json`

The built-in linear-profile fallback is rejected. The generator uses 401
vertical points, perturbation amplitude `1e-4`, and
`saturation_width=1e-8`, matching the DNS switch
`0.5*(1+tanh(1e8*(q-qs)))`.

After generation, submit normally. `simexec` reads `drizzle_init.dat`,
verifies its parameters and boundary values against the current `bou.in`,
sets `q=q_drizzle`, and initializes
`b=b_drizzle+1e-4*sin(pi*z/H)*N(0,1)`. Therefore the perturbation is added
to the drizzle buoyancy rather than to a linear buoyancy profile.

The case builder only installs the generator and the newest `simexec`; it
does not silently generate or submit a case. This keeps the two required
steps explicit: generate drizzle, then run/submit.
