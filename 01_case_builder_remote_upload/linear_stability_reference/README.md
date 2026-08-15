# Drizzle base-state solver reference

`stability_solver.py` is an unchanged copy of the user-supplied file:

`transfer_linear_stability_Ra8e6_20260728/stability_solver.py`

SHA-256:

`a52049b108e12205c08a7d062de2ceba51ae49aa9d3e6135964ba1d3da342781`

The case builder calls `moist_base_state` from this file. A failed nonlinear
solve is treated as an error; its built-in linear-profile fallback is never
accepted as a DNS initial condition.
