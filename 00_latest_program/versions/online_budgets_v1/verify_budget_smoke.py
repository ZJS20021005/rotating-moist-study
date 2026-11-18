"""Summarize complete DNS smoke runs in baseline, new1 and new2 subfolders."""
import csv
import json
from pathlib import Path
import re
import sys

import h5py
import numpy as np


def read_budget(run, name):
    with (run / "diagnostics" / "budget" / f"{name}_budget.csv").open() as stream:
        reader = csv.reader(stream)
        columns = next(reader)
        rows = list(reader)
    if any(len(row) != len(columns) for row in rows):
        raise ValueError(f"Ragged CSV: {name} in {run}")
    return dict(zip(columns, np.asarray(rows, dtype=float).T))


def checkpoint_error(before, after, field):
    path = f"continua_{field}.h5"
    with h5py.File(before / path) as old, h5py.File(after / path) as new:
        return max(float(np.max(np.abs(old[key][:] - new[key][:]))) for key in old)


root = Path(sys.argv[1])
baseline, new1, new2 = (root / run for run in ("baseline", "new1", "new2"))
results = {}
for name, rates in {
    "lmse": ("d_sigma_m2_dt", "LMSE_rhs"),
    "ke": ("dK_dt", "KE_rhs"),
    "ke2d": ("dK_2D_dt", "KE2D_rhs"),
}.items():
    a = read_budget(new1, name)
    b = read_budget(new2, name)
    assert len(b["time"]) > 2 and np.all(np.diff(b["time"]) > 0)
    assert np.array_equal(a["time"], b["time"])
    np.testing.assert_allclose(b["time"], np.loadtxt(new2 / "data" / "avgvar.out")[:, 0], atol=5e-5, rtol=0)
    assert all(np.isnan(b[key][0]) and np.isfinite(b[key][1:]).all()
               for key in (rates[0], "rhs_interval", "closure"))
    for key in b:
        if key not in (rates[0], "rhs_interval", "closure"):
            assert np.isfinite(b[key]).all(), (name, key)
    peak_closure = float(np.max(np.abs(b["closure"][1:])))
    rhs_peak = float(np.max(np.abs(b["rhs_interval"][1:])))
    results[name] = dict(
        samples=len(b["time"]),
        peak_interval_closure=peak_closure,
        peak_interval_rhs=rhs_peak,
        peak_closure_over_peak_rhs=peak_closure / rhs_peak if rhs_peak else None,
        two_rank_difference={key: float(np.nanmax(np.abs(a[key] - b[key]))) for key in b},
    )
    if name == "lmse":
        results[name]["peak_production_residual"] = float(np.max(np.abs(b["production_residual"])))
        results[name]["peak_P_grad"] = float(np.max(np.abs(b["P_grad"])))
        legacy = np.loadtxt(new2 / "diagnostics" / "thermo" / "mprime_square.dat")
        matching = [b["sigma_m2"][np.argmin(np.abs(b["time"] - t))] for t in legacy[:, 0]]
        np.testing.assert_allclose(matching, legacy[:, 4], rtol=1e-8, atol=1e-18)
    if name == "ke2d":
        ratio = b["E2D_over_Etot_legacy"]
        assert np.all((ratio >= 0) & (ratio <= 1 + 1e-12))
        np.testing.assert_allclose(b["K_2D_legacy"], ratio*b["K_legacy"], rtol=1e-12, atol=1e-30)
        movie = new2 / "movie"
        xmf = sorted(movie.glob("field*.xmf"))[-1]
        time = float(re.search(r'<Time\s+Value="([^"]+)"', xmf.read_text()).group(1))
        with h5py.File(movie / "cordin_info.h5") as coordinates:
            z = coordinates["z"][:]
        faces = np.zeros(z.size+1)
        for k, height in enumerate(z):
            faces[k+1] = 2*height-faces[k]
        weights = np.diff(faces)/faces[-1]
        index = xmf.stem[len("field"):]
        with h5py.File(movie / f"horizontal_velocity{index}.h5") as horizontal:
            u, v = horizontal["VX_me"][:], horizontal["VY_me"][:]
        with h5py.File(xmf.with_suffix(".h5")) as field:
            w = field["VZ_me"][:]
        e2 = 0.5*np.mean(np.sum(weights[:, None, None]*u, axis=0)**2
                         + np.sum(weights[:, None, None]*v, axis=0)**2)
        etot = 0.5*np.dot(weights, np.mean(u*u+v*v+w*w, axis=(1, 2)))
        i = np.flatnonzero(np.isclose(b["time"], time, rtol=0, atol=1e-10)).item()
        np.testing.assert_allclose([b["K_2D_legacy"][i], b["K_legacy"][i]], [e2, etot], rtol=1e-12)
        results[name]["legacy_movie_ratio_abs_error"] = float(abs(ratio[i]-e2/etot))

assert all(not list((new2 / "diagnostics").glob(f"{group}/*")) for group in ("vortex", "vortex_q"))
results["two_rank_baseline_checkpoint_max_abs"] = {
    field: checkpoint_error(baseline, new2, field)
    for field in ("q1", "q2", "q3", "dsal", "qvap")
}
results["elapsed_seconds"] = {
    run.name: float((run / "wall_seconds.txt").read_text()) for run in (baseline, new1, new2)
}
results["budget_timing"] = {
    run.name: next(line for line in (run / "run.log").read_text().splitlines()
                   if line.startswith("Online budgets:"))
    for run in (new1, new2)
}
print(json.dumps(results, indent=2))
