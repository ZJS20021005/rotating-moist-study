"""Check GMSE smoke runs that started from the same DNS checkpoint."""

import argparse
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np


COLUMNS = (
    "time", "GMSE_variance", "GMSE_PH", "GMSE_DH", "GMSE_SB", "GMSE_rhs",
    "d_GMSE_variance_dt", "rhs_interval", "closure", "closure_abs", "closure_rel",
    "GMSE_PH_ibp", "GMSE_PH_residual",
)
INTERVAL = ("d_GMSE_variance_dt", "rhs_interval", "closure", "closure_abs", "closure_rel")
MPI_DIAGNOSTICS = ("GMSE_variance", "GMSE_PH", "GMSE_DH", "GMSE_SB", "GMSE_rhs")


def peak(values):
    return float(np.max(np.abs(values)))


def check_close(actual, expected, scale=None):
    if scale is None:
        scale = peak(expected)
    np.testing.assert_allclose(actual, expected, rtol=4e-12, atol=4e-12 * scale + 1e-30)


def read_run(run, closure_rtol, ibp_rtol):
    path = run / "diagnostics" / "budget" / "gmse_budget.csv"
    budget = np.genfromtxt(path, names=True, delimiter=",", ndmin=1, encoding="utf-8-sig")
    if budget.dtype.names != COLUMNS or budget.size < 3:
        raise ValueError(f"{path} must contain the GMSE columns and at least three samples")
    for name in COLUMNS:
        values = budget[name]
        if name in INTERVAL:
            if not np.isnan(values[0]):
                raise ValueError(f"{path}: expected NaN in the first row of {name}")
            values = values[1:]
        if not np.isfinite(values).all():
            raise ValueError(f"{path}: unexpected non-finite values in {name}")

    time = budget["time"]
    dt = np.diff(time)
    if not np.all(dt > 0):
        raise ValueError(f"{path}: output times must increase")
    avgvar_time = np.loadtxt(run / "data" / "avgvar.out", usecols=0, ndmin=1)
    np.testing.assert_allclose(time, avgvar_time, atol=5e-5, rtol=0)
    lmse_time = np.loadtxt(run / "diagnostics" / "budget" / "lmse_budget.csv",
                          delimiter=",", skiprows=1, usecols=0, ndmin=1)
    np.testing.assert_array_equal(time, lmse_time)

    variance = budget["GMSE_variance"]
    if not np.all(variance >= 0) or peak(variance) == 0:
        raise ValueError(f"{path}: expected nonzero, nonnegative GMSE variance")
    ph, dh, sb = (budget[name] for name in ("GMSE_PH", "GMSE_DH", "GMSE_SB"))
    rhs = 2 * (ph - dh + sb)
    check_close(budget["GMSE_rhs"], rhs, 2 * (peak(ph) + peak(dh) + peak(sb)))
    trapezoid = (budget["GMSE_rhs"][1:] + budget["GMSE_rhs"][:-1]) / 2
    rate = np.diff(variance) / dt
    closure = rate - trapezoid
    check_close(budget["d_GMSE_variance_dt"][1:], rate)
    check_close(budget["rhs_interval"][1:], trapezoid)
    check_close(budget["closure"][1:], closure, max(peak(rate), peak(trapezoid)))
    check_close(budget["closure_abs"][1:], np.abs(closure), max(peak(rate), peak(trapezoid)))
    denom = np.maximum(np.maximum(np.abs(rate), np.abs(trapezoid)), np.finfo(float).tiny)
    check_close(budget["closure_rel"][1:], np.abs(closure) / denom, 1)

    ph_ibp = budget["GMSE_PH_ibp"]
    ibp = ph - ph_ibp
    ph_scale = max(peak(ph), peak(ph_ibp))
    check_close(budget["GMSE_PH_residual"], ibp, ph_scale)
    if ph_scale == 0 or peak(dh) == 0 or peak(sb) == 0:
        raise ValueError(f"{path}: initial condition must exercise all production and diffusion terms")
    ibp_error = peak(ibp) / ph_scale
    if ibp_error > ibp_rtol:
        raise ValueError(f"{path}: integration-by-parts error {ibp_error:.3g} > {ibp_rtol:.3g}")

    rhs_scale = max(peak(rate), peak(trapezoid))
    if rhs_scale == 0:
        raise ValueError(f"{path}: the budget time series is stationary")
    closure_error = peak(closure) / rhs_scale
    if closure_error > closure_rtol:
        raise ValueError(f"{path}: normalized closure error {closure_error:.3g} > {closure_rtol:.3g}")
    result = {
        "directory": str(run.resolve()),
        "samples": int(budget.size),
        "time_start": float(time[0]),
        "time_end": float(time[-1]),
        "peak_abs": {name: peak(budget[name]) for name in MPI_DIAGNOSTICS},
        "peak_interval_rate": peak(rate),
        "peak_interval_rhs": peak(trapezoid),
        "peak_absolute_closure": peak(closure),
        "peak_absolute_closure_over_peak_rate_or_rhs": closure_error,
        "peak_pointwise_relative_closure": peak(budget["closure_rel"][1:]),
        "signed_integrated_closure": float(np.dot(closure, dt)),
        "peak_horizontal_integration_by_parts_residual": peak(ibp),
        "peak_horizontal_integration_by_parts_relative_error": ibp_error,
    }
    wall = run / "wall_seconds.txt"
    if wall.is_file():
        result["wall_seconds"] = float(wall.read_text().strip())
    log = run / "run.log"
    if log.is_file():
        for line in log.read_text(errors="replace").splitlines():
            if line.startswith("Online budgets:"):
                result["budget_timing"] = line
    return budget, result


def verify_movie(run, budget):
    import h5py

    # The movie contains the unaveraged m' and exact scalar cell centres.
    movie = run / "movie"
    with h5py.File(movie / "cordin_info.h5") as coordinates:
        z = coordinates["z"][:]
    faces = np.zeros(z.size + 1)
    for k, height in enumerate(z):
        faces[k + 1] = 2 * height - faces[k]
    weights = np.diff(faces)
    if not np.all(weights > 0):
        raise ValueError(f"{movie}: invalid cell weights from movie coordinates")

    errors = []
    for xmf in sorted(movie.glob("mprime*.xmf")):
        time = float(ET.parse(xmf).find(".//Time").attrib["Value"])
        rows = np.flatnonzero(np.isclose(budget["time"], time, atol=5e-9, rtol=0))
        if rows.size != 1:
            continue
        with h5py.File(xmf.with_suffix(".h5")) as frame:
            mprime = frame["MPRIME_me"][:]
        column = np.tensordot(weights, mprime, axes=1)
        variance = float(np.mean(column**2))
        np.testing.assert_allclose(budget["GMSE_variance"][rows[0]], variance, rtol=1e-10, atol=1e-24)
        errors.append(abs(float(budget["GMSE_variance"][rows[0]]) - variance))
    if not errors:
        raise ValueError(f"{movie}: no m' movie frame matches a GMSE output time")
    return {"matching_movie_frames": len(errors), "peak_variance_abs_error": max(errors)}


def verify_baseline(baseline, run):
    import h5py

    # Compare whole solver states from equal-rank runs of the old and new executable.
    errors = {}
    for field in ("q1", "q2", "q3", "dsal", "qvap"):
        path = f"continua_{field}.h5"
        with h5py.File(baseline / path) as old, h5py.File(run / path) as new:
            if list(old.keys()) != list(new.keys()):
                raise ValueError(f"{path}: the old and new runs produced different datasets")
            for name in old.keys():
                before, after = old[name][:], new[name][:]
                np.testing.assert_array_equal(after, before)
                errors[field] = peak(after - before)
    for path in ("data/avgvar.out", "diagnostics/budget/lmse_budget.csv",
                 "diagnostics/budget/ke_budget.csv", "diagnostics/budget/ke2d_budget.csv"):
        if (baseline / path).read_bytes() != (run / path).read_bytes():
            raise ValueError(f"{path}: an existing diagnostic changed")
    result = {"two_rank_checkpoint_max_abs_error": errors}
    wall = baseline / "wall_seconds.txt"
    if wall.is_file():
        result["wall_seconds"] = float(wall.read_text().strip())
    return result


def verify_append(original, continued, budget):
    for name in ("gmse", "lmse", "ke", "ke2d"):
        path = Path("diagnostics") / "budget" / f"{name}_budget.csv"
        before = (original / path).read_bytes()
        after = (continued / path).read_bytes()
        header = before.splitlines()[0]
        if not after.startswith(before) or after.splitlines().count(header) != 1:
            raise ValueError(f"{path}: continuation altered existing output or repeated the header")

    path = continued / "diagnostics" / "budget" / "gmse_budget.csv"
    rows = np.genfromtxt(path, names=True, delimiter=",", ndmin=1)[budget.size:]
    if rows.size < 2 or not np.all(np.diff(np.r_[budget["time"][-1], rows["time"]]) > 0):
        raise ValueError(f"{path}: continuation must provide strictly later samples")
    for name in INTERVAL:
        if not np.isnan(rows[name][0]) or not np.isfinite(rows[name][1:]).all():
            raise ValueError(f"{path}: invalid continuation interval values for {name}")
    expected = np.diff(rows["GMSE_variance"]) / np.diff(rows["time"])
    check_close(rows["d_GMSE_variance_dt"][1:], expected)
    continuation_times = np.loadtxt(continued / "diagnostics" / "budget" / "lmse_budget.csv",
                                   delimiter=",", skiprows=1, usecols=0, ndmin=1)
    np.testing.assert_array_equal(rows["time"], continuation_times[budget.size:])
    return {"existing_samples": int(budget.size), "appended_samples": int(rows.size),
            "time_start": float(rows["time"][0]), "time_end": float(rows["time"][-1]),
            "single_header_preserved": True, "first_continuation_interval_reset": True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--one-rank", type=Path, required=True)
    parser.add_argument("--two-rank", type=Path, required=True)
    parser.add_argument("--baseline", type=Path, help="Optional equal-configuration v1 two-rank run")
    parser.add_argument("--append-run", type=Path, help="Optional continuation of the two-rank run")
    parser.add_argument("--offline", action="store_true", help="Check variance from matching HDF5 frames")
    parser.add_argument("--closure-rtol", type=float, default=5e-3)
    parser.add_argument("--ibp-rtol", type=float, default=1e-8)
    parser.add_argument("--mpi-rtol", type=float, default=5e-7)
    args = parser.parse_args()
    one, one_summary = read_run(args.one_rank, args.closure_rtol, args.ibp_rtol)
    two, two_summary = read_run(args.two_rank, args.closure_rtol, args.ibp_rtol)
    np.testing.assert_array_equal(one["time"], two["time"])
    result = {"one_rank": one_summary, "two_ranks": two_summary, "mpi_rank_comparison": {}}
    for name in MPI_DIAGNOSTICS:
        difference = peak(one[name] - two[name])
        scale = max(peak(one[name]), peak(two[name]), np.finfo(float).tiny)
        error = difference / scale
        if error > args.mpi_rtol:
            raise ValueError(f"{name}: one/two-rank error {error:.3g} > {args.mpi_rtol:.3g}")
        result["mpi_rank_comparison"][name] = {"peak_abs": difference, "peak_normalized": error}
    if args.offline:
        result["offline_variance"] = verify_movie(args.two_rank, two)
    if args.baseline:
        result["baseline"] = verify_baseline(args.baseline, args.two_rank)
    if args.append_run:
        result["append"] = verify_append(args.two_rank, args.append_run, two)
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
