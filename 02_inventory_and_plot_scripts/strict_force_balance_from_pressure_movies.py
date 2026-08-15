#!/usr/bin/env python3
"""Strict rotating moist-convection force-balance diagnostics from movie fields.

This script is designed for the Rainy-Benard rotating_case runs after the DNS
has started writing `movie/pressureNNNNN.h5:/PR_me`.

For every available synchronized snapshot it computes, at each height z,

    F_I = -(u · grad) u
    F_C = (invRo*v, -invRo*u, 0)
    F_P = -grad p
    F_V = sqrt(Pr/Ra) laplacian(u)
    F_B = b' zhat, with b'=b-<b>_xy(z,t)

The component-wise horizontal mean is removed from every force before RMS
profiles are formed.  It also computes the horizontal geostrophic residual,
vertical momentum residual, time-derivative term, and a momentum-budget
closure residual.

Outputs:
  strict_force_balance_zt.npz
  strict_force_balance_bulk_timeseries.csv
  strict_force_balance_timeseries_long.csv
  strict_force_balance_profiles_long.csv
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np


TIME_RE = re.compile(r'<Time\s+Value="([+\-0-9.eEdD]+)"')
FRAME_RE = re.compile(r"field(\d+)\.h5$")


def ffloat(text: str) -> float:
    return float(text.replace("D", "e").replace("d", "e"))


def clean_values(line: str) -> list[str]:
    return line.split("!")[0].strip().split()


def parse_bou(run: Path) -> dict[str, float]:
    """Parse the current bou.in rows used by this project."""
    lines = (run / "bou.in").read_text(errors="ignore").splitlines()
    vals: dict[str, float] = {}
    # Fixed layout used by the maintained rotating_case builder.
    grid = clean_values(lines[2])
    geom = clean_values(lines[11])
    control = clean_values(lines[17])
    cfl = clean_values(lines[23])
    vals.update({
        "N1": ffloat(grid[0]),
        "N2": ffloat(grid[1]),
        "N3": ffloat(grid[2]),
        "ALX3D": ffloat(geom[0]),
        "REXT1": ffloat(geom[1]),
        "REXT2": ffloat(geom[2]),
        "Ra": ffloat(control[0]),
        "Pr": ffloat(control[1]),
        "invRo": ffloat(control[2]),
        "alpha": ffloat(control[3]),
        "gamma": ffloat(control[4]),
        "Sm": ffloat(control[5]),
        "alphaqs": ffloat(control[6]),
        "betaqs": ffloat(control[7]),
        "tau_cond": ffloat(control[8]),
        "dt": ffloat(cfl[3]),
    })
    return vals


def read_time(xmf: Path, frame: int) -> float:
    if xmf.exists():
        text = xmf.read_text(errors="ignore")
        m = TIME_RE.search(text)
        if m:
            return ffloat(m.group(1))
    return 10.0 * frame


def snapshots(run: Path) -> list[dict[str, object]]:
    movie = run / "movie"
    out = []
    if not movie.is_dir():
        return out
    for field in sorted(movie.glob("field[0-9][0-9][0-9][0-9][0-9].h5")):
        m = FRAME_RE.search(field.name)
        if not m:
            continue
        frame = int(m.group(1))
        tag = f"{frame:05d}"
        hvel = movie / f"horizontal_velocity{tag}.h5"
        pressure = movie / f"pressure{tag}.h5"
        coord = movie / "cordin_info.h5"
        xmf = movie / f"field{tag}.xmf"
        if field.exists() and hvel.exists() and pressure.exists() and coord.exists():
            out.append({
                "frame": frame,
                "time": read_time(xmf, frame),
                "field": field,
                "hvel": hvel,
                "pressure": pressure,
                "coord": coord,
            })
    out.sort(key=lambda r: (float(r["time"]), int(r["frame"])))
    return out


def dump_dataset(h5_path: Path, dataset: str, output_path: Path) -> np.ndarray:
    subprocess.check_call(
        ["h5dump", "-d", "/" + dataset, "-b", "LE", "-o", str(output_path), str(h5_path)],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    return np.fromfile(str(output_path), dtype="<f8")


def read_snapshot(snap: dict[str, object], work: Path) -> tuple[np.ndarray, ...]:
    coord = Path(snap["coord"])
    field = Path(snap["field"])
    hvel = Path(snap["hvel"])
    pressure = Path(snap["pressure"])
    tag = f"{int(snap['frame']):05d}"
    x = dump_dataset(coord, "x", work / f"x_{tag}.bin")
    y = dump_dataset(coord, "y", work / f"y_{tag}.bin")
    z = dump_dataset(coord, "z", work / f"z_{tag}.bin")
    shape = (len(z), len(y), len(x))
    u = dump_dataset(hvel, "VX_me", work / f"u_{tag}.bin").reshape(shape)
    v = dump_dataset(hvel, "VY_me", work / f"v_{tag}.bin").reshape(shape)
    w = dump_dataset(field, "VZ_me", work / f"w_{tag}.bin").reshape(shape)
    b = dump_dataset(field, "DSAL_me", work / f"b_{tag}.bin").reshape(shape)
    p = dump_dataset(pressure, "PR_me", work / f"p_{tag}.bin").reshape(shape)
    return x, y, z, u, v, w, b, p


def periodic_first(a: np.ndarray, dx: float, axis: int) -> np.ndarray:
    return (np.roll(a, -1, axis=axis) - np.roll(a, 1, axis=axis)) / (2.0 * dx)


def periodic_second(a: np.ndarray, dx: float, axis: int) -> np.ndarray:
    return (np.roll(a, -1, axis=axis) - 2.0 * a + np.roll(a, 1, axis=axis)) / (dx * dx)


def vertical_first(a: np.ndarray, z: np.ndarray) -> np.ndarray:
    return np.gradient(a, z, axis=0, edge_order=2)


def vertical_second(a: np.ndarray, z: np.ndarray) -> np.ndarray:
    return np.gradient(np.gradient(a, z, axis=0, edge_order=2), z, axis=0, edge_order=2)


def demean_components(comps: tuple[np.ndarray, np.ndarray, np.ndarray]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    return tuple(c - c.mean(axis=(1, 2), keepdims=True) for c in comps)  # type: ignore[return-value]


def profile_total(comps: tuple[np.ndarray, np.ndarray, np.ndarray]) -> np.ndarray:
    return np.sqrt(np.mean(comps[0] ** 2 + comps[1] ** 2 + comps[2] ** 2, axis=(1, 2)))


def profile_h(comps: tuple[np.ndarray, np.ndarray, np.ndarray]) -> np.ndarray:
    return np.sqrt(np.mean(comps[0] ** 2 + comps[1] ** 2, axis=(1, 2)))


def profile_z(comps: tuple[np.ndarray, np.ndarray, np.ndarray]) -> np.ndarray:
    return np.sqrt(np.mean(comps[2] ** 2, axis=(1, 2)))


def z_bulk(profile: np.ndarray, z: np.ndarray, zmin: float, zmax: float) -> float:
    mask = (z >= zmin) & (z <= zmax)
    if np.count_nonzero(mask) < 2:
        return float(np.mean(profile[mask]))
    zz = z[mask]
    pp = profile[mask]
    return float(np.trapz(pp, zz) / (zz[-1] - zz[0]))


def compute_forces(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    u: np.ndarray,
    v: np.ndarray,
    w: np.ndarray,
    b: np.ndarray,
    p: np.ndarray,
    u_prev: tuple[np.ndarray, np.ndarray, np.ndarray] | None,
    t_prev: float | None,
    u_next: tuple[np.ndarray, np.ndarray, np.ndarray] | None,
    t_next: float | None,
    params: dict[str, float],
) -> dict[str, np.ndarray]:
    dx = float(np.mean(np.diff(x)))
    dy = float(np.mean(np.diff(y)))
    nu = math.sqrt(params["Pr"] / params["Ra"])
    inv_ro = params["invRo"]

    fields = {"u": u, "v": v, "w": w}
    deriv: dict[tuple[str, str], np.ndarray] = {}
    for name, arr in fields.items():
        deriv[(name, "x")] = periodic_first(arr, dx, 2)
        deriv[(name, "y")] = periodic_first(arr, dy, 1)
        deriv[(name, "z")] = vertical_first(arr, z)

    inertia = tuple(
        -(u * deriv[(name, "x")] + v * deriv[(name, "y")] + w * deriv[(name, "z")])
        for name in ("u", "v", "w")
    )
    coriolis = (inv_ro * v, -inv_ro * u, np.zeros_like(w))
    pressure = (-periodic_first(p, dx, 2), -periodic_first(p, dy, 1), -vertical_first(p, z))
    viscous = tuple(
        nu * (periodic_second(arr, dx, 2) + periodic_second(arr, dy, 1) + vertical_second(arr, z))
        for arr in (u, v, w)
    )
    b_anom = b - b.mean(axis=(1, 2), keepdims=True)
    buoyancy = (np.zeros_like(w), np.zeros_like(w), b_anom)

    if u_prev is not None and u_next is not None and t_prev is not None and t_next is not None and t_next != t_prev:
        transient = tuple((un - up) / (t_next - t_prev) for up, un in zip(u_prev, u_next))
    else:
        transient = (np.full_like(u, np.nan), np.full_like(v, np.nan), np.full_like(w, np.nan))

    raw = {
        "I": inertia,
        "C": coriolis,
        "P": pressure,
        "V": viscous,
        "B": buoyancy,
        "T": transient,
    }
    primed = {name: demean_components(comps) for name, comps in raw.items()}

    out: dict[str, np.ndarray] = {}
    for name in ("I", "C", "P", "V", "B", "T"):
        out[f"F_{name}"] = profile_total(primed[name])
        out[f"F_{name}_h"] = profile_h(primed[name])
        out[f"F_{name}_z"] = profile_z(primed[name])

    cp_x = primed["C"][0] + primed["P"][0]
    cp_y = primed["C"][1] + primed["P"][1]
    out["R_G"] = np.sqrt(np.mean(cp_x ** 2 + cp_y ** 2, axis=(1, 2)))

    rz = primed["I"][2] + primed["P"][2] + primed["V"][2] + primed["B"][2]
    out["R_z"] = np.sqrt(np.mean(rz ** 2, axis=(1, 2)))

    closure = tuple(
        primed["T"][i] - (primed["I"][i] + primed["C"][i] + primed["P"][i] + primed["V"][i] + primed["B"][i])
        for i in range(3)
    )
    out["R_momentum"] = profile_total(closure)
    out["R_momentum_over_FT"] = out["R_momentum"] / np.maximum(out["F_T"], 1.0e-300)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--output-dir", type=Path, default=None)
    ap.add_argument("--zmin", type=float, default=0.1)
    ap.add_argument("--zmax", type=float, default=0.9)
    args = ap.parse_args()

    run = args.run_dir.resolve()
    outdir = args.output_dir or (run / "diagnostics" / "force_balance_strict")
    outdir.mkdir(parents=True, exist_ok=True)

    params = parse_bou(run)
    snaps = snapshots(run)
    if len(snaps) < 1:
        raise SystemExit(f"No synchronized field/horizontal_velocity/pressure movie snapshots in {run}/movie")

    work = Path(tempfile.mkdtemp(prefix="strict_force_", dir=str(outdir)))
    try:
        cached: dict[int, tuple[np.ndarray, ...]] = {}
        times = np.array([float(s["time"]) for s in snaps])
        frames = np.array([int(s["frame"]) for s in snaps])
        z_ref = None
        series: dict[str, list[np.ndarray]] = {}
        bulk_rows = []
        profile_rows = []

        for idx, snap in enumerate(snaps):
            data = read_snapshot(snap, work)
            cached[idx] = data
            x, y, z, u, v, w, b, p = data
            z_ref = z
            if idx > 0:
                prev = cached.get(idx - 1) or read_snapshot(snaps[idx - 1], work)
                u_prev = (prev[3], prev[4], prev[5])
                t_prev = float(snaps[idx - 1]["time"])
            else:
                u_prev = None
                t_prev = None
            if idx + 1 < len(snaps):
                nxt = cached.get(idx + 1) or read_snapshot(snaps[idx + 1], work)
                cached[idx + 1] = nxt
                u_next = (nxt[3], nxt[4], nxt[5])
                t_next = float(snaps[idx + 1]["time"])
            else:
                u_next = None
                t_next = None

            profiles = compute_forces(x, y, z, u, v, w, b, p, u_prev, t_prev, u_next, t_next, params)
            for key, value in profiles.items():
                series.setdefault(key, []).append(value)

            bulk = {
                "time": float(snap["time"]),
                "frame": int(snap["frame"]),
                "zmin": args.zmin,
                "zmax": args.zmax,
            }
            for key, value in profiles.items():
                bulk[key] = z_bulk(value, z, args.zmin, args.zmax)
            bulk_rows.append(bulk)

            for k, zz in enumerate(z):
                row = {"time": float(snap["time"]), "frame": int(snap["frame"]), "z": float(zz)}
                for key, value in profiles.items():
                    row[key] = float(value[k])
                profile_rows.append(row)

        arrays = {
            "time": times,
            "frame": frames,
            "z": z_ref,
            "zmin": np.array([args.zmin]),
            "zmax": np.array([args.zmax]),
            "Ra": np.array([params["Ra"]]),
            "Pr": np.array([params["Pr"]]),
            "invRo": np.array([params["invRo"]]),
        }
        arrays.update({key: np.vstack(value) for key, value in series.items()})
        # Keep an explicit bulk time-series copy in the NPZ as well.  The
        # z-t arrays above are the primary diagnostic; these are the
        # 0.1H--0.9H column averages used for case-to-case time-series plots.
        for key in series:
            arrays[f"bulk_{key}"] = np.asarray([row[key] for row in bulk_rows], dtype=float)
        np.savez_compressed(outdir / "strict_force_balance_zt.npz", **arrays)

        with (outdir / "strict_force_balance_bulk_timeseries.csv").open("w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(bulk_rows[0].keys()))
            writer.writeheader()
            writer.writerows(bulk_rows)

        # Long-form output makes every force, horizontal/vertical component,
        # residual, and closure quantity an explicit time series.  This is
        # intentionally separate from the wide CSV so plotting code can select
        # quantities without relying on column-name conventions.
        force_keys = [key for key in series if key.startswith(("F_", "R_"))]
        with (outdir / "strict_force_balance_timeseries_long.csv").open("w", newline="") as fh:
            writer = csv.DictWriter(
                fh,
                fieldnames=["time", "frame", "zmin", "zmax", "quantity", "value"],
            )
            writer.writeheader()
            for row in bulk_rows:
                for key in force_keys:
                    writer.writerow({
                        "time": row["time"],
                        "frame": row["frame"],
                        "zmin": row["zmin"],
                        "zmax": row["zmax"],
                        "quantity": key,
                        "value": row[key],
                    })

        with (outdir / "strict_force_balance_profiles_long.csv").open("w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(profile_rows[0].keys()))
            writer.writeheader()
            writer.writerows(profile_rows)

        print(f"wrote {outdir}")
        print(f"frames={len(snaps)} time=[{times.min()}, {times.max()}]")
    finally:
        shutil.rmtree(str(work), ignore_errors=True)


if __name__ == "__main__":
    main()
