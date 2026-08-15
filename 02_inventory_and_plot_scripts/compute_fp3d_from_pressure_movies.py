"""Compute the strict three-dimensional pressure-force RMS from pressure movies.

For each snapshot and height:
  F'_P = -grad(p) - < -grad(p) >_xy
  F_P(z,t) = sqrt(<F'_Px^2 + F'_Py^2 + F'_Pz^2>_xy)

The reported bulk value is the z-average of F_P(z,t) over zmin <= z <= zmax.
The script uses h5dump because the cluster Python environment has no h5py.
"""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np


FRAME_RE = re.compile(r"pressure(\d+)\.h5$")
TIME_RE = re.compile(r'<Time\s+Value="([+\-0-9.eEdD]+)"')


def ffloat(text: str) -> float:
    return float(text.replace("D", "e").replace("d", "e"))


def dump_dataset(path: Path, dataset: str, output: Path) -> np.ndarray:
    subprocess.run(
        ["h5dump", "-d", "/" + dataset, "-b", "LE", "-o", str(output), str(path)],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    return np.fromfile(output, dtype="<f8")


def read_time(xmf: Path, frame: int) -> float:
    if xmf.exists():
        match = TIME_RE.search(xmf.read_text(errors="ignore"))
        if match:
            return ffloat(match.group(1))
    return 10.0 * frame


def compute_one(run: Path, outdir: Path, zmin: float, zmax: float) -> None:
    movie = run / "movie"
    pressure_paths = []
    for path in movie.glob("pressure[0-9][0-9][0-9][0-9][0-9].h5"):
        if FRAME_RE.search(path.name):
            pressure_paths.append(path)
    pressure_paths.sort(key=lambda p: int(FRAME_RE.search(p.name).group(1)))
    if not pressure_paths:
        return

    outdir.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix="fp3d_", dir=str(outdir)))
    rows = []
    profiles = []
    try:
        coord_path = movie / "cordin_info.h5"
        x = dump_dataset(coord_path, "x", work / "x.bin")
        y = dump_dataset(coord_path, "y", work / "y.bin")
        z = dump_dataset(coord_path, "z", work / "z.bin")
        shape = (len(z), len(y), len(x))
        dx = float(np.mean(np.diff(x)))
        dy = float(np.mean(np.diff(y)))
        zmask = (z >= zmin) & (z <= zmax)
        if np.count_nonzero(zmask) < 2:
            raise RuntimeError(f"not enough z points in [{zmin}, {zmax}]")

        for index, pressure_path in enumerate(pressure_paths, start=1):
            frame = int(FRAME_RE.search(pressure_path.name).group(1))
            p = dump_dataset(pressure_path, "PR_me", work / "p.bin").reshape(shape)
            dpdx = (np.roll(p, -1, axis=2) - np.roll(p, 1, axis=2)) / (2.0 * dx)
            dpdy = (np.roll(p, -1, axis=1) - np.roll(p, 1, axis=1)) / (2.0 * dy)
            dpdz = np.gradient(p, z, axis=0, edge_order=2)
            fpx = -dpdx
            fpy = -dpdy
            fpz = -dpdz
            fpx -= fpx.mean(axis=(1, 2), keepdims=True)
            fpy -= fpy.mean(axis=(1, 2), keepdims=True)
            fpz -= fpz.mean(axis=(1, 2), keepdims=True)
            profile = np.sqrt(np.mean(fpx * fpx + fpy * fpy + fpz * fpz, axis=(1, 2)))
            zsel = z[zmask]
            psel = profile[zmask]
            bulk = float(np.trapz(psel, zsel) / (zsel[-1] - zsel[0]))
            time = read_time(movie / f"field{frame:05d}.xmf", frame)
            rows.append({"time": time, "frame": frame, "F_P_3d_bulk": bulk})
            profiles.append(profile)
            if index % 10 == 0 or index == len(pressure_paths):
                print(f"{run}: {index}/{len(pressure_paths)}", flush=True)

        with (outdir / "fp3d_timeseries.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["time", "frame", "F_P_3d_bulk"])
            writer.writeheader()
            writer.writerows(rows)
        np.savez_compressed(outdir / "fp3d_profiles.npz", time=np.array([r["time"] for r in rows]), frame=np.array([r["frame"] for r in rows]), z=z, F_P_3d=np.asarray(profiles))
    finally:
        shutil.rmtree(str(work), ignore_errors=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--zmin", type=float, default=0.1)
    parser.add_argument("--zmax", type=float, default=0.9)
    args = parser.parse_args()
    compute_one(args.run.resolve(), args.output_dir.resolve(), args.zmin, args.zmax)


if __name__ == "__main__":
    main()
