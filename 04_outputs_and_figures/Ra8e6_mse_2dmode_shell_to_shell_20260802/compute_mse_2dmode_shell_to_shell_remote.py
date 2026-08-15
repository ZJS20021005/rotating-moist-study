#!/usr/bin/env python3
"""Fast vertically averaged MSE-mode shell-to-shell transfer for all cases."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
from pathlib import Path

import numpy as np


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def write_csv(path: Path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def process_case(label, anchor, output_dir, helpers, shell, nframes, matrix_dk, kmax):
    frames = helpers.select_latest_complete_frames(anchor, nframes)
    params = helpers.read_bou(Path(frames[-1]["field"]).parent.parent / "bou.in")
    shape = helpers.dataset_shape(Path(frames[-1]["mprime"]), "MPRIME_me")
    nz, ny, nx = shape
    coordinate = Path(frames[-1]["field"]).parent / "cordin_info.h5"
    z = helpers.h5dump_array(coordinate, "z")
    zweights = helpers.vertical_weights(z, 0.0, float(params["Lz"]))
    kx, ky, kh, bin_id, centers, edges, masks, complete = shell.geometry(
        nx,
        ny,
        float(params["Lx"]),
        float(params["Ly"]),
        matrix_dk,
        kmax,
    )

    results = []
    for iframe, frame in enumerate(frames, start=1):
        mprime = helpers.h5dump_array(Path(frame["mprime"]), "MPRIME_me", shape)
        mprime -= mprime.mean(axis=(-2, -1), keepdims=True)
        u = helpers.h5dump_array(Path(frame["horizontal"]), "VX_me", shape)
        v = helpers.h5dump_array(Path(frame["horizontal"]), "VY_me", shape)
        result = shell.vertical_mode_matrix(
            mprime, u, v, zweights, kx, ky, bin_id, masks, matrix_dk
        )
        result["time"] = float(frame["time"])
        results.append(result)
        print(f"{label}: t={result['time']:.6g} ({iframe}/{len(frames)})", flush=True)

    names = [
        "T_2D",
        "T_2D_row_sum",
        "T_2D_exchange",
        "T_2D_exchange_row_sum",
        "T_2D_symmetric",
        "T_2D_exact_shell",
        "F_3D_to_2D_shell",
    ]
    scalar_names = [
        "M_2D_variance",
        "T_2D_antisymmetry_relative_error",
        "UV_horizontal_divergence_rms_over_velocity_rms",
    ]
    arrays = {
        name + "_frames": np.stack([np.asarray(item[name]) for item in results])
        for name in names
    }
    scalars = {
        name + "_frames": np.array([float(item[name]) for item in results])
        for name in scalar_names
    }
    times = np.array([item["time"] for item in results], dtype=float)
    case_dir = output_dir / label
    case_dir.mkdir(parents=True, exist_ok=True)
    npz_path = case_dir / "mse_2dmode_shell_to_shell_latest10.npz"
    np.savez_compressed(
        npz_path,
        times=times,
        k_centers=centers,
        k_edges=edges,
        matrix_dk=matrix_dk,
        z=z,
        zweights=zweights,
        **arrays,
        **scalars,
    )
    manifest = {
        "case": label,
        "anchor": str(anchor),
        "parameters": params,
        "movie_shape_zyx": shape,
        "times": times.tolist(),
        "matrix_dk": matrix_dk,
        "matrix_kmax_edge": float(edges[-1]),
        "complete_isotropic_kmax": complete,
        "n_bins": len(centers),
        "vertical_average": "M=sum_z(w_z*mprime), U=sum_z(w_z*u), V=sum_z(w_z*v); w_z are stretched-grid cell-thickness weights",
        "mprime": "m-<m>_xy independently at each z; stored MPRIME_me used",
        "T_2D": "-<M_K (U partial_x+V partial_y) M_P>_xy/dk; P donor, K receiver",
        "plotted_exchange": "0.5*(T_2D-T_2D.T)",
        "F_3D_to_2D": "shell input from exact -div_h<uh*mprime>_z minus T_2D row sum",
        "averaging": "matrix computed in every frame, then latest-ten arithmetic mean",
        "antisymmetry_relative_error_max": float(
            np.max(scalars["T_2D_antisymmetry_relative_error_frames"])
        ),
        "npz": str(npz_path),
    }
    (case_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    return manifest


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--helpers", type=Path, required=True)
    parser.add_argument("--shell-module", type=Path, required=True)
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--nframes", type=int, default=10)
    parser.add_argument("--matrix-dk", type=float, default=np.pi / 2.0)
    parser.add_argument("--kmax", type=float, default=60.0)
    parser.add_argument("--cases", nargs="*", default=None)
    args = parser.parse_args()

    helpers = load_module(args.helpers, "mse_helpers")
    shell = load_module(args.shell_module, "mse_shell")
    labels = list(helpers.CASE_RELATIVE_PATHS)
    if args.cases:
        requested = set(args.cases)
        labels = [label for label in labels if label in requested]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifests = []
    for label in labels:
        manifests.append(
            process_case(
                label,
                args.base / helpers.CASE_RELATIVE_PATHS[label],
                args.output_dir,
                helpers,
                shell,
                args.nframes,
                args.matrix_dk,
                args.kmax,
            )
        )

    rows = []
    for item in manifests:
        rows.append(
            {
                "case": item["case"],
                "Ek": item["parameters"]["Ek"],
                "Lx": item["parameters"]["Lx"],
                "movie_shape_zyx": str(tuple(item["movie_shape_zyx"])),
                "t_first": item["times"][0],
                "t_last": item["times"][-1],
                "nframes": len(item["times"]),
                "matrix_dk": item["matrix_dk"],
                "matrix_kmax_edge": item["matrix_kmax_edge"],
                "n_bins": item["n_bins"],
                "antisymmetry_relative_error_max": item[
                    "antisymmetry_relative_error_max"
                ],
            }
        )
    write_csv(args.output_dir / "mse_2dmode_shell_to_shell_summary.csv", rows)
    (args.output_dir / "manifest_all_cases.json").write_text(
        json.dumps(manifests, indent=2), encoding="utf-8"
    )
    print(f"saved {args.output_dir}", flush=True)


if __name__ == "__main__":
    main()
