#!/usr/bin/env python3
"""Shell-to-shell transfer of moist-static-energy variance.

Two related diagnostics are produced from every snapshot.

1. Full-volume horizontal-shell transfer

   T_full(K,P) = -< m'_K u_i partial_i m'_P >_V / Delta k.

   This is the conservative nonlinear redistribution of the complete 3-D
   MSE-anomaly variance among horizontal Fourier shells.

2. Vertically coherent (2-D-mode) transfer

   M=<m'>_z, U=<u>_z, V=<v>_z,
   T_2D(K,P) = -< M_K (U partial_x + V partial_y) M_P >_xy / Delta k.

   It measures self-transfer within the vertically coherent scalar and
   horizontal-velocity modes.  The exact vertically averaged scalar flux
   also contains <u'_h m''>_z; its residual forcing spectrum is saved as
   F_3D_to_2D.

The donor shell is P (matrix column) and the receiver shell is K (row).
Positive T(K,P) means that P supplies MSE variance to K.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
from pathlib import Path

import numpy as np


def load_helpers(path: Path):
    spec = importlib.util.spec_from_file_location("mse_gaussian_helpers", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def fft2f(field: np.ndarray) -> np.ndarray:
    return np.fft.fft2(field, axes=(-2, -1), norm="forward")


def ifft2f(field: np.ndarray) -> np.ndarray:
    return np.fft.ifft2(field, axes=(-2, -1), norm="forward").real


def bin_sum(values: np.ndarray, bin_id: np.ndarray, nbins: int, dk: float) -> np.ndarray:
    valid = bin_id >= 0
    return (
        np.bincount(
            bin_id[valid].ravel(),
            weights=np.asarray(values)[valid].ravel(),
            minlength=nbins,
        )[:nbins]
        / dk
    )


def geometry(nx: int, ny: int, lx: float, ly: float, matrix_dk: float, kmax: float):
    dx, dy = lx / nx, ly / ny
    kx_1d = 2.0 * np.pi * np.fft.fftfreq(nx, d=dx)
    ky_1d = 2.0 * np.pi * np.fft.fftfreq(ny, d=dy)
    kx, ky = np.meshgrid(kx_1d, ky_1d, indexing="xy")
    kh = np.sqrt(kx * kx + ky * ky)
    complete = 0.90 * min(np.pi / dx, np.pi / dy)
    upper = min(float(kmax), complete)
    nbins = int(math.floor(upper / matrix_dk))
    if nbins < 2:
        raise ValueError("Too few transfer bins")
    upper = nbins * matrix_dk
    bin_id = np.floor(kh / matrix_dk).astype(int)
    bin_id[(kh <= 0.0) | (kh >= upper)] = -1
    centers = (np.arange(nbins, dtype=float) + 0.5) * matrix_dk
    edges = np.arange(nbins + 1, dtype=float) * matrix_dk
    masks = [bin_id == index for index in range(nbins)]
    return kx, ky, kh, bin_id, centers, edges, masks, complete


def matrix_antisymmetry_error(matrix: np.ndarray) -> float:
    scale = float(np.max(np.abs(matrix)))
    if scale == 0.0:
        return 0.0
    return float(np.max(np.abs(matrix + matrix.T)) / scale)


def full_volume_matrix(
    mprime: np.ndarray,
    u: np.ndarray,
    v: np.ndarray,
    w: np.ndarray,
    z: np.ndarray,
    zweights: np.ndarray,
    kx: np.ndarray,
    ky: np.ndarray,
    bin_id: np.ndarray,
    masks: list[np.ndarray],
    dk: float,
):
    mhat = fft2f(mprime)
    nbins = len(masks)
    horizontal = np.zeros((nbins, nbins), dtype=float)
    vertical = np.zeros_like(horizontal)

    # Direct complete nonlinear transfer, used to quantify omitted high-k
    # donors/receivers when the matrix is truncated at the common kmax.
    mx = ifft2f(1j * kx[None, :, :] * mhat)
    my = ifft2f(1j * ky[None, :, :] * mhat)
    mz = np.gradient(mprime, z, axis=0, edge_order=2)
    nhat_direct = fft2f(-(u * mx + v * my + w * mz))
    direct_mode = np.tensordot(
        zweights, np.real(np.conj(mhat) * nhat_direct), axes=(0, 0)
    )
    direct_shell = bin_sum(direct_mode, bin_id, nbins, dk)
    del mx, my, mz, nhat_direct, direct_mode

    for donor_index, mask in enumerate(masks):
        donor_hat = mhat * mask[None, :, :]
        donor_x = ifft2f(1j * kx[None, :, :] * donor_hat)
        donor_y = ifft2f(1j * ky[None, :, :] * donor_hat)
        donor_field = ifft2f(donor_hat)
        donor_z = np.gradient(donor_field, z, axis=0, edge_order=2)

        nh_hat = fft2f(-(u * donor_x + v * donor_y))
        nv_hat = fft2f(-(w * donor_z))
        mode_h = np.tensordot(
            zweights, np.real(np.conj(mhat) * nh_hat), axes=(0, 0)
        )
        mode_v = np.tensordot(
            zweights, np.real(np.conj(mhat) * nv_hat), axes=(0, 0)
        )
        horizontal[:, donor_index] = bin_sum(mode_h, bin_id, nbins, dk)
        vertical[:, donor_index] = bin_sum(mode_v, bin_id, nbins, dk)
        del donor_hat, donor_x, donor_y, donor_field, donor_z, nh_hat, nv_hat, mode_h, mode_v

    total = horizontal + vertical
    row_sum = total.sum(axis=1)
    exchange = 0.5 * (total - total.T)
    symmetric = 0.5 * (total + total.T)
    return {
        "T_full_horizontal": horizontal,
        "T_full_vertical": vertical,
        "T_full_total": total,
        "T_full_row_sum": row_sum,
        "T_full_exchange": exchange,
        "T_full_exchange_row_sum": exchange.sum(axis=1),
        "T_full_symmetric": symmetric,
        "T_full_direct_shell": direct_shell,
        "T_full_truncation_residual": direct_shell - row_sum,
        "full_antisymmetry_relative_error": matrix_antisymmetry_error(total),
    }


def vertical_mode_matrix(
    mprime: np.ndarray,
    u: np.ndarray,
    v: np.ndarray,
    zweights: np.ndarray,
    kx: np.ndarray,
    ky: np.ndarray,
    bin_id: np.ndarray,
    masks: list[np.ndarray],
    dk: float,
):
    # Stretched-grid volume averages.  M is not an equal-weight average over
    # the stored z planes.
    M = np.tensordot(zweights, mprime, axes=(0, 0))
    U = np.tensordot(zweights, u, axes=(0, 0))
    V = np.tensordot(zweights, v, axes=(0, 0))
    M -= M.mean()
    Mhat = fft2f(M)
    nbins = len(masks)
    matrix = np.zeros((nbins, nbins), dtype=float)

    for donor, mask in enumerate(masks):
        donor_hat = Mhat * mask
        donor_x = ifft2f(1j * kx * donor_hat)
        donor_y = ifft2f(1j * ky * donor_hat)
        nonlinear_hat = fft2f(-(U * donor_x + V * donor_y))
        transfer_mode = np.real(np.conj(Mhat) * nonlinear_hat)
        matrix[:, donor] = bin_sum(transfer_mode, bin_id, nbins, dk)

    row_sum = matrix.sum(axis=1)
    exchange = 0.5 * (matrix - matrix.T)
    symmetric = 0.5 * (matrix + matrix.T)

    # Exact nonlinear forcing of the vertically averaged scalar:
    # -div_h <u_h m'>_z = -div_h(U M) - div_h<u''_h m''>_z.
    flux_x = np.tensordot(zweights, u * mprime, axes=(0, 0))
    flux_y = np.tensordot(zweights, v * mprime, axes=(0, 0))
    forcing_hat = -(1j * kx * fft2f(flux_x) + 1j * ky * fft2f(flux_y))
    exact_mode = np.real(np.conj(Mhat) * forcing_hat)
    exact_shell = bin_sum(exact_mode, bin_id, nbins, dk)
    eddy_forcing = exact_shell - row_sum

    divergence = ifft2f(1j * kx * fft2f(U) + 1j * ky * fft2f(V))
    uv_scale = max(float(np.sqrt(np.mean(U * U + V * V))), 1.0e-300)
    return {
        "M_2D_variance": float(0.5 * np.mean(M * M)),
        "T_2D": matrix,
        "T_2D_row_sum": row_sum,
        "T_2D_exchange": exchange,
        "T_2D_exchange_row_sum": exchange.sum(axis=1),
        "T_2D_symmetric": symmetric,
        "T_2D_exact_shell": exact_shell,
        "F_3D_to_2D_shell": eddy_forcing,
        "T_2D_antisymmetry_relative_error": matrix_antisymmetry_error(matrix),
        "UV_horizontal_divergence_rms_over_velocity_rms": float(
            np.sqrt(np.mean(divergence * divergence)) / uv_scale
        ),
    }


def write_csv(path: Path, rows: list[dict[str, object]]):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def process_case(label, anchor, output_dir, helpers, nframes, matrix_dk, kmax):
    frames = helpers.select_latest_complete_frames(anchor, nframes)
    params = helpers.read_bou(Path(frames[-1]["field"]).parent.parent / "bou.in")
    shape = helpers.dataset_shape(Path(frames[-1]["mprime"]), "MPRIME_me")
    nz, ny, nx = shape
    coordinate = Path(frames[-1]["field"]).parent / "cordin_info.h5"
    z = helpers.h5dump_array(coordinate, "z")
    zweights = helpers.vertical_weights(z, 0.0, float(params["Lz"]))
    kx, ky, kh, bin_id, centers, edges, masks, complete = geometry(
        nx,
        ny,
        float(params["Lx"]),
        float(params["Ly"]),
        matrix_dk,
        kmax,
    )

    frame_results = []
    for iframe, frame in enumerate(frames, start=1):
        mprime = helpers.h5dump_array(Path(frame["mprime"]), "MPRIME_me", shape)
        # Remove only the tiny stored horizontal-mean roundoff at each z.
        mprime -= mprime.mean(axis=(-2, -1), keepdims=True)
        u = helpers.h5dump_array(Path(frame["horizontal"]), "VX_me", shape)
        v = helpers.h5dump_array(Path(frame["horizontal"]), "VY_me", shape)
        w = helpers.h5dump_array(Path(frame["field"]), "VZ_me", shape)

        result = {"time": float(frame["time"])}
        result.update(
            full_volume_matrix(
                mprime, u, v, w, z, zweights, kx, ky, bin_id, masks, matrix_dk
            )
        )
        result.update(
            vertical_mode_matrix(
                mprime, u, v, zweights, kx, ky, bin_id, masks, matrix_dk
            )
        )
        frame_results.append(result)
        print(f"{label}: t={result['time']:.6g} ({iframe}/{len(frames)})", flush=True)

    times = np.array([item["time"] for item in frame_results], dtype=float)
    array_names = [
        "T_full_horizontal",
        "T_full_vertical",
        "T_full_total",
        "T_full_row_sum",
        "T_full_exchange",
        "T_full_exchange_row_sum",
        "T_full_symmetric",
        "T_full_direct_shell",
        "T_full_truncation_residual",
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
        "full_antisymmetry_relative_error",
        "T_2D_antisymmetry_relative_error",
        "UV_horizontal_divergence_rms_over_velocity_rms",
    ]
    arrays = {
        name + "_frames": np.stack([np.asarray(item[name]) for item in frame_results])
        for name in array_names
    }
    scalars = {
        name + "_frames": np.array([float(item[name]) for item in frame_results])
        for name in scalar_names
    }

    case_dir = output_dir / label
    case_dir.mkdir(parents=True, exist_ok=True)
    npz_path = case_dir / "mse_variance_shell_to_shell_latest_frames.npz"
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
        "vertical_average": "sum_z cell_thickness_weight(z)*field(z); weights sum to one",
        "mprime": "m-<m>_xy at every z; stored MPRIME_me used",
        "T_full": "-<mprime_K u_i partial_i mprime_P>_V/dk; P donor, K receiver",
        "T_2D": "-<M_K (U partial_x+V partial_y) M_P>_xy/dk, M=<mprime>_z, U=<u>_z, V=<v>_z",
        "F_3D_to_2D": "exact -div_h<uh*mprime>_z shell input minus T_2D row sum",
        "averaging": "matrix first for each frame, then arithmetic mean",
        "full_antisymmetry_relative_error_max": float(
            np.max(scalars["full_antisymmetry_relative_error_frames"])
        ),
        "T_2D_antisymmetry_relative_error_max": float(
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
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--nframes", type=int, default=10)
    parser.add_argument("--matrix-dk", type=float, default=math.pi / 2.0)
    parser.add_argument("--kmax", type=float, default=60.0)
    parser.add_argument("--cases", nargs="*", default=None)
    args = parser.parse_args()

    helpers = load_helpers(args.helpers)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    labels = list(helpers.CASE_RELATIVE_PATHS)
    if args.cases:
        requested = set(args.cases)
        labels = [label for label in labels if label in requested]
    manifests = []
    for label in labels:
        anchor = args.base / helpers.CASE_RELATIVE_PATHS[label]
        manifests.append(
            process_case(
                label,
                anchor,
                args.output_dir,
                helpers,
                args.nframes,
                args.matrix_dk,
                args.kmax,
            )
        )
    write_csv(
        args.output_dir / "mse_shell_to_shell_case_summary.csv",
        [
            {
                "case": item["case"],
                "Ek": item["parameters"]["Ek"],
                "Lx": item["parameters"]["Lx"],
                "Ly": item["parameters"]["Ly"],
                "movie_shape_zyx": str(tuple(item["movie_shape_zyx"])),
                "t_first": item["times"][0],
                "t_last": item["times"][-1],
                "nframes": len(item["times"]),
                "matrix_dk": item["matrix_dk"],
                "matrix_kmax_edge": item["matrix_kmax_edge"],
                "n_bins": item["n_bins"],
                "full_antisymmetry_relative_error_max": item[
                    "full_antisymmetry_relative_error_max"
                ],
                "T_2D_antisymmetry_relative_error_max": item[
                    "T_2D_antisymmetry_relative_error_max"
                ],
            }
            for item in manifests
        ],
    )
    (args.output_dir / "manifest_all_cases.json").write_text(
        json.dumps(manifests, indent=2), encoding="utf-8"
    )
    print(f"saved {args.output_dir}", flush=True)


if __name__ == "__main__":
    main()
