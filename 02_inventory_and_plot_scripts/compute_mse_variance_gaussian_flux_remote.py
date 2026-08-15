#!/usr/bin/env python3
"""MSE-variance spectra, budgets, and Gaussian coarse-grained SGS flux.

This script is designed to run on the remote cluster beside the HDF5 movie
files.  It reads only local remote files and writes compact NPZ/CSV products.

Definitions
-----------
m = b + gamma*q, q = RH*exp(alphaqs*(b-betaqs*z))
m' = m - <m>_xy(z)
G_ell(k_h) = exp(-ell**2*k_h**2/2), k_c = 1/ell
Pi_m(ell) = -< tau_i(u_i,m') * d_i m'_ell >_V

For an incompressible, impermeable volume the resolved-advection part
integrates to zero, so the volume-mean flux is evaluated efficiently as

Pi_m = -< (u_i*m')_ell * d_i m'_ell >_V.

Positive Pi_m is forward/downscale MSE-variance transfer; negative Pi_m is
inverse/upscale transfer.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import subprocess
import tempfile
from pathlib import Path

import numpy as np


CASE_RELATIVE_PATHS = {
    "Ek1p5e-4": "Ek1p5e-4/AR4/Beta1p02/qbot0p5_qtop0p004978/N257x257x65",
    "Ek2e-4": "Ek2e-4/AR4/Beta1p02/qbot0p5_qtop0p004978/N257x257x65",
    "Ek3e-3": "Ek3e-3/AR16/Beta1p02/qbot0p5_qtop0p004978/N257x257x65",
    "Ek5e-3": "Ek5e-3/AR16/Beta1p02/qbot0p5_qtop0p004978/N257x257x65",
    "Ek7e-3": "Ek7e-3/AR16/Beta1p02/qbot0p5_qtop0p004978/N257x257x65",
    "Ek1e-2": "Ek1e-2/AR16/Beta1p02/qbot0p5_qtop0p004978/N257x257x65",
    "Ek3e-2": "Ek3e-2/AR16/Beta1p02/qbot0p5_qtop0p004978/N385x385x65",
}


def as_float(token: str) -> float:
    return float(token.replace("D", "e").replace("d", "e"))


def row_after(lines: list[str], marker: str) -> list[str]:
    for index, line in enumerate(lines):
        if marker in line:
            for later in lines[index + 1 :]:
                clean = later.split("!")[0].strip()
                if clean:
                    return clean.split()
    return []


def read_bou(path: Path) -> dict[str, float | int]:
    lines = path.read_text(errors="ignore").splitlines()
    grid = row_after(lines, "N1      N2")
    geometry = row_after(lines, "ALX3D")
    control = row_after(lines, "Ra       Pr")
    if len(grid) < 3 or len(geometry) < 3 or len(control) < 9:
        raise ValueError(f"Cannot parse {path}")
    ra, pr, invro = map(as_float, control[:3])
    ek = math.sqrt(pr / ra) / invro if abs(invro) > 1.0e-15 else math.inf
    return {
        "Nx_bou": int(as_float(grid[0])),
        "Ny_bou": int(as_float(grid[1])),
        "Nz_bou": int(as_float(grid[2])),
        "Lz": as_float(geometry[0]),
        "Lx": as_float(geometry[1]),
        "Ly": as_float(geometry[2]),
        "Ra": ra,
        "Pr": pr,
        "Ek": ek,
        "gamma": as_float(control[4]),
        "Sm": as_float(control[5]),
        "alphaqs": as_float(control[6]),
        "betaqs": as_float(control[7]),
        "tau_cond": as_float(control[8]),
    }


def dataset_shape(path: Path, dataset: str) -> tuple[int, ...]:
    output = subprocess.check_output(
        ["h5dump", "-H", "-d", "/" + dataset, str(path)], text=True
    )
    match = re.search(r"DATASPACE\s+SIMPLE\s*\{\s*\(\s*([^)]*)\)", output)
    if not match:
        raise ValueError(f"Cannot parse shape of {dataset} in {path}")
    return tuple(int(value.strip()) for value in match.group(1).split(","))


def h5dump_array(path: Path, dataset: str, shape: tuple[int, ...] | None = None) -> np.ndarray:
    if shape is None:
        shape = dataset_shape(path, dataset)
    descriptor, temporary = tempfile.mkstemp(suffix=".bin")
    os.close(descriptor)
    try:
        subprocess.run(
            ["h5dump", "-d", "/" + dataset, "-b", "LE", "-o", temporary, str(path)],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=True,
        )
        values = np.fromfile(temporary, dtype="<f8")
    finally:
        try:
            os.remove(temporary)
        except OSError:
            pass
    expected = int(np.prod(shape))
    if values.size != expected:
        raise ValueError(f"{path}:{dataset} has {values.size}, expected {expected}")
    return values.reshape(shape)


def frame_time(xmf: Path) -> float:
    text = xmf.read_text(errors="ignore")
    match = re.search(r'<Time\s+Value="\s*([^\"]+)', text)
    if not match:
        raise ValueError(f"No time in {xmf}")
    return float(match.group(1))


def continuation_rank(path: Path) -> tuple[int, int, str]:
    values = [int(value) for value in re.findall(r"conti(?:nuation)?(\d+)", str(path), re.I)]
    return (sum(values), len(values), str(path))


def select_latest_complete_frames(anchor: Path, nframes: int) -> list[dict[str, object]]:
    merged: dict[float, dict[str, object]] = {}
    for mprime in anchor.rglob("movie/mprime*.h5"):
        suffix = mprime.stem.replace("mprime", "")
        movie = mprime.parent
        field = movie / f"field{suffix}.h5"
        horizontal = movie / f"horizontal_velocity{suffix}.h5"
        xmf = movie / f"mprime{suffix}.xmf"
        if not (field.exists() and horizontal.exists() and xmf.exists()):
            continue
        time = round(frame_time(xmf), 9)
        item = {
            "time": time,
            "field": field,
            "horizontal": horizontal,
            "mprime": mprime,
            "xmf": xmf,
            "rank": continuation_rank(mprime),
        }
        if time not in merged or item["rank"] > merged[time]["rank"]:
            merged[time] = item
    frames = [merged[key] for key in sorted(merged)]
    if len(frames) < nframes:
        raise RuntimeError(f"{anchor}: only {len(frames)} complete fields")
    return frames[-nframes:]


def vertical_weights(z: np.ndarray, lower: float, upper: float) -> np.ndarray:
    edges = np.empty(z.size + 1, dtype=float)
    edges[0] = lower
    edges[-1] = upper
    edges[1:-1] = 0.5 * (z[:-1] + z[1:])
    weights = np.diff(edges)
    if np.any(weights <= 0.0):
        raise ValueError("Non-monotone vertical grid")
    return weights / weights.sum()


def fft2f(field: np.ndarray) -> np.ndarray:
    return np.fft.fft2(field, axes=(-2, -1), norm="forward")


def ifft2f(field: np.ndarray) -> np.ndarray:
    return np.fft.ifft2(field, axes=(-2, -1), norm="forward").real


def shell_sum(mode_values: np.ndarray, shell_id: np.ndarray, n_shell: int, dk: float) -> np.ndarray:
    return np.bincount(
        shell_id.ravel(), weights=mode_values.ravel(), minlength=n_shell
    )[:n_shell] / dk


def reconstruct_mprime(
    b: np.ndarray,
    rh: np.ndarray,
    z: np.ndarray,
    gamma: float,
    alphaqs: float,
    betaqs: float,
) -> tuple[np.ndarray, np.ndarray]:
    temperature = b - betaqs * z[:, None, None]
    qsat = np.exp(alphaqs * temperature)
    qvap = rh * qsat
    moist = b + gamma * qvap
    mean_profile = moist.mean(axis=(-2, -1))
    mprime = moist - mean_profile[:, None, None]
    return mprime, mean_profile


def process_frame(
    frame: dict[str, object],
    params: dict[str, float | int],
    z: np.ndarray,
    zweights: np.ndarray,
    shape: tuple[int, int, int],
    kx: np.ndarray,
    ky: np.ndarray,
    kh2: np.ndarray,
    shell_id: np.ndarray,
    n_shell: int,
    dk: float,
    kc: np.ndarray,
    validate: bool,
) -> dict[str, np.ndarray | float]:
    field = Path(frame["field"])
    horizontal = Path(frame["horizontal"])
    b = h5dump_array(field, "DSAL_me", shape)
    w = h5dump_array(field, "VZ_me", shape)
    rh = h5dump_array(field, "RH_me", shape)
    u = h5dump_array(horizontal, "VX_me", shape)
    v = h5dump_array(horizontal, "VY_me", shape)

    mprime, mean_profile = reconstruct_mprime(
        b,
        rh,
        z,
        float(params["gamma"]),
        float(params["alphaqs"]),
        float(params["betaqs"]),
    )
    mprime_validation_max = math.nan
    if validate:
        stored = h5dump_array(Path(frame["mprime"]), "MPRIME_me", shape)
        mprime_validation_max = float(np.max(np.abs(stored - mprime)))
        del stored

    mhat = fft2f(mprime)
    what = fft2f(w)
    mzhat = np.gradient(mhat, z, axis=0, edge_order=2)
    mean_gradient = np.gradient(mean_profile, z, edge_order=2)

    # Physical nonlinear term N_m = -(u.grad)m'.
    mx = ifft2f(1j * kx[None, :, :] * mhat)
    nonlinear = -u * mx
    del mx
    my = ifft2f(1j * ky[None, :, :] * mhat)
    nonlinear -= v * my
    del my
    mz = np.gradient(mprime, z, axis=0, edge_order=2)
    nonlinear -= w * mz
    nhat = fft2f(nonlinear)
    del nonlinear, mz

    # Volume-integrated horizontal shell spectra/budget densities.
    energy_mode = np.tensordot(zweights, 0.5 * np.abs(mhat) ** 2, axes=(0, 0))
    transfer_mode = np.tensordot(
        zweights, np.real(np.conj(mhat) * nhat), axes=(0, 0)
    )
    production_mode = np.tensordot(
        zweights,
        -mean_gradient[:, None, None] * np.real(np.conj(mhat) * what),
        axes=(0, 0),
    )
    kappa = 1.0 / math.sqrt(float(params["Ra"]) * float(params["Pr"]))
    dissipation_mode = np.tensordot(
        zweights,
        -kappa * (kh2[None, :, :] * np.abs(mhat) ** 2 + np.abs(mzhat) ** 2),
        axes=(0, 0),
    )
    spectrum = shell_sum(energy_mode, shell_id, n_shell, dk)
    transfer = shell_sum(transfer_mode, shell_id, n_shell, dk)
    production = shell_sum(production_mode, shell_id, n_shell, dk)
    dissipation = shell_sum(dissipation_mode, shell_id, n_shell, dk)
    del nhat, energy_mode, transfer_mode, production_mode, dissipation_mode

    # Direct Gaussian SGS MSE-variance flux.  For a closed incompressible
    # volume, <u_l m_l grad m_l>=0, leaving the product-filter contribution.
    cross_h = np.zeros(kh2.shape, dtype=float)
    umhat = fft2f(u * mprime)
    cross_h += np.tensordot(
        zweights,
        np.real(umhat * np.conj(1j * kx[None, :, :] * mhat)),
        axes=(0, 0),
    )
    del umhat
    vmhat = fft2f(v * mprime)
    cross_h += np.tensordot(
        zweights,
        np.real(vmhat * np.conj(1j * ky[None, :, :] * mhat)),
        axes=(0, 0),
    )
    del vmhat
    wmhat = fft2f(w * mprime)
    cross_v = np.tensordot(
        zweights, np.real(wmhat * np.conj(mzhat)), axes=(0, 0)
    )
    del wmhat
    cross_total = cross_h + cross_v

    production_cumulative_mode = np.tensordot(
        zweights,
        -mean_gradient[:, None, None] * np.real(np.conj(mhat) * what),
        axes=(0, 0),
    )
    dissipation_cumulative_mode = np.tensordot(
        zweights,
        kappa * (kh2[None, :, :] * np.abs(mhat) ** 2 + np.abs(mzhat) ** 2),
        axes=(0, 0),
    )
    variance_cumulative_mode = np.tensordot(
        zweights, 0.5 * np.abs(mhat) ** 2, axes=(0, 0)
    )

    pi_gaussian = np.empty_like(kc)
    p_gaussian = np.empty_like(kc)
    chi_gaussian = np.empty_like(kc)
    variance_gaussian = np.empty_like(kc)
    for index, cutoff in enumerate(kc):
        ell = 1.0 / cutoff
        gaussian_squared = np.exp(-ell * ell * kh2)
        pi_gaussian[index] = -float(np.sum(gaussian_squared * cross_total))
        p_gaussian[index] = float(np.sum(gaussian_squared * production_cumulative_mode))
        chi_gaussian[index] = float(np.sum(gaussian_squared * dissipation_cumulative_mode))
        variance_gaussian[index] = float(np.sum(gaussian_squared * variance_cumulative_mode))

    output = {
        "time": float(frame["time"]),
        "spectrum": spectrum,
        "transfer": transfer,
        "production": production,
        "dissipation": dissipation,
        "pi_gaussian": pi_gaussian,
        "production_gaussian": p_gaussian,
        "dissipation_gaussian": chi_gaussian,
        "variance_gaussian": variance_gaussian,
        "mprime_variance": float(np.sum(variance_cumulative_mode)),
        "mprime_validation_max_abs": mprime_validation_max,
        "transfer_sum": float(np.sum(transfer) * dk),
        "production_sum": float(np.sum(production) * dk),
        "dissipation_sum": float(np.sum(dissipation) * dk),
    }
    return output


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def process_case(label: str, anchor: Path, output_dir: Path, nframes: int, nkc: int) -> dict[str, object]:
    frames = select_latest_complete_frames(anchor, nframes)
    params = read_bou(Path(frames[-1]["field"]).parent.parent / "bou.in")
    shape = dataset_shape(Path(frames[-1]["field"]), "DSAL_me")
    if len(shape) != 3:
        raise ValueError(f"Expected 3D movie field, got {shape}")
    nz, ny, nx = shape
    coordinate = Path(frames[-1]["field"]).parent / "cordin_info.h5"
    z = h5dump_array(coordinate, "z")
    if z.size != nz:
        raise ValueError(f"z size {z.size}, field nz {nz}")
    zweights = vertical_weights(z, 0.0, float(params["Lz"]))

    lx, ly = float(params["Lx"]), float(params["Ly"])
    dx, dy = lx / nx, ly / ny
    kx_1d = 2.0 * np.pi * np.fft.fftfreq(nx, d=dx)
    ky_1d = 2.0 * np.pi * np.fft.fftfreq(ny, d=dy)
    kx, ky = np.meshgrid(kx_1d, ky_1d, indexing="xy")
    kh2 = kx * kx + ky * ky
    kh = np.sqrt(kh2)
    dk = min(2.0 * np.pi / lx, 2.0 * np.pi / ly)
    shell_id = np.floor(kh / dk + 1.0e-12).astype(int)
    complete_kmax = 0.90 * min(np.pi / dx, np.pi / dy)
    n_shell = int(math.floor(complete_kmax / dk)) + 1
    k_shell = dk * np.arange(n_shell)
    kc = np.logspace(math.log10(dk), math.log10(complete_kmax), nkc)

    frame_results: list[dict[str, np.ndarray | float]] = []
    for index, frame in enumerate(frames):
        result = process_frame(
            frame,
            params,
            z,
            zweights,
            shape,
            kx,
            ky,
            kh2,
            shell_id,
            n_shell,
            dk,
            kc,
            validate=index == 0,
        )
        frame_results.append(result)
        print(f"{label}: t={result['time']:.6g} ({index + 1}/{len(frames)})", flush=True)

    def stack(name: str) -> np.ndarray:
        return np.stack([np.asarray(result[name]) for result in frame_results], axis=0)

    times = np.array([float(result["time"]) for result in frame_results])
    spectrum_frames = stack("spectrum")
    transfer_frames = stack("transfer")
    production_frames = stack("production")
    dissipation_frames = stack("dissipation")
    pi_frames = stack("pi_gaussian")
    pg_frames = stack("production_gaussian")
    chig_frames = stack("dissipation_gaussian")
    vg_frames = stack("variance_gaussian")
    dt = times[-1] - times[0]
    dvariance_dt = (vg_frames[-1] - vg_frames[0]) / dt
    gaussian_closure = pg_frames.mean(axis=0) - chig_frames.mean(axis=0) - pi_frames.mean(axis=0) - dvariance_dt

    spectral_dedt = (spectrum_frames[-1] - spectrum_frames[0]) / dt
    spectral_closure = (
        transfer_frames.mean(axis=0)
        + production_frames.mean(axis=0)
        + dissipation_frames.mean(axis=0)
        - spectral_dedt
    )

    case_dir = output_dir / label
    case_dir.mkdir(parents=True, exist_ok=True)
    npz_path = case_dir / "mse_variance_spectrum_budget_gaussian_flux.npz"
    np.savez_compressed(
        npz_path,
        times=times,
        k_shell=k_shell,
        dk=dk,
        spectrum_frames=spectrum_frames,
        transfer_frames=transfer_frames,
        production_frames=production_frames,
        dissipation_frames=dissipation_frames,
        spectral_dedt=spectral_dedt,
        spectral_closure=spectral_closure,
        k_c=kc,
        ell=1.0 / kc,
        pi_gaussian_frames=pi_frames,
        production_gaussian_frames=pg_frames,
        dissipation_gaussian_frames=chig_frames,
        variance_gaussian_frames=vg_frames,
        dvariance_dt=dvariance_dt,
        gaussian_closure=gaussian_closure,
        z=z,
        zweights=zweights,
    )

    spectral_rows = []
    for index, wave in enumerate(k_shell):
        spectral_rows.append(
            {
                "case": label,
                "Ek": params["Ek"],
                "k_h": wave,
                "E_m": spectrum_frames[:, index].mean(),
                "T_m": transfer_frames[:, index].mean(),
                "P_m": production_frames[:, index].mean(),
                "D_m": dissipation_frames[:, index].mean(),
                "dE_m_dt": spectral_dedt[index],
                "closure_residual": spectral_closure[index],
            }
        )
    write_csv(case_dir / "mse_variance_spectral_budget_mean.csv", spectral_rows)

    gaussian_rows = []
    for index, cutoff in enumerate(kc):
        gaussian_rows.append(
            {
                "case": label,
                "Ek": params["Ek"],
                "k_c": cutoff,
                "ell": 1.0 / cutoff,
                "Pi_m": pi_frames[:, index].mean(),
                "P_m_resolved": pg_frames[:, index].mean(),
                "chi_m_resolved": chig_frames[:, index].mean(),
                "V_m_resolved": vg_frames[:, index].mean(),
                "dV_m_dt": dvariance_dt[index],
                "closure_residual": gaussian_closure[index],
            }
        )
    write_csv(case_dir / "mse_variance_gaussian_budget_mean.csv", gaussian_rows)

    validation = float(frame_results[0]["mprime_validation_max_abs"])
    transfer_totals = [float(item["transfer_sum"]) for item in frame_results]
    manifest = {
        "case": label,
        "anchor": str(anchor),
        "parameters": params,
        "movie_shape_zyx": shape,
        "frames": [
            {
                "time": float(frame["time"]),
                "field": str(frame["field"]),
                "horizontal_velocity": str(frame["horizontal"]),
            }
            for frame in frames
        ],
        "mprime_reconstruction_max_abs_error_first_frame": validation,
        "filter": "G_ell(k_h)=exp(-ell^2*k_h^2/2), k_c=1/ell",
        "flux": "Pi_m=-<tau_i(u_i,mprime)*partial_i(mprime_ell)>_V",
        "flux_sign": "Pi_m>0 forward/downscale; Pi_m<0 inverse/upscale",
        "averaging": "compute every frame first, then average latest frames",
        "transfer_sum_per_frame": transfer_totals,
        "npz": str(npz_path),
    }
    (case_dir / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--nframes", type=int, default=10)
    parser.add_argument("--nkc", type=int, default=50)
    parser.add_argument("--cases", nargs="*", default=list(CASE_RELATIVE_PATHS))
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifests = []
    for label in args.cases:
        if label not in CASE_RELATIVE_PATHS:
            raise KeyError(f"Unknown case {label}")
        anchor = args.root / CASE_RELATIVE_PATHS[label]
        manifests.append(process_case(label, anchor, args.output_dir, args.nframes, args.nkc))
    (args.output_dir / "all_cases_manifest.json").write_text(
        json.dumps(manifests, indent=2), encoding="utf-8"
    )
    print(f"Saved compact outputs to {args.output_dir}", flush=True)


if __name__ == "__main__":
    main()
