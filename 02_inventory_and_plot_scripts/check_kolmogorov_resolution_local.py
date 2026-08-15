from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path

import h5py
import numpy as np


VELOCITY_FILES = {
    "u": ("continua_q1.h5", "Vth"),
    "v": ("continua_q2.h5", "Vr"),
    "w": ("continua_q3.h5", "Vz"),
}


def parse_bou(path: Path) -> dict[str, float]:
    def numbers(line: str) -> list[float]:
        return [float(token.replace("d", "e").replace("D", "e")) for token in line.split()]

    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    values: dict[str, float] = {}
    for index, line in enumerate(lines):
        key = line.strip().lower()
        if key.startswith("n1") and "n2" in key and "n3" in key:
            values["N1"], values["N2"], values["N3"] = numbers(lines[index + 1])[:3]
        elif key.startswith("alx3d"):
            _, values["Lx"], values["Ly"] = numbers(lines[index + 1])[:3]
        elif key.startswith("ra") and "invro" in key:
            row = numbers(lines[index + 1])
            values.update(Ra=row[0], Pr=row[1], invRo=row[2], Sm=row[5])
        elif key.startswith("dtmax") and "cflmax" in key:
            values["dt"] = numbers(lines[index + 1])[3]
    required = {"N1", "N2", "N3", "Lx", "Ly", "Ra", "Pr", "invRo", "Sm", "dt"}
    missing = required - values.keys()
    if missing:
        raise ValueError(f"missing bou.in values: {sorted(missing)}")
    return values


def parse_ek(case_name: str, inv_ro: float, ra: float, pr: float) -> float:
    if case_name.lower() in {"norotating", "nonrotating"} or inv_ro == 0.0:
        return math.nan
    match = re.fullmatch(r"Ek(?P<mant>[0-9]+(?:p[0-9]+)?)e-(?P<exp>[0-9]+)", case_name)
    if match:
        mantissa = float(match.group("mant").replace("p", "."))
        return mantissa * 10.0 ** (-int(match.group("exp")))
    return 1.0 / (inv_ro * math.sqrt(ra / pr))


def z_cell_widths(z: np.ndarray) -> np.ndarray:
    edges = np.empty(z.size + 1)
    edges[1:-1] = 0.5 * (z[:-1] + z[1:])
    edges[0] = 0.0
    edges[-1] = 1.0
    return np.diff(edges)


def find_compatible_grid(run: Path, params: dict[str, float]) -> Path:
    direct = run / "field_gridc.h5"
    if direct.exists() and direct.stat().st_size > 0:
        return direct
    for candidate in sorted(run.parent.parent.glob("*/run/field_gridc.h5")):
        if candidate.stat().st_size == 0:
            continue
        other = parse_bou(candidate.parent / "bou.in")
        if all(other[key] == params[key] for key in ("N1", "N2", "N3", "Lx", "Ly")):
            return candidate
    raise FileNotFoundError(f"no compatible nonempty field_gridc.h5 for {run}")


def load_centered_velocity(
    run: Path, params: dict[str, float]
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, Path]:
    arrays = {}
    for component, (filename, dataset) in VELOCITY_FILES.items():
        with h5py.File(run / filename, "r") as handle:
            arrays[component] = handle[dataset][:].astype(np.float64)

    grid_path = find_compatible_grid(run, params)
    with h5py.File(grid_path, "r") as handle:
        zc = handle["zc"][:].astype(np.float64)
        zf = handle["zf"][:].astype(np.float64)

    # Restart arrays retain the duplicated periodic endpoint in x and y.
    ny = arrays["u"].shape[1] - 1
    nx = arrays["u"].shape[2] - 1
    u = 0.5 * (arrays["u"][:, :ny, :nx] + arrays["u"][:, :ny, 1 : nx + 1])
    v = 0.5 * (arrays["v"][:, :ny, :nx] + arrays["v"][:, 1 : ny + 1, :nx])

    w_face = arrays["w"][:, :ny, :nx]
    w_face = np.concatenate((w_face, np.zeros((1, ny, nx), dtype=w_face.dtype)), axis=0)
    zf_full = np.concatenate((zf, [1.0]))
    lower_weight = (zf_full[1:] - zc) / np.diff(zf_full)
    upper_weight = (zc - zf_full[:-1]) / np.diff(zf_full)
    w = lower_weight[:, None, None] * w_face[:-1] + upper_weight[:, None, None] * w_face[1:]
    return u, v, w, zc, grid_path


def horizontal_derivatives(field: np.ndarray, length: float, axis: int) -> np.ndarray:
    n = field.shape[axis]
    spacing = length / n
    wave = 2.0 * np.pi * np.fft.fftfreq(n, d=spacing)
    shape = [1] * field.ndim
    shape[axis] = n
    transformed = np.fft.fft(field, axis=axis)
    return np.fft.ifft(1j * wave.reshape(shape) * transformed, axis=axis).real


def velocity_gradients(
    u: np.ndarray, v: np.ndarray, w: np.ndarray, z: np.ndarray, lx: float, ly: float
) -> tuple[np.ndarray, ...]:
    derivatives = []
    for field in (u, v, w):
        derivatives.extend(
            (
                horizontal_derivatives(field, lx, axis=2),
                horizontal_derivatives(field, ly, axis=1),
                np.gradient(field, z, axis=0, edge_order=2),
            )
        )
    return tuple(derivatives)


def profile_mean(field: np.ndarray) -> np.ndarray:
    return np.mean(field, axis=(1, 2), dtype=np.float64)


def compute_case(case_dir: Path) -> dict[str, object]:
    run = case_dir / "run"
    params = parse_bou(run / "bou.in")
    u, v, w, z, grid_path = load_centered_velocity(run, params)
    dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz = velocity_gradients(
        u, v, w, z, params["Lx"], params["Ly"]
    )
    del u, v, w

    nu = math.sqrt(params["Pr"] / params["Ra"])
    strain_square = dudx**2 + dvdy**2 + dwdz**2
    strain_square += 0.5 * (dudy + dvdx) ** 2
    strain_square += 0.5 * (dudz + dwdx) ** 2
    strain_square += 0.5 * (dvdz + dwdy) ** 2
    epsilon_strain = 2.0 * nu * profile_mean(strain_square)
    del strain_square

    epsilon_gradient = nu * sum(
        profile_mean(derivative**2)
        for derivative in (dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz)
    )
    del dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz

    dz = z_cell_widths(z)
    weights = dz / dz.sum()
    epsilon_global = float(np.sum(epsilon_strain * weights))
    epsilon_gradient_global = float(np.sum(epsilon_gradient * weights))
    eta = (nu**3 / np.maximum(epsilon_strain, 1.0e-300)) ** 0.25
    eta_b = eta / math.sqrt(params["Pr"])

    nx = int(params["N1"] - 1)
    ny = int(params["N2"] - 1)
    dx = params["Lx"] / nx
    dy = params["Ly"] / ny
    eta_min = float(eta.min())
    eta_b_min = float(eta_b.min())
    wall = (z <= 0.1) | (z >= 0.9)
    ratio_z = dz / eta
    ratio_z_b = dz / eta_b
    kmax = min((nx / 3.0) * 2.0 * np.pi / params["Lx"], (ny / 3.0) * 2.0 * np.pi / params["Ly"])

    ek = parse_ek(case_dir.name, params["invRo"], params["Ra"], params["Pr"])
    delta_e = math.sqrt(ek) if math.isfinite(ek) else math.nan
    if math.isfinite(delta_e):
        ek_mask = (z <= delta_e) | (z >= 1.0 - delta_e)
        n_ek_bottom = int(np.count_nonzero(z <= delta_e))
        n_ek_top = int(np.count_nonzero(z >= 1.0 - delta_e))
    else:
        ek_mask = np.zeros_like(z, dtype=bool)
        n_ek_bottom = n_ek_top = 0

    kmax_eta = kmax * eta_min
    wall_z_eta = float(np.max(ratio_z[wall]))
    tau_eta_min = float(np.sqrt(nu / np.max(epsilon_strain)))
    dt_over_tau_eta_min = params["dt"] / tau_eta_min
    if kmax_eta >= 1.0 and wall_z_eta <= 1.0:
        spectral_verdict = "resolved"
    elif kmax_eta >= 0.7 and wall_z_eta <= 2.0:
        spectral_verdict = "marginal"
    else:
        spectral_verdict = "under-resolved"

    max_spacing_eta = max(dx / eta_min, dy / eta_min, float(ratio_z.max()))
    max_spacing_eta_b = max(dx / eta_b_min, dy / eta_b_min, float(ratio_z_b.max()))
    zhang_velocity_pass = max_spacing_eta <= 0.57
    zhang_scalar_pass = max_spacing_eta_b <= 0.48
    ekman_points_pass = n_ek_bottom >= 8 and n_ek_top >= 8 if math.isfinite(ek) else False
    time_resolution_pass = dt_over_tau_eta_min < 0.01
    if not math.isfinite(ek):
        combined_status = "not_confirmed_nonrotating_BL"
    elif spectral_verdict == "resolved" and ekman_points_pass and time_resolution_pass:
        combined_status = "resolved"
    elif spectral_verdict in {"resolved", "marginal"} and n_ek_bottom >= 5 and n_ek_top >= 5:
        combined_status = "marginal"
    else:
        combined_status = "under-resolved"
    return {
        "case": case_dir.name,
        "source_run": str(run),
        "grid_source": str(grid_path),
        "Ra": params["Ra"],
        "Pr": params["Pr"],
        "Ek": ek,
        "AR": params["Lx"],
        "grid_effective": f"{nx}x{ny}x{z.size}",
        "nu": nu,
        "epsilon_strain_global": epsilon_global,
        "epsilon_gradient_global": epsilon_gradient_global,
        "epsilon_form_relative_difference": abs(epsilon_global - epsilon_gradient_global)
        / max(epsilon_global, 1.0e-300),
        "eta_K_global": float((nu**3 / epsilon_global) ** 0.25),
        "eta_K_min": eta_min,
        "z_eta_K_min": float(z[np.argmin(eta)]),
        "eta_B_min_Pr": eta_b_min,
        "dt": params["dt"],
        "tau_eta_min": tau_eta_min,
        "dt_over_tau_eta_min": dt_over_tau_eta_min,
        "dx_over_eta_K_min": dx / eta_min,
        "dy_over_eta_K_min": dy / eta_min,
        "max_dz_over_eta_K_all": float(ratio_z.max()),
        "max_dz_over_eta_K_wall01": wall_z_eta,
        "max_dz_over_eta_B_all": float(ratio_z_b.max()),
        "kmax_dealias_eta_K_min": kmax_eta,
        "delta_E_sqrtEk": delta_e,
        "N_E_bottom": n_ek_bottom,
        "N_E_top": n_ek_top,
        "max_dz_over_eta_K_EkBL": float(np.max(ratio_z[ek_mask])) if np.any(ek_mask) else math.nan,
        "spectral_screening": spectral_verdict,
        "ekman_points_ge_8": ekman_points_pass,
        "time_resolution_dt_over_tau_lt_0p01": time_resolution_pass,
        "combined_resolution_status": combined_status,
        "zhang2017_velocity_spacing_pass": zhang_velocity_pass,
        "zhang2017_scalar_spacing_pass_Pr": zhang_scalar_pass,
        "zhang2017_max_spacing_over_eta_K": max_spacing_eta,
        "zhang2017_max_spacing_over_eta_B": max_spacing_eta_b,
    }


def write_markdown(rows: list[dict[str, object]], path: Path) -> None:
    columns = [
        "case",
        "Ek",
        "AR",
        "grid_effective",
        "eta_K_min",
        "kmax_dealias_eta_K_min",
        "max_dz_over_eta_K_wall01",
        "N_E_bottom",
        "spectral_screening",
        "combined_resolution_status",
        "zhang2017_velocity_spacing_pass",
    ]
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
    for row in rows:
        values = []
        for column in columns:
            value = row[column]
            values.append(f"{value:.5g}" if isinstance(value, float) else str(value))
        lines.append("| " + " | ".join(values) + " |")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Audit restart-field resolution using local Kolmogorov scales.")
    parser.add_argument("--cases-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for case_dir in sorted(path for path in args.cases_root.iterdir() if path.is_dir()):
        print(f"Computing {case_dir.name}", flush=True)
        try:
            rows.append(compute_case(case_dir))
        except Exception as error:
            rows.append({"case": case_dir.name, "source_run": str(case_dir / "run"), "error": repr(error)})

    fieldnames = sorted({key for row in rows for key in row})
    csv_path = args.output_dir / "kolmogorov_resolution_current_cases_restart.csv"
    with csv_path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    valid_rows = [row for row in rows if "error" not in row]
    write_markdown(valid_rows, args.output_dir / "kolmogorov_resolution_current_cases_restart.md")
    print(csv_path)


if __name__ == "__main__":
    main()
