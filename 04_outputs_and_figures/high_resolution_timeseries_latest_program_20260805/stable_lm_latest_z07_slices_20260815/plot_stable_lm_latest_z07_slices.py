from __future__ import annotations

import csv
import math
import re
from pathlib import Path

import h5py
import matplotlib as mpl
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize, TwoSlopeNorm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


OUT = Path(__file__).resolve().parent
Z_TARGET = 0.7
GAMMA = 1.1
STABLE_CASES = {
    "Ek1e-3": (1.0e-3, Path(r"G:\moist convection\Ek1e-3")),
    "Ek3e-3": (3.0e-3, Path(r"G:\moist convection\Ek3e-3")),
    "Ek5e-3": (5.0e-3, Path(r"G:\moist convection\Ek5e-3")),
    "Ek7e-3": (7.0e-3, Path(r"G:\moist convection\Ek7e-3")),
    "Ek1e-2": (1.0e-2, Path(r"G:\moist convection\Ek1e-2")),
}


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "axes.linewidth": 4.5,
            "axes.labelsize": 24,
            "axes.titlesize": 22,
            "xtick.labelsize": 13,
            "ytick.labelsize": 13,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


def valid(path: Path) -> bool:
    text = str(path).lower()
    return "qbot0p5_qtop0p004978" in text and "\\run\\run\\" not in text


def read_numeric(path: Path, minimum_columns: int) -> np.ndarray:
    rows = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            try:
                values = [float(item.replace("D", "E").replace("d", "e")) for item in line.split()]
            except ValueError:
                continue
            if len(values) >= minimum_columns and np.all(np.isfinite(values[:minimum_columns])):
                rows.append(values)
    if not rows:
        raise ValueError(f"No numeric rows in {path}")
    return np.asarray(rows, dtype=float)


def interpolate(dataset: h5py.Dataset, z: np.ndarray, target: float) -> np.ndarray:
    hi = int(np.searchsorted(z, target))
    hi = min(max(hi, 1), len(z) - 1)
    lo = hi - 1
    weight = float((target - z[lo]) / (z[hi] - z[lo]))
    return (1.0 - weight) * np.asarray(dataset[lo], dtype=float) + weight * np.asarray(dataset[hi], dtype=float)


def frame_index(path: Path) -> int:
    match = re.search(r"(\d{5})\.h5$", path.name)
    return int(match.group(1)) if match else -1


def h5_has(path: Path, keys: tuple[str, ...]) -> bool:
    try:
        with h5py.File(path, "r") as handle:
            return all(key in handle for key in keys)
    except OSError:
        return False


def latest_movie(root: Path) -> tuple[Path, Path, Path, Path, float] | None:
    candidates = []
    for field in root.rglob("field[0-9][0-9][0-9][0-9][0-9].h5"):
        if not valid(field):
            continue
        movie = field.parent
        index = frame_index(field)
        mprime = movie / f"mprime{index:05d}.h5"
        coordinate = movie / "cordin_info.h5"
        run = movie.parent
        if not (mprime.exists() and coordinate.exists()):
            continue
        if not h5_has(field, ("DSAL_me", "VZ_me")) or not h5_has(mprime, ("MPRIME_me",)):
            continue
        candidates.append((index, field.stat().st_mtime, field, mprime, coordinate, run, index * 10.0))
    if not candidates:
        return None
    _, _, field, mprime, coordinate, run, time = max(candidates, key=lambda item: (item[0], item[1]))
    return field, mprime, coordinate, run, time


def latest_continuation(root: Path) -> tuple[Path, Path, Path, Path, Path, float] | None:
    candidates = []
    for dsal in root.rglob("continua_dsal.h5"):
        if not valid(dsal):
            continue
        run = dsal.parent
        qvap = run / "continua_qvap.h5"
        q3 = run / "continua_q3.h5"
        grid = run / "field_gridc.h5"
        convective = run / "diagnostics" / "scale" / "convective_scale.dat"
        moist = run / "diagnostics" / "scale" / "moist_integral_scale.dat"
        if not all(path.exists() for path in (qvap, q3, grid, convective, moist)):
            continue
        time = min(read_numeric(convective, 3)[-1, 0], read_numeric(moist, 4)[-1, 0])
        candidates.append((time, dsal.stat().st_mtime, dsal, qvap, q3, grid, run))
    if not candidates:
        return None
    time, _, dsal, qvap, q3, grid, run = max(candidates, key=lambda item: (item[0], item[1]))
    return dsal, qvap, q3, grid, run, float(time)


def scales_at(run: Path, time: float) -> tuple[float, float, float, float]:
    moist_path = run / "diagnostics" / "scale" / "moist_integral_scale.dat"
    conv_path = run / "diagnostics" / "scale" / "convective_scale.dat"
    moist = read_numeric(moist_path, 4)
    conv = read_numeric(conv_path, 3)
    im = int(np.argmin(np.abs(moist[:, 0] - time)))
    ic = int(np.argmin(np.abs(conv[:, 0] - time)))
    lambda_m = 2.0 * math.pi * float(moist[im, 1])
    lambda_c = 2.0 * math.pi * float(conv[ic, 1])
    return lambda_m, lambda_c, float(moist[im, 0]), float(conv[ic, 0])


def load_case(label: str, ek: float, root: Path) -> dict[str, object]:
    movie = latest_movie(root)
    if movie is not None:
        field_path, mprime_path, coordinate_path, run, time = movie
        with h5py.File(coordinate_path, "r") as coordinates:
            x = np.asarray(coordinates["x"], dtype=float)
            y = np.asarray(coordinates["y"], dtype=float)
            z = np.asarray(coordinates["z"], dtype=float)
        with h5py.File(field_path, "r") as field, h5py.File(mprime_path, "r") as moist:
            b = interpolate(field["DSAL_me"], z, Z_TARGET)
            w = interpolate(field["VZ_me"], z, Z_TARGET)
            mprime = interpolate(moist["MPRIME_me"], z, Z_TARGET)
        source_type = "latest movie frame"
        sources = f"{field_path}; {mprime_path}"
    else:
        continuation = latest_continuation(root)
        if continuation is None:
            raise FileNotFoundError(f"No usable latest field under {root}")
        dsal_path, qvap_path, q3_path, grid_path, run, time = continuation
        with h5py.File(grid_path, "r") as coordinates:
            x = np.asarray(coordinates["xc"], dtype=float)[:-1]
            y = np.asarray(coordinates["yc"], dtype=float)[:-1]
            zc = np.asarray(coordinates["zc"], dtype=float)
            zf = np.asarray(coordinates["zf"], dtype=float)
        with h5py.File(dsal_path, "r") as dsal, h5py.File(qvap_path, "r") as qvap, h5py.File(q3_path, "r") as q3:
            b = interpolate(dsal["dsal"], zc, Z_TARGET)[:-1, :-1]
            q = interpolate(qvap["qvap"], zc, Z_TARGET)[:-1, :-1]
            w = interpolate(q3["Vz"], zf, Z_TARGET)[:-1, :-1]
        m = b + GAMMA * q
        mprime = m - np.mean(m)
        source_type = "latest continuation snapshot"
        sources = f"{dsal_path}; {qvap_path}; {q3_path}"

    lambda_m, lambda_c, lm_time, lc_time = scales_at(run, time)
    return {
        "label": label,
        "Ek": ek,
        "time": time,
        "x": x,
        "y": y,
        "mprime": np.asarray(mprime, dtype=float),
        "b": np.asarray(b, dtype=float),
        "w": np.asarray(w, dtype=float),
        "lambda_m": lambda_m,
        "lambda_c": lambda_c,
        "lm_time": lm_time,
        "lc_time": lc_time,
        "source_type": source_type,
        "sources": sources,
        "run": str(run),
    }


def robust_limits(cases: list[dict[str, object]]) -> dict[str, tuple[float, float]]:
    m_abs = np.concatenate([np.abs(case["mprime"]).ravel() for case in cases])
    w_abs = np.concatenate([np.abs(case["w"]).ravel() for case in cases])
    b_all = np.concatenate([case["b"].ravel() for case in cases])
    return {
        "mprime": (-float(np.quantile(m_abs, 0.995)), float(np.quantile(m_abs, 0.995))),
        "w": (-float(np.quantile(w_abs, 0.995)), float(np.quantile(w_abs, 0.995))),
        "b": (float(np.quantile(b_all, 0.005)), float(np.quantile(b_all, 0.995))),
    }


def style_map_axis(ax: plt.Axes, show_y: bool = True) -> None:
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(0.0, 16.0)
    ax.set_ylim(0.0, 16.0)
    ax.set_xticks([0, 4, 8, 12, 16])
    ax.set_yticks([0, 4, 8, 12, 16])
    ax.set_xlabel(r"$x/H$")
    if show_y:
        ax.set_ylabel(r"$y/H$")
    else:
        ax.set_yticklabels([])
    ax.tick_params(which="major", direction="in", top=True, right=True, length=10, width=1.2)
    ax.tick_params(which="minor", direction="in", top=True, right=True, length=5, width=1.0)
    ax.minorticks_on()
    for spine in ax.spines.values():
        spine.set_linewidth(4.5)


def add_scale_bars(ax: plt.Axes, lambda_m: float, lambda_c: float, compact: bool = False) -> None:
    x_end = 15.35
    y_m, y_c = (1.95, 0.85) if not compact else (1.60, 0.65)
    line_m = ax.plot([x_end - lambda_m, x_end], [y_m, y_m], color="#FFCC00", lw=5.0, solid_capstyle="butt", zorder=8)[0]
    line_c = ax.plot([x_end - lambda_c, x_end], [y_c, y_c], color="#00D7FF", lw=5.0, solid_capstyle="butt", zorder=8)[0]
    for line in (line_m, line_c):
        line.set_path_effects([pe.Stroke(linewidth=8.0, foreground="black"), pe.Normal()])
    fontsize = 12 if compact else 14
    text_effect = [pe.Stroke(linewidth=3.0, foreground="white"), pe.Normal()]
    ax.text(x_end, y_m + 0.22, rf"$2\pi L_m={lambda_m:.2f}H$", ha="right", va="bottom", fontsize=fontsize, color="black", path_effects=text_effect, zorder=9)
    ax.text(x_end, y_c + 0.22, rf"$2\pi l_c={lambda_c:.2f}H$", ha="right", va="bottom", fontsize=fontsize, color="black", path_effects=text_effect, zorder=9)


def draw_case(case: dict[str, object], limits: dict[str, tuple[float, float]], output: bool = True) -> plt.Figure:
    fig, axes = plt.subplots(1, 3, figsize=(19.5, 6.65), facecolor="white")
    fig.subplots_adjust(left=0.055, right=0.985, bottom=0.12, top=0.78, wspace=0.12)
    variables = [
        ("mprime", r"$m'(z=0.7)$", "RdBu_r", TwoSlopeNorm(vmin=limits["mprime"][0], vcenter=0.0, vmax=limits["mprime"][1])),
        ("b", r"$b(z=0.7)$", "viridis", Normalize(vmin=limits["b"][0], vmax=limits["b"][1])),
        ("w", r"$w(z=0.7)$", "RdBu_r", TwoSlopeNorm(vmin=limits["w"][0], vcenter=0.0, vmax=limits["w"][1])),
    ]
    for index, (ax, (key, title, cmap, norm)) in enumerate(zip(axes, variables)):
        image = ax.imshow(case[key], origin="lower", extent=(0, 16, 0, 16), interpolation="bilinear", cmap=cmap, norm=norm, rasterized=True)
        style_map_axis(ax, show_y=index == 0)
        # Keep the variable name above the colorbar so neither is obscured.
        ax.text(0.5, 1.155, title, transform=ax.transAxes, ha="center", va="bottom", fontsize=22, clip_on=False)
        cax = inset_axes(ax, width="78%", height="3.2%", loc="upper center", bbox_to_anchor=(0.0, 0.10, 1.0, 1.0), bbox_transform=ax.transAxes, borderpad=0)
        cb = fig.colorbar(image, cax=cax, orientation="horizontal")
        cb.ax.tick_params(labelsize=10, direction="in", length=4, width=0.8)
        cb.outline.set_linewidth(1.2)
    axes[0].text(0.45, 15.45, rf"$Ek={case['Ek']:.0e}$, $t={case['time']:.1f}$", ha="left", va="top", fontsize=17, color="black", bbox=dict(facecolor="white", edgecolor="none", alpha=0.78, pad=2.0))
    add_scale_bars(axes[2], float(case["lambda_m"]), float(case["lambda_c"]))
    if output:
        stem = OUT / f"Ra8e6_{case['label']}_latest_z07_mprime_b_w_with_Lm_Lc"
        fig.savefig(stem.with_suffix(".png"), dpi=300)
        fig.savefig(stem.with_suffix(".pdf"), dpi=300)
    return fig


def draw_overview(cases: list[dict[str, object]], limits: dict[str, tuple[float, float]]) -> None:
    fig, axes = plt.subplots(len(cases), 3, figsize=(15.6, 25.5), facecolor="white")
    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.035, top=0.90, hspace=0.18, wspace=0.08)
    specs = [
        ("mprime", r"$m'(z=0.7)$", "RdBu_r", TwoSlopeNorm(vmin=limits["mprime"][0], vcenter=0.0, vmax=limits["mprime"][1])),
        ("b", r"$b(z=0.7)$", "viridis", Normalize(vmin=limits["b"][0], vmax=limits["b"][1])),
        ("w", r"$w(z=0.7)$", "RdBu_r", TwoSlopeNorm(vmin=limits["w"][0], vcenter=0.0, vmax=limits["w"][1])),
    ]
    first_images = []
    for row, case in enumerate(cases):
        for col, (key, title, cmap, norm) in enumerate(specs):
            ax = axes[row, col]
            image = ax.imshow(case[key], origin="lower", extent=(0, 16, 0, 16), interpolation="bilinear", cmap=cmap, norm=norm, rasterized=True)
            if row == 0:
                ax.text(0.5, 1.14, title, transform=ax.transAxes, ha="center", va="bottom", fontsize=22, clip_on=False)
                first_images.append(image)
            style_map_axis(ax, show_y=col == 0)
            if row < len(cases) - 1:
                ax.set_xticklabels([])
                ax.set_xlabel("")
            if col == 0:
                ax.text(
                    0.35,
                    15.45,
                    rf"$Ek={case['Ek']:.0e}$",
                    ha="left",
                    va="top",
                    fontsize=15,
                    bbox=dict(facecolor="white", edgecolor="none", alpha=0.78, pad=1.8),
                )
        add_scale_bars(axes[row, 2], float(case["lambda_m"]), float(case["lambda_c"]), compact=True)
    for col, image in enumerate(first_images):
        cax = inset_axes(axes[0, col], width="82%", height="3.0%", loc="upper center", bbox_to_anchor=(0.0, 0.095, 1.0, 1.0), bbox_transform=axes[0, col].transAxes, borderpad=0)
        cb = fig.colorbar(image, cax=cax, orientation="horizontal")
        cb.ax.tick_params(labelsize=9, direction="in", length=3)
        cb.outline.set_linewidth(1.0)
    fig.savefig(OUT / "Ra8e6_stable_Lm_cases_latest_z07_mprime_b_w_overview.png", dpi=240)
    fig.savefig(OUT / "Ra8e6_stable_Lm_cases_latest_z07_mprime_b_w_overview.pdf", dpi=300)
    plt.close(fig)


def save_summary(cases: list[dict[str, object]]) -> None:
    rows = []
    for case in cases:
        rows.append(
            {
                "case": case["label"],
                "Ek": case["Ek"],
                "field_time": case["time"],
                "z": Z_TARGET,
                "two_pi_Lm_H": case["lambda_m"],
                "two_pi_lc_z05_H": case["lambda_c"],
                "Lm_scale_time": case["lm_time"],
                "lc_scale_time": case["lc_time"],
                "field_source_type": case["source_type"],
                "run": case["run"],
                "field_sources": case["sources"],
            }
        )
    with (OUT / "Ra8e6_stable_Lm_cases_latest_z07_slice_scale_summary.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    configure_style()
    cases = [load_case(label, ek, root) for label, (ek, root) in STABLE_CASES.items()]
    limits = robust_limits(cases)
    for case in cases:
        figure = draw_case(case, limits)
        plt.close(figure)
    draw_overview(cases, limits)
    save_summary(cases)
    print(OUT)


if __name__ == "__main__":
    main()
