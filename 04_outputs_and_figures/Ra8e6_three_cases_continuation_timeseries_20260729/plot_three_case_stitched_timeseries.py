from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
RAW = HERE / "raw_segments"

CASES = [
    {
        "key": "Ek7e-3_AR16",
        "label": r"$Ek=7\times10^{-3},\ AR=16$",
        "color": (0.74, 0.14, 0.18),
    },
    {
        "key": "Ek1p5e-4_AR4",
        "label": r"$Ek=1.5\times10^{-4},\ AR=4$",
        "color": (0.11, 0.44, 0.71),
    },
]

FILE_TYPES = {
    "avgvar": "avgvar.out",
    "convective": "convective_scale.dat",
    "moist": "moist_integral_scale.dat",
    "mprime": "mprime_square.dat",
}

FIGSIZE = (6.5, 5.84)
AX_RECT = [0.186, 0.260, 0.792, 0.705]
LINE_WIDTH = 3.5
TWO_PI = 2.0 * np.pi
LINEAR_INTERVALS = {"Ek7e-3_AR16": (8.5, 28.6)}


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Times New Roman",
            "mathtext.fontset": "stix",
            "axes.linewidth": 4.5,
            "axes.labelsize": 24,
            "xtick.labelsize": 13,
            "ytick.labelsize": 13,
            "legend.fontsize": 10,
            "axes.formatter.use_mathtext": True,
            "figure.dpi": 180,
            "savefig.dpi": 300,
        }
    )


def make_axes() -> tuple[plt.Figure, plt.Axes]:
    fig = plt.figure(figsize=FIGSIZE, facecolor="white")
    ax = fig.add_axes(AX_RECT)
    ax.set_box_aspect(5.2 / 6.5)
    for spine in ax.spines.values():
        spine.set_linewidth(4.5)
    ax.tick_params(
        which="major",
        direction="in",
        length=12,
        width=1.2,
        top=True,
        right=True,
    )
    ax.minorticks_on()
    ax.tick_params(
        which="minor",
        direction="in",
        length=6,
        width=1.0,
        top=True,
        right=True,
    )
    return fig, ax


def clean_rows(data: np.ndarray) -> np.ndarray:
    data = np.atleast_2d(data)
    mask = np.all(np.isfinite(data), axis=1)
    data = data[mask]
    if data.size == 0:
        raise ValueError("No finite rows were found.")
    return data


def segment_paths(case_key: str, suffix: str) -> list[Path]:
    case_dir = RAW / case_key
    paths = [case_dir / f"run_{suffix}"]
    paths.extend(sorted(case_dir.glob(f"conti*_run_{suffix}")))
    return [path for path in paths if path.exists()]


def stitch(case_key: str, suffix: str) -> tuple[np.ndarray, dict[str, object]]:
    paths = segment_paths(case_key, suffix)
    if not paths:
        raise FileNotFoundError(f"No segments found for {case_key}: {suffix}")

    rows: dict[float, tuple[int, np.ndarray]] = {}
    segment_info: list[dict[str, object]] = []
    for priority, path in enumerate(paths):
        data = clean_rows(np.loadtxt(path))
        if data.shape[1] < 2:
            raise ValueError(f"{path} contains fewer than two columns.")
        segment_info.append(
            {
                "file": str(path),
                "rows": int(data.shape[0]),
                "time_start": float(data[0, 0]),
                "time_end": float(data[-1, 0]),
            }
        )
        for row in data:
            key = round(float(row[0]), 8)
            rows[key] = (priority, row.copy())

    stitched = np.vstack([rows[key][1] for key in sorted(rows)])
    if np.any(np.diff(stitched[:, 0]) <= 0.0):
        raise ValueError(f"Non-increasing time after stitching {case_key}: {suffix}")

    metadata = {
        "segments": segment_info,
        "stitched_rows": int(stitched.shape[0]),
        "time_start": float(stitched[0, 0]),
        "time_end": float(stitched[-1, 0]),
        "duplicate_times_removed": int(
            sum(item["rows"] for item in segment_info) - stitched.shape[0]
        ),
    }
    return stitched, metadata


def save_figure(fig: plt.Figure, stem: str) -> None:
    fig.savefig(HERE / f"{stem}.png", facecolor="white")
    fig.savefig(HERE / f"{stem}.pdf", facecolor="white")
    plt.close(fig)


def save_case_figure(fig: plt.Figure, case_key: str, stem: str) -> None:
    output = HERE / "separate_figures" / case_key
    output.mkdir(parents=True, exist_ok=True)
    fig.savefig(output / f"{stem}.png", facecolor="white")
    fig.savefig(output / f"{stem}.pdf", facecolor="white")
    plt.close(fig)


def style_legend(ax: plt.Axes, loc: str = "best", ncol: int = 1) -> None:
    ax.legend(frameon=False, loc=loc, ncol=ncol, handlelength=2.8)


def plot_kinetic(
    datasets: dict[str, dict[str, np.ndarray]], *, log_y: bool
) -> None:
    fig, ax = make_axes()
    for case in CASES:
        data = datasets[case["key"]]["avgvar"]
        kinetic = 0.5 * data[:, 3]
        plot = ax.semilogy if log_y else ax.plot
        plot(
            data[:, 0],
            kinetic,
            color=case["color"],
            lw=LINE_WIDTH,
            label=case["label"],
        )
    ax.set_xlim(left=0.0)
    if not log_y:
        ax.set_ylim(bottom=0.0)
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$K$")
    style_legend(ax, loc="upper right")
    stem = "kinetic_energy_stitched_logy" if log_y else "kinetic_energy_stitched"
    save_figure(fig, stem)


def plot_convective_scale(datasets: dict[str, dict[str, np.ndarray]]) -> None:
    fig, ax = make_axes()
    for case in CASES:
        data = datasets[case["key"]]["convective"]
        ax.plot(
            data[:, 0],
            TWO_PI * data[:, 1],
            color=case["color"],
            lw=LINE_WIDTH,
            label=case["label"],
        )
        ax.plot(
            data[:, 0],
            TWO_PI * data[:, 2],
            color=case["color"],
            lw=LINE_WIDTH,
            ls="--",
            alpha=0.82,
        )
    ax.set_xlim(left=0.0)
    ax.set_ylim(bottom=0.0)
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"Convective wavelength, $2\pi l_c$")
    style_legend(ax, loc="upper right")
    ax.text(
        0.96,
        0.06,
        r"solid: $z\simeq0.5$" "\n" r"dashed: $z$-weighted",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=10,
    )
    save_figure(fig, "convective_scale_stitched")


def plot_moist_integral_scale(
    datasets: dict[str, dict[str, np.ndarray]]
) -> None:
    fig, ax = make_axes()
    for case in CASES:
        data = datasets[case["key"]]["moist"]
        ax.plot(
            data[:, 0],
            TWO_PI * data[:, 1],
            color=case["color"],
            lw=LINE_WIDTH,
            label=case["label"],
        )
    ax.set_xlim(left=0.0)
    ax.set_ylim(bottom=0.0)
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"MSE integral scale, $2\pi L_m$")
    style_legend(ax, loc="upper right")
    save_figure(fig, "moist_integral_scale_stitched")


def plot_individual_cases(
    datasets: dict[str, dict[str, np.ndarray]]
) -> None:
    for case in CASES:
        case_key = case["key"]
        color = case["color"]
        avg = datasets[case_key]["avgvar"]
        convective = datasets[case_key]["convective"]
        moist = datasets[case_key]["moist"]
        mprime = datasets[case_key]["mprime"]
        kinetic = 0.5 * avg[:, 3]

        fig, ax = make_axes()
        ax.plot(avg[:, 0], kinetic, color=color, lw=LINE_WIDTH)
        ax.set_xlim(float(avg[0, 0]), float(avg[-1, 0]))
        ax.set_ylim(bottom=0.0)
        ax.set_xlabel(r"$t$")
        ax.set_ylabel(r"$K$")
        save_case_figure(fig, case_key, "kinetic_energy_stitched")

        fig, ax = make_axes()
        if case_key in LINEAR_INTERVALS:
            ax.axvspan(
                *LINEAR_INTERVALS[case_key],
                color=(0.72, 0.72, 0.72),
                alpha=0.28,
                linewidth=0,
                zorder=0,
            )
        ax.semilogy(avg[:, 0], kinetic, color=color, lw=LINE_WIDTH)
        ax.set_xlim(float(avg[0, 0]), float(avg[-1, 0]))
        ax.set_xlabel(r"$t$")
        ax.set_ylabel(r"$K$")
        save_case_figure(fig, case_key, "kinetic_energy_stitched_logy")

        fig, ax = make_axes()
        if case_key in LINEAR_INTERVALS:
            ax.axvspan(
                *LINEAR_INTERVALS[case_key],
                color=(0.72, 0.72, 0.72),
                alpha=0.28,
                linewidth=0,
                zorder=0,
            )
        ax.plot(
            convective[:, 0],
            TWO_PI * convective[:, 1],
            color=color,
            lw=LINE_WIDTH,
            label=r"$z\simeq0.5$",
        )
        ax.plot(
            convective[:, 0],
            TWO_PI * convective[:, 2],
            color=color,
            lw=LINE_WIDTH,
            ls="--",
            alpha=0.82,
            label=r"$z$-weighted",
        )
        ax.set_xlim(float(convective[0, 0]), float(convective[-1, 0]))
        ax.set_ylim(bottom=0.0)
        ax.set_xlabel(r"$t$")
        ax.set_ylabel(r"Convective wavelength, $2\pi l_c$")
        style_legend(ax, loc="best")
        save_case_figure(fig, case_key, "convective_scale_stitched")

        fig, ax = make_axes()
        if case_key in LINEAR_INTERVALS:
            ax.axvspan(
                *LINEAR_INTERVALS[case_key],
                color=(0.72, 0.72, 0.72),
                alpha=0.28,
                linewidth=0,
                zorder=0,
            )
        ax.plot(
            moist[:, 0],
            TWO_PI * moist[:, 1],
            color=color,
            lw=LINE_WIDTH,
        )
        ax.set_xlim(float(moist[0, 0]), float(moist[-1, 0]))
        ax.set_ylim(bottom=0.0)
        ax.set_xlabel(r"$t$")
        ax.set_ylabel(r"MSE integral scale, $2\pi L_m$")
        save_case_figure(fig, case_key, "moist_integral_scale_stitched")

        fig, ax = make_axes()
        if case_key in LINEAR_INTERVALS:
            ax.axvspan(
                *LINEAR_INTERVALS[case_key],
                color=(0.72, 0.72, 0.72),
                alpha=0.28,
                linewidth=0,
                zorder=0,
            )
        ax.plot(
            mprime[:, 0],
            mprime[:, 4],
            color=color,
            lw=LINE_WIDTH,
        )
        ax.set_xlim(float(mprime[0, 0]), float(mprime[-1, 0]))
        ax.set_ylim(bottom=0.0)
        ax.set_xlabel(r"$t$")
        ax.set_ylabel(
            r"$\left\langle\left\langle m'^2\right\rangle_{xy}\right\rangle_z$"
        )
        save_case_figure(fig, case_key, "mprime_variance_stitched")


def write_reduced_csv(datasets: dict[str, dict[str, np.ndarray]]) -> None:
    output = HERE / "stitched_timeseries.csv"
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "case",
                "time",
                "kinetic_energy",
                "lc_w_z05_raw",
                "lc_w_z_weighted_raw",
                "convective_wavelength_z05_2pi_lc",
                "convective_wavelength_z_weighted_2pi_lc",
                "Lm_raw",
                "MSE_integral_wavelength_2pi_Lm",
                "mprime2_z025",
                "mprime2_z050",
                "mprime2_z075",
                "mprime2_z_weighted",
            ]
        )
        for case in CASES:
            case_data = datasets[case["key"]]
            avg = {round(float(row[0]), 8): row for row in case_data["avgvar"]}
            conv = {
                round(float(row[0]), 8): row for row in case_data["convective"]
            }
            moist = {round(float(row[0]), 8): row for row in case_data["moist"]}
            mprime = {
                round(float(row[0]), 8): row for row in case_data["mprime"]
            }
            times = sorted(set(avg) | set(conv) | set(moist) | set(mprime))
            for time in times:
                avg_row = avg.get(time)
                conv_row = conv.get(time)
                moist_row = moist.get(time)
                mprime_row = mprime.get(time)
                writer.writerow(
                    [
                        case["key"],
                        time,
                        "" if avg_row is None else 0.5 * avg_row[3],
                        "" if conv_row is None else conv_row[1],
                        "" if conv_row is None else conv_row[2],
                        "" if conv_row is None else TWO_PI * conv_row[1],
                        "" if conv_row is None else TWO_PI * conv_row[2],
                        "" if moist_row is None else moist_row[1],
                        "" if moist_row is None else TWO_PI * moist_row[1],
                        "" if mprime_row is None else mprime_row[1],
                        "" if mprime_row is None else mprime_row[2],
                        "" if mprime_row is None else mprime_row[3],
                        "" if mprime_row is None else mprime_row[4],
                    ]
                )


def main() -> None:
    configure_style()
    datasets: dict[str, dict[str, np.ndarray]] = {}
    metadata: dict[str, object] = {
        "definitions": {
            "kinetic_energy": "K = 0.5 * avgvar one-based column 4",
            "convective_scale_file_columns": [
                "time",
                "lc_w_z_approximately_0.5",
                "z_weighted_lc_w",
            ],
            "lc_w_raw": "sum(E_w) / sum(k_h E_w), excluding k_h=0",
            "plotted_convective_wavelength": "2*pi*lc_w",
            "moist_integral_scale_file_columns": [
                "time",
                "Lm",
                "Lq",
                "Lvorticity",
            ],
            "Lm_raw": "sum(E_m) / sum(k_h E_m), excluding k_h=0",
            "plotted_moist_integral_wavelength": "2*pi*Lm",
            "mprime_square_file_columns": [
                "time",
                "mean_xy(mprime^2)_z025",
                "mean_xy(mprime^2)_z050",
                "mean_xy(mprime^2)_z075",
                "stretched_grid_weighted_vertical_mean_xy(mprime^2)",
            ],
            "plotted_mprime_variance": (
                "stretched-grid-weighted vertical mean of the per-height "
                "horizontal variance of mprime, with "
                "mprime=m-mean_xy(m) and m=b+gamma*q"
            ),
            "stitch_rule": (
                "Sort by physical time; at duplicate times retain the later "
                "continuation segment. No time offset is applied."
            ),
        },
        "cases": {},
    }

    for case in CASES:
        case_sets: dict[str, np.ndarray] = {}
        case_meta: dict[str, object] = {}
        for data_type, suffix in FILE_TYPES.items():
            data, info = stitch(case["key"], suffix)
            case_sets[data_type] = data
            case_meta[data_type] = info
        case_meta["has_continuation"] = any(
            len(case_meta[data_type]["segments"]) > 1 for data_type in FILE_TYPES
        )
        datasets[case["key"]] = case_sets
        metadata["cases"][case["key"]] = case_meta

    plot_individual_cases(datasets)
    write_reduced_csv(datasets)
    with (HERE / "stitch_metadata.json").open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2, ensure_ascii=False)


if __name__ == "__main__":
    main()
