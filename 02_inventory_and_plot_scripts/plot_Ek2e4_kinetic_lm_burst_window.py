from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch


RUN = Path(
    r"G:\moist convection\Ek2e-4\AR4\Beta1p02\qbot0p5_qtop0p004978"
    r"\N257x257x65\conti1\conti_strict_force_500\run"
)
OUTPUT = Path(
    r"E:\moist RB\rotating_case_inventory\04_outputs_and_figures"
    r"\high_resolution_timeseries_latest_program_20260805"
    r"\force_num_vs_Ek_20260812"
)

TIME_MIN = 1210.0
TIME_MAX = 1700.0

FRAME_WIDTH = 4.5
LINE_WIDTH = 3.5
FIGSIZE = (7.8, 5.84)
AX_RECT = [0.160, 0.260, 0.660, 0.705]
BOX_ASPECT = 5.2 / 6.5

KINETIC_RED = (0.72, 0.08, 0.32)
LM_BLUE = (0.00, 0.25, 0.90)
BURST_GRAY = (0.55, 0.55, 0.55)


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Times New Roman",
            "mathtext.fontset": "stix",
            "axes.linewidth": FRAME_WIDTH,
            "axes.labelsize": 24,
            "xtick.labelsize": 13,
            "ytick.labelsize": 13,
            "legend.fontsize": 10,
            "figure.dpi": 180,
            "savefig.dpi": 300,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def style_axis(axis: plt.Axes) -> None:
    for spine in axis.spines.values():
        spine.set_linewidth(FRAME_WIDTH)
        spine.set_color("black")
    axis.tick_params(
        which="major",
        direction="in",
        length=12,
        width=1.2,
        top=True,
        right=True,
        colors="black",
    )
    axis.minorticks_on()
    axis.tick_params(
        which="minor",
        direction="in",
        length=6,
        width=1.0,
        top=True,
        right=True,
        colors="black",
    )
    axis.set_box_aspect(BOX_ASPECT)


def load_series() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    avg = np.loadtxt(RUN / "data" / "avgvar.out")
    moist = np.loadtxt(RUN / "diagnostics" / "scale" / "moist_integral_scale.dat")

    kinetic_mask = (avg[:, 0] >= TIME_MIN) & (avg[:, 0] <= TIME_MAX)
    lm_mask = (moist[:, 0] >= TIME_MIN) & (moist[:, 0] <= TIME_MAX)

    kinetic_time = avg[kinetic_mask, 0]
    # Project convention: K = 0.5 * one-based column 4 of avgvar.out.
    kinetic = 0.5 * avg[kinetic_mask, 3]
    lm_time = moist[lm_mask, 0]
    # Project convention: the plotted aggregation scale is 2*pi*L_m.
    two_pi_lm = 2.0 * np.pi * moist[lm_mask, 1]
    return kinetic_time, kinetic, lm_time, two_pi_lm


def two_means(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    centers = np.quantile(values, (0.25, 0.75))
    for _ in range(200):
        labels = np.argmin(np.abs(values[:, None] - centers[None, :]), axis=1)
        updated = np.asarray(
            [np.mean(values[labels == index]) for index in range(2)], dtype=float
        )
        if np.allclose(updated, centers, rtol=1.0e-10, atol=1.0e-12):
            centers = updated
            break
        centers = updated
    return centers, labels


def identify_flat_low_lm_windows(
    time: np.ndarray, values: np.ndarray
) -> tuple[list[tuple[float, float]], np.ndarray]:
    dt = float(np.median(np.diff(time)))
    smoothing_points = max(3, int(round(10.0 / dt)))
    smooth = np.convolve(
        values,
        np.ones(smoothing_points, dtype=float) / smoothing_points,
        mode="same",
    )
    indices = np.arange(len(time))
    valid = (indices >= smoothing_points // 2) & (
        indices < len(time) - smoothing_points // 2
    )
    valid_time = time[valid]
    valid_smooth = smooth[valid]

    centers, labels = two_means(valid_smooth)
    low_mask = labels == int(np.argmin(centers))
    slope = np.gradient(valid_smooth, valid_time)
    flat_low = low_mask & (np.abs(slope) <= 3.0e-3)

    edges = np.flatnonzero(np.diff(np.r_[False, flat_low, False])).reshape(-1, 2)
    windows: list[tuple[float, float]] = []
    for start, stop in edges:
        left = float(valid_time[start])
        right = float(valid_time[stop - 1])
        if right - left >= 20.0:
            windows.append((left, right))
    return windows, smooth


def main() -> None:
    configure_style()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    kinetic_time, kinetic, lm_time, two_pi_lm = load_series()
    windows, smoothed_lm = identify_flat_low_lm_windows(lm_time, two_pi_lm)

    figure = plt.figure(figsize=FIGSIZE, facecolor="white")
    left = figure.add_axes(AX_RECT)
    right = left.twinx()
    style_axis(left)
    style_axis(right)
    left.yaxis.set_offset_position("left")
    right.yaxis.set_offset_position("right")

    for start, stop in windows:
        left.axvspan(
            start,
            stop,
            facecolor=BURST_GRAY,
            alpha=0.18,
            edgecolor="none",
            zorder=0,
        )

    kinetic_line = left.plot(
        kinetic_time,
        kinetic,
        color=KINETIC_RED,
        lw=LINE_WIDTH,
        label=r"$K$",
        zorder=3,
    )[0]
    lm_line = right.plot(
        lm_time,
        two_pi_lm,
        color=LM_BLUE,
        lw=LINE_WIDTH,
        label=r"$2\pi L_m$",
        zorder=4,
    )[0]

    left.set_xlim(TIME_MIN, TIME_MAX)
    left.set_ylim(bottom=0.0)
    right.set_ylim(bottom=0.0)
    left.set_xlabel(r"$t$")
    left.set_ylabel(r"$K$")
    right.set_ylabel(r"$2\pi L_m$")
    left.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2), useMathText=True)
    left.yaxis.get_offset_text().set_fontsize(13)

    stem = OUTPUT / "Ra8e6_Ek2e-4_kinetic_lm_dual_axis_burst_windows"
    figure.savefig(stem.with_suffix(".png"), facecolor="white")
    figure.savefig(stem.with_suffix(".pdf"), facecolor="white")
    plt.close(figure)

    legend_figure = plt.figure(figsize=(6.5, 1.0), facecolor="white")
    legend_axis = legend_figure.add_axes([0.0, 0.0, 1.0, 1.0])
    legend_axis.axis("off")
    burst_patch = Patch(
        facecolor=BURST_GRAY,
        alpha=0.18,
        edgecolor="none",
        label="plume burst",
    )
    legend_axis.legend(
        [kinetic_line, lm_line, burst_patch],
        [r"$K$", r"$2\pi L_m$", "plume burst"],
        frameon=False,
        loc="center",
        ncol=3,
        handlelength=2.6,
        columnspacing=1.8,
        fontsize=16,
    )
    legend_figure.savefig(
        stem.with_name(stem.name + "_legend").with_suffix(".png"),
        facecolor="white",
    )
    legend_figure.savefig(
        stem.with_name(stem.name + "_legend").with_suffix(".pdf"),
        facecolor="white",
    )
    plt.close(legend_figure)

    merged = pd.DataFrame(
        {
            "time": kinetic_time,
            "kinetic_energy": kinetic,
        }
    ).merge(
        pd.DataFrame({"time": lm_time, "two_pi_Lm": two_pi_lm}),
        on="time",
        how="outer",
    )
    merged.to_csv(stem.with_name(stem.name + "_plot_data.csv"), index=False)
    pd.DataFrame(windows, columns=["burst_start", "burst_end"]).to_csv(
        stem.with_name(stem.name + "_burst_windows.csv"), index=False
    )

    # Save the smoothed scale used only for auditing the plateau selection.
    pd.DataFrame({"time": lm_time, "two_pi_Lm_smoothed": smoothed_lm}).to_csv(
        stem.with_name(stem.name + "_lm_plateau_detection.csv"), index=False
    )


if __name__ == "__main__":
    main()
