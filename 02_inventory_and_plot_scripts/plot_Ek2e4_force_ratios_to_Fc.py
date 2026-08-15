from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


SOURCE = Path(
    r"G:\moist convection\Ek2e-4\AR4\Beta1p02\qbot0p5_qtop0p004978"
    r"\N257x257x65\conti1\conti_strict_force_500\run\diagnostics"
    r"\force_balance_strict_fp3d_20260806"
    r"\strict_force_balance_bulk_timeseries.csv"
)
OUTPUT = Path(
    r"E:\moist RB\rotating_case_inventory\04_outputs_and_figures"
    r"\high_resolution_timeseries_latest_program_20260805"
    r"\force_num_vs_Ek_20260812"
)

TIME_MIN = 1210.0
TIME_MAX = 1700.0
BURST_WINDOWS = [(1219.4, 1321.8), (1523.1, 1597.9)]

FRAME_WIDTH = 4.5
LINE_WIDTH = 3.5
FIGSIZE = (6.5, 5.84)
AX_RECT = [0.186, 0.260, 0.792, 0.705]
BOX_ASPECT = 5.2 / 6.5
BURST_GRAY = (0.55, 0.55, 0.55)

COLORS = {
    "F_I": (0.74, 0.14, 0.18),
    "F_C": (0.11, 0.44, 0.71),
    "F_P": (0.58, 0.05, 0.90),
    "F_V": (0.96, 0.58, 0.19),
    "F_B": (0.88, 0.05, 0.65),
}
LABELS = {
    "F_I": r"$F_I/F_C$",
    "F_C": r"$F_C/F_C$",
    "F_P": r"$F_P/F_C$",
    "F_V": r"$F_V/F_C$",
    "F_B": r"$F_B/F_C$",
}
GEOSTROPHIC_COLOR = (0.00, 0.65, 0.12)
GEOSTROPHIC_LABEL = r"$|\mathbf{F}_{C,h}+\mathbf{F}_{P,h}|/F_C$"


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
    axis.tick_params(
        which="major", direction="in", length=12, width=1.2, top=True, right=True
    )
    axis.minorticks_on()
    axis.tick_params(
        which="minor", direction="in", length=6, width=1.0, top=True, right=True
    )
    axis.set_box_aspect(BOX_ASPECT)


def main() -> None:
    configure_style()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(SOURCE).sort_values("time").drop_duplicates("time", keep="last")
    data = data[(data["time"] >= TIME_MIN) & (data["time"] <= TIME_MAX)].copy()
    denominator = pd.to_numeric(data["F_C"], errors="coerce")

    ratio_data = pd.DataFrame({"time": data["time"]})
    for field in COLORS:
        numerator = pd.to_numeric(data[field], errors="coerce")
        ratio_data[f"{field}_over_F_C"] = numerator / denominator
    ratio_data["R_G_over_F_C"] = (
        pd.to_numeric(data["R_G"], errors="coerce") / denominator
    )

    figure = plt.figure(figsize=FIGSIZE, facecolor="white")
    axis = figure.add_axes(AX_RECT)
    style_axis(axis)
    for start, stop in BURST_WINDOWS:
        axis.axvspan(
            start,
            stop,
            facecolor=BURST_GRAY,
            alpha=0.18,
            edgecolor="none",
            zorder=0,
        )

    handles: list[Line2D] = []
    for field, color in COLORS.items():
        values = ratio_data[f"{field}_over_F_C"].to_numpy(float)
        time = ratio_data["time"].to_numpy(float)
        valid = np.isfinite(time) & np.isfinite(values) & (values > 0.0)
        linestyle = "--" if field == "F_C" else "-"
        zorder = 2 if field == "F_C" else 3
        line = axis.semilogy(
            time[valid],
            values[valid],
            color=color,
            lw=LINE_WIDTH,
            ls=linestyle,
            zorder=zorder,
        )[0]
        handles.append(line)

    geostrophic_ratio = ratio_data["R_G_over_F_C"].to_numpy(float)
    time = ratio_data["time"].to_numpy(float)
    valid = np.isfinite(time) & np.isfinite(geostrophic_ratio) & (geostrophic_ratio > 0.0)
    geostrophic_line = axis.semilogy(
        time[valid],
        geostrophic_ratio[valid],
        color=GEOSTROPHIC_COLOR,
        lw=LINE_WIDTH,
        ls="-.",
        zorder=4,
    )[0]
    handles.append(geostrophic_line)

    axis.set_xlim(TIME_MIN, TIME_MAX)
    axis.set_xlabel(r"$t$")
    axis.set_ylabel(r"$F_\alpha/F_C$")

    stem = OUTPUT / "Ra8e6_Ek2e-4_force_ratios_to_Fc_timeseries"
    figure.savefig(stem.with_suffix(".png"), facecolor="white")
    figure.savefig(stem.with_suffix(".pdf"), facecolor="white")
    plt.close(figure)

    legend_figure = plt.figure(figsize=(7.6, 1.35), facecolor="white")
    legend_axis = legend_figure.add_axes([0.0, 0.0, 1.0, 1.0])
    legend_axis.axis("off")
    burst_patch = Patch(
        facecolor=BURST_GRAY,
        alpha=0.18,
        edgecolor="none",
        label="plume burst",
    )
    legend_axis.legend(
        [*handles, burst_patch],
        [*[LABELS[field] for field in COLORS], GEOSTROPHIC_LABEL, "plume burst"],
        frameon=False,
        loc="center",
        ncol=4,
        handlelength=2.5,
        columnspacing=1.4,
        fontsize=15,
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

    ratio_data.to_csv(stem.with_name(stem.name + "_plot_data.csv"), index=False)


if __name__ == "__main__":
    main()
