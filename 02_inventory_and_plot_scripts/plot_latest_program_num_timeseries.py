from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D


PROJECT = Path(r"E:\moist RB\rotating_case_inventory")
SOURCE = PROJECT / "04_outputs_and_figures" / "high_resolution_timeseries_latest_program_20260805"
INPUT = SOURCE / "high_resolution_timeseries_long.csv"
OUTPUT = SOURCE / "num_timeseries"

FRAME_WIDTH = 4.5
LINE_WIDTH = 3.5
AXIS_LABEL_SIZE = 24
TICK_LABEL_SIZE = 13
LEGEND_SIZE = 12
BOX_ASPECT = 5.2 / 6.5
BLOCK_WIDTH = 5.0


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Times New Roman",
            "mathtext.fontset": "stix",
            "axes.linewidth": FRAME_WIDTH,
            "axes.labelsize": AXIS_LABEL_SIZE,
            "xtick.labelsize": TICK_LABEL_SIZE,
            "ytick.labelsize": TICK_LABEL_SIZE,
            "legend.fontsize": LEGEND_SIZE,
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


def make_axis() -> tuple[plt.Figure, plt.Axes]:
    figure = plt.figure(figsize=(6.5, 5.84), facecolor="white")
    axis = figure.add_axes([0.186, 0.260, 0.792, 0.705])
    style_axis(axis)
    axis.set_xlabel(r"$t$")
    axis.set_ylabel(r"$Nu_m(t)$")
    return figure, axis


def ek_label(ek: float) -> str:
    if not math.isfinite(ek):
        return "Nonrotating"
    scientific = f"{ek:.8e}"
    coefficient_text, exponent_text = scientific.split("e")
    coefficient = float(coefficient_text)
    exponent = int(exponent_text)
    if abs(coefficient - 1.0) < 1.0e-8:
        return rf"$10^{{{exponent}}}$"
    return rf"${coefficient:g}\times10^{{{exponent}}}$"


def split_at_gaps(frame: pd.DataFrame) -> list[pd.DataFrame]:
    frame = frame.sort_values("time").drop_duplicates("time", keep="last")
    time = frame["time"].to_numpy(float)
    if len(time) < 2:
        return [frame]
    positive = np.diff(time)
    positive = positive[positive > 0.0]
    nominal = float(np.median(positive)) if len(positive) else 0.1
    breaks = np.flatnonzero(np.diff(time) > max(5.0 * nominal, 2.0 * BLOCK_WIDTH)) + 1
    return [piece for piece in np.split(frame, breaks) if len(piece)]


def block_average(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.copy()
    origin = float(frame["time"].iloc[0])
    frame["block"] = np.floor((frame["time"] - origin) / BLOCK_WIDTH).astype(int)
    return frame.groupby("block", as_index=False).agg(time=("time", "mean"), value=("value", "mean"))


def save(figure: plt.Figure, stem: Path) -> None:
    figure.savefig(stem.with_suffix(".png"), facecolor="white")
    figure.savefig(stem.with_suffix(".pdf"), facecolor="white")
    plt.close(figure)


def case_key(frame: pd.DataFrame) -> tuple:
    row = frame.iloc[0]
    ek = float(row["Ek"]) if pd.notna(row["Ek"]) else math.nan
    return (
        1 if not math.isfinite(ek) else 0,
        ek if math.isfinite(ek) else math.inf,
        int(row["Nx"]),
    )


def main() -> None:
    configure_style()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(INPUT)
    data = data[data["metric"].eq("num")].copy()
    group_columns = ["Ra", "Ek", "AR", "beta", "qbot", "Nx", "Ny", "Nz"]
    grouped = sorted(
        (group for _, group in data.groupby(group_columns, dropna=False)),
        key=case_key,
    )

    finite_count = sum(pd.notna(group.iloc[0]["Ek"]) for group in grouped)
    colors = iter(plt.cm.turbo(np.linspace(0.04, 0.94, max(finite_count, 1))))
    case_styles = []
    for group in grouped:
        row = group.iloc[0]
        ek = float(row["Ek"]) if pd.notna(row["Ek"]) else math.nan
        if math.isfinite(ek):
            color, linestyle = tuple(next(colors)), "-"
        else:
            color, linestyle = (0.35, 0.35, 0.35, 1.0), "--"
        case_styles.append((group, ek, color, linestyle))

    figure, axis = make_axis()
    handles = []
    for group, ek, color, linestyle in case_styles:
        for segment in split_at_gaps(group):
            reduced = block_average(segment)
            axis.plot(
                reduced["time"],
                reduced["value"],
                color=color,
                linestyle=linestyle,
                linewidth=LINE_WIDTH,
            )
        handles.append(
            Line2D([], [], color=color, linestyle=linestyle, linewidth=LINE_WIDTH, label=ek_label(ek))
        )
    axis.set_xlim(left=0.0)
    all_stem = OUTPUT / "Ra8e6_latest_program_all_cases_Num_timeseries"
    save(figure, all_stem)

    legend_figure = plt.figure(figsize=(6.5, 2.5), facecolor="white")
    legend_axis = legend_figure.add_axes([0.02, 0.02, 0.96, 0.96])
    legend_axis.axis("off")
    legend_axis.legend(
        handles=handles,
        title=r"$Ek$",
        frameon=False,
        ncol=4,
        loc="center",
        fontsize=LEGEND_SIZE,
        title_fontsize=LEGEND_SIZE + 1,
    )
    save(legend_figure, OUTPUT / "Ra8e6_latest_program_all_cases_Num_legend")

    index_rows = []
    for group, ek, color, linestyle in case_styles:
        row = group.iloc[0]
        figure, axis = make_axis()
        for segment in split_at_gaps(group):
            reduced = block_average(segment)
            axis.plot(
                reduced["time"], reduced["value"], color=color,
                linestyle=linestyle, linewidth=LINE_WIDTH,
            )
        axis.set_xlim(left=0.0)
        token = "NR" if not math.isfinite(ek) else f"{ek:.8g}".replace(".", "p")
        stem = OUTPUT / "individual_cases" / f"Ra8e6_Ek{token}_Num_timeseries"
        stem.parent.mkdir(parents=True, exist_ok=True)
        save(figure, stem)
        index_rows.append(
            {
                "Ek": "" if not math.isfinite(ek) else ek,
                "AR": row["AR"],
                "grid": f"{int(row['Nx'])}x{int(row['Ny'])}x{int(row['Nz'])}",
                "time_min": group["time"].min(),
                "time_max": group["time"].max(),
                "rows": len(group),
                "png": str(stem.with_suffix(".png")),
            }
        )

    pd.DataFrame(index_rows).to_csv(
        OUTPUT / "Ra8e6_latest_program_Num_case_index.csv", index=False,
        encoding="utf-8-sig",
    )
    data.to_csv(
        OUTPUT / "Ra8e6_latest_program_Num_timeseries_normalized.csv",
        index=False, encoding="utf-8-sig",
    )
    print(all_stem.with_suffix(".png"))
    print(f"cases={len(grouped)}")


if __name__ == "__main__":
    main()
