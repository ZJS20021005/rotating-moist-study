#!/usr/bin/env python3
"""Plot current remote-reduced Re and Nu_m diagnostics."""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D


ROOT = Path(r"E:\moist RB\rotating_case_inventory")
DATA_DIR = ROOT / "04_outputs_and_figures" / "re_num_20260805"
SERIES = DATA_DIR / "re_num_timeseries.csv"
SUMMARY = DATA_DIR / "re_num_stability_summary.csv"

FRAME = 4.5
LABEL = 24
TICK = 13
LINE = 3.5
BOX = 5.2 / 6.5
AX = [0.186, 0.260, 0.792, 0.705]


def style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Times New Roman",
            "mathtext.fontset": "stix",
            "axes.linewidth": FRAME,
            "axes.labelsize": LABEL,
            "xtick.labelsize": TICK,
            "ytick.labelsize": TICK,
            "legend.fontsize": 10,
            "figure.dpi": 180,
            "savefig.dpi": 300,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def axis(fig: plt.Figure) -> plt.Axes:
    ax = fig.add_axes(AX)
    for spine in ax.spines.values():
        spine.set_linewidth(FRAME)
    ax.tick_params(which="major", direction="in", length=12, width=1.2, top=True, right=True)
    ax.minorticks_on()
    ax.tick_params(which="minor", direction="in", length=6, width=1.0, top=True, right=True)
    ax.set_box_aspect(BOX)
    return ax


def ek_label(value: float) -> str:
    return f"$Ek={value:.2g}$"


def save(fig: plt.Figure, stem: Path) -> None:
    fig.savefig(stem.with_suffix(".png"), facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
    plt.close(fig)


def main() -> None:
    style()
    data = pd.read_csv(SERIES)
    summary = pd.read_csv(SUMMARY)
    data["Ek"] = pd.to_numeric(data["Ek"], errors="coerce")
    summary["Ek"] = pd.to_numeric(summary["Ek"], errors="coerce")
    finite = summary[np.isfinite(summary["Ek"])].sort_values("Ek")
    baseline = summary[~np.isfinite(summary["Ek"])].copy()
    colors = plt.cm.turbo(np.linspace(0.06, 0.94, len(finite)))
    cmap = {float(row.Ek): color for row, color in zip(finite.itertuples(), colors)}

    for metric, ylabel, stem_name in (
        ("Re", r"$Re(t)$", "Ra8e6_Re_timeseries_all_cases"),
        ("Num", r"$Nu_m(t)$", "Ra8e6_Num_timeseries_all_cases"),
    ):
        fig = plt.figure(figsize=(6.5, 5.84), facecolor="white")
        ax = axis(fig)
        for ek, group in data[data["metric"].eq(metric)].groupby("Ek"):
            group = group.sort_values("time")
            color = cmap.get(float(ek), (0.2, 0.2, 0.2, 1.0))
            ax.plot(group["time"], group["value"], color=color, lw=LINE)
        base_group = data[data["metric"].eq(metric) & ~np.isfinite(data["Ek"])]
        if not base_group.empty:
            base_group = base_group.sort_values("time")
            ax.plot(
                base_group["time"],
                base_group["value"],
                color="0.35",
                lw=LINE,
                ls="--",
            )
        ax.set_xlabel(r"$t$")
        ax.set_ylabel(ylabel)
        ax.set_xlim(left=0.0)
        ax.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2), useMathText=True)
        ax.yaxis.get_offset_text().set_fontsize(TICK)
        save(fig, DATA_DIR / stem_name)

    legend_fig = plt.figure(figsize=(6.5, 2.5), facecolor="white")
    legend_ax = legend_fig.add_axes([0.02, 0.02, 0.96, 0.96])
    legend_ax.axis("off")
    handles = [
        Line2D([], [], color=cmap[float(row.Ek)], lw=LINE, label=ek_label(float(row.Ek)))
        for row in finite.itertuples()
    ]
    legend_ax.legend(
        handles=handles,
        title=r"$Ek$",
        frameon=False,
        ncol=4,
        loc="center",
        fontsize=10,
        title_fontsize=11,
    )
    save(legend_fig, DATA_DIR / "Ra8e6_Re_Num_timeseries_legend")

    stable = finite[finite["stable_for_Ek_scaling"].astype(str).str.lower().eq("true")].copy()
    stable = stable.sort_values("Ek")
    for metric, ylabel, stem_name, color in (
        ("Re", r"$Re$", "Ra8e6_Re_vs_Ek_stable", (0.00, 0.18, 0.75)),
        ("Num", r"$Nu_m$", "Ra8e6_Num_vs_Ek_stable", (0.82, 0.18, 0.05)),
    ):
        column = "Re_mean_200" if metric == "Re" else "Num_mean_200"
        fig = plt.figure(figsize=(6.5, 5.84), facecolor="white")
        ax = axis(fig)
        ax.loglog(
            stable["Ek"].to_numpy(float),
            stable[column].to_numpy(float),
            marker="o",
            ms=9,
            mfc=color,
            mec="black",
            mew=1.5,
            linestyle="None",
        )
        base_value = baseline.loc[
            baseline["stable_for_Ek_scaling"].astype(str).str.lower().eq("true"),
            "Re_mean_200" if metric == "Re" else "Num_mean_200",
        ]
        if not base_value.empty:
            ax.axhline(float(base_value.iloc[0]), color="0.35", ls="--", lw=2.5)
        ax.set_xlabel(r"$Ek$")
        ax.set_ylabel(ylabel)
        ax.grid(False)
        save(fig, DATA_DIR / stem_name)

    stable.to_csv(DATA_DIR / "stable_cases_used_for_ek_scaling.csv", index=False, encoding="utf-8-sig")
    print(f"saved figures and tables under {DATA_DIR}")
    print("stable Ek:", ", ".join(f"{x:.6g}" for x in stable["Ek"]))


if __name__ == "__main__":
    main()
