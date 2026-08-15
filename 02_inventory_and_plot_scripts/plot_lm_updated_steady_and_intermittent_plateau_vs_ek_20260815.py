from __future__ import annotations

from pathlib import Path

import json

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FuncFormatter, LogFormatterMathtext, LogLocator
import numpy as np
import pandas as pd


BASE = Path(r"E:\moist RB\rotating_case_inventory\04_outputs_and_figures\high_resolution_timeseries_latest_program_20260805")
INTERMITTENT_SOURCE = BASE / "lm_low_valley_scaling" / "Ra8e6_lm_low_valley_and_stable_time_averages.csv"
STEADY_SOURCE = BASE / "stable_lm_latest_z07_slices_20260815" / "lc_lm_skill_redraw_large_Ek_20260815" / "Ra8e6_Lm_late500_summary.csv"
OUT = BASE / "stable_lm_latest_z07_slices_20260815" / "Lm_updated_steady_and_intermittent_plateau_vs_Ek_20260815"

INTERMITTENT_COLOR = (1.00, 0.72, 0.00)
STEADY_COLOR = (0.00, 0.32, 0.90)
FRAME = 4.5
LINE = 3.5
MARKER_SIZE = 12
EDGE = 1.8


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Times New Roman",
            "mathtext.fontset": "stix",
            "axes.linewidth": FRAME,
            "axes.labelsize": 24,
            "xtick.labelsize": 13,
            "ytick.labelsize": 13,
            "legend.fontsize": 15,
            "figure.dpi": 180,
            "savefig.dpi": 300,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def make_axis() -> tuple[plt.Figure, plt.Axes]:
    fig = plt.figure(figsize=(6.5, 5.84), facecolor="white")
    ax = fig.add_axes([0.186, 0.260, 0.792, 0.705])
    ax.set_box_aspect(5.2 / 6.5)
    for spine in ax.spines.values():
        spine.set_linewidth(FRAME)
    ax.tick_params(which="major", direction="in", length=12, width=1.2, top=True, right=True)
    ax.tick_params(which="minor", direction="in", length=6, width=1.0, top=True, right=True)
    ax.minorticks_on()
    return fig, ax


def fit_power(x: np.ndarray, y: np.ndarray) -> dict[str, float]:
    power, log_prefactor = np.polyfit(np.log(x), np.log(y), 1)
    prefactor = float(np.exp(log_prefactor))
    prediction = prefactor * x**power
    residual = float(np.sum((np.log(y) - np.log(prediction)) ** 2))
    total = float(np.sum((np.log(y) - np.mean(np.log(y))) ** 2))
    return {
        "power": float(power),
        "prefactor": prefactor,
        "R2_log_space": float(1.0 - residual / total),
        "point_count": int(len(x)),
    }


def load_data() -> pd.DataFrame:
    intermittent = pd.read_csv(INTERMITTENT_SOURCE)
    intermittent = intermittent[intermittent["selection"] == "intermittent low valley"].copy()
    intermittent = intermittent.rename(
        columns={
            "two_pi_Lm_mean": "mean",
            "two_pi_Lm_std": "std",
            "time_window": "averaging_window",
            "sample_count": "n_samples",
        }
    )
    intermittent["case"] = intermittent["Ek"].map(lambda value: f"Ek{value:g}")
    intermittent["group"] = "intermittent plateau"

    steady = pd.read_csv(STEADY_SOURCE).copy()
    steady["averaging_window"] = steady.apply(
        lambda row: f"{row['time_start']:.1f}-{row['time_end']:.1f}", axis=1
    )
    steady["group"] = "late-time steady"

    columns = ["case", "Ek", "group", "averaging_window", "n_samples", "mean", "std"]
    return pd.concat([intermittent[columns], steady[columns]], ignore_index=True).sort_values("Ek")


def main() -> None:
    configure_style()
    OUT.mkdir(parents=True, exist_ok=True)
    data = load_data()
    intermittent = data[data["group"] == "intermittent plateau"]
    steady = data[data["group"] == "late-time steady"]

    fits = {
        "combined": fit_power(data["Ek"].to_numpy(float), data["mean"].to_numpy(float)),
        "intermittent_plateau": fit_power(intermittent["Ek"].to_numpy(float), intermittent["mean"].to_numpy(float)),
        "late_time_steady": fit_power(steady["Ek"].to_numpy(float), steady["mean"].to_numpy(float)),
    }
    data.to_csv(OUT / "Ra8e6_Lm_updated_steady_and_intermittent_plateau_means.csv", index=False, encoding="utf-8-sig")
    (OUT / "Ra8e6_Lm_updated_steady_and_intermittent_powerlaw_fits.json").write_text(
        json.dumps(fits, indent=2), encoding="utf-8"
    )

    fig, ax = make_axis()
    ax.set_xscale("log")
    ax.set_yscale("log")

    fit = fits["combined"]
    xline = np.geomspace(1.75e-4, 1.15e-2, 300)
    ax.plot(
        xline,
        fit["prefactor"] * xline ** fit["power"],
        color="black",
        linestyle="--",
        linewidth=LINE,
        zorder=1,
    )
    for frame, color in ((intermittent, INTERMITTENT_COLOR), (steady, STEADY_COLOR)):
        ax.plot(
            frame["Ek"],
            frame["mean"],
            linestyle="none",
            marker="o",
            markersize=MARKER_SIZE,
            markerfacecolor=color,
            markeredgecolor="black",
            markeredgewidth=EDGE,
            zorder=5,
        )

    ax.set_xlim(1.45e-4, 1.25e-2)
    ax.set_ylim(0.27, 4.15)
    ax.set_xlabel(r"$Ek$")
    ax.set_ylabel(r"$\left\langle2\pi L_m\right\rangle_t/H$")
    ax.xaxis.set_major_locator(LogLocator(base=10.0, numticks=4))
    ax.xaxis.set_major_formatter(LogFormatterMathtext(base=10.0))
    ax.xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1))
    ax.yaxis.set_major_locator(FixedLocator([0.3, 0.5, 1.0, 2.0, 3.0, 4.0]))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}"))
    ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1))

    ax.text(
        0.965,
        0.075,
        rf"$2\pi L_m\propto Ek^{{{fit['power']:.3f}}}$" + "\n" + rf"$R^2={fit['R2_log_space']:.3f}$",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=15,
    )
    ax.legend(
        handles=[
            Line2D([], [], linestyle="none", marker="o", markersize=9, markerfacecolor=INTERMITTENT_COLOR, markeredgecolor="black", markeredgewidth=1.5, label="Intermittent plateau"),
            Line2D([], [], linestyle="none", marker="o", markersize=9, markerfacecolor=STEADY_COLOR, markeredgecolor="black", markeredgewidth=1.5, label="Late-time steady"),
        ],
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(0.035, 0.955),
        borderaxespad=0.0,
        fontsize=15,
    )

    stem = OUT / "Ra8e6_Lm_updated_steady_and_intermittent_plateau_vs_Ek_combined_fit"
    fig.savefig(stem.with_suffix(".png"), facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
    plt.close(fig)

    print(stem.with_suffix(".png"))
    print(json.dumps(fits, indent=2))
    print(data.to_string(index=False))


if __name__ == "__main__":
    main()
