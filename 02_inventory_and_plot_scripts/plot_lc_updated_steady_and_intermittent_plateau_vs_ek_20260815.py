from __future__ import annotations

from pathlib import Path
import json

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FuncFormatter, LogFormatterMathtext, LogLocator, NullFormatter
import numpy as np
import pandas as pd


BASE = Path(r"E:\moist RB\rotating_case_inventory\04_outputs_and_figures\high_resolution_timeseries_latest_program_20260805")
SERIES = BASE / "high_resolution_timeseries_long.csv"
INTERMITTENT_SOURCE = BASE / "lm_low_valley_scaling" / "Ra8e6_lm_low_valley_and_stable_time_averages.csv"
STEADY_SOURCE = BASE / "stable_lm_latest_z07_slices_20260815" / "lc_lm_skill_redraw_large_Ek_20260815" / "Ra8e6_lc_late500_summary.csv"
OUT = BASE / "stable_lm_latest_z07_slices_20260815" / "lc_updated_steady_and_intermittent_plateau_vs_Ek_20260815"

INTERMITTENT_COLOR = (1.00, 0.72, 0.00)
STEADY_COLOR = (0.00, 0.32, 0.90)
FRAME = 4.5
LINE = 3.5


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


def fit_power(frame: pd.DataFrame) -> dict[str, float]:
    x = frame["Ek"].to_numpy(float)
    y = frame["mean"].to_numpy(float)
    power, log_prefactor = np.polyfit(np.log(x), np.log(y), 1)
    prefactor = float(np.exp(log_prefactor))
    prediction = prefactor * x**power
    residual = float(np.sum((np.log(y) - np.log(prediction)) ** 2))
    total = float(np.sum((np.log(y) - np.mean(np.log(y))) ** 2))
    return {
        "power": float(power),
        "prefactor": prefactor,
        "R2_log_space": float(1.0 - residual / total),
        "point_count": int(len(frame)),
    }


def parse_windows(text: str) -> list[tuple[float, float]]:
    windows = []
    for item in text.split(";"):
        start, end = item.strip().split("-")
        windows.append((float(start), float(end)))
    return windows


def load_intermittent() -> pd.DataFrame:
    selection = pd.read_csv(INTERMITTENT_SOURCE)
    selection = selection[selection["selection"] == "intermittent low valley"].copy()
    targets = selection["Ek"].to_numpy(float)

    pieces = []
    for chunk in pd.read_csv(SERIES, chunksize=300_000):
        ek = chunk["Ek"].to_numpy(float)
        target_mask = np.isclose(ek[:, None], targets[None, :], rtol=0.0, atol=1.0e-10).any(axis=1)
        mask = (chunk["Ra"] == 8.0e6) & (chunk["metric"] == "convective_mid") & chunk["Ek"].notna() & target_mask
        if mask.any():
            pieces.append(chunk.loc[mask, ["Ek", "time", "value"]])
    history = pd.concat(pieces, ignore_index=True).drop_duplicates(["Ek", "time"], keep="last")

    rows = []
    for item in selection.itertuples(index=False):
        case = history[np.isclose(history["Ek"], item.Ek, rtol=0.0, atol=1.0e-10)]
        windows = parse_windows(item.time_window)
        selected = pd.concat([case[(case["time"] >= start) & (case["time"] <= end)] for start, end in windows])
        selected = selected.drop_duplicates("time", keep="last")
        rows.append(
            {
                "case": f"Ek{item.Ek:g}",
                "Ek": float(item.Ek),
                "group": "intermittent plateau",
                "averaging_window": item.time_window,
                "n_samples": int(len(selected)),
                "mean": float(selected["value"].mean()),
                "std": float(selected["value"].std(ddof=1)),
            }
        )
    return pd.DataFrame(rows)


def load_steady() -> pd.DataFrame:
    steady = pd.read_csv(STEADY_SOURCE)
    steady["group"] = "late-time steady"
    steady["averaging_window"] = steady.apply(lambda row: f"{row.time_start:.1f}-{row.time_end:.1f}", axis=1)
    return steady[["case", "Ek", "group", "averaging_window", "n_samples", "mean", "std"]]


def main() -> None:
    configure_style()
    OUT.mkdir(parents=True, exist_ok=True)
    intermittent = load_intermittent()
    steady = load_steady()
    data = pd.concat([intermittent, steady], ignore_index=True).sort_values("Ek")

    fits = {
        "combined": fit_power(data),
        "intermittent_plateau": fit_power(intermittent),
        "late_time_steady": fit_power(steady),
    }
    data.to_csv(OUT / "Ra8e6_lc_updated_steady_and_intermittent_plateau_means.csv", index=False, encoding="utf-8-sig")
    (OUT / "Ra8e6_lc_updated_steady_and_intermittent_powerlaw_fits.json").write_text(json.dumps(fits, indent=2), encoding="utf-8")

    fig, ax = make_axis()
    ax.set_xscale("log")
    ax.set_yscale("log")
    fit = fits["combined"]
    xline = np.geomspace(1.75e-4, 1.15e-2, 300)
    ax.plot(xline, fit["prefactor"] * xline ** fit["power"], color="black", ls="--", lw=LINE, zorder=1)
    for frame, color in ((intermittent, INTERMITTENT_COLOR), (steady, STEADY_COLOR)):
        ax.plot(
            frame["Ek"], frame["mean"], linestyle="none", marker="o", markersize=12,
            markerfacecolor=color, markeredgecolor="black", markeredgewidth=1.8, zorder=5,
        )

    ax.set_xlim(1.45e-4, 1.25e-2)
    ax.set_ylim(0.245, 1.02)
    ax.set_xlabel(r"$Ek$")
    ax.set_ylabel(r"$\left\langle2\pi l_c(z\simeq0.5)\right\rangle_t/H$")
    ax.xaxis.set_major_locator(LogLocator(base=10.0, numticks=4))
    ax.xaxis.set_major_formatter(LogFormatterMathtext(base=10.0))
    ax.xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1))
    ax.yaxis.set_major_locator(FixedLocator([0.25, 0.3, 0.4, 0.5, 0.7, 1.0]))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}"))
    ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1))
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.text(
        0.965, 0.075,
        rf"$2\pi l_c\propto Ek^{{{fit['power']:.3f}}}$" + "\n" + rf"$R^2={fit['R2_log_space']:.3f}$",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=15,
    )
    ax.legend(
        handles=[
            Line2D([], [], linestyle="none", marker="o", markersize=9, markerfacecolor=INTERMITTENT_COLOR, markeredgecolor="black", markeredgewidth=1.5, label="Intermittent plateau"),
            Line2D([], [], linestyle="none", marker="o", markersize=9, markerfacecolor=STEADY_COLOR, markeredgecolor="black", markeredgewidth=1.5, label="Late-time steady"),
        ],
        frameon=False, loc="upper left", bbox_to_anchor=(0.035, 0.955), borderaxespad=0.0, fontsize=15,
    )

    stem = OUT / "Ra8e6_lc_updated_steady_and_intermittent_plateau_vs_Ek_combined_fit"
    fig.savefig(stem.with_suffix(".png"), facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
    plt.close(fig)

    print(stem.with_suffix(".png"))
    print(json.dumps(fits, indent=2))
    print(data.to_string(index=False))


if __name__ == "__main__":
    main()
