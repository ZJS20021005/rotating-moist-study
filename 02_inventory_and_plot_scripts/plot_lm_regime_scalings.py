from __future__ import annotations

import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FuncFormatter, LogFormatterMathtext, LogLocator


PROJECT_ROOT = Path(r"E:\moist RB\rotating_case_inventory")
SERIES_PARENT = PROJECT_ROOT / "04_outputs_and_figures"
SERIES_DIR = max(
    (
        path
        for path in SERIES_PARENT.glob("high_resolution_timeseries_latest_program_*")
        if (path / "lm_low_valley_scaling" / "Ra8e6_lm_low_valley_and_stable_time_averages.csv").exists()
    ),
    key=lambda path: path.name,
)
INPUT = (
    SERIES_DIR
    / "lm_low_valley_scaling"
    / "Ra8e6_lm_low_valley_and_stable_time_averages.csv"
)
OUTPUT_DIR = SERIES_DIR / "lm_regime_scalings"

INTERMITTENT_COLOR = (1.00, 0.72, 0.00)
FUNNEL_COLOR = (0.00, 0.32, 0.90)
FRAME_WIDTH = 4.5
LINE_WIDTH = 3.5
AXIS_LABEL_SIZE = 24
TICK_LABEL_SIZE = 13
BOX_ASPECT = 5.2 / 6.5
MARKER_SIZE = 12
MARKER_EDGE_WIDTH = 1.8
FUNNEL_FIXED_POWER = 0.69


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Times New Roman",
            "mathtext.fontset": "stix",
            "axes.linewidth": FRAME_WIDTH,
            "axes.labelsize": AXIS_LABEL_SIZE,
            "xtick.labelsize": TICK_LABEL_SIZE,
            "ytick.labelsize": TICK_LABEL_SIZE,
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
    axis.tick_params(
        which="minor", direction="in", length=6, width=1.0, top=True, right=True
    )
    axis.minorticks_on()
    axis.set_box_aspect(BOX_ASPECT)


def make_axis() -> tuple[plt.Figure, plt.Axes]:
    figure = plt.figure(figsize=(6.5, 5.84), facecolor="white")
    axis = figure.add_axes([0.186, 0.260, 0.792, 0.705])
    style_axis(axis)
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlabel(r"$Ek$")
    axis.set_ylabel(r"$\left\langle 2\pi L_m\right\rangle_t$")
    return figure, axis


def save(figure: plt.Figure, stem: Path) -> None:
    figure.savefig(stem.with_suffix(".png"), facecolor="white")
    figure.savefig(stem.with_suffix(".pdf"), facecolor="white")
    plt.close(figure)


def add_shared_multiplier(
    axis: plt.Axes,
    axis_name: str,
    exponent: int,
    ticks: np.ndarray,
) -> None:
    scale = 10.0**exponent
    target = axis.xaxis if axis_name == "x" else axis.yaxis
    target.set_major_locator(FixedLocator(ticks))
    target.set_major_formatter(FuncFormatter(lambda value, _: f"{value / scale:g}"))
    if axis_name == "x":
        axis.text(
            0.985,
            -0.125,
            rf"$\times10^{{{exponent}}}$",
            transform=axis.transAxes,
            ha="right",
            va="top",
            fontsize=TICK_LABEL_SIZE,
        )
    else:
        axis.text(
            0.005,
            0.985,
            rf"$\times10^{{{exponent}}}$",
            transform=axis.transAxes,
            ha="left",
            va="top",
            fontsize=TICK_LABEL_SIZE,
        )


def marker(axis: plt.Axes, x: np.ndarray, y: np.ndarray, color: tuple) -> None:
    axis.plot(
        x,
        y,
        linestyle="none",
        marker="o",
        markersize=MARKER_SIZE,
        markerfacecolor=color,
        markeredgecolor="black",
        markeredgewidth=MARKER_EDGE_WIDTH,
        zorder=5,
    )


def load_rows() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    intermittent = []
    funnel = []
    with INPUT.open("r", newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            item = (float(row["Ek"]), float(row["two_pi_Lm_mean"]))
            if row["selection"] == "intermittent low valley":
                intermittent.append(item)
            elif row["selection"] == "late statistical steady":
                funnel.append(item)
    intermittent = np.asarray(sorted(intermittent), dtype=float)
    funnel = np.asarray(sorted(funnel), dtype=float)
    return intermittent[:, 0], intermittent[:, 1], funnel[:, 0], funnel[:, 1]


def free_fit(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float]:
    power, log_prefactor = np.polyfit(np.log(x), np.log(y), 1)
    prefactor = float(np.exp(log_prefactor))
    prediction = prefactor * x**power
    residual = np.sum((np.log(y) - np.log(prediction)) ** 2)
    total = np.sum((np.log(y) - np.mean(np.log(y))) ** 2)
    return float(power), prefactor, float(1.0 - residual / total)


def fixed_fit(
    x: np.ndarray, y: np.ndarray, power: float
) -> tuple[float, float]:
    prefactor = float(np.exp(np.mean(np.log(y) - power * np.log(x))))
    prediction = prefactor * x**power
    residual = np.sum((np.log(y) - np.log(prediction)) ** 2)
    total = np.sum((np.log(y) - np.mean(np.log(y))) ** 2)
    return prefactor, float(1.0 - residual / total)


def main() -> None:
    configure_style()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    intermittent_ek, intermittent_lm, funnel_ek, funnel_lm = load_rows()
    intermittent_power, intermittent_prefactor, intermittent_r2 = free_fit(
        intermittent_ek, intermittent_lm
    )
    funnel_prefactor, funnel_r2 = fixed_fit(
        funnel_ek, funnel_lm, FUNNEL_FIXED_POWER
    )

    # Figure 1: intermittent plume-burst scaling.
    figure, axis = make_axis()
    xline = np.geomspace(1.85e-4, 7.6e-4, 200)
    axis.plot(
        xline,
        intermittent_prefactor * xline**intermittent_power,
        color="black",
        linestyle="--",
        linewidth=LINE_WIDTH,
        zorder=1,
    )
    marker(axis, intermittent_ek, intermittent_lm, INTERMITTENT_COLOR)
    axis.set_xlim(1.65e-4, 8.2e-4)
    axis.set_ylim(0.28, 0.86)
    add_shared_multiplier(axis, "x", -4, np.arange(2.0, 8.0) * 1.0e-4)
    add_shared_multiplier(axis, "y", -1, np.arange(3.0, 9.0) * 1.0e-1)
    axis.xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1))
    axis.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1))
    axis.text(
        1.78e-4,
        0.74,
        rf"$2\pi L_m\propto Ek^{{{intermittent_power:.2f}}}$" + "\n"
        + rf"$R^2={intermittent_r2:.3f}$",
        ha="left",
        va="top",
        fontsize=13,
    )
    intermittent_stem = OUTPUT_DIR / "Ra8e6_intermittent_lm_scaling"
    save(figure, intermittent_stem)

    # Figure 2: funnel-regime fixed 0.69 scaling.
    figure, axis = make_axis()
    xline = np.geomspace(8.8e-4, 1.12e-2, 220)
    axis.plot(
        xline,
        funnel_prefactor * xline**FUNNEL_FIXED_POWER,
        color="black",
        linestyle="--",
        linewidth=LINE_WIDTH,
        zorder=1,
    )
    marker(axis, funnel_ek, funnel_lm, FUNNEL_COLOR)
    axis.set_xlim(7.5e-4, 1.25e-2)
    axis.set_ylim(0.70, 4.2)
    add_shared_multiplier(
        axis, "x", -3, np.asarray([1, 2, 3, 5, 7, 10], dtype=float) * 1.0e-3
    )
    axis.yaxis.set_major_locator(FixedLocator(np.asarray([1, 2, 3, 4], dtype=float)))
    axis.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}"))
    axis.xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1))
    axis.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1))
    axis.text(
        9.0e-4,
        3.45,
        rf"$2\pi L_m\propto Ek^{{{FUNNEL_FIXED_POWER:.2f}}}$" + "\n"
        + rf"$R^2={funnel_r2:.3f}$",
        ha="left",
        va="top",
        fontsize=13,
    )
    funnel_stem = OUTPUT_DIR / "Ra8e6_funnel_lm_scaling"
    save(figure, funnel_stem)

    # Figure 3: both regimes, no single fit. Draw each regime-local scaling
    # segment through its points and extend it just beyond the outer markers.
    figure, axis = make_axis()
    marker(axis, intermittent_ek, intermittent_lm, INTERMITTENT_COLOR)
    # The combined regime figure intentionally stops below Ek=1e-2.
    combined_funnel = funnel_ek < 1.0e-2
    marker(
        axis,
        funnel_ek[combined_funnel],
        funnel_lm[combined_funnel],
        FUNNEL_COLOR,
    )
    ix = np.geomspace(1.75e-4, 7.8e-4, 150)
    fx = np.geomspace(8.5e-4, 8.0e-3, 180)
    axis.plot(
        ix,
        intermittent_prefactor * ix**intermittent_power,
        color="black",
        linestyle="--",
        linewidth=LINE_WIDTH,
        zorder=1,
    )
    axis.plot(
        fx,
        funnel_prefactor * fx**FUNNEL_FIXED_POWER,
        color="black",
        linestyle="--",
        linewidth=LINE_WIDTH,
        zorder=1,
    )
    axis.set_xlim(1.0e-4, 8.5e-3)
    axis.set_ylim(0.28, 5.2)
    axis.xaxis.set_major_locator(LogLocator(base=10.0, numticks=4))
    axis.xaxis.set_major_formatter(LogFormatterMathtext(base=10.0))
    axis.yaxis.set_major_locator(FixedLocator(np.asarray([0.3, 0.5, 1, 2, 3, 5])))
    axis.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}"))
    axis.xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1))
    axis.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1))
    axis.text(2.0e-4, 0.50, rf"$Ek^{{{intermittent_power:.2f}}}$", fontsize=15)
    axis.text(1.65e-3, 2.65, rf"$Ek^{{{FUNNEL_FIXED_POWER:.2f}}}$", fontsize=15)
    axis.legend(
        handles=(
            Line2D(
                [], [], linestyle="none", marker="o", markersize=9,
                markerfacecolor=INTERMITTENT_COLOR, markeredgecolor="black",
                markeredgewidth=1.5, label="Intermittent plume burst",
            ),
            Line2D(
                [], [], linestyle="none", marker="o", markersize=9,
                markerfacecolor=FUNNEL_COLOR, markeredgecolor="black",
                markeredgewidth=1.5, label="Funnel",
            ),
        ),
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(0.035, 0.925),
        borderaxespad=0.0,
        fontsize=15,
    )
    combined_stem = OUTPUT_DIR / "Ra8e6_intermittent_and_funnel_lm_scalings"
    save(figure, combined_stem)

    summary = OUTPUT_DIR / "Ra8e6_lm_regime_scaling_summary.csv"
    with summary.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("regime", "power", "prefactor", "r2", "point_count"),
        )
        writer.writeheader()
        writer.writerow(
            {
                "regime": "intermittent plume burst",
                "power": intermittent_power,
                "prefactor": intermittent_prefactor,
                "r2": intermittent_r2,
                "point_count": len(intermittent_ek),
            }
        )
        writer.writerow(
            {
                "regime": "funnel",
                "power": FUNNEL_FIXED_POWER,
                "prefactor": funnel_prefactor,
                "r2": funnel_r2,
                "point_count": len(funnel_ek),
            }
        )

    print(intermittent_stem.with_suffix(".png"))
    print(funnel_stem.with_suffix(".png"))
    print(combined_stem.with_suffix(".png"))
    print(summary)
    print(
        f"intermittent: {intermittent_prefactor:.8g} Ek^{intermittent_power:.8g}, "
        f"R2={intermittent_r2:.8g}"
    )
    print(
        f"funnel: {funnel_prefactor:.8g} Ek^{FUNNEL_FIXED_POWER:.8g}, "
        f"R2={funnel_r2:.8g}"
    )


if __name__ == "__main__":
    main()
