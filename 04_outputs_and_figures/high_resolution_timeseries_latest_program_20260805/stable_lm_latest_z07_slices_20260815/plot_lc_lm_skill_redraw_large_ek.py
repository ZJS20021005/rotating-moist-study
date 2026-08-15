from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter, LogFormatterMathtext, LogLocator
import numpy as np
import pandas as pd


BASE = Path(r"E:\moist RB\rotating_case_inventory\04_outputs_and_figures\high_resolution_timeseries_latest_program_20260805")
SOURCE = BASE / "high_resolution_timeseries_long.csv"
OUT = BASE / "stable_lm_latest_z07_slices_20260815" / "lc_lm_skill_redraw_large_Ek_20260815"
OUT.mkdir(parents=True, exist_ok=True)

CASES = {
    "Ek1e-3": (1.0e-3, Path(r"G:\moist convection\Ek1e-3")),
    "Ek3e-3": (3.0e-3, Path(r"G:\moist convection\Ek3e-3")),
    "Ek5e-3": (5.0e-3, Path(r"G:\moist convection\Ek5e-3")),
    "Ek7e-3": (7.0e-3, Path(r"G:\moist convection\Ek7e-3")),
    "Ek1e-2": (1.0e-2, Path(r"G:\moist convection\Ek1e-2")),
    "Ek3e-2": (3.0e-2, Path(r"G:\moist convection\Ek3e-2")),
    "Ek5e-2": (5.0e-2, Path(r"G:\moist convection\Ek5e-2")),
    "Ek1e-1": (1.0e-1, Path(r"H:\rotating_case\Pr0p7\Ra8e6\Ek1e-1")),
    "norotating": (None, Path(r"G:\moist convection\norotating")),
}
STABLE_LABELS = ("Ek1e-3", "Ek3e-3", "Ek5e-3", "Ek7e-3", "Ek1e-2")
COLORS = {
    "Ek1e-3": (0.74, 0.14, 0.18),
    "Ek3e-3": (0.93, 0.32, 0.23),
    "Ek5e-3": (0.96, 0.58, 0.19),
    "Ek7e-3": (0.16, 0.56, 0.80),
    "Ek1e-2": (0.11, 0.44, 0.71),
    "Ek3e-2": (0.00, 0.65, 0.12),
    "Ek5e-2": (0.58, 0.05, 0.90),
    "Ek1e-1": (0.90, 0.00, 0.50),
    "norotating": (0.42, 0.42, 0.42),
}
WINDOW = 500.0


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "axes.linewidth": 4.5,
            "axes.labelsize": 24,
            "xtick.labelsize": 13,
            "ytick.labelsize": 13,
            "legend.fontsize": 10,
            "figure.dpi": 180,
            "savefig.dpi": 300,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


def valid(path: Path) -> bool:
    text = str(path).lower()
    return "qbot0p5_qtop0p004978" in text and "\\run\\run\\" not in text


def read_raw(path: Path, minimum_columns: int = 3) -> list[tuple[float, float]]:
    rows = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            try:
                values = [float(item.replace("D", "E").replace("d", "e")) for item in line.split()]
            except ValueError:
                continue
            if len(values) >= minimum_columns and np.all(np.isfinite(values[:minimum_columns])) and values[1] > 0.0:
                rows.append((values[0], 2.0 * math.pi * values[1]))
    return rows


def merge_metric(old: pd.DataFrame, label: str, ek: float | None, root: Path, metric: str, filename: str) -> np.ndarray:
    values: dict[float, float] = {}
    if ek is None:
        historical = old[(old["metric"] == metric) & old["Ek"].isna()]
    else:
        historical = old[(old["metric"] == metric) & np.isclose(old["Ek"].to_numpy(dtype=float), ek, rtol=0.0, atol=1.0e-12)]
    for row in historical.itertuples(index=False):
        values[round(float(row.time), 8)] = float(row.value)
    try:
        raw_files = [path for path in root.rglob(filename) if valid(path)]
    except OSError:
        raw_files = []
    raw_files.sort(key=lambda path: (path.stat().st_mtime, len(path.parts), str(path)))
    for path in raw_files:
        try:
            for time, value in read_raw(path):
                values[round(time, 8)] = value
        except OSError:
            continue
    if not values:
        raise RuntimeError(f"No {metric} data for {label}")
    return np.asarray([[time, values[time]] for time in sorted(values)], dtype=float)


def load_series() -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    old = pd.read_csv(SOURCE)
    old = old[old["Ra"] == 8.0e6]
    lc: dict[str, np.ndarray] = {}
    lm: dict[str, np.ndarray] = {}
    for label, (ek, root) in CASES.items():
        lc[label] = merge_metric(old, label, ek, root, "convective_mid", "convective_scale.dat")
        lm[label] = merge_metric(old, label, ek, root, "mse_spectral", "moist_integral_scale.dat")
    return lc, lm


def ek_label(ek: float) -> str:
    exponent = int(math.floor(math.log10(ek)))
    coefficient = ek / 10.0**exponent
    if math.isclose(coefficient, 1.0, rel_tol=0.0, abs_tol=1.0e-10):
        return rf"$Ek=10^{{{exponent}}}$"
    return rf"$Ek={coefficient:g}\times10^{{{exponent}}}$"


def style_axis(ax: plt.Axes) -> None:
    ax.set_box_aspect(5.2 / 6.5)
    for spine in ax.spines.values():
        spine.set_linewidth(4.5)
    ax.tick_params(which="major", direction="in", length=12, width=1.2, top=True, right=True)
    ax.tick_params(which="minor", direction="in", length=6, width=1.0, top=True, right=True)
    ax.minorticks_on()


def pair_figure() -> tuple[plt.Figure, tuple[plt.Axes, plt.Axes]]:
    fig = plt.figure(figsize=(13.0, 5.84), facecolor="white")
    axes = (
        fig.add_axes([0.093, 0.260, 0.396, 0.705]),
        fig.add_axes([0.593, 0.260, 0.396, 0.705]),
    )
    for ax in axes:
        style_axis(ax)
    return fig, axes


def single_figure() -> tuple[plt.Figure, plt.Axes]:
    fig = plt.figure(figsize=(6.5, 5.84), facecolor="white")
    ax = fig.add_axes([0.186, 0.260, 0.792, 0.705])
    style_axis(ax)
    return fig, ax


def late_summary(series: dict[str, np.ndarray], labels: tuple[str, ...]) -> pd.DataFrame:
    rows = []
    for label in labels:
        data = series[label]
        end = float(data[-1, 0])
        selected = data[data[:, 0] >= end - WINDOW]
        rows.append(
            {
                "case": label,
                "Ek": CASES[label][0],
                "time_start": float(selected[0, 0]),
                "time_end": end,
                "n_samples": len(selected),
                "mean": float(np.mean(selected[:, 1])),
                "std": float(np.std(selected[:, 1], ddof=1)),
            }
        )
    return pd.DataFrame(rows)


def power_fit(summary: pd.DataFrame) -> dict[str, float]:
    log_x = np.log(summary["Ek"].to_numpy(dtype=float))
    log_y = np.log(summary["mean"].to_numpy(dtype=float))
    alpha, log_c = np.polyfit(log_x, log_y, 1)
    prediction = alpha * log_x + log_c
    r2 = 1.0 - np.sum((log_y - prediction) ** 2) / np.sum((log_y - np.mean(log_y)) ** 2)
    return {"alpha": float(alpha), "C": float(np.exp(log_c)), "R2_log_space": float(r2)}


def configure_log_axis(ax: plt.Axes, quantity: str) -> None:
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$Ek$")
    ax.xaxis.set_major_locator(FixedLocator([1.0e-3, 1.0e-2]))
    ax.xaxis.set_major_formatter(LogFormatterMathtext(10))
    ax.xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1))
    if quantity == "lc":
        ax.set_ylabel(r"$\left\langle2\pi l_c(z\simeq0.5)\right\rangle_t/H$")
        ticks = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
        ax.set_ylim(0.43, 1.03)
    else:
        ax.set_ylabel(r"$\left\langle2\pi L_m\right\rangle_t/H$")
        ticks = [1.0, 2.0, 3.0]
        ax.set_ylim(0.78, 3.75)
    ax.yaxis.set_major_locator(FixedLocator(ticks))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}"))
    ax.set_xlim(7.5e-4, 1.3e-2)


def draw_scaling_panel(ax: plt.Axes, summary: pd.DataFrame, fit: dict[str, float], quantity: str) -> None:
    configure_log_axis(ax, quantity)
    x = summary["Ek"].to_numpy(dtype=float)
    y = summary["mean"].to_numpy(dtype=float)
    xline = np.logspace(np.log10(0.85 * x.min()), np.log10(1.15 * x.max()), 300)
    ax.plot(xline, fit["C"] * xline ** fit["alpha"], color="black", ls="--", lw=3.5, zorder=1, label=rf"$Ek^{{{fit['alpha']:.3f}}}$")
    ax.plot(x, y, linestyle="none", marker="o", markersize=12, markerfacecolor=(0.11, 0.44, 0.71), markeredgecolor="black", markeredgewidth=2.0, zorder=3, label="late-time mean")
    ax.legend(frameon=False, loc="lower right", handlelength=2.2, fontsize=10)


def draw_scaling(lc: dict[str, np.ndarray], lm: dict[str, np.ndarray]) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, float], dict[str, float]]:
    lc_summary = late_summary(lc, STABLE_LABELS)
    lm_summary = late_summary(lm, STABLE_LABELS)
    lc_fit = power_fit(lc_summary)
    lm_fit = power_fit(lm_summary)
    fig, axes = pair_figure()
    draw_scaling_panel(axes[0], lc_summary, lc_fit, "lc")
    draw_scaling_panel(axes[1], lm_summary, lm_fit, "lm")
    fig.savefig(OUT / "Ra8e6_lc_Lm_vs_Ek_loglog_skill_redraw.png")
    fig.savefig(OUT / "Ra8e6_lc_Lm_vs_Ek_loglog_skill_redraw.pdf")
    plt.close(fig)
    for quantity, summary, fit in (("lc", lc_summary, lc_fit), ("Lm", lm_summary, lm_fit)):
        fig, ax = single_figure()
        draw_scaling_panel(ax, summary, fit, "lc" if quantity == "lc" else "lm")
        fig.savefig(OUT / f"Ra8e6_{quantity}_vs_Ek_loglog_skill_redraw.png")
        fig.savefig(OUT / f"Ra8e6_{quantity}_vs_Ek_loglog_skill_redraw.pdf")
        plt.close(fig)
    return lc_summary, lm_summary, lc_fit, lm_fit


def with_gap_breaks(data: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    if len(data) < 2:
        return data[:, 0], data[:, 1]
    differences = np.diff(data[:, 0])
    positive = differences[differences > 0.0]
    typical_dt = float(np.median(positive))
    x_values: list[float] = [float(data[0, 0])]
    y_values: list[float] = [float(data[0, 1])]
    for index in range(1, len(data)):
        if data[index, 0] - data[index - 1, 0] > 1.51 * typical_dt:
            x_values.append(math.nan)
            y_values.append(math.nan)
        x_values.append(float(data[index, 0]))
        y_values.append(float(data[index, 1]))
    return np.asarray(x_values), np.asarray(y_values)


def draw_timeseries_panel(ax: plt.Axes, series: dict[str, np.ndarray], quantity: str, show_legend: bool = True) -> None:
    for label, (ek, _) in CASES.items():
        data = series[label]
        x_values, y_values = with_gap_breaks(data)
        if ek is None:
            ax.plot(x_values, y_values, color=COLORS[label], lw=3.5, ls="--", label="nonrotating")
        else:
            ax.plot(x_values, y_values, color=COLORS[label], lw=3.5, label=ek_label(ek))
    ax.set_xlabel(r"$t$")
    ax.set_xlim(left=0.0)
    ax.set_ylim(bottom=0.0)
    if quantity == "lc":
        ax.set_ylabel(r"$2\pi l_c(z\simeq0.5)/H$")
    else:
        ax.set_ylabel(r"$2\pi L_m/H$")
    if show_legend:
        ax.legend(frameon=False, ncol=2, loc="lower right", handlelength=2.0, columnspacing=0.9, fontsize=9)


def draw_timeseries(lc: dict[str, np.ndarray], lm: dict[str, np.ndarray]) -> None:
    fig, axes = pair_figure()
    draw_timeseries_panel(axes[0], lc, "lc", show_legend=False)
    draw_timeseries_panel(axes[1], lm, "lm", show_legend=False)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, ncol=5, loc="lower center", bbox_to_anchor=(0.5, 0.012), handlelength=2.0, columnspacing=1.0, fontsize=9)
    fig.savefig(OUT / "Ra8e6_lc_Lm_timeseries_Ek1e-3_to_1e-1_with_nonrotating_skill_redraw.png")
    fig.savefig(OUT / "Ra8e6_lc_Lm_timeseries_Ek1e-3_to_1e-1_with_nonrotating_skill_redraw.pdf")
    plt.close(fig)
    for quantity, series in (("lc", lc), ("Lm", lm)):
        fig, ax = single_figure()
        draw_timeseries_panel(ax, series, "lc" if quantity == "lc" else "lm")
        fig.savefig(OUT / f"Ra8e6_{quantity}_timeseries_Ek1e-3_to_1e-1_with_nonrotating_skill_redraw.png")
        fig.savefig(OUT / f"Ra8e6_{quantity}_timeseries_Ek1e-3_to_1e-1_with_nonrotating_skill_redraw.pdf")
        plt.close(fig)


def save_series(lc: dict[str, np.ndarray], lm: dict[str, np.ndarray]) -> None:
    with (OUT / "Ra8e6_lc_Lm_timeseries_Ek1e-3_to_1e-1_with_nonrotating.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(["case", "Ek", "metric", "time", "value"])
        for metric, series in (("two_pi_lc_z05_H", lc), ("two_pi_Lm_H", lm)):
            for label, data in series.items():
                for time, value in data:
                    writer.writerow([label, CASES[label][0], metric, time, value])


def main() -> None:
    configure_style()
    lc, lm = load_series()
    lc_summary, lm_summary, lc_fit, lm_fit = draw_scaling(lc, lm)
    draw_timeseries(lc, lm)
    save_series(lc, lm)
    lc_summary.to_csv(OUT / "Ra8e6_lc_late500_summary.csv", index=False, encoding="utf-8-sig")
    lm_summary.to_csv(OUT / "Ra8e6_Lm_late500_summary.csv", index=False, encoding="utf-8-sig")
    (OUT / "Ra8e6_lc_Lm_powerlaw_fits.json").write_text(json.dumps({"lc": lc_fit, "Lm": lm_fit}, indent=2), encoding="utf-8")
    print(json.dumps({"lc": lc_fit, "Lm": lm_fit}, indent=2))
    for label in CASES:
        print(label, lc[label][0, 0], lc[label][-1, 0], lm[label][0, 0], lm[label][-1, 0])
    print(OUT)


if __name__ == "__main__":
    main()
