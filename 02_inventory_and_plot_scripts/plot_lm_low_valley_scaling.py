from __future__ import annotations

import csv
import math
from collections import defaultdict
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
        if (path / "high_resolution_timeseries_long.csv").exists()
    ),
    key=lambda path: path.name,
)
SOURCE = SERIES_DIR / "high_resolution_timeseries_long.csv"
OUTPUT_DIR = SERIES_DIR / "lm_low_valley_scaling"

TARGET_RA = 8.0e6
INTERMITTENT_CASES = (2.0e-4, 5.0e-4, 7.0e-4)
STABLE_CASES = {
    1.0e-3: 500.0,
    3.0e-3: 500.0,
    5.0e-3: 500.0,
    7.0e-3: 500.0,
    1.0e-2: 1000.0,
}
SMOOTHING_TIME = 10.0
MIN_LOW_INTERVAL = 20.0
ANALYSIS_START = 100.0
FLAT_SLOPE_LIMIT = 3.0e-3
LARGE_VORTEX_POWER = 0.69

BLUE = (0.00, 0.25, 0.90)
LIGHT_BLUE = (0.15, 0.55, 1.00)
FRAME_WIDTH = 4.5
LINE_WIDTH = 3.5
AXIS_LABEL_SIZE = 24
TICK_LABEL_SIZE = 13
LEGEND_SIZE = 10
BOX_ASPECT = 5.2 / 6.5


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Times New Roman",
            "mathtext.fontset": "stix",
            "axes.linewidth": FRAME_WIDTH,
            "axes.labelsize": AXIS_LABEL_SIZE,
            "xtick.labelsize": TICK_LABEL_SIZE,
            "ytick.labelsize": TICK_LABEL_SIZE,
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


def close(value: float, target: float) -> bool:
    return math.isclose(value, target, rel_tol=1.0e-8, abs_tol=1.0e-12)


def load_histories() -> dict[float, np.ndarray]:
    histories: dict[float, list[tuple[float, float]]] = defaultdict(list)
    with SOURCE.open("r", newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            if row.get("metric") != "mse_spectral":
                continue
            if not close(float(row["Ra"]), TARGET_RA) or not row.get("Ek"):
                continue
            histories[float(row["Ek"])].append(
                (float(row["time"]), float(row["value"]))
            )
    return {
        ek: np.asarray(sorted(values), dtype=float) for ek, values in histories.items()
    }


def matching_history(histories: dict[float, np.ndarray], target: float) -> np.ndarray:
    key = min(histories, key=lambda value: abs(value - target))
    if not close(key, target):
        raise RuntimeError(f"No history found for Ek={target:g}")
    return histories[key]


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


def identify_low_valleys(history: np.ndarray) -> tuple[dict, np.ndarray, np.ndarray]:
    working = history[history[:, 0] >= ANALYSIS_START]
    dt = float(np.median(np.diff(working[:, 0])))
    window = max(3, int(round(SMOOTHING_TIME / dt)))
    smoothed_full = np.convolve(
        working[:, 1], np.ones(window, dtype=float) / window, mode="same"
    )
    index = np.arange(len(working))
    valid = (index >= window // 2) & (index < len(working) - window // 2)
    time = working[valid, 0]
    raw = working[valid, 1]
    smoothed = smoothed_full[valid]

    centers, labels = two_means(smoothed)
    low_cluster = int(np.argmin(centers))
    low_mask = labels == low_cluster

    # The initially tiny scale belongs to linear onset, not to a plume valley.
    high_indices = np.flatnonzero(~low_mask)
    if len(high_indices) == 0:
        raise RuntimeError("No separated high-Lm state was detected")
    low_mask[: int(high_indices[0])] = False

    # A low value alone is insufficient: reject the falling and regrowing
    # shoulders and retain only the nearly horizontal plume-scale plateau.
    slope = np.gradient(smoothed, time)
    flat_low_mask = low_mask & (np.abs(slope) <= FLAT_SLOPE_LIMIT)
    accepted = np.zeros_like(flat_low_mask, dtype=bool)
    edges = np.flatnonzero(
        np.diff(np.r_[False, flat_low_mask, False])
    ).reshape(-1, 2)
    intervals = []
    for start, stop in edges:
        duration = float(time[stop - 1] - time[start])
        if duration < MIN_LOW_INTERVAL:
            continue
        accepted[start:stop] = True
        intervals.append(
            {
                "time_min": float(time[start]),
                "time_max": float(time[stop - 1]),
                "duration": duration,
                "sample_count": int(stop - start),
            }
        )
    if not np.any(accepted):
        raise RuntimeError("No persistent low-Lm interval survived the duration test")

    selected = raw[accepted]
    summary = {
        "low_center": float(np.min(centers)),
        "high_center": float(np.max(centers)),
        "cluster_separation_ratio": float(np.max(centers) / np.min(centers)),
        "selection_threshold": float(np.mean(centers)),
        "flat_slope_limit": FLAT_SLOPE_LIMIT,
        "two_pi_Lm_mean": float(np.mean(selected)),
        "two_pi_Lm_std": float(np.std(selected)),
        "sample_count": int(len(selected)),
        "interval_count": int(len(intervals)),
        "intervals": intervals,
    }
    diagnostic = np.column_stack((time, raw, smoothed, accepted.astype(float)))
    return summary, diagnostic, accepted


def late_mean(history: np.ndarray, window: float) -> dict:
    time_max = float(history[-1, 0])
    selected = history[history[:, 0] >= time_max - window]
    return {
        "time_min": float(selected[0, 0]),
        "time_max": time_max,
        "sample_count": int(len(selected)),
        "two_pi_Lm_mean": float(np.mean(selected[:, 1])),
        "two_pi_Lm_std": float(np.std(selected[:, 1])),
    }


def main() -> None:
    configure_style()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    histories = load_histories()
    rows = []
    diagnostics = {}

    for ek in INTERMITTENT_CASES:
        history = matching_history(histories, ek)
        result, diagnostic, _ = identify_low_valleys(history)
        diagnostics[ek] = diagnostic
        rows.append(
            {
                "Ra": TARGET_RA,
                "Ek": ek,
                "selection": "intermittent low valley",
                "time_window": "; ".join(
                    f"{item['time_min']:.1f}-{item['time_max']:.1f}"
                    for item in result["intervals"]
                ),
                "sample_count": result["sample_count"],
                "two_pi_Lm_mean": result["two_pi_Lm_mean"],
                "two_pi_Lm_std": result["two_pi_Lm_std"],
                "low_center": result["low_center"],
                "high_center": result["high_center"],
                "cluster_separation_ratio": result["cluster_separation_ratio"],
            }
        )

    for ek, window in STABLE_CASES.items():
        result = late_mean(matching_history(histories, ek), window)
        rows.append(
            {
                "Ra": TARGET_RA,
                "Ek": ek,
                "selection": "late statistical steady",
                "time_window": f"{result['time_min']:.1f}-{result['time_max']:.1f}",
                "sample_count": result["sample_count"],
                "two_pi_Lm_mean": result["two_pi_Lm_mean"],
                "two_pi_Lm_std": result["two_pi_Lm_std"],
                "low_center": "",
                "high_center": "",
                "cluster_separation_ratio": "",
            }
        )
    rows.sort(key=lambda row: row["Ek"])

    table = OUTPUT_DIR / "Ra8e6_lm_low_valley_and_stable_time_averages.csv"
    with table.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    ek = np.asarray([row["Ek"] for row in rows], dtype=float)
    lm = np.asarray([row["two_pi_Lm_mean"] for row in rows], dtype=float)
    large_rows = [
        row for row in rows if row["selection"] == "late statistical steady"
    ]
    large_ek = np.asarray([row["Ek"] for row in large_rows], dtype=float)
    large_lm = np.asarray(
        [row["two_pi_Lm_mean"] for row in large_rows], dtype=float
    )
    # Only the prefactor is fitted. The exponent is the previously selected
    # large-vortex value and is not inferred from the intermittent points.
    fit_prefactor = float(
        np.exp(np.mean(np.log(large_lm) - LARGE_VORTEX_POWER * np.log(large_ek)))
    )
    fitted = fit_prefactor * large_ek**LARGE_VORTEX_POWER
    residual = float(np.sum((np.log(large_lm) - np.log(fitted)) ** 2))
    total = float(np.sum((np.log(large_lm) - np.mean(np.log(large_lm))) ** 2))
    fit_r2 = 1.0 - residual / total

    stale_fit = OUTPUT_DIR / "Ra8e6_lm_low_valley_and_stable_powerlaw_fit.csv"
    if stale_fit.exists():
        stale_fit.unlink()
    fit_table = OUTPUT_DIR / "Ra8e6_large_vortex_fixed_069_fit.csv"
    with fit_table.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("fixed_power", "fit_prefactor", "fit_r2", "point_count"),
        )
        writer.writeheader()
        writer.writerow(
            {
                "fixed_power": LARGE_VORTEX_POWER,
                "fit_prefactor": fit_prefactor,
                "fit_r2": fit_r2,
                "point_count": len(rows),
            }
        )

    intermittent_rows = [
        row for row in rows if row["selection"] == "intermittent low valley"
    ]
    intermittent_ek = np.asarray(
        [row["Ek"] for row in intermittent_rows], dtype=float
    )
    intermittent_lm = np.asarray(
        [row["two_pi_Lm_mean"] for row in intermittent_rows], dtype=float
    )
    intermittent_power, intermittent_log_prefactor = np.polyfit(
        np.log(intermittent_ek), np.log(intermittent_lm), 1
    )
    intermittent_prefactor = float(np.exp(intermittent_log_prefactor))
    intermittent_fitted = intermittent_prefactor * intermittent_ek**intermittent_power
    intermittent_residual = float(
        np.sum((np.log(intermittent_lm) - np.log(intermittent_fitted)) ** 2)
    )
    intermittent_total = float(
        np.sum((np.log(intermittent_lm) - np.mean(np.log(intermittent_lm))) ** 2)
    )
    intermittent_r2 = 1.0 - intermittent_residual / intermittent_total
    intermittent_fit_table = OUTPUT_DIR / "Ra8e6_intermittent_flat_valley_fit.csv"
    with intermittent_fit_table.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("fit_power", "fit_prefactor", "fit_r2", "point_count"),
        )
        writer.writeheader()
        writer.writerow(
            {
                "fit_power": intermittent_power,
                "fit_prefactor": intermittent_prefactor,
                "fit_r2": intermittent_r2,
                "point_count": len(intermittent_rows),
            }
        )

    figure = plt.figure(figsize=(6.5, 5.84), facecolor="white")
    axis = figure.add_axes([0.186, 0.260, 0.792, 0.705])
    style_axis(axis)
    axis.set_xscale("log")
    axis.set_yscale("log")
    # Draw short, parallel reference segments beside the points. They are not
    # connecting lines and do not visually pass through the markers.
    early_line_ek = np.geomspace(1.9e-4, 7.4e-4, 120)
    large_line_ek = np.geomspace(1.15e-3, 9.0e-3, 160)
    axis.plot(
        early_line_ek,
        intermittent_prefactor * early_line_ek**intermittent_power,
        color="black",
        linestyle="--",
        linewidth=LINE_WIDTH,
        zorder=1,
    )
    axis.plot(
        large_line_ek,
        1.28 * fit_prefactor * large_line_ek**LARGE_VORTEX_POWER,
        color="black",
        linestyle="--",
        linewidth=LINE_WIDTH,
        zorder=1,
    )
    for row in rows:
        intermittent = row["selection"] == "intermittent low valley"
        axis.plot(
            row["Ek"],
            row["two_pi_Lm_mean"],
            linestyle="none",
            marker="^" if intermittent else "s",
            markersize=12,
            markerfacecolor=BLUE,
            markeredgecolor=BLUE,
            color=BLUE,
            zorder=3,
        )
    axis.set_xlabel(r"$Ek$")
    axis.set_ylabel(r"$\left\langle 2\pi L_m\right\rangle_t$")
    axis.set_xlim(1.7e-4, 1.35e-2)
    axis.set_ylim(0.25, 5.2)
    axis.xaxis.set_major_locator(LogLocator(base=10.0))
    axis.xaxis.set_major_formatter(LogFormatterMathtext(base=10.0))
    axis.xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1))
    axis.yaxis.set_major_locator(LogLocator(base=10.0))
    axis.yaxis.set_major_formatter(LogFormatterMathtext(base=10.0))
    axis.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1))
    axis.text(
        2.25e-4, 0.43,
        rf"$Ek^{{{intermittent_power:.2f}}}$",
        ha="left", va="bottom", fontsize=15,
    )
    axis.text(1.75e-3, 2.45, r"$Ek^{0.69}$", ha="left", va="bottom", fontsize=15)
    axis.legend(
        handles=(
            Line2D(
                [], [], linestyle="none", marker="^", markersize=10,
                markerfacecolor=BLUE, markeredgecolor=BLUE,
                label="Intermittent low valley",
            ),
            Line2D(
                [], [], linestyle="none", marker="s", markersize=10,
                markerfacecolor=BLUE, markeredgecolor=BLUE,
                label="Statistical steady",
            ),
        ),
        frameon=False,
        loc="lower right",
        fontsize=LEGEND_SIZE,
    )
    stem = OUTPUT_DIR / "Ra8e6_lm_low_valley_and_stable_vs_ek_loglog"
    figure.savefig(stem.with_suffix(".png"), facecolor="white")
    figure.savefig(stem.with_suffix(".pdf"), facecolor="white")
    plt.close(figure)

    # Standalone fit requested for the intermittent plume-valley regime.
    intermittent_figure = plt.figure(figsize=(6.5, 5.84), facecolor="white")
    intermittent_axis = intermittent_figure.add_axes([0.186, 0.260, 0.792, 0.705])
    style_axis(intermittent_axis)
    intermittent_axis.set_xscale("log")
    intermittent_axis.set_yscale("log")
    intermittent_line_ek = np.geomspace(1.85e-4, 7.6e-4, 200)
    intermittent_axis.plot(
        intermittent_line_ek,
        intermittent_prefactor * intermittent_line_ek**intermittent_power,
        color="black",
        linestyle="--",
        linewidth=LINE_WIDTH,
        zorder=1,
    )
    intermittent_axis.plot(
        intermittent_ek,
        intermittent_lm,
        linestyle="none",
        marker="^",
        markersize=13,
        markerfacecolor=BLUE,
        markeredgecolor=BLUE,
        color=BLUE,
        zorder=3,
    )
    intermittent_axis.set_xlabel(r"$Ek$")
    intermittent_axis.set_ylabel(r"$\left\langle 2\pi L_m\right\rangle_t$")
    intermittent_axis.set_xlim(1.65e-4, 8.2e-4)
    intermittent_axis.set_ylim(0.28, 0.86)
    # Narrow log ranges use shared scientific multipliers per the project
    # style: plain mantissas on ticks, x exponent at lower right, y exponent
    # at upper left. Do not repeat 10^n in every tick label.
    intermittent_axis.xaxis.set_major_locator(
        FixedLocator(np.arange(2.0, 8.0) * 1.0e-4)
    )
    intermittent_axis.xaxis.set_major_formatter(
        FuncFormatter(lambda value, _: f"{value / 1.0e-4:g}")
    )
    intermittent_axis.xaxis.set_minor_locator(
        LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1)
    )
    intermittent_axis.yaxis.set_major_locator(
        FixedLocator(np.arange(3.0, 9.0) * 1.0e-1)
    )
    intermittent_axis.yaxis.set_major_formatter(
        FuncFormatter(lambda value, _: f"{value / 1.0e-1:g}")
    )
    intermittent_axis.yaxis.set_minor_locator(
        LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1)
    )
    intermittent_axis.text(
        0.985,
        -0.125,
        r"$\times10^{-4}$",
        transform=intermittent_axis.transAxes,
        ha="right",
        va="top",
        fontsize=TICK_LABEL_SIZE,
    )
    intermittent_axis.text(
        0.005,
        0.985,
        r"$\times10^{-1}$",
        transform=intermittent_axis.transAxes,
        ha="left",
        va="top",
        fontsize=TICK_LABEL_SIZE,
    )
    intermittent_axis.text(
        1.78e-4,
        0.74,
        rf"$2\pi L_m\propto Ek^{{{intermittent_power:.2f}}}$" + "\n"
        + rf"$R^2={intermittent_r2:.3f}$",
        ha="left",
        va="top",
        fontsize=13,
    )
    intermittent_stem = (
        OUTPUT_DIR / "Ra8e6_intermittent_flat_valley_lm_vs_ek_loglog_fit"
    )
    intermittent_figure.savefig(
        intermittent_stem.with_suffix(".png"), facecolor="white"
    )
    intermittent_figure.savefig(
        intermittent_stem.with_suffix(".pdf"), facecolor="white"
    )
    plt.close(intermittent_figure)

    # Selection audit: show exactly which time intervals entered each low-valley mean.
    audit = plt.figure(figsize=(6.5, 9.0), facecolor="white")
    axes = audit.subplots(3, 1)
    for axis, ek in zip(axes, INTERMITTENT_CASES):
        style_axis(axis)
        data = diagnostics[ek]
        axis.plot(data[:, 0], data[:, 1], color=LIGHT_BLUE, linewidth=2.0)
        axis.plot(data[:, 0], data[:, 2], color=BLUE, linewidth=LINE_WIDTH)
        chosen = data[:, 3] > 0.5
        axis.fill_between(
            data[:, 0], 0.0, data[:, 2], where=chosen,
            color=BLUE, alpha=0.18, linewidth=0.0,
        )
        axis.set_ylabel(r"$2\pi L_m$")
        axis.text(0.96, 0.90, rf"$Ek={ek:g}$", transform=axis.transAxes,
                  ha="right", va="top", fontsize=13)
    axes[-1].set_xlabel(r"$t$")
    audit.subplots_adjust(left=0.19, right=0.96, bottom=0.09, top=0.98, hspace=0.16)
    audit_stem = OUTPUT_DIR / "Ra8e6_intermittent_lm_low_valley_selection"
    audit.savefig(audit_stem.with_suffix(".png"), facecolor="white")
    audit.savefig(audit_stem.with_suffix(".pdf"), facecolor="white")
    plt.close(audit)

    print(table)
    print(fit_table)
    print(intermittent_fit_table)
    print(stem.with_suffix(".png"))
    print(intermittent_stem.with_suffix(".png"))
    print(audit_stem.with_suffix(".png"))
    print(
        f"large-vortex fixed-slope comparison: "
        f"2piLm={fit_prefactor:.8g}*Ek^{LARGE_VORTEX_POWER:.8g}, "
        f"R2={fit_r2:.8g}"
    )
    print(
        f"intermittent flat-valley fit: "
        f"2piLm={intermittent_prefactor:.8g}*Ek^{intermittent_power:.8g}, "
        f"R2={intermittent_r2:.8g}"
    )
    for row in rows:
        print(
            f"Ek={row['Ek']:.4g}, selection={row['selection']}, "
            f"window={row['time_window']}, mean={row['two_pi_Lm_mean']:.8g}"
        )


if __name__ == "__main__":
    main()
