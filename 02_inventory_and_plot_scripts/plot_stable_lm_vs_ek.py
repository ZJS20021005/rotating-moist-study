from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import LogFormatterMathtext, LogLocator


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
OUTPUT_DIR = SERIES_DIR / "lm_stable_scaling"

TARGET_RA = 8.0e6
# Current stability audit (2026-08-03 13:16): include only established,
# statistically stationary, and newly stable cases.  The near-steady Ek=3e-2
# point is deliberately excluded because its 2*pi*Lm still drifts by about
# 12% over the latest 500 time units.  Statistical steady cases use a longer
# averaging window to sample their low-frequency fluctuations.
TARGET_CASES = {
    1.0e-3: {"status": "newly stable", "late_window": 500.0},
    3.0e-3: {"status": "stable", "late_window": 500.0},
    5.0e-3: {"status": "stable", "late_window": 500.0},
    7.0e-3: {"status": "stable", "late_window": 500.0},
    1.0e-2: {"status": "statistical steady", "late_window": 1000.0},
}
BLUE = (0.00, 0.25, 0.90)

FRAME_WIDTH = 4.5
AXIS_LABEL_SIZE = 24
TICK_LABEL_SIZE = 13
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


def load_late_means() -> list[dict]:
    histories: dict[float, list[tuple[float, float]]] = defaultdict(list)
    with SOURCE.open("r", newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            if row.get("metric") != "mse_spectral":
                continue
            if not close(float(row["Ra"]), TARGET_RA):
                continue
            if not row.get("Ek"):
                continue
            ek = float(row["Ek"])
            match = next((target for target in TARGET_CASES if close(ek, target)), None)
            if match is not None:
                histories[match].append((float(row["time"]), float(row["value"])))

    summary = []
    for ek, selection in TARGET_CASES.items():
        values = np.asarray(sorted(histories[ek]), dtype=float)
        if values.size == 0:
            raise RuntimeError(f"No mse_spectral history for Ek={ek:g}")
        time_max = float(values[-1, 0])
        late_window = float(selection["late_window"])
        late = values[values[:, 0] >= time_max - late_window]
        summary.append(
            {
                "Ra": TARGET_RA,
                "Ek": ek,
                "time_min": float(late[0, 0]),
                "time_max": time_max,
                "sample_count": int(len(late)),
                "late_window": late_window,
                "two_pi_Lm_mean": float(np.mean(late[:, 1])),
                "stability_category": selection["status"],
                "included_in_fit": True,
            }
        )
    return summary


def main() -> None:
    configure_style()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    # Remove superseded selections so the active result folder contains only
    # the current strict-stability version.
    for stale_name in (
        "Ra8e6_stable_cases_lm_vs_ek_last200.csv",
        "Ra8e6_stable_cases_lm_vs_ek_powerlaw_fit.csv",
        "Ra8e6_stable_cases_lm_vs_ek_loglog.png",
        "Ra8e6_stable_cases_lm_vs_ek_loglog.pdf",
        "Ra8e6_stable_and_near_stable_lm_vs_ek_last200.csv",
        "Ra8e6_stable_and_near_stable_lm_vs_ek_powerlaw_fit.csv",
        "Ra8e6_stable_and_near_stable_lm_vs_ek_loglog.png",
        "Ra8e6_stable_and_near_stable_lm_vs_ek_loglog.pdf",
    ):
        stale_path = OUTPUT_DIR / stale_name
        if stale_path.exists():
            stale_path.unlink()
    summary = load_late_means()

    table = OUTPUT_DIR / "Ra8e6_stable_only_lm_vs_ek_time_averaged.csv"
    with table.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary[0]))
        writer.writeheader()
        writer.writerows(summary)

    ek = np.asarray([row["Ek"] for row in summary])
    lm = np.asarray([row["two_pi_Lm_mean"] for row in summary])
    fit_power, fit_log_prefactor = np.polyfit(np.log(ek), np.log(lm), 1)
    fit_prefactor = float(np.exp(fit_log_prefactor))
    fitted_at_points = fit_prefactor * ek**fit_power
    residual = np.sum((np.log(lm) - np.log(fitted_at_points)) ** 2)
    total = np.sum((np.log(lm) - np.mean(np.log(lm))) ** 2)
    fit_r2 = float(1.0 - residual / total)

    fit_table = OUTPUT_DIR / "Ra8e6_stable_only_lm_vs_ek_powerlaw_fit.csv"
    with fit_table.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("fit_power", "fit_prefactor", "fit_r2", "point_count"),
        )
        writer.writeheader()
        writer.writerow(
            {
                "fit_power": fit_power,
                "fit_prefactor": fit_prefactor,
                "fit_r2": fit_r2,
                "point_count": len(ek),
            }
        )

    figure = plt.figure(figsize=(6.5, 5.84), facecolor="white")
    axis = figure.add_axes([0.186, 0.260, 0.792, 0.705])
    style_axis(axis)
    axis.set_xscale("log")
    axis.set_yscale("log")
    fit_ek = np.geomspace(float(np.min(ek)), float(np.max(ek)), 200)
    axis.plot(
        fit_ek,
        fit_prefactor * fit_ek**fit_power,
        color="black",
        linestyle="--",
        linewidth=3.5,
        zorder=1,
    )
    marker_by_status = {
        "stable": ("s", "Stable"),
        "statistical steady": ("D", "Statistical steady"),
        "newly stable": ("o", "Newly stable"),
    }
    plotted_labels = set()
    for row in summary:
        marker, label = marker_by_status[row["stability_category"]]
        axis.plot(
            row["Ek"],
            row["two_pi_Lm_mean"],
            linestyle="none",
            marker=marker,
            markersize=12,
            markerfacecolor=BLUE,
            markeredgecolor=BLUE,
            color=BLUE,
            label=label if label not in plotted_labels else None,
            zorder=2,
        )
        plotted_labels.add(label)
    axis.set_xlabel(r"$Ek$")
    axis.set_ylabel(r"$\left\langle 2\pi L_m\right\rangle_t$")
    axis.set_xlim(7.0e-4, 1.35e-2)
    axis.set_ylim(0.7, 4.2)
    # Logarithmic axes use explicit math-scientific notation (10^n).  This is
    # the log-axis counterpart of the shared scientific multiplier used by
    # the project style on linear axes.
    axis.xaxis.set_major_locator(LogLocator(base=10.0))
    axis.xaxis.set_major_formatter(LogFormatterMathtext(base=10.0))
    axis.xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1))
    axis.yaxis.set_major_locator(LogLocator(base=10.0))
    axis.yaxis.set_major_formatter(LogFormatterMathtext(base=10.0))
    axis.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1))
    axis.text(
        0.07,
        0.91,
        rf"$2\pi L_m\propto Ek^{{{fit_power:.2f}}}$" + "\n" + rf"$R^2={fit_r2:.4f}$",
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=13,
    )
    axis.legend(frameon=False, loc="lower right", fontsize=10)

    stem = OUTPUT_DIR / "Ra8e6_stable_only_lm_vs_ek_loglog"
    figure.savefig(stem.with_suffix(".png"), facecolor="white")
    figure.savefig(stem.with_suffix(".pdf"), facecolor="white")
    plt.close(figure)

    print(table)
    print(fit_table)
    print(stem.with_suffix(".png"))
    print(f"fit: 2piLm={fit_prefactor:.8g}*Ek^{fit_power:.8g}, R2={fit_r2:.8g}")
    for row in summary:
        print(f"Ek={row['Ek']:.4g}, <2piLm>={row['two_pi_Lm_mean']:.8g}")


if __name__ == "__main__":
    main()
