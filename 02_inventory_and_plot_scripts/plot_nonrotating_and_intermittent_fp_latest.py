#!/usr/bin/env python3
"""Plot the stable nonrotating forces and intermittent pressure-force series."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(r"E:\moist RB\rotating_case_inventory")
FORCE_ROOT = (
    ROOT
    / "04_outputs_and_figures"
    / "high_resolution_timeseries_latest_program_20260805"
    / "force_balance"
)
OUTPUT = FORCE_ROOT / "nonrotating_and_intermittent_fp_20260806"
REDUCED = OUTPUT / "reduced_strict_fp3d"
NONROTATING = (
    ROOT
    / "04_outputs_and_figures"
    / "re_num_20260805"
    / "norotating_force_balance_bulk_summary.csv"
)
PRESSURE_2D = FORCE_ROOT / "online_pressure_force"

FRAME = 4.5
LINE = 3.5
BOX = 5.2 / 6.5
PRESSURE = (0.58, 0.05, 0.90)
COLORS = {
    "F_I": (0.74, 0.14, 0.18),
    "F_C": (0.11, 0.44, 0.71),
    "F_P": PRESSURE,
    "F_V": (0.96, 0.58, 0.19),
    "F_B": (0.00, 0.65, 0.12),
    "F_T": (0.15, 0.15, 0.15),
}
CASES = {
    2.0e-4: "Ek2e-4",
    5.0e-4: "Ek5e-4",
    7.0e-4: "Ek7e-4",
}


def configure() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Times New Roman",
            "mathtext.fontset": "stix",
            "axes.linewidth": FRAME,
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
        spine.set_linewidth(FRAME)
    axis.tick_params(
        which="major", direction="in", length=12, width=1.2, top=True, right=True
    )
    axis.minorticks_on()
    axis.tick_params(
        which="minor", direction="in", length=6, width=1.0, top=True, right=True
    )
    axis.set_box_aspect(BOX)


def single_axis() -> tuple[plt.Figure, plt.Axes]:
    figure = plt.figure(figsize=(6.5, 5.84), facecolor="white")
    axis = figure.add_axes([0.186, 0.260, 0.792, 0.705])
    style_axis(axis)
    return figure, axis


def three_axes() -> tuple[plt.Figure, list[plt.Axes]]:
    figure = plt.figure(figsize=(19.5, 5.84), facecolor="white")
    axes = [
        figure.add_axes([0.055, 0.22, 0.255, 0.72]),
        figure.add_axes([0.373, 0.22, 0.255, 0.72]),
        figure.add_axes([0.691, 0.22, 0.255, 0.72]),
    ]
    for axis in axes:
        style_axis(axis)
    return figure, axes


def save(figure: plt.Figure, stem: Path) -> None:
    figure.savefig(stem.with_suffix(".png"), facecolor="white")
    figure.savefig(stem.with_suffix(".pdf"), facecolor="white")
    plt.close(figure)


def plot_nonrotating() -> None:
    data = pd.read_csv(NONROTATING).sort_values("time")
    end = float(data["time"].max())
    window = data[data["time"] >= end - 200.0].copy()
    fields = [
        ("inertia_bulk", "F_I"),
        ("viscous_bulk", "F_V"),
        ("buoyancy_bulk", "F_B"),
    ]
    rows = []
    for field, quantity in fields:
        values = window[field].dropna().to_numpy(float)
        rows.append(
            {
                "quantity": quantity,
                "mean": float(np.mean(values)),
                "standard_deviation": float(np.std(values, ddof=1)),
                "time_start": float(window["time"].min()),
                "time_end": end,
                "n_samples": int(values.size),
            }
        )
    summary = pd.DataFrame(rows)
    summary.to_csv(
        OUTPUT / "Ra8e6_nonrotating_stable_force_magnitudes.csv",
        index=False,
        encoding="utf-8-sig",
    )

    figure, axis = single_axis()
    x = np.arange(len(summary), dtype=float)
    for index, row in summary.iterrows():
        axis.scatter(
            x[index],
            row["mean"],
            s=150,
            marker="o",
            facecolor=COLORS[row["quantity"]],
            edgecolor="black",
            linewidth=1.8,
            zorder=3,
        )
    axis.set_yscale("log")
    axis.set_xticks(x, [rf"${name}$" for name in summary["quantity"]])
    axis.set_xlim(-0.55, len(summary) - 0.45)
    axis.set_ylabel(r"$F_\alpha$")
    save(figure, OUTPUT / "Ra8e6_nonrotating_stable_force_magnitudes")


def read_pressure_2d(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep=r"\s+", names=["time", "F_P_h"], engine="python")
    frame = frame.replace([np.inf, -np.inf], np.nan).dropna().sort_values("time")
    return frame.drop_duplicates("time", keep="last")


def plot_pressure_2d() -> None:
    frames = []
    figure, axes = three_axes()
    for axis, (ek, token) in zip(axes, CASES.items()):
        frame = read_pressure_2d(PRESSURE_2D / f"{token}_pressure_force.out")
        frame.insert(0, "Ek", ek)
        frames.append(frame)
        axis.semilogy(frame["time"], frame["F_P_h"], color=PRESSURE, lw=LINE)
        axis.set_xlabel(r"$t$")
        axis.text(
            0.95,
            0.94,
            rf"$Ek={ek:.1g}$",
            transform=axis.transAxes,
            ha="right",
            va="top",
            fontsize=16,
        )
    axes[0].set_ylabel(r"$F_{P,h}$")
    pd.concat(frames, ignore_index=True).to_csv(
        OUTPUT / "Ra8e6_intermittent_Fp_horizontal_timeseries.csv",
        index=False,
        encoding="utf-8-sig",
    )
    save(figure, OUTPUT / "Ra8e6_intermittent_Fp_horizontal_timeseries")


def plot_strict_3d_if_available() -> int:
    available: list[tuple[float, str, pd.DataFrame]] = []
    for ek, token in CASES.items():
        path = REDUCED / token / "strict_force_balance_bulk_timeseries.csv"
        if path.exists():
            available.append((ek, token, pd.read_csv(path).sort_values("time")))
    if not available:
        return 0

    long_rows = []
    for ek, _, frame in available:
        for _, row in frame.iterrows():
            long_rows.append(
                {
                    "Ek": ek,
                    "time": row["time"],
                    "frame": row["frame"],
                    "F_P": row["F_P"],
                    "F_P_h": row["F_P_h"],
                    "F_P_z": row["F_P_z"],
                }
            )
    pd.DataFrame(long_rows).to_csv(
        OUTPUT / "Ra8e6_intermittent_Fp_strict_3d_timeseries.csv",
        index=False,
        encoding="utf-8-sig",
    )

    figure, axes = three_axes()
    axis_by_token = {token: axis for axis, token in zip(axes, CASES.values())}
    for ek, token, frame in available:
        axis = axis_by_token[token]
        axis.semilogy(frame["time"], frame["F_P"], color=PRESSURE, lw=LINE, label=r"$F_P$")
        axis.semilogy(
            frame["time"],
            frame["F_P_h"],
            color=PRESSURE,
            lw=LINE,
            ls="--",
            label=r"$F_{P,h}$",
        )
        axis.set_xlabel(r"$t$")
        axis.text(
            0.95, 0.94, rf"$Ek={ek:.1g}$", transform=axis.transAxes,
            ha="right", va="top", fontsize=16,
        )
    axes[0].set_ylabel(r"pressure-force magnitude")
    axes[0].legend(frameon=False, loc="best")
    save(figure, OUTPUT / "Ra8e6_intermittent_Fp_strict_3d_and_horizontal_timeseries")

    for ek, token, frame in available:
        figure, axis = single_axis()
        for quantity in ("F_I", "F_C", "F_P", "F_V", "F_B", "F_T"):
            values = frame[quantity].to_numpy(float)
            valid = np.isfinite(values) & (values > 0.0)
            axis.semilogy(
                frame.loc[valid, "time"],
                values[valid],
                color=COLORS[quantity],
                lw=LINE,
                label=rf"${quantity}$",
            )
        axis.set_xlabel(r"$t$")
        axis.set_ylabel(r"$F_\alpha$")
        axis.legend(frameon=False, loc="best")
        save(figure, OUTPUT / f"Ra8e6_{token}_strict_force_timeseries_with_Fp")
    return len(available)


def main() -> None:
    configure()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    REDUCED.mkdir(parents=True, exist_ok=True)
    plot_nonrotating()
    plot_pressure_2d()
    count = plot_strict_3d_if_available()
    print(f"output={OUTPUT}")
    print(f"strict_3d_cases={count}")


if __name__ == "__main__":
    main()
