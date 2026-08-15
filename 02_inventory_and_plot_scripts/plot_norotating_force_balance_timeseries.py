#!/usr/bin/env python3
"""Plot the nonrotating force-balance baseline without pressure."""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(r"E:\moist RB\rotating_case_inventory")
INPUT = ROOT / "04_outputs_and_figures" / "re_num_20260805" / "norotating_force_balance_bulk_summary.csv"
OUTPUT = ROOT / "04_outputs_and_figures" / "re_num_20260805" / "norotating_force_balance_plots"

FRAME = 4.5
LINE = 3.5
BOX = 5.2 / 6.5
AX_RECT = [0.186, 0.260, 0.792, 0.705]


def style() -> None:
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


def make_axis():
    fig = plt.figure(figsize=(6.5, 5.84), facecolor="white")
    ax = fig.add_axes(AX_RECT)
    for spine in ax.spines.values():
        spine.set_linewidth(FRAME)
    ax.tick_params(which="major", direction="in", length=12, width=1.2, top=True, right=True)
    ax.minorticks_on()
    ax.tick_params(which="minor", direction="in", length=6, width=1.0, top=True, right=True)
    ax.set_box_aspect(BOX)
    return fig, ax


def stable_window(frame: pd.DataFrame, window: float = 200.0) -> pd.DataFrame:
    end = float(frame["time"].max())
    start = max(float(frame["time"].min()), end - window)
    return frame[frame["time"] >= start].copy()


def main() -> None:
    style()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(INPUT)
    data = data.sort_values("time")
    if data.empty:
        raise SystemExit(f"no data in {INPUT}")

    components = [
        ("inertia_bulk", r"$F_I$", (0.74, 0.14, 0.18)),
        ("viscous_bulk", r"$F_V$", (0.96, 0.58, 0.19)),
        ("buoyancy_bulk", r"$F_B$", (0.11, 0.44, 0.71)),
    ]
    window = stable_window(data, 200.0)
    means = {field: float(window[field].mean()) for field, _, _ in components}

    fig, ax = make_axis()
    for field, label, color in components:
        ax.semilogy(data["time"], data[field], color=color, lw=LINE, label=label)
        ax.hlines(
            max(means[field], 1e-30),
            xmin=float(window["time"].min()),
            xmax=float(window["time"].max()),
            colors="0.35",
            linestyles="--",
            linewidth=2.2,
            zorder=1,
    )
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"bulk force magnitude")
    ax.set_xlim(left=0.0)
    ax.legend(frameon=False, loc="best")
    fig.savefig(OUTPUT / "norotating_force_balance_timeseries.png", facecolor="white")
    fig.savefig(OUTPUT / "norotating_force_balance_timeseries.pdf", facecolor="white")
    plt.close(fig)

    summary = pd.DataFrame([{
        "time_start": float(window["time"].min()),
        "time_end": float(window["time"].max()),
        "n_points": int(len(window)),
        **{f"{field}_mean": means[field] for field, _, _ in components},
    }])
    summary.to_csv(OUTPUT / "norotating_force_balance_stable_window.csv", index=False, encoding="utf-8-sig")
    print(f"saved to {OUTPUT}")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
