#!/usr/bin/env python3
"""Late-time force balance for the genuine Ra8e6 AR16 nonrotating case."""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(r"E:\moist RB\rotating_case_inventory")
SOURCE = ROOT / "04_outputs_and_figures" / "re_num_20260805" / "norotating_force_balance_bulk_summary.csv"
ENERGY_SOURCE = ROOT / "04_outputs_and_figures" / "high_resolution_timeseries_latest_program_20260805" / "high_resolution_timeseries_long.csv"
OUTPUT = ROOT / "04_outputs_and_figures" / "high_resolution_timeseries_latest_program_20260805" / "force_balance" / "norotating_late_balance_20260815"

FRAME = 4.5
LINE = 3.5
BOX = 5.2 / 6.5
AX_RECT = [0.186, 0.260, 0.792, 0.705]
LATE_START = 1200.0

FORCES = [
    ("inertia_bulk", r"$F_I$", (0.74, 0.14, 0.18)),
    ("viscous_bulk", r"$F_V$", (0.96, 0.58, 0.19)),
    ("buoyancy_bulk", r"$F_B$", (0.86, 0.05, 0.58)),
]


def set_style() -> None:
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
    ax.set_box_aspect(BOX)
    for spine in ax.spines.values():
        spine.set_linewidth(FRAME)
    ax.tick_params(which="major", direction="in", length=12, width=1.2, top=True, right=True)
    ax.minorticks_on()
    ax.tick_params(which="minor", direction="in", length=6, width=1.0, top=True, right=True)
    return fig, ax


def load_energy() -> pd.DataFrame:
    parts = []
    for chunk in pd.read_csv(ENERGY_SOURCE, chunksize=300_000):
        mask = (
            (chunk["Ra"] == 8.0e6)
            & chunk["Ek"].isna()
            & (chunk["AR"] == 16.0)
            & (chunk["beta"] == 1.02)
            & (chunk["qbot"] == 0.5)
            & (chunk["metric"] == "kinetic")
        )
        if mask.any():
            parts.append(chunk.loc[mask, ["time", "value"]])
    if not parts:
        raise RuntimeError("No genuine nonrotating kinetic-energy series was found")
    return pd.concat(parts).drop_duplicates("time", keep="last").sort_values("time")


def main() -> None:
    set_style()
    OUTPUT.mkdir(parents=True, exist_ok=True)

    data = pd.read_csv(SOURCE).drop_duplicates("time", keep="last").sort_values("time")
    late = data[data["time"] >= LATE_START].copy()
    energy = load_energy()
    energy_late = energy[energy["time"] >= LATE_START].copy()

    rows = []
    for field, label, _ in FORCES:
        values = late[field].dropna().to_numpy(float)
        rows.append(
            {
                "quantity": label.replace("$", ""),
                "source_column": field,
                "time_start": float(late["time"].min()),
                "time_end": float(late["time"].max()),
                "n_snapshots": int(len(values)),
                "mean": float(np.mean(values)),
                "standard_deviation": float(np.std(values, ddof=1)),
                "coefficient_of_variation": float(np.std(values, ddof=1) / np.mean(values)),
            }
        )

    t = energy_late["time"].to_numpy(float)
    e = energy_late["value"].to_numpy(float)
    slope = float(np.polyfit(t, e, 1)[0])
    pd.DataFrame(rows).to_csv(OUTPUT / "norotating_late_force_means.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(
        [
            {
                "time_start": float(t[0]),
                "time_end": float(t[-1]),
                "n_points": int(len(t)),
                "kinetic_energy_mean": float(np.mean(e)),
                "kinetic_energy_std": float(np.std(e, ddof=1)),
                "coefficient_of_variation": float(np.std(e, ddof=1) / np.mean(e)),
                "linear_slope": slope,
                "relative_change_from_linear_fit": float(slope * (t[-1] - t[0]) / np.mean(e)),
            }
        ]
    ).to_csv(OUTPUT / "norotating_late_energy_stationarity.csv", index=False, encoding="utf-8-sig")

    fig, ax = make_axis()
    ax.axvspan(LATE_START, float(data["time"].max()), color="0.88", alpha=0.55, lw=0, zorder=0)
    for field, label, color in FORCES:
        ax.semilogy(data["time"], data[field], color=color, lw=LINE, label=label, zorder=3)
    ax.set_xlim(float(data["time"].min()), float(data["time"].max()))
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"Force magnitude")
    ax.legend(frameon=False, loc="lower right", ncol=1)
    fig.savefig(OUTPUT / "norotating_force_balance_timeseries.png", facecolor="white")
    fig.savefig(OUTPUT / "norotating_force_balance_timeseries.pdf", facecolor="white")
    plt.close(fig)

    fig, ax = make_axis()
    for field, label, color in FORCES:
        ax.semilogy(late["time"], late[field], color=color, lw=LINE, label=label, zorder=3)
        ax.hlines(
            float(late[field].mean()),
            float(late["time"].min()),
            float(late["time"].max()),
            color=color,
            lw=2.0,
            linestyle="--",
            zorder=2,
        )
    ax.set_xlim(float(late["time"].min()), float(late["time"].max()))
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"Force magnitude")
    ax.legend(frameon=False, loc="lower right", ncol=1)
    fig.savefig(OUTPUT / "norotating_force_balance_late_timeseries.png", facecolor="white")
    fig.savefig(OUTPUT / "norotating_force_balance_late_timeseries.pdf", facecolor="white")
    plt.close(fig)

    fig, ax = make_axis()
    x = np.arange(len(rows))
    means = np.array([r["mean"] for r in rows])
    colors = [c for _, _, c in FORCES]
    ax.scatter(x, means, s=150, c=colors, edgecolors="black", linewidths=1.8, zorder=4)
    ax.set_yscale("log")
    ax.set_xlim(-0.65, len(rows) - 0.35)
    ax.set_xticks(x)
    ax.set_xticklabels([label for _, label, _ in FORCES])
    ax.set_ylabel(r"Late-time mean force")
    ax.set_xlabel(r"Force term")
    fig.savefig(OUTPUT / "norotating_force_balance_late_mean.png", facecolor="white")
    fig.savefig(OUTPUT / "norotating_force_balance_late_mean.pdf", facecolor="white")
    plt.close(fig)

    provenance = (
        "Genuine nonrotating source paths embedded in norotating_force_balance_bulk_summary.csv; "
        "Ra=8e6, Pr=0.7, AR=16, beta=1.02, qbot=0.5, N385x385x65. "
        "Late statistics use t=1200..1640. F_C=0. Full three-dimensional F_P is not included "
        "because the local genuine nonrotating snapshots do not contain pressure movies. "
        "G:/moist convection/.../conti_strict_force_500 is excluded because bou.in has invRo=1.47901994577 "
        "and CONTINUATION_INFO identifies Ek2e-4/AR4 as its source."
    )
    (OUTPUT / "README.txt").write_text(provenance, encoding="utf-8")

    print(OUTPUT)
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    main()
