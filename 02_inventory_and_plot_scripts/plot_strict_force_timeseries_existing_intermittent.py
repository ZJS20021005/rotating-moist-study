from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


INPUT = Path(
    r"G:\moist convection\Ek2e-4\AR4\Beta1p02\qbot0p5_qtop0p004978"
    r"\N257x257x65\conti1\conti_strict_force_500\run\diagnostics"
    r"\force_balance_strict_fp3d_20260806\strict_force_balance_bulk_timeseries.csv"
)
OUTPUT = Path(
    r"E:\moist RB\rotating_case_inventory\04_outputs_and_figures"
    r"\high_resolution_timeseries_latest_program_20260805\force_balance"
    r"\Ek2e-4_strict_3d_force_timeseries"
)

FRAME = 4.5
LINE = 3.5
FIGSIZE = (6.5, 5.84)
AX_RECT = [0.186, 0.260, 0.792, 0.705]

COLORS = {
    "I": (0.74, 0.14, 0.18),
    "C": (0.11, 0.44, 0.71),
    "P": (0.58, 0.05, 0.90),
    "V": (0.96, 0.58, 0.19),
    "B": (0.00, 0.65, 0.12),
    "T": (0.20, 0.20, 0.20),
}


def configure_style() -> None:
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


def style_axis(ax: plt.Axes) -> None:
    for spine in ax.spines.values():
        spine.set_linewidth(FRAME)
    ax.tick_params(which="major", direction="in", length=12, width=1.2, top=True, right=True)
    ax.minorticks_on()
    ax.tick_params(which="minor", direction="in", length=6, width=1.0, top=True, right=True)
    ax.set_box_aspect(5.2 / 6.5)


def make_axis() -> tuple[plt.Figure, plt.Axes]:
    fig = plt.figure(figsize=FIGSIZE, facecolor="white")
    ax = fig.add_axes(AX_RECT)
    style_axis(ax)
    return fig, ax


def plot_series(ax: plt.Axes, data: pd.DataFrame, specs: list[tuple[str, str, str]]) -> None:
    time = data["time"].to_numpy(float)
    for field, label, key in specs:
        values = data[field].to_numpy(float)
        valid = np.isfinite(time) & np.isfinite(values) & (values > 0.0)
        ax.semilogy(time[valid], values[valid], color=COLORS[key], lw=LINE, label=label)
    ax.set_xlim(float(np.nanmin(time)), float(np.nanmax(time)))
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"Force magnitude")
    ax.legend(frameon=False, ncol=2, loc="best", handlelength=2.6, columnspacing=1.2)


def save_single(data: pd.DataFrame, specs: list[tuple[str, str, str]], stem: str) -> None:
    fig, ax = make_axis()
    plot_series(ax, data, specs)
    fig.savefig(OUTPUT / f"{stem}.png", facecolor="white")
    fig.savefig(OUTPUT / f"{stem}.pdf", facecolor="white")
    plt.close(fig)


def save_overview(
    data: pd.DataFrame,
    groups: list[tuple[str, list[tuple[str, str, str]]]],
) -> None:
    fig = plt.figure(figsize=(19.5, 5.84), facecolor="white")
    axes = [
        fig.add_axes([0.055, 0.260, 0.254, 0.705]),
        fig.add_axes([0.376, 0.260, 0.254, 0.705]),
        fig.add_axes([0.697, 0.260, 0.254, 0.705]),
    ]
    for ax, (_, specs) in zip(axes, groups):
        style_axis(ax)
        plot_series(ax, data, specs)
    fig.savefig(OUTPUT / "Ra8e6_Ek2e-4_strict_force_timeseries_overview.png", facecolor="white")
    fig.savefig(OUTPUT / "Ra8e6_Ek2e-4_strict_force_timeseries_overview.pdf", facecolor="white")
    plt.close(fig)


def main() -> None:
    configure_style()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(INPUT).sort_values("time").drop_duplicates("time", keep="last")

    total = [
        ("F_I", r"$F_I$", "I"),
        ("F_C", r"$F_C$", "C"),
        ("F_P", r"$F_P$", "P"),
        ("F_V", r"$F_V$", "V"),
        ("F_B", r"$F_B$", "B"),
        ("F_T", r"$F_T$", "T"),
    ]
    horizontal = [
        ("F_I_h", r"$F_{I,h}$", "I"),
        ("F_C_h", r"$F_{C,h}$", "C"),
        ("F_P_h", r"$F_{P,h}$", "P"),
        ("F_V_h", r"$F_{V,h}$", "V"),
        ("F_T_h", r"$F_{T,h}$", "T"),
    ]
    vertical = [
        ("F_I_z", r"$F_{I,z}$", "I"),
        ("F_P_z", r"$F_{P,z}$", "P"),
        ("F_V_z", r"$F_{V,z}$", "V"),
        ("F_B_z", r"$F_{B,z}$", "B"),
        ("F_T_z", r"$F_{T,z}$", "T"),
    ]

    groups = [("total", total), ("horizontal", horizontal), ("vertical", vertical)]
    for name, specs in groups:
        save_single(data, specs, f"Ra8e6_Ek2e-4_strict_force_timeseries_{name}")
    save_overview(data, groups)

    selected = ["time"] + [item[0] for _, specs in groups for item in specs]
    selected = list(dict.fromkeys(selected))
    data[selected].to_csv(OUTPUT / "Ra8e6_Ek2e-4_strict_force_timeseries_plot_data.csv", index=False)


if __name__ == "__main__":
    main()
