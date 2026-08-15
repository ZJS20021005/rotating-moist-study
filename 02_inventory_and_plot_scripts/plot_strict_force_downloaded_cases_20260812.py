from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(r"E:\moist RB\rotating_case_inventory")
SOURCE = ROOT / "04_outputs_and_figures" / "high_resolution_timeseries_latest_program_20260805"
OUT = SOURCE / "force_num_vs_Ek_20260812"
DOWNLOADED = OUT / "downloaded_force_balance"

EK2E3 = DOWNLOADED / "Ek2e-3" / "force_balance.out"
EK2E3_Z = DOWNLOADED / "Ek2e-3" / "force_balance_z.out"
EK2E4 = Path(
    r"G:\moist convection\Ek2e-4\AR4\Beta1p02\qbot0p5_qtop0p004978"
    r"\N257x257x65\conti1\conti_strict_force_500\run\diagnostics"
    r"\force_balance_strict_fp3d_20260806\strict_force_balance_bulk_timeseries.csv"
)
EK2E4_Z = EK2E4.with_name("strict_force_balance_profiles_long.csv")

FRAME = 4.5
LINE = 3.5
FIGSIZE = (6.5, 5.84)
AX_RECT = [0.186, 0.260, 0.792, 0.705]

COLORS = {
    "F_I": (0.74, 0.14, 0.18),
    "F_C": (0.11, 0.44, 0.71),
    "F_P": (0.58, 0.05, 0.90),
    "F_V": (0.96, 0.58, 0.19),
    "F_B": (0.88, 0.05, 0.65),
    "F_CP": (0.00, 0.65, 0.12),
}

LINE_STYLES = {
    "F_I": "-",
    "F_C": (0, (8, 3)),
    "F_P": "-",
    "F_V": (0, (4, 2)),
    "F_B": "-.",
    "F_CP": "-",
}

LABELS = {
    "F_I": r"$F_I$",
    "F_C": r"$F_C$",
    "F_P": r"$F_P$",
    "F_V": r"$F_V$",
    "F_B": r"$F_B$",
    "F_CP": r"$|\mathbf{F}_{C}+\mathbf{F}_{P}|$",
}

BURST_WINDOWS = {
    2.0e-4: [(1219.4, 1321.8), (1523.1, 1597.9)],
}


def configure_style() -> None:
    mpl.rcParams.update({
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
    })


def make_axis() -> tuple[plt.Figure, plt.Axes]:
    fig = plt.figure(figsize=FIGSIZE, facecolor="white")
    ax = fig.add_axes(AX_RECT)
    for spine in ax.spines.values():
        spine.set_linewidth(FRAME)
    ax.tick_params(which="major", direction="in", length=12, width=1.2, top=True, right=True)
    ax.minorticks_on()
    ax.tick_params(which="minor", direction="in", length=6, width=1.0, top=True, right=True)
    ax.set_box_aspect(5.2 / 6.5)
    return fig, ax


def load_online(path: Path) -> pd.DataFrame:
    names = [
        "time",
        "F_I", "F_C", "F_P", "F_V", "F_B", "F_T",
        "F_I_h", "F_C_h", "F_P_h", "F_V_h", "F_B_h", "F_T_h",
        "F_I_z", "F_C_z", "F_P_z", "F_V_z", "F_B_z", "F_T_z",
        "R_G", "R_z", "R_momentum",
        "F_I_x", "F_C_x", "F_P_x", "F_V_x", "F_B_x", "F_T_x",
        "F_I_y", "F_C_y", "F_P_y", "F_V_y", "F_B_y", "F_T_y",
    ]
    return pd.read_csv(path, sep=r"\s+", header=None, names=names).sort_values("time").drop_duplicates("time", keep="last")


def add_total_coriolis_pressure_from_online_z(
    bulk: pd.DataFrame, profile_path: Path
) -> pd.DataFrame:
    names = [
        "time", "z",
        "F_I", "F_C", "F_P", "F_V", "F_B", "F_T",
        "F_I_h", "F_C_h", "F_P_h", "F_V_h", "F_B_h", "F_T_h",
        "F_I_z", "F_C_z", "F_P_z", "F_V_z", "F_B_z", "F_T_z",
        "R_G", "R_z", "R_momentum",
        "F_I_x", "F_C_x", "F_P_x", "F_V_x", "F_B_x", "F_T_x",
        "F_I_y", "F_C_y", "F_P_y", "F_V_y", "F_B_y", "F_T_y",
    ]
    profile = pd.read_csv(profile_path, sep=r"\s+", header=None, names=names)
    for field in ("time", "z", "R_G", "F_P_z"):
        profile[field] = pd.to_numeric(profile[field], errors="coerce")
    profile = profile[
        profile["time"].notna()
        & profile["z"].between(0.1, 0.9)
        & profile["R_G"].notna()
        & profile["F_P_z"].notna()
    ].copy()
    profile["F_CP"] = np.sqrt(profile["R_G"] ** 2 + profile["F_P_z"] ** 2)
    cp = profile.groupby("time", as_index=False)["F_CP"].mean()
    return bulk.merge(cp, on="time", how="left")


def load_offline(path: Path, profile_path: Path) -> pd.DataFrame:
    bulk = pd.read_csv(path).sort_values("time").drop_duplicates("time", keep="last")
    profile = pd.read_csv(profile_path)
    for field in ("time", "z", "R_G", "F_P_z"):
        profile[field] = pd.to_numeric(profile[field], errors="coerce")
    profile = profile[
        profile["time"].notna()
        & profile["z"].between(0.1, 0.9)
        & profile["R_G"].notna()
        & profile["F_P_z"].notna()
    ].copy()
    # F_C,z=0, so the complete three-dimensional vector residual is the
    # plane RMS sqrt(R_G^2 + F_P,z^2), followed by the configured z average.
    profile["F_CP"] = np.sqrt(profile["R_G"] ** 2 + profile["F_P_z"] ** 2)
    cp = profile.groupby("time", as_index=False)["F_CP"].mean()
    return bulk.merge(cp, on="time", how="left")


def plot_case(data: pd.DataFrame, case: str, ek: float) -> None:
    fig, ax = make_axis()
    time = data["time"].to_numpy(float)
    for field, color in COLORS.items():
        values = pd.to_numeric(data[field], errors="coerce").to_numpy(float)
        valid = np.isfinite(time) & np.isfinite(values) & (values > 0.0)
        ax.semilogy(
            time[valid],
            values[valid],
            color=color,
            lw=LINE + 0.3 if field in ("F_P", "F_CP") else LINE,
            ls=LINE_STYLES[field],
            solid_capstyle="round",
            dash_capstyle="round",
            label=LABELS[field],
            zorder=4 if field == "F_CP" else 3,
        )
    for start, stop in BURST_WINDOWS.get(ek, []):
        ax.axvspan(start, stop, color=(1.00, 0.72, 0.00), alpha=0.11, lw=0, zorder=0)
    ax.set_xlim(float(np.nanmin(time)), float(np.nanmax(time)))
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"Force magnitude")
    ax.legend(frameon=False, ncol=2, loc="best", handlelength=3.4, columnspacing=1.2)
    fig.savefig(OUT / f"Ra8e6_{case}_strict_force_timeseries.png", facecolor="white")
    fig.savefig(OUT / f"Ra8e6_{case}_strict_force_timeseries.pdf", facecolor="white")
    plt.close(fig)

    keep = ["time", *COLORS]
    data[keep].to_csv(OUT / f"Ra8e6_{case}_strict_force_timeseries_plot_data.csv", index=False)


def main() -> None:
    configure_style()
    OUT.mkdir(parents=True, exist_ok=True)
    if EK2E3.exists() and EK2E3_Z.exists():
        plot_case(
            add_total_coriolis_pressure_from_online_z(load_online(EK2E3), EK2E3_Z),
            "Ek2e-3",
            2.0e-3,
        )
    if EK2E4.exists() and EK2E4_Z.exists():
        plot_case(load_offline(EK2E4, EK2E4_Z), "Ek2e-4", 2.0e-4)


if __name__ == "__main__":
    main()
