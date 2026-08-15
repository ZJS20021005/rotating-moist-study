"""Plot Ra=8e6 bulk force magnitudes with strict 3-D pressure movies.

The pressure series is produced by compute_fp3d_from_pressure_movies.py and
includes the horizontal and vertical pressure gradients after horizontal
demeaning at each height.
"""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(r"E:\moist RB\rotating_case_inventory\04_outputs_and_figures\high_resolution_timeseries_latest_program_20260805")
FORCE_SUMMARY = ROOT / "force_balance" / "bulk_no_pressure" / "Ra8e6_stable_bulk_force_summary_no_pressure.csv"
PRESSURE_DIR = ROOT / "force_balance" / "fp3d_remote"
OUT = ROOT / "force_balance" / "bulk_with_fp3d"

FRAME = 4.5
LINE = 3.5
AX_RECT = [0.186, 0.260, 0.792, 0.705]
INTERMITTENT_SPAN = (2.0e-4, 7.0e-4)
YELLOW = (1.0, 0.82, 0.12)

SERIES = [
    ("inertia_bulk", r"$F_I$", (0.55, 0.55, 0.55), "*", 15),
    ("coriolis_bulk", r"$F_C$", (0.00, 0.42, 0.85), "o", 9),
    ("pressure_3d_bulk", r"$F_P$", (0.20, 0.20, 0.20), "^", 8),
    ("viscous_bulk", r"$F_V$", (0.00, 0.00, 0.00), "*", 10),
    ("buoyancy_bulk", r"$F_B$", (0.85, 0.00, 0.85), "s", 7),
]


def configure():
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


def style_axis(ax):
    for spine in ax.spines.values():
        spine.set_linewidth(FRAME)
    ax.tick_params(which="major", direction="in", length=12, width=1.2, top=True, right=True)
    ax.minorticks_on()
    ax.tick_params(which="minor", direction="in", length=6, width=1.0, top=True, right=True)
    ax.set_box_aspect(5.2 / 6.5)


def ek_from_name(path):
    token = path.parent.name
    if token.startswith("Ek"):
        token = token[2:]
    return float(token.replace("p", "."))


def strict_pressure_summary():
    rows = []
    for case_dir in sorted(PRESSURE_DIR.glob("Ek*")):
        path = case_dir / "fp3d_timeseries.csv"
        if not path.exists():
            continue
        values = pd.read_csv(path)
        values = values[["time", "F_P_3d_bulk"]].apply(pd.to_numeric, errors="coerce").dropna()
        values = values.apply(pd.to_numeric, errors="coerce").dropna()
        if values.empty:
            continue
        tmax = float(values["time"].max())
        late = values[values["time"] >= tmax - 500.0]
        time = late["time"].to_numpy(float)
        force = late["F_P_3d_bulk"].to_numpy(float)
        average = float(np.trapz(force, time) / (time[-1] - time[0])) if len(time) > 1 and time[-1] > time[0] else float(force.mean())
        rows.append({
            "Ek": ek_from_name(path),
            "pressure_3d_bulk": average,
            "pressure_time_min": float(late["time"].min()),
            "pressure_time_max": float(late["time"].max()),
            "pressure_n": int(len(late)),
            "pressure_source": str(path),
        })
    return pd.DataFrame(rows)


def make_axis():
    fig = plt.figure(figsize=(6.5, 5.84), facecolor="white")
    ax = fig.add_axes(AX_RECT)
    style_axis(ax)
    return fig, ax


def plot(summary, xfield, xlabel, stem):
    fig, ax = make_axis()
    if xfield == "Ek":
        ax.axvspan(*INTERMITTENT_SPAN, color=YELLOW, alpha=0.14, lw=0, zorder=0)
        ax.axvline(INTERMITTENT_SPAN[0], color="0.35", ls="--", lw=1.6, zorder=1)
        ax.axvline(INTERMITTENT_SPAN[1], color="0.35", ls="--", lw=1.6, zorder=1)
        ax.text(3.7e-4, 6.0e-2, "intermittent", ha="center", va="center", fontsize=13, color="0.25")
    for field, label, color, marker, size in SERIES:
        valid = np.isfinite(summary[field].to_numpy(float)) & (summary[field].to_numpy(float) > 0)
        ax.plot(
            summary.loc[valid, xfield], summary.loc[valid, field], ls="None", marker=marker,
            ms=size, color=color, markerfacecolor=color, markeredgecolor="black",
            markeredgewidth=0.8, label=label, zorder=3,
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"Forces")
    ax.legend(frameon=True, fancybox=False, edgecolor="black", loc="lower right", ncol=2)
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{stem}.png", facecolor="white")
    fig.savefig(OUT / f"{stem}.pdf", facecolor="white")
    plt.close(fig)


def main():
    configure()
    OUT.mkdir(parents=True, exist_ok=True)
    summary = pd.read_csv(FORCE_SUMMARY)
    pressure = strict_pressure_summary()
    summary["Ek_key"] = summary["Ek"].round(12)
    pressure["Ek_key"] = pressure["Ek"].round(12)
    pressure = pressure.drop(columns="Ek")
    summary = summary.merge(pressure, on="Ek_key", how="left").drop(columns="Ek_key")
    summary["fp3d_status"] = np.where(summary["pressure_3d_bulk"].notna(), "available", "missing_pressure_movie")
    summary.to_csv(OUT / "Ra8e6_stable_bulk_force_summary_with_fp3d.csv", index=False)
    plot_summary = summary[summary["pressure_3d_bulk"].notna()].copy()
    plot(plot_summary, "Ek", r"$Ek$", "Ra8e6_stable_bulk_forces_vs_Ek_with_fp3d")
    plot(plot_summary, "Ro_c", r"$Ro_c$", "Ra8e6_stable_bulk_forces_vs_Roc_with_fp3d")
    print(summary[["Ek", "pressure_3d_bulk", "pressure_time_min", "pressure_time_max", "pressure_n"]].to_string(index=False))
    print(f"saved to {OUT}")


if __name__ == "__main__":
    main()
