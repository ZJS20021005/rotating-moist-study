from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(r"E:\moist RB\rotating_case_inventory\04_outputs_and_figures\high_resolution_timeseries_latest_program_20260805")
OUT = ROOT / "force_balance" / "bulk_no_pressure"
OUT.mkdir(parents=True, exist_ok=True)

CELL_AND_BURST_SOURCE = ROOT / "force_balance" / "pressure_timeseries" / "Ek1p5e4_Ek2e4_force_balance_pressure_timeseries.csv"
INTERMITTENT_SOURCE = ROOT / "force_balance" / "intermittent_timeseries" / "force_balance_timeseries.csv"
NONINTERMITTENT_SOURCE = ROOT / "force_balance" / "nonintermittent_timeseries" / "force_balance_timeseries.csv"

FRAME = 4.5
LINE = 3.5
BURST_WINDOWS = {2.0e-4: [(353.3, 432.6), (622.9, 732.1), (939.4, 1029.8)]}
YELLOW = (1.0, 0.82, 0.12)

FORCE_STYLE = {
    "inertia_bulk": (r"$F_I$", (0.55, 0.55, 0.55), "*", 15),
    "coriolis_bulk": (r"$F_C$", (0.00, 0.42, 0.85), "o", 9),
    "viscous_bulk": (r"$F_V$", (0.00, 0.00, 0.00), "*", 10),
    "buoyancy_bulk": (r"$F_B$", (0.85, 0.00, 0.85), "s", 7),
}
SERIES_COLORS = {
    "inertia_bulk": (0.74, 0.14, 0.18),
    "coriolis_bulk": (0.11, 0.44, 0.71),
    "viscous_bulk": (0.96, 0.58, 0.19),
    "buoyancy_bulk": (0.00, 0.65, 0.12),
}
STABLE_EKS = [
    1.5e-4,
    1.0e-3,
    2.0e-3,
    3.0e-3,
    5.0e-3,
    7.0e-3,
    1.0e-2,
    3.0e-2,
    5.0e-2,
    1.0e-1,
]
INTERMITTENT_EKS = [2.0e-4, 5.0e-4, 7.0e-4]
INTERMITTENT_SPAN = (2.0e-4, 7.0e-4)


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


def style(ax):
    for spine in ax.spines.values():
        spine.set_linewidth(FRAME)
    ax.tick_params(which="major", direction="in", length=12, width=1.2, top=True, right=True)
    ax.minorticks_on()
    ax.tick_params(which="minor", direction="in", length=6, width=1.0, top=True, right=True)
    ax.set_box_aspect(5.2 / 6.5)


def ek_token(ek):
    return f"{ek:.2g}".replace(".", "p").replace("-", "m")


def ek_label(ek):
    if np.isclose(ek, 1.5e-4):
        return r"$Ek=1.5\times10^{-4}$"
    if np.isclose(ek, 2.0e-4):
        return r"$Ek=2\times10^{-4}$"
    return fr"$Ek={ek:.0e}$"


def load_bulk():
    frames = []
    for path in (CELL_AND_BURST_SOURCE, INTERMITTENT_SOURCE, NONINTERMITTENT_SOURCE):
        if path.exists():
            frames.append(pd.read_csv(path))
    data = pd.concat(frames, ignore_index=True)
    data = data.sort_values(["Ek", "time"])
    return data.drop_duplicates(subset=["Ek", "time", "frame"], keep="first")


def shade_bursts(ax, ek):
    for start, stop in BURST_WINDOWS.get(round(float(ek), 10), []):
        ax.axvspan(start, stop, color=YELLOW, alpha=0.20, lw=0, zorder=0)


def plot_case_timeseries(data, ek):
    case = data[np.isclose(data["Ek"], ek, rtol=0, atol=1e-10)].sort_values("time")
    if case.empty:
        return []
    fig = plt.figure(figsize=(13.0, 5.84), facecolor="white")
    left = fig.add_axes([0.08, 0.22, 0.38, 0.72])
    right = fig.add_axes([0.58, 0.22, 0.38, 0.72])
    style(left)
    style(right)
    shade_bursts(left, ek)
    shade_bursts(right, ek)
    time = case["time"].to_numpy(float)
    for field, (label, _, _, _) in FORCE_STYLE.items():
        value = case[field].to_numpy(float)
        valid = np.isfinite(value) & (value > 0)
        left.semilogy(time[valid], value[valid], color=SERIES_COLORS[field], lw=LINE, label=label)
    coriolis = np.maximum(case["coriolis_bulk"].to_numpy(float), 1e-30)
    for field, (label, _, _, _) in FORCE_STYLE.items():
        if field == "coriolis_bulk":
            continue
        ratio = case[field].to_numpy(float) / coriolis
        valid = np.isfinite(ratio) & (ratio > 0)
        right.semilogy(time[valid], ratio[valid], color=SERIES_COLORS[field], lw=LINE, label=label + r"/$F_C$")
    right.axhline(1.0, color="0.35", ls="--", lw=2.5)
    left.set_xlabel(r"$t$")
    right.set_xlabel(r"$t$")
    left.set_ylabel(r"bulk force magnitude")
    right.set_ylabel(r"bulk force ratio")
    left.legend(frameon=False, loc="best")
    right.legend(frameon=False, loc="best")
    left.text(0.95, 0.94, ek_label(ek), transform=left.transAxes, ha="right", va="top", fontsize=16)
    stem = OUT / f"Ra8e6_Ek{ek_token(ek)}_bulk_force_timeseries_no_pressure"
    fig.savefig(stem.with_suffix(".png"), facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
    plt.close(fig)
    return [stem.with_suffix(".png"), stem.with_suffix(".pdf")]


def late_mean(case):
    tmax = float(case["time"].max())
    window = case[case["time"] >= tmax - 500.0]
    if len(window) < 3:
        window = case.tail(min(5, len(case)))
    values = {
        "time_min": float(window["time"].min()),
        "time_max": float(window["time"].max()),
        "n_samples": int(len(window)),
    }
    for field in FORCE_STYLE:
        values[field] = float(window[field].mean())
    return values


def make_summary(data):
    rows = []
    for ek in STABLE_EKS:
        if any(np.isclose(ek, bad, rtol=0, atol=1e-10) for bad in INTERMITTENT_EKS):
            continue
        case = data[np.isclose(data["Ek"], ek, rtol=0, atol=1e-10)].sort_values("time")
        if case.empty:
            continue
        values = late_mean(case)
        values.update({
            "Ra": float(case["Ra"].iloc[0]),
            "Pr": float(case["Pr"].iloc[0]),
            "Ek": float(case["Ek"].iloc[0]),
            "AR": float(case["AR"].iloc[0]),
        })
        values["Ro_c"] = np.sqrt(values["Ra"] * values["Ek"] ** 2 / values["Pr"])
        qbot = 0.5
        qtop = 0.004978
        gamma = 1.1
        delta_m = gamma * (qbot - qtop)
        values["Ra_m_Ek_4_3"] = values["Ra"] * delta_m * values["Ek"] ** (4.0 / 3.0)
        values["supercritical_m_over_8p7"] = values["Ra_m_Ek_4_3"] / 8.7
        rows.append(values)
    summary = pd.DataFrame(rows).sort_values("Ro_c")
    summary.to_csv(OUT / "Ra8e6_stable_bulk_force_summary_no_pressure.csv", index=False)
    return summary


def plot_summary(summary, xfield, xlabel, stem_name):
    fig = plt.figure(figsize=(6.5, 5.84), facecolor="white")
    ax = fig.add_axes([0.186, 0.260, 0.792, 0.705])
    style(ax)
    if xfield == "Ek":
        ax.axvspan(INTERMITTENT_SPAN[0], INTERMITTENT_SPAN[1], color=YELLOW, alpha=0.14, lw=0, zorder=0)
        ax.axvline(INTERMITTENT_SPAN[0], color="0.35", ls="--", lw=1.6, zorder=1)
        ax.axvline(INTERMITTENT_SPAN[1], color="0.35", ls="--", lw=1.6, zorder=1)
    for field, (label, color, marker, size) in FORCE_STYLE.items():
        ax.plot(
            summary[xfield],
            summary[field],
            ls="None",
            marker=marker,
            ms=size,
            color=color,
            markerfacecolor=color,
            markeredgecolor="black",
            markeredgewidth=0.8,
            label=label,
            zorder=3,
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"Forces")
    if xfield == "Ek":
        ax.text(
            3.7e-4, 6.0e-2, "intermittent",
            ha="center", va="center", fontsize=13, color="0.25",
            rotation=0, zorder=2,
        )
    ax.legend(frameon=True, fancybox=False, edgecolor="black", loc="lower right", ncol=2)
    stem = OUT / stem_name
    fig.savefig(stem.with_suffix(".png"), facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
    plt.close(fig)
    return [stem.with_suffix(".png"), stem.with_suffix(".pdf")]


def plot_ro_l(summary):
    fig = plt.figure(figsize=(6.5, 5.84), facecolor="white")
    ax = fig.add_axes([0.186, 0.260, 0.792, 0.705])
    style(ax)
    ax.axvspan(INTERMITTENT_SPAN[0], INTERMITTENT_SPAN[1], color=YELLOW, alpha=0.14, lw=0, zorder=0)
    ax.axvline(INTERMITTENT_SPAN[0], color="0.35", ls="--", lw=1.6, zorder=1)
    ax.axvline(INTERMITTENT_SPAN[1], color="0.35", ls="--", lw=1.6, zorder=1)
    ro_l = summary["inertia_bulk"] / summary["coriolis_bulk"]
    ax.plot(summary["Ek"], ro_l, color=(0.74, 0.14, 0.18), lw=3.5, marker="o", ms=8, markeredgecolor="black", markerfacecolor=(0.74, 0.14, 0.18))
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$Ek$")
    ax.set_ylabel(r"$Ro_l^{(F)}=F_I/F_C$")
    ax.text(3.7e-4, 2.0e-1, "intermittent", ha="center", va="center", fontsize=13, color="0.25")
    stem = OUT / "Ra8e6_Rol_forcebalance_vs_Ek_no_pressure"
    fig.savefig(stem.with_suffix(".png"), facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
    plt.close(fig)
    return [stem.with_suffix(".png"), stem.with_suffix(".pdf")]


def main():
    configure()
    data = load_bulk()
    outputs = []
    for ek in [1.5e-4, 1.0e-3, 2.0e-3, 3.0e-3, 5.0e-3, 7.0e-3, 1.0e-2]:
        outputs.extend(plot_case_timeseries(data, ek))
    for ek in INTERMITTENT_EKS:
        outputs.extend(plot_case_timeseries(data, ek))
    summary = make_summary(data)
    outputs.extend(plot_summary(summary, "Ek", r"$Ek$", "Ra8e6_stable_bulk_forces_vs_Ek_no_pressure"))
    outputs.extend(plot_summary(summary, "Ro_c", r"$Ro_c$", "Ra8e6_stable_bulk_forces_vs_Roc_no_pressure"))
    outputs.extend(plot_summary(summary, "supercritical_m_over_8p7", r"$Ra_m Ek^{4/3}/8.7$", "Ra8e6_stable_bulk_forces_vs_supercriticality_no_pressure"))
    outputs.extend(plot_ro_l(summary))
    for path in outputs:
        print(path)


if __name__ == "__main__":
    main()
