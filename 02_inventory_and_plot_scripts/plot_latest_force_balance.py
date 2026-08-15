from pathlib import Path
import math

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(r"E:\moist RB\rotating_case_inventory\04_outputs_and_figures\high_resolution_timeseries_latest_program_20260805\force_balance")
FRAME, LINE = 4.5, 3.5


def configure():
    mpl.rcParams.update({"font.family": "Times New Roman", "mathtext.fontset": "stix",
        "axes.linewidth": FRAME, "axes.labelsize": 24, "xtick.labelsize": 13,
        "ytick.labelsize": 13, "legend.fontsize": 10, "figure.dpi": 180,
        "savefig.dpi": 300, "pdf.fonttype": 42})


def axis(xlabel, ylabel):
    fig = plt.figure(figsize=(6.5, 5.84), facecolor="white")
    ax = fig.add_axes([0.186, 0.260, 0.792, 0.705])
    for spine in ax.spines.values(): spine.set_linewidth(FRAME)
    ax.tick_params(which="major", direction="in", length=12, width=1.2, top=True, right=True)
    ax.minorticks_on(); ax.tick_params(which="minor", direction="in", length=6, width=1.0, top=True, right=True)
    ax.set_box_aspect(5.2 / 6.5); ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    return fig, ax


def save(fig, name):
    fig.savefig(ROOT / (name + ".png"), facecolor="white")
    fig.savefig(ROOT / (name + ".pdf"), facecolor="white")
    plt.close(fig)


def regime(row):
    if not np.isfinite(row.Ek):
        return "Nonrotating: I-B-V", "high"
    ro = row.local_Ro_force
    if ro >= 3:
        return "Inertial / weakly rotation-affected", "medium"
    if ro >= 1:
        return "Rotation-affected transitional", "medium"
    # Rotational constraint is established, but pressure was not saved.
    if row.inertia_bulk >= 0.5 * row.buoyancy_bulk and row.inertia_bulk >= 2 * row.viscous_bulk:
        return "Geostrophic-compatible; CIA-like residual", "medium"
    if row.viscous_bulk >= row.inertia_bulk:
        return "Geostrophic-compatible; VAC/B mixed residual", "medium"
    return "Geostrophic-compatible; mixed residual", "medium"


def main():
    configure()
    bulk = pd.read_csv(ROOT / "force_balance_bulk_summary.csv")
    profiles = pd.read_csv(ROOT / "force_balance_profiles.csv")
    finite = bulk[np.isfinite(bulk.Ek)].sort_values("Ek")

    fig, ax = axis(r"$Ek$", r"bulk force ratio")
    for column, label, color, marker in [
        ("local_Ro_force", r"$F_I/F_C$", (0.74, 0.14, 0.18), "o"),
        ("viscous_bulk", r"$F_V/F_C$", (0.96, 0.58, 0.19), "s"),
        ("buoyancy_bulk", r"$F_B/F_C$", (0.11, 0.44, 0.71), "D")]:
        value = finite[column].to_numpy(float)
        if column != "local_Ro_force": value = value / finite.coriolis_bulk.to_numpy(float)
        ax.loglog(finite.Ek, value, "--", color=color, lw=LINE, marker=marker,
                  ms=8, mec="black", mew=1.2, mfc=color, label=label)
    ax.axhline(1, color="0.35", ls=":", lw=2.5)
    ax.legend(frameon=False, loc="upper left")
    save(fig, "Ra8e6_bulk_force_ratios_vs_Ek")

    selected = finite.iloc[np.linspace(0, len(finite)-1, 6).round().astype(int)]
    fig, ax = axis(r"$F/F_C$", r"$z$")
    colors = plt.cm.turbo(np.linspace(0.06, 0.94, len(selected)))
    for (_, row), color in zip(selected.iterrows(), colors):
        frame = profiles[np.isclose(profiles.Ek, row.Ek)].sort_values("z")
        ratio = frame.inertia / np.maximum(frame.coriolis, 1e-30)
        ax.semilogx(ratio, frame.z, color=color, lw=LINE, label=fr"${row.Ek:.0e}$")
    ax.axvline(1, color="0.35", ls=":", lw=2.5)
    ax.set_ylim(0, 1); ax.legend(title=r"$Ek$", frameon=False, loc="best")
    save(fig, "Ra8e6_force_Rossby_profiles_selected_Ek")

    classifications = []
    for _, row in bulk.sort_values("Ek", na_position="last").iterrows():
        label, confidence = regime(row)
        classifications.append({"Ek": row.Ek, "frame": row.frame,
            "F_I/F_C": row.local_Ro_force, "F_V/F_C": row.viscous_bulk / row.coriolis_bulk if row.coriolis_bulk else np.nan,
            "F_B/F_C": row.buoyancy_bulk / row.coriolis_bulk if row.coriolis_bulk else np.nan,
            "buoyancy_power": row.buoyancy_power_bulk, "viscous_dissipation": row.viscous_dissipation_bulk,
            "classification": label, "confidence": confidence,
            "pressure_caveat": "pressure not saved; leading C-P cancellation not directly tested" if np.isfinite(row.Ek) else "not applicable"})
    pd.DataFrame(classifications).to_csv(ROOT / "force_balance_case_classification.csv", index=False, encoding="utf-8-sig")


if __name__ == "__main__": main()
