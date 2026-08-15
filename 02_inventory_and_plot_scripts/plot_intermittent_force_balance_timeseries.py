from pathlib import Path
import argparse

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

DEFAULT_ROOT = Path(r"E:\moist RB\rotating_case_inventory\04_outputs_and_figures\high_resolution_timeseries_latest_program_20260805\force_balance\intermittent_timeseries")
FRAME, LINE = 4.5, 3.5
COLORS = {
    "inertia_bulk": (0.74, 0.14, 0.18),
    "coriolis_bulk": (0.11, 0.44, 0.71),
    "viscous_bulk": (0.96, 0.58, 0.19),
    "buoyancy_bulk": (0.00, 0.65, 0.12),
}
PRESSURE_COLOR = (0.58, 0.05, 0.90)
BURST_WINDOWS = {
    2.0e-4: [(353.3, 432.6), (622.9, 732.1), (939.4, 1029.8)],
    5.0e-4: [(390.5, 429.0), (451.8, 488.5)],
    7.0e-4: [(447.9, 469.6), (491.1, 523.1), (580.6, 684.6)],
}
LABELS = {"inertia_bulk": r"$F_I$", "coriolis_bulk": r"$F_C$",
          "viscous_bulk": r"$F_V$", "buoyancy_bulk": r"$F_B$"}


def configure():
    mpl.rcParams.update({"font.family": "Times New Roman", "mathtext.fontset": "stix",
        "axes.linewidth": FRAME, "axes.labelsize": 24, "xtick.labelsize": 13,
        "ytick.labelsize": 13, "legend.fontsize": 11, "figure.dpi": 180,
        "savefig.dpi": 300, "pdf.fonttype": 42})


def style(ax):
    for spine in ax.spines.values(): spine.set_linewidth(FRAME)
    ax.tick_params(which="major", direction="in", length=12, width=1.2, top=True, right=True)
    ax.minorticks_on(); ax.tick_params(which="minor", direction="in", length=6, width=1.0, top=True, right=True)
    ax.set_box_aspect(5.2 / 6.5)


def token(ek): return f"{ek:.2g}".replace("-", "m").replace(".", "p")


def shade_bursts(axis, ek):
    for target, windows in BURST_WINDOWS.items():
        if np.isclose(ek, target, rtol=1e-8, atol=1e-12):
            for start, stop in windows:
                axis.axvspan(start, stop, color=(1.0, 0.82, 0.12), alpha=0.20, lw=0, zorder=0)
            break


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default=str(DEFAULT_ROOT / "force_balance_timeseries.csv"))
    parser.add_argument("--output", default=str(DEFAULT_ROOT))
    args = parser.parse_args()
    output = Path(args.output); output.mkdir(parents=True, exist_ok=True)
    configure(); data = pd.read_csv(args.input).sort_values(["Ek", "time"], na_position="last")
    for ek, case in data.groupby("Ek", dropna=False):
        if not np.isfinite(ek):
            fig = plt.figure(figsize=(6.5, 5.84), facecolor="white")
            left = fig.add_axes([0.186, 0.260, 0.792, 0.705]); style(left)
            for field in ("inertia_bulk", "viscous_bulk", "buoyancy_bulk"):
                value = case[field].to_numpy(float); valid = value > 0
                left.semilogy(case.time.to_numpy(float)[valid], value[valid], color=COLORS[field], lw=LINE, label=LABELS[field])
            left.set_xlabel(r"$t$"); left.set_ylabel(r"bulk force magnitude")
            left.legend(frameon=False, loc="best")
            stem = output / "Ra8e6_nonrotating_force_balance_timeseries"
            fig.savefig(stem.with_suffix(".png"), facecolor="white"); fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
            plt.close(fig); continue
        has_pressure = "pressure_plane" in case and case.pressure_plane.notna().any()
        if has_pressure:
            fig = plt.figure(figsize=(19.5, 5.84), facecolor="white")
            left = fig.add_axes([0.045, 0.22, 0.255, 0.72])
            middle = fig.add_axes([0.372, 0.22, 0.255, 0.72])
            right = fig.add_axes([0.699, 0.22, 0.255, 0.72])
            style(left); style(middle); style(right)
        else:
            fig = plt.figure(figsize=(13.0, 5.84), facecolor="white")
            left = fig.add_axes([0.09, 0.22, 0.38, 0.72]); right = fig.add_axes([0.58, 0.22, 0.38, 0.72])
            style(left); style(right); middle = None
        for field in COLORS:
            value = case[field].to_numpy(float)
            valid = value > 0
            left.semilogy(case.time.to_numpy(float)[valid], value[valid], color=COLORS[field], lw=LINE, label=LABELS[field])
        coriolis = (case.coriolis_plane if has_pressure else case.coriolis_bulk).to_numpy(float)
        valid_c = coriolis > 1e-15
        ratio_fields = (
            (("inertia_plane", r"$F_{I,h}/F_{C,h}$", COLORS["inertia_bulk"]),
             ("viscous_plane", r"$F_{V,h}/F_{C,h}$", COLORS["viscous_bulk"]))
            if has_pressure else
            (("inertia_bulk", r"$F_I/F_C$", COLORS["inertia_bulk"]),
             ("viscous_bulk", r"$F_V/F_C$", COLORS["viscous_bulk"]),
             ("buoyancy_bulk", r"$F_B/F_C$", COLORS["buoyancy_bulk"]))
        )
        for field, ratio_label, ratio_color in ratio_fields:
            ratio = case[field].to_numpy(float) / np.maximum(coriolis, 1e-30)
            valid = valid_c & np.isfinite(ratio) & (ratio > 0)
            right.semilogy(case.time.to_numpy(float)[valid], ratio[valid], color=ratio_color, lw=LINE,
                           label=ratio_label)
        right.axhline(1, color="0.35", ls="--", lw=2.5)
        if has_pressure:
            plane_fields = [
                ("inertia_plane", r"$F_{I,h}$", COLORS["inertia_bulk"]),
                ("coriolis_plane", r"$F_{C,h}$", COLORS["coriolis_bulk"]),
                ("viscous_plane", r"$F_{V,h}$", COLORS["viscous_bulk"]),
                ("pressure_plane", r"$F_{P,h}$", PRESSURE_COLOR),
            ]
            for field, label, color in plane_fields:
                value = case[field].to_numpy(float); valid = np.isfinite(value) & (value > 0)
                middle.semilogy(case.time.to_numpy(float)[valid], value[valid], color=color, lw=LINE, label=label)
            fp_fc = case.pressure_plane.to_numpy(float) / np.maximum(case.coriolis_plane.to_numpy(float), 1e-30)
            valid = np.isfinite(fp_fc) & (fp_fc > 0)
            right.semilogy(case.time.to_numpy(float)[valid], fp_fc[valid], color=PRESSURE_COLOR, lw=LINE, label=r"$F_{P,h}/F_{C,h}$")
            residual = case.cp_residual_over_coriolis.to_numpy(float)
            valid = np.isfinite(residual) & (residual > 0)
            right.semilogy(case.time.to_numpy(float)[valid], residual[valid], color="black", lw=LINE, ls=":", label=r"$|\mathbf{F}_{P,h}+\mathbf{F}_{C,h}|/F_{C,h}$")
            middle.set_xlabel(r"$t$"); middle.set_ylabel(r"horizontal force at $z\simeq0.048$")
            middle.legend(frameon=False, loc="best")
        left.set_xlabel(r"$t$"); right.set_xlabel(r"$t$")
        left.set_ylabel(r"bulk force magnitude")
        right.set_ylabel(r"horizontal force ratio" if has_pressure else r"bulk force ratio")
        left.legend(frameon=False, loc="best"); right.legend(frameon=False, loc="best")
        shade_bursts(left, ek); shade_bursts(right, ek)
        if middle is not None: shade_bursts(middle, ek)
        left.text(0.96, 0.95, fr"$Ek={ek:.2g}$", transform=left.transAxes, ha="right", va="top", fontsize=16)
        stem = output / f"Ra8e6_Ek{token(ek)}_force_balance_timeseries"
        fig.savefig(stem.with_suffix(".png"), facecolor="white")
        fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
        plt.close(fig)


if __name__ == "__main__": main()
