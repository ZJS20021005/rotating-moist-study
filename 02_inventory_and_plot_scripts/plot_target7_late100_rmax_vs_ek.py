from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


INPUT = Path(
    r"E:\moist RB\rotating_case_inventory\04_outputs_and_figures"
    r"\target7_latest_program_timeseries_20260802"
    r"\maximum_mse_vortex_radius_timeseries.csv"
)
OUTDIR = INPUT.parent
EXCLUDED_EK = 2.0e-4
LATE_WINDOW = 100.0
LIMIT_EK = 1.5e-4
LIMIT_R = 0.16  # 2*vortex_lc in this case; zero output means no component passed it.


def main() -> None:
    data = pd.read_csv(INPUT)
    rows: list[dict[str, float | int | str]] = []

    for ek, group in data.groupby("Ek", sort=True):
        ek = float(ek)
        if np.isclose(ek, EXCLUDED_EK, rtol=0.0, atol=1e-12):
            continue
        group = group.sort_values("time")
        t_end = float(group["time"].max())
        t_start = t_end - LATE_WINDOW
        window = group[group["time"] >= t_start]
        rows.append(
            {
                "Ek": ek,
                "AR": int(window["AR"].iloc[-1]),
                "Nx": int(window["Nx"].iloc[-1]),
                "Ny": int(window["Ny"].iloc[-1]),
                "Nz": int(window["Nz"].iloc[-1]),
                "t_start": t_start,
                "t_end": t_end,
                "n_samples": int(len(window)),
                "mean_Rmax_positive": float(window["Rmax_positive"].mean()),
                "mean_Rmax_core": float(window["Rmax_core"].mean()),
                "status": (
                    "upper_limit_below_0.16"
                    if np.isclose(ek, LIMIT_EK, rtol=0.0, atol=1e-12)
                    else "measured"
                ),
            }
        )

    summary = pd.DataFrame(rows).sort_values("Ek")
    csv_path = OUTDIR / "late100_mean_maximum_mse_vortex_radius_vs_Ek.csv"
    summary.to_csv(csv_path, index=False)

    mpl.rcParams.update(
        {
            "font.family": "Times New Roman",
            "mathtext.fontset": "stix",
            "axes.linewidth": 4.5,
            "axes.labelsize": 24,
            "xtick.labelsize": 13,
            "ytick.labelsize": 13,
            "legend.fontsize": 10,
            "figure.dpi": 180,
            "savefig.dpi": 300,
        }
    )

    fig = plt.figure(figsize=(6.5, 5.84), facecolor="white")
    ax = fig.add_axes([0.186, 0.260, 0.792, 0.705])
    for spine in ax.spines.values():
        spine.set_linewidth(4.5)
    ax.tick_params(direction="in", length=12, width=1.2, top=True, right=True)
    ax.minorticks_on()
    ax.tick_params(
        which="minor", direction="in", length=6, width=1.0, top=True, right=True
    )
    ax.set_box_aspect(5.2 / 6.5)
    ax.set_xscale("log")
    ax.set_yscale("log")

    measured = summary[summary["status"] == "measured"]
    blue = (0.00, 0.00, 1.00)
    ax.plot(
        measured["Ek"],
        measured["mean_Rmax_positive"],
        color=blue,
        lw=3.5,
        ls="-",
        marker="o",
        ms=10,
        mfc=blue,
        mec=blue,
        label=r"$m'_{2D}>0$",
        zorder=3,
    )
    ax.plot(
        measured["Ek"],
        measured["mean_Rmax_core"],
        color=blue,
        lw=3.5,
        ls="--",
        marker="s",
        ms=9,
        mfc=blue,
        mec=blue,
        label=r"$m'_{2D}>\sigma_m$",
        zorder=3,
    )

    # A zero value cannot be shown on logarithmic axes.  It is a censored
    # diagnostic value: no retained component exceeded R=2*vortex_lc=0.16.
    ax.scatter(
        [LIMIT_EK],
        [LIMIT_R],
        s=115,
        marker="v",
        facecolor=blue,
        edgecolor=blue,
        linewidth=1.2,
        zorder=5,
    )
    ax.annotate(
        "",
        xy=(LIMIT_EK, LIMIT_R / 1.38),
        xytext=(LIMIT_EK, LIMIT_R),
        arrowprops=dict(arrowstyle="-|>", color=blue, lw=2.2),
        annotation_clip=False,
    )
    ax.text(
        LIMIT_EK * 1.14,
        LIMIT_R * 1.08,
        r"$R_{\max}<0.16$",
        color=blue,
        fontsize=12,
        ha="left",
        va="bottom",
    )

    ax.set_xlabel(r"$Ek$")
    ax.set_ylabel(r"$\left\langle R_{\max}\right\rangle_{t_f-100:t_f}$")
    ax.set_xlim(1.0e-4, 5.0e-2)
    ax.set_ylim(0.10, 10.0)
    ax.legend(frameon=False, loc="lower right", handlelength=3.0)

    png = OUTDIR / "Ra8e6_late100_mean_maximum_mse_vortex_radius_vs_Ek_loglog.png"
    pdf = OUTDIR / "Ra8e6_late100_mean_maximum_mse_vortex_radius_vs_Ek_loglog.pdf"
    fig.savefig(png, facecolor="white")
    fig.savefig(pdf, facecolor="white")
    plt.close(fig)

    print(summary.to_string(index=False))
    print(f"CSV: {csv_path}")
    print(f"PNG: {png}")
    print(f"PDF: {pdf}")


if __name__ == "__main__":
    main()
