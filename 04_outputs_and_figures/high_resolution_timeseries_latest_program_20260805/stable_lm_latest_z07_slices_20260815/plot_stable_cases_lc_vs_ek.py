from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "lc_timeseries" / "Ra8e6_stable_Lm_cases_two_pi_lc_z05_timeseries.csv"
OUT = ROOT / "lc_vs_ek_scaling"
OUT.mkdir(parents=True, exist_ok=True)
WINDOW = 500.0


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "axes.linewidth": 4.5,
            "axes.labelsize": 24,
            "xtick.labelsize": 13,
            "ytick.labelsize": 13,
            "legend.fontsize": 12,
            "figure.dpi": 180,
            "savefig.dpi": 300,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


def late_means(data: pd.DataFrame, window: float) -> pd.DataFrame:
    rows = []
    for ek, group in data.groupby("Ek"):
        group = group.sort_values("time")
        time_max = float(group["time"].max())
        selected = group[group["time"] >= time_max - window]
        rows.append(
            {
                "Ek": float(ek),
                "time_start": float(selected["time"].min()),
                "time_end": time_max,
                "window": window,
                "n_samples": int(len(selected)),
                "mean_two_pi_lc_z05_H": float(selected["two_pi_lc_z05_H"].mean()),
                "std_two_pi_lc_z05_H": float(selected["two_pi_lc_z05_H"].std(ddof=1)),
            }
        )
    return pd.DataFrame(rows).sort_values("Ek").reset_index(drop=True)


def fit_power_law(summary: pd.DataFrame) -> dict[str, float]:
    x = np.log(summary["Ek"].to_numpy(dtype=float))
    y = np.log(summary["mean_two_pi_lc_z05_H"].to_numpy(dtype=float))
    coefficients, covariance = np.polyfit(x, y, 1, cov=True)
    alpha, log_c = map(float, coefficients)
    prediction = alpha * x + log_c
    residual = y - prediction
    r2 = 1.0 - float(np.sum(residual**2) / np.sum((y - np.mean(y)) ** 2))
    stderr_alpha = float(np.sqrt(covariance[0, 0]))
    # Student-t 97.5% quantile for n-2=3 degrees of freedom.
    ci_half_width = 3.182446305 * stderr_alpha
    return {
        "alpha": alpha,
        "C": float(np.exp(log_c)),
        "R2_log_space": r2,
        "stderr_alpha": stderr_alpha,
        "alpha_95ci_low": alpha - ci_half_width,
        "alpha_95ci_high": alpha + ci_half_width,
        "n_cases": int(len(summary)),
        "late_window": WINDOW,
    }


def window_sensitivity(data: pd.DataFrame) -> list[dict[str, float]]:
    results = []
    for window in (300.0, 500.0, 800.0, 1000.0):
        summary = late_means(data, window)
        fit = fit_power_law(summary)
        fit["late_window"] = window
        results.append(fit)
    return results


def draw(summary: pd.DataFrame, fit: dict[str, float]) -> None:
    x = summary["Ek"].to_numpy(dtype=float)
    y = summary["mean_two_pi_lc_z05_H"].to_numpy(dtype=float)
    xline = np.logspace(np.log10(0.85 * x.min()), np.log10(1.15 * x.max()), 300)
    yline = fit["C"] * xline ** fit["alpha"]

    fig = plt.figure(figsize=(6.5, 5.84), facecolor="white")
    ax = fig.add_axes([0.186, 0.260, 0.792, 0.705])
    ax.set_box_aspect(5.2 / 6.5)
    for spine in ax.spines.values():
        spine.set_linewidth(4.5)
    ax.tick_params(which="major", direction="in", length=12, width=1.2, top=True, right=True)
    ax.tick_params(which="minor", direction="in", length=6, width=1.0, top=True, right=True)
    ax.minorticks_on()
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$Ek$")
    ax.set_ylabel(r"$\left\langle2\pi l_c(z\simeq0.5)\right\rangle_{t}/H$")

    ax.plot(
        xline,
        yline,
        color="black",
        lw=3.5,
        ls="--",
        zorder=1,
        label=rf"$Ek^{{{fit['alpha']:.3f}}}$",
    )
    ax.plot(
        x,
        y,
        linestyle="none",
        marker="o",
        markersize=12,
        markerfacecolor=(0.11, 0.44, 0.71),
        markeredgecolor="black",
        markeredgewidth=2.0,
        zorder=3,
        label="late-time mean",
    )
    ax.legend(frameon=False, loc="lower right", handlelength=2.2)
    ax.set_xlim(0.75 * x.min(), 1.3 * x.max())
    ax.set_ylim(0.88 * y.min(), 1.13 * y.max())

    stem = OUT / "Ra8e6_stable_cases_two_pi_lc_z05_vs_Ek_loglog_powerlaw_fit"
    fig.savefig(stem.with_suffix(".png"))
    fig.savefig(stem.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    configure_style()
    data = pd.read_csv(SOURCE)
    summary = late_means(data, WINDOW)
    fit = fit_power_law(summary)
    sensitivity = window_sensitivity(data)
    summary.to_csv(OUT / "Ra8e6_stable_cases_two_pi_lc_z05_late500_means.csv", index=False, encoding="utf-8-sig")
    with (OUT / "Ra8e6_stable_cases_two_pi_lc_z05_powerlaw_fit.json").open("w", encoding="utf-8") as handle:
        json.dump(fit, handle, indent=2)
    with (OUT / "Ra8e6_stable_cases_two_pi_lc_z05_fit_window_sensitivity.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(sensitivity[0]))
        writer.writeheader()
        writer.writerows(sensitivity)
    draw(summary, fit)
    print(summary.to_string(index=False))
    print(json.dumps(fit, indent=2))
    print(OUT)


if __name__ == "__main__":
    main()
