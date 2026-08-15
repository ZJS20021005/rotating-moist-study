#!/usr/bin/env python3
"""Publication plots for vertically averaged MSE-mode shell transfer."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import TwoSlopeNorm


ROOT = Path(
    r"E:\moist RB\rotating_case_inventory\04_outputs_and_figures"
    r"\Ra8e6_mse_2dmode_shell_to_shell_20260802"
)
DATA = ROOT / "final_2dmode_latest10"
OUT = ROOT / "figures"

CASE_ORDER = ["Ek1p5e-4", "Ek2e-4", "Ek3e-3", "Ek5e-3", "Ek7e-3", "Ek1e-2", "Ek3e-2"]
CASE_EK = {
    "Ek1p5e-4": 1.5e-4,
    "Ek2e-4": 2.0e-4,
    "Ek3e-3": 3.0e-3,
    "Ek5e-3": 5.0e-3,
    "Ek7e-3": 7.0e-3,
    "Ek1e-2": 1.0e-2,
    "Ek3e-2": 3.0e-2,
}

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
        "axes.unicode_minus": True,
    }
)


def ek_text(value: float) -> str:
    exponent = int(np.floor(np.log10(value)))
    mantissa = value / 10.0**exponent
    if abs(mantissa - 1.0) < 1e-10:
        return rf"10^{{{exponent}}}"
    return rf"{mantissa:g}\times10^{{{exponent}}}"


def style_axis(axis):
    for spine in axis.spines.values():
        spine.set_linewidth(4.5)
    axis.tick_params(direction="in", length=12, width=1.2, top=True, right=True)
    axis.minorticks_on()
    axis.tick_params(
        which="minor", direction="in", length=6, width=1.0, top=True, right=True
    )


def save(figure, stem):
    OUT.mkdir(parents=True, exist_ok=True)
    png = OUT / f"{stem}.png"
    pdf = OUT / f"{stem}.pdf"
    figure.savefig(png)
    figure.savefig(pdf)
    plt.close(figure)
    return png, pdf


def load_case(label):
    return np.load(DATA / label / "mse_2dmode_shell_to_shell_latest10.npz")


def per_case_heatmap(label, plot_limit=None):
    data = load_case(label)
    matrix = data["T_2D_exchange_frames"].mean(axis=0)
    row_sum = data["T_2D_exchange_row_sum_frames"].mean(axis=0)
    edges = data["k_edges"]
    centers = data["k_centers"]
    if plot_limit is not None:
        count = int(np.searchsorted(edges, plot_limit, side="right") - 1)
        count = max(2, min(count, centers.size))
        matrix = matrix[:count, :count]
        row_sum = row_sum[:count]
        centers = centers[:count]
        edges = edges[: count + 1]
    vmax = float(np.max(np.abs(matrix)))
    if vmax == 0.0:
        vmax = 1.0

    figure = plt.figure(figsize=(9.2, 6.2), facecolor="white")
    axis = figure.add_axes([0.11, 0.12, 0.54, 0.80])
    marginal = figure.add_axes([0.68, 0.12, 0.18, 0.80], sharey=axis)
    colorbar_axis = figure.add_axes([0.89, 0.12, 0.025, 0.80])

    mesh = axis.pcolormesh(
        edges,
        edges,
        matrix,
        cmap="RdBu_r",
        norm=TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax),
        shading="flat",
        rasterized=True,
    )
    axis.plot([edges[0], edges[-1]], [edges[0], edges[-1]], color="black", lw=1.5)
    axis.set_xlim(edges[0], edges[-1])
    axis.set_ylim(edges[0], edges[-1])
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlabel(r"Donor wavenumber $P$")
    axis.set_ylabel(r"Receiver wavenumber $K$")
    style_axis(axis)

    marginal.axvline(0.0, color="0.35", lw=2.0, ls="--", zorder=1)
    marginal.plot(row_sum, centers, color=(0.00, 0.38, 0.95), lw=3.5, zorder=2)
    marginal.set_ylim(edges[0], edges[-1])
    marginal.set_xlabel(r"$\sum_P T_{m,\mathrm{ex}}^{2D}$", fontsize=19)
    marginal.tick_params(labelleft=False)
    marginal.ticklabel_format(axis="x", style="sci", scilimits=(0, 0), useMathText=True)
    style_axis(marginal)

    colorbar = figure.colorbar(mesh, cax=colorbar_axis)
    colorbar.outline.set_linewidth(4.5)
    colorbar.ax.tick_params(direction="in", length=8, width=1.2, labelsize=13)
    colorbar.set_label(r"$T_{m,\mathrm{ex}}^{2D}(K,P)$", fontsize=20)
    colorbar.formatter.set_powerlimits((0, 0))
    colorbar.update_ticks()

    figure.text(
        0.38,
        0.965,
        rf"$Ek={ek_text(CASE_EK[label])}$",
        ha="center",
        va="top",
        fontsize=21,
    )
    suffix = "" if plot_limit is None else f"_k{int(plot_limit)}"
    return save(figure, f"{label}_MSE_2Dmode_TKP_exchange{suffix}_latest10")


def normalized_gallery():
    figure, axes = plt.subplots(2, 4, figsize=(13.0, 11.68), facecolor="white")
    axes = axes.ravel()
    for axis, label in zip(axes, CASE_ORDER):
        data = load_case(label)
        matrix = data["T_2D_exchange_frames"].mean(axis=0)
        edges = data["k_edges"]
        scale = float(np.max(np.abs(matrix)))
        normalized = matrix / scale if scale > 0 else matrix
        axis.pcolormesh(
            edges,
            edges,
            normalized,
            cmap="RdBu_r",
            vmin=-1.0,
            vmax=1.0,
            shading="flat",
            rasterized=True,
        )
        axis.plot([edges[0], edges[-1]], [edges[0], edges[-1]], color="black", lw=1.2)
        axis.set_aspect("equal", adjustable="box")
        axis.set_xlim(edges[0], edges[-1])
        axis.set_ylim(edges[0], edges[-1])
        axis.text(
            0.50,
            0.92,
            rf"$Ek={ek_text(CASE_EK[label])}$",
            transform=axis.transAxes,
            ha="center",
            va="top",
            fontsize=16,
        )
        style_axis(axis)
    axes[-1].axis("off")
    for axis in axes[:4]:
        axis.set_xticklabels([])
    for axis in [axes[1], axes[2], axes[3], axes[5], axes[6]]:
        axis.set_yticklabels([])
    figure.text(0.53, 0.035, r"Donor wavenumber $P$", ha="center", fontsize=24)
    figure.text(0.025, 0.53, r"Receiver wavenumber $K$", va="center", rotation=90, fontsize=24)
    figure.subplots_adjust(left=0.08, right=0.98, bottom=0.10, top=0.98, wspace=0.08, hspace=0.12)
    return save(figure, "Ra8e6_MSE_2Dmode_TKP_exchange_normalized_gallery_latest10")


def forcing_decomposition(label):
    data = load_case(label)
    k = data["k_centers"]
    self_budget = data["T_2D_row_sum_frames"].mean(axis=0)
    eddy = data["F_3D_to_2D_shell_frames"].mean(axis=0)
    exact = data["T_2D_exact_shell_frames"].mean(axis=0)

    figure = plt.figure(figsize=(6.5, 5.84), facecolor="white")
    axis = figure.add_axes([0.186, 0.260, 0.792, 0.705])
    axis.axhline(0.0, color="0.35", lw=2.0, ls="--")
    axis.plot(k, self_budget, color=(0.00, 0.38, 0.95), lw=3.5, label=r"$T_{2D\rightarrow2D}$")
    axis.plot(k, eddy, color=(0.95, 0.15, 0.08), lw=3.5, ls="--", label=r"$F_{3D\rightarrow2D}$")
    axis.plot(k, exact, color=(0.10, 0.70, 0.18), lw=3.5, ls="-.", label=r"$N_M$")
    axis.set_xlim(k[0], k[-1])
    axis.set_xlabel(r"$k_h$")
    axis.set_ylabel(r"MSE variance input")
    axis.ticklabel_format(axis="y", style="sci", scilimits=(0, 0), useMathText=True)
    style_axis(axis)
    axis.set_box_aspect(5.2 / 6.5)
    axis.legend(frameon=False, loc="upper right")
    return save(figure, f"{label}_MSE_2Dmode_forcing_decomposition_latest10")


def summary_rows():
    rows = []
    for label in CASE_ORDER:
        data = load_case(label)
        matrix = data["T_2D_exchange_frames"].mean(axis=0)
        raw = data["T_2D_frames"].mean(axis=0)
        k = data["k_centers"]
        dk = float(data["matrix_dk"])
        receiver = k[:, None]
        donor = k[None, :]
        upscale = float(
            np.sum(np.where(receiver < donor, np.clip(matrix, 0.0, None), 0.0)) * dk
        )
        downscale = float(
            np.sum(np.where(receiver > donor, np.clip(matrix, 0.0, None), 0.0)) * dk
        )
        raw_scale = max(float(np.max(np.abs(raw))), 1.0e-300)
        anti_error = float(np.max(np.abs(raw + raw.T)) / raw_scale)
        eddy = data["F_3D_to_2D_shell_frames"].mean(axis=0)
        self_net = data["T_2D_row_sum_frames"].mean(axis=0)
        rows.append(
            {
                "case": label,
                "Ek": CASE_EK[label],
                "t_first": float(data["times"][0]),
                "t_last": float(data["times"][-1]),
                "nframes": int(data["times"].size),
                "matrix_dk": dk,
                "matrix_kmax_edge": float(data["k_edges"][-1]),
                "max_abs_T_exchange": float(np.max(np.abs(matrix))),
                "positive_upscale_transfer": upscale,
                "positive_downscale_transfer": downscale,
                "upscale_to_downscale_positive_ratio": upscale / max(downscale, 1.0e-300),
                "mean_matrix_antisymmetry_relative_error": anti_error,
                "max_abs_2D_self_net": float(np.max(np.abs(self_net))),
                "max_abs_3D_to_2D_forcing": float(np.max(np.abs(eddy))),
                "3D_forcing_to_2D_self_net_ratio": float(
                    np.max(np.abs(eddy)) / max(np.max(np.abs(self_net)), 1.0e-300)
                ),
            }
        )
    return rows


def main():
    outputs = []
    for label in CASE_ORDER:
        outputs.extend(per_case_heatmap(label))
        outputs.extend(per_case_heatmap(label, plot_limit=30.0))
        outputs.extend(forcing_decomposition(label))
    outputs.extend(normalized_gallery())
    rows = summary_rows()
    csv_path = OUT / "Ra8e6_MSE_2Dmode_TKP_summary_latest10.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print("\n".join(str(path) for path in outputs))
    print(csv_path)


if __name__ == "__main__":
    main()
