#!/usr/bin/env python3
"""Plot latest-10-frame mean periodic MSE Voronoi diagrams."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection


STYLE = {
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


def style_axis(ax):
    for spine in ax.spines.values():
        spine.set_linewidth(4.5)
    ax.tick_params(direction="in", length=12, width=1.2, top=True, right=True)
    ax.minorticks_on()
    ax.tick_params(which="minor", direction="in", length=6, width=1.0, top=True, right=True)


def case_title(row):
    ra = float(row["Ra"])
    ra_exp = int(round(math.log10(ra)))
    if row["Ek"] == "norotating":
        return rf"$Ra=10^{{{ra_exp}}}$, no rotation"
    ek = float(row["Ek"])
    exponent = int(math.floor(math.log10(ek)))
    coefficient = ek / 10**exponent
    if abs(coefficient - 1.0) < 1e-10:
        ek_text = rf"10^{{{exponent}}}"
    else:
        ek_text = rf"{coefficient:g}\times10^{{{exponent}}}"
    return rf"$Ra=10^{{{ra_exp}}}$, $Ek={ek_text}$"


def draw_voronoi(ax, data):
    field = data["normalized_m2d"]
    lx = float(data["lx"])
    ly = float(data["ly"])
    image = ax.imshow(
        field,
        origin="lower",
        extent=(0.0, lx, 0.0, ly),
        cmap="RdBu_r",
        vmin=-3.0,
        vmax=3.0,
        interpolation="bilinear",
        aspect="equal",
        rasterized=True,
    )
    segments = data["segments"]
    if len(segments):
        collection = LineCollection(segments, colors="black", linewidths=0.8)
        collection.set_clip_on(True)
        ax.add_collection(collection)
    points = data["points_xy"]
    if len(points):
        ax.scatter(points[:, 0], points[:, 1], s=20, c="black", edgecolors="none", zorder=4)
    ax.set_xlim(0.0, lx)
    ax.set_ylim(0.0, ly)
    ax.set_aspect("equal", adjustable="box")
    style_axis(ax)
    return image


def save(fig, output, stem):
    fig.savefig(output / f"{stem}.png", dpi=300)
    fig.savefig(output / f"{stem}.pdf")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    source = args.input.resolve()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    mpl.rcParams.update(STYLE)

    with (source / "voronoi_latest10_metadata.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    rows.sort(key=lambda row: (float(row["Ra"]), row["Ek"] == "norotating", -float(row["Ek"]) if row["Ek"] != "norotating" else 0.0))

    loaded = []
    for row in rows:
        path = source / f"{row['case_id']}_latest10_mean_voronoi.npz"
        loaded.append((row, np.load(path)))
        fig = plt.figure(figsize=(7.4, 6.5), facecolor="white")
        ax = fig.add_axes([0.145, 0.145, 0.70, 0.76])
        image = draw_voronoi(ax, loaded[-1][1])
        ax.set_xlabel(r"$x/H$")
        ax.set_ylabel(r"$y/H$")
        ax.set_title(case_title(row), fontsize=20, pad=10)
        cax = fig.add_axes([0.875, 0.19, 0.035, 0.67])
        cbar = fig.colorbar(image, cax=cax)
        cbar.set_label(r"$\overline{m'_{2D}}^{,10}/\sigma_{\overline{m}}$", fontsize=20)
        cbar.ax.tick_params(labelsize=13, direction="in", length=7, width=1.2)
        for spine in cbar.ax.spines.values():
            spine.set_linewidth(2.0)
        stem = f"最新10帧平均周期维诺图_Latest10_mean_periodic_Voronoi_{row['case_id']}"
        save(fig, output, stem)

    ncols = 3
    nrows = int(math.ceil(len(loaded) / ncols))
    fig = plt.figure(figsize=(18.0, 17.5), facecolor="white")
    grid = fig.add_gridspec(nrows, ncols, left=0.06, right=0.91, bottom=0.06, top=0.96, wspace=0.28, hspace=0.30)
    image = None
    for index, (row, data) in enumerate(loaded):
        ax = fig.add_subplot(grid[index // ncols, index % ncols])
        image = draw_voronoi(ax, data)
        ax.set_title(case_title(row), fontsize=18, pad=8)
        ax.set_xlabel(r"$x/H$")
        ax.set_ylabel(r"$y/H$")
    for index in range(len(loaded), nrows * ncols):
        ax = fig.add_subplot(grid[index // ncols, index % ncols])
        ax.axis("off")
    cax = fig.add_axes([0.935, 0.19, 0.018, 0.64])
    cbar = fig.colorbar(image, cax=cax)
    cbar.set_label(r"$\overline{m'_{2D}}^{,10}/\sigma_{\overline{m}}$", fontsize=20)
    cbar.ax.tick_params(labelsize=13, direction="in", length=7, width=1.2)
    save(fig, output, "所有成熟AR10算例最新10帧平均周期维诺图_All_mature_AR10_latest10_mean_periodic_Voronoi")


if __name__ == "__main__":
    main()
