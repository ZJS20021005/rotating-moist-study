from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
METRICS = ROOT / "aggregation_metrics"
ONLINE = ROOT / "run" / "data"

BLUE = (0.0, 0.0, 1.0)
RED = (1.0, 0.0, 0.0)
GREEN = (0.0, 0.65, 0.12)
ORANGE = (0.96, 0.45, 0.05)

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


def style_axis(ax, box_aspect=True):
    for spine in ax.spines.values():
        spine.set_linewidth(4.5)
    ax.tick_params(direction="in", length=12, width=1.2, top=True, right=True)
    ax.minorticks_on()
    ax.tick_params(which="minor", direction="in", length=6, width=1.0, top=True, right=True)
    if box_aspect:
        ax.set_box_aspect(5.2 / 6.5)


def single_axis():
    fig = plt.figure(figsize=(6.5, 5.84), facecolor="white")
    ax = fig.add_axes([0.186, 0.260, 0.792, 0.705])
    style_axis(ax)
    return fig, ax


def save(fig, stem):
    fig.savefig(ROOT / f"{stem}.png", dpi=300)
    fig.savefig(ROOT / f"{stem}.pdf")
    plt.close(fig)


def load_online_mse():
    values = np.loadtxt(ONLINE / "mse_aggregation.out")
    if values.ndim == 1:
        values = values[None, :]
    values = values[np.all(np.isfinite(values[:, :3]), axis=1)]
    order = np.argsort(values[:, 0])
    values = values[order]
    _, reverse_index = np.unique(values[::-1, 0], return_index=True)
    keep = np.sort(len(values) - 1 - reverse_index)
    return values[keep]


def load_online_profile():
    values = np.loadtxt(ONLINE / "mse_variance_profile.out")
    if values.ndim == 1:
        values = values[None, :]
    values = values[np.all(np.isfinite(values[:, :4]), axis=1)]
    # Retain only complete 64-level time groups if the running solver was
    # writing the file while it was copied.
    times, counts = np.unique(values[:, 0], return_counts=True)
    complete = times[counts == 64]
    values = values[np.isin(values[:, 0], complete)]
    values = values[np.lexsort((values[:, 1], values[:, 0]))]
    time = np.unique(values[:, 0])
    z = np.unique(values[:, 1])
    if values.shape[0] != time.size * z.size:
        raise RuntimeError("Online profile is not a complete time-height grid")
    variance = values[:, 3].reshape(time.size, z.size)
    return time, z, variance


def voronoi_scale_by_frame(times):
    records = np.load(METRICS / "voronoi_area.npy")
    mean = np.full(times.size, np.nan)
    std = np.full(times.size, np.nan)
    for frame in range(times.size):
        values = records[records["frame"] == frame]["sqrt_area"]
        values = values[np.isfinite(values)]
        if values.size:
            mean[frame] = np.mean(values)
            std[frame] = np.std(values)
    return mean, std


def main():
    online = load_online_mse()
    online_time, online_am, online_rms = online[:, 0], online[:, 1], online[:, 2]
    frame_time = np.load(METRICS / "time.npy")
    field_am = np.load(METRICS / "MSE_variance.npy")
    l_peak = np.load(METRICS / "L_peak.npy")
    l_integral = np.load(METRICS / "L_integral.npy")
    cluster0 = np.load(METRICS / "cluster_summary_threshold0.npy")
    cluster_sigma = np.load(METRICS / "cluster_summary_threshold_sigma.npy")
    vor_mean, vor_std = voronoi_scale_by_frame(frame_time)

    fig, ax = single_axis()
    ax.plot(online_time, online_am, color=BLUE, lw=3.5, label="online")
    ax.plot(frame_time, field_am, ls="none", marker="o", ms=7.5, color=RED, label="full-field check")
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$A_m$")
    ax.set_xlim(0.0, online_time.max())
    ax.set_ylim(bottom=0.0)
    ax.legend(frameon=False, loc="upper right")
    save(fig, "Ek3e-3_MSE_horizontal_variance_timeseries")

    profile_time, z, variance = load_online_profile()
    fig = plt.figure(figsize=(7.8, 5.84), facecolor="white")
    ax = fig.add_axes([0.15, 0.20, 0.68, 0.72])
    image = ax.pcolormesh(profile_time, z, variance.T, shading="auto", cmap="magma")
    style_axis(ax)
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$z/H$")
    cax = fig.add_axes([0.86, 0.20, 0.035, 0.72])
    colorbar = fig.colorbar(image, cax=cax)
    colorbar.set_label(r"${\rm Var}_{xy}(m)$", fontsize=20)
    colorbar.ax.tick_params(labelsize=13, width=1.2, length=7, direction="in")
    save(fig, "Ek3e-3_MSE_variance_height_time")

    fig, ax = single_axis()
    ax.plot(frame_time, l_peak, color=RED, lw=3.5, marker="o", ms=7.5, label=r"$L_{\rm peak}$")
    ax.plot(frame_time, l_integral, color=BLUE, lw=3.5, marker="s", ms=6.5, ls="--", label=r"$L_{\rm int}$")
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"aggregation scale")
    ax.set_xlim(0.0, online_time.max())
    ax.set_ylim(bottom=0.0)
    ax.legend(frameon=False, loc="upper left")
    save(fig, "Ek3e-3_MSE_spectral_scales_timeseries")

    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.84), facecolor="white")
    for ax, data, title in (
        (axes[0], cluster0, r"$m'_{2D}>0$"),
        (axes[1], cluster_sigma, r"$m'_{2D}>\sigma_m$"),
    ):
        ax.plot(frame_time, data[:, 1], color=BLUE, lw=3.5, marker="o", ms=6.5, label="mean radius")
        ax.plot(frame_time, data[:, 2], color=RED, lw=3.5, marker="s", ms=6.0, ls="--", label="maximum radius")
        style_axis(ax)
        ax.set_xlabel(r"$t$")
        ax.set_ylabel(r"cluster radius")
        ax.set_xlim(0.0, online_time.max())
        ax.set_ylim(bottom=0.0)
        ax.set_title(title, fontsize=20)
        ax.legend(frameon=False, loc="best")
    fig.subplots_adjust(left=0.09, right=0.97, bottom=0.19, top=0.90, wspace=0.32)
    save(fig, "Ek3e-3_cluster_scales_timeseries")

    fig, ax = single_axis()
    ax.plot(frame_time, vor_mean, color=GREEN, lw=3.5, marker="o", ms=7.0)
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$L_{\rm Vor}=\langle\sqrt{A_V}\rangle$")
    ax.set_xlim(0.0, online_time.max())
    ax.set_ylim(bottom=0.0)
    save(fig, "Ek3e-3_Voronoi_scale_timeseries")

    np.savetxt(
        ROOT / "Ek3e-3_online_MSE_variance_timeseries.csv",
        online,
        delimiter=",",
        header="time,A_m,sqrt_A_m",
        comments="",
    )
    columns = np.c_[
        frame_time,
        field_am,
        l_peak,
        l_integral,
        cluster0,
        cluster_sigma,
        vor_mean,
        vor_std,
    ]
    np.savetxt(
        ROOT / "Ek3e-3_field_aggregation_scales_timeseries.csv",
        columns,
        delimiter=",",
        header=(
            "time,A_m_field,L_peak,L_integral,"
            "cluster_number_threshold0,cluster_mean_radius_threshold0,cluster_max_radius_threshold0,"
            "cluster_number_threshold_sigma,cluster_mean_radius_threshold_sigma,cluster_max_radius_threshold_sigma,"
            "Voronoi_mean_sqrt_area,Voronoi_std_sqrt_area"
        ),
        comments="",
    )

    post_min = int(np.argmin(online_am))
    late = online_time >= 40.0
    late_indices = np.flatnonzero(late)
    late_max = late_indices[int(np.argmax(online_am[late]))]
    summary = {
        "case": "Ra1e8Pr07Ek3e-3AR10",
        "gamma": 1.1,
        "online_last_time": float(online_time[-1]),
        "online_records": int(online_time.size),
        "field_frames": int(frame_time.size),
        "A_m_min": float(online_am[post_min]),
        "A_m_min_time": float(online_time[post_min]),
        "A_m_late_max": float(online_am[late_max]),
        "A_m_late_max_time": float(online_time[late_max]),
        "A_m_latest": float(online_am[-1]),
        "L_integral_first": float(l_integral[0]),
        "L_integral_latest": float(l_integral[-1]),
        "L_peak_latest": float(l_peak[-1]),
        "cluster_sigma_mean_radius_first": float(cluster_sigma[0, 1]),
        "cluster_sigma_mean_radius_latest": float(cluster_sigma[-1, 1]),
        "cluster_sigma_count_first": int(cluster_sigma[0, 0]),
        "cluster_sigma_count_latest": int(cluster_sigma[-1, 0]),
        "Voronoi_scale_first": float(vor_mean[0]),
        "Voronoi_scale_latest": float(vor_mean[-1]),
        "remote_source": (
            "/share/org/SHUTUANL/shu_zhangjs/rainy model/ns/transition_study/beta1/"
            "lowrestest/aspect_ratio_study/Ra1e8/AR10/Ra1e8Pr07Ek3e-3AR10"
        ),
    }
    (ROOT / "Ek3e-3_aggregation_current_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
