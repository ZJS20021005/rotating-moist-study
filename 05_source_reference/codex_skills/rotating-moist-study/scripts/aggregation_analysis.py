#!/usr/bin/env python3
"""Moist-convective self-aggregation diagnostics for one DNS case.

The analysis is deliberately based on horizontal MSE anomalies

    m'(x,y,z,t) = m - <m>_xy,   m = b + gamma*q,

so the imposed vertical background profile is not counted as aggregation.
Every diagnostic is computed on each field snapshot first; the resulting time
series can subsequently be averaged or compared between cases.
"""

from __future__ import annotations

import argparse
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path

import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from scipy import ndimage
from scipy.spatial import Voronoi


FIELD_RE = re.compile(r"field0*(\d+)\.h5$", re.I)
TIME_RE = re.compile(r'<Time\s+Value="\s*([^\"]+)"', re.I)
CONTI_RE = re.compile(r"conti(\d*)", re.I)

MPL_STYLE = {
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


@dataclass(frozen=True)
class FrameRecord:
    path: Path
    time: float
    number: int
    continuation: tuple[int, ...]


def continuation_rank(path: Path) -> tuple[int, ...]:
    values = []
    for token in CONTI_RE.findall(str(path)):
        values.append(int(token) if token else 1)
    return tuple(values)


def read_time(path: Path) -> float:
    xmf = path.with_suffix(".xmf")
    if xmf.is_file():
        try:
            match = TIME_RE.search(xmf.read_text(errors="replace"))
            if match:
                return float(match.group(1).replace("D", "E").replace("d", "e"))
        except OSError:
            pass
    return math.nan


def discover_frames(case_root: Path, stride: int, last: int | None) -> list[FrameRecord]:
    records = []
    for path in case_root.rglob("field*.h5"):
        match = FIELD_RE.fullmatch(path.name)
        if not match or path.parent.name != "movie":
            continue
        records.append(
            FrameRecord(
                path=path,
                time=read_time(path),
                number=int(match.group(1)),
                continuation=continuation_rank(path),
            )
        )
    if not records:
        raise FileNotFoundError(f"No movie/fieldNNNNN.h5 below {case_root}")

    # Continuations can repeat field numbers.  Physical time is the preferred
    # key; if a time occurs twice, retain the file from the later continuation.
    if any(np.isfinite(record.time) for record in records):
        records.sort(
            key=lambda record: (
                record.time if np.isfinite(record.time) else -np.inf,
                record.continuation,
                record.number,
            )
        )
        unique: dict[float | tuple, FrameRecord] = {}
        for record in records:
            key = round(record.time, 10) if np.isfinite(record.time) else (
                record.continuation,
                record.number,
            )
            unique[key] = record
        records = sorted(
            unique.values(),
            key=lambda record: (
                record.time if np.isfinite(record.time) else -np.inf,
                record.continuation,
                record.number,
            ),
        )
    else:
        records.sort(key=lambda record: (record.continuation, record.number, str(record.path)))

    records = records[:: max(1, stride)]
    if last is not None and last > 0:
        records = records[-last:]
    return records


def find_grid(frame: Path, case_root: Path) -> Path:
    for parent in frame.parents:
        candidate = parent / "field_gridc.h5"
        if candidate.is_file():
            return candidate
        if parent == case_root:
            break
    grids = sorted(case_root.rglob("field_gridc.h5"), key=lambda p: (continuation_rank(p), str(p)))
    if not grids:
        raise FileNotFoundError(f"No field_gridc.h5 below {case_root}")
    return grids[-1]


def read_grid(grid_path: Path, shape: tuple[int, int, int]):
    nz, ny, nx = shape
    with h5py.File(grid_path, "r") as source:
        z = np.asarray(source["zc"][:nz], dtype=float)
        x = np.asarray(source["xc"][:nx], dtype=float)
        y = np.asarray(source["yc"][:ny], dtype=float)
        zf = np.asarray(source["zf"][:nz], dtype=float) if "zf" in source else None
        xf = np.asarray(source["xf"][: nx + 1], dtype=float) if "xf" in source else None
        yf = np.asarray(source["yf"][: ny + 1], dtype=float) if "yf" in source else None

    if zf is not None and zf.size == nz:
        top = 2.0 * z[-1] - zf[-1]
        z_faces = np.r_[zf, top]
    else:
        z_faces = np.empty(nz + 1)
        z_faces[1:-1] = 0.5 * (z[:-1] + z[1:])
        z_faces[0] = z[0] - 0.5 * (z[1] - z[0])
        z_faces[-1] = z[-1] + 0.5 * (z[-1] - z[-2])
    dz = np.diff(z_faces)

    dx = float(np.median(np.diff(x)))
    dy = float(np.median(np.diff(y)))
    lx = float(xf[nx] - xf[0]) if xf is not None and xf.size > nx else nx * dx
    ly = float(yf[ny] - yf[0]) if yf is not None and yf.size > ny else ny * dy
    return x, y, z, dx, dy, dz, lx, ly


def read_m(frame: Path, gamma: float) -> np.ndarray:
    with h5py.File(frame, "r") as source:
        if "DSAL_me" not in source or "QVAP_me" not in source:
            raise KeyError(f"{frame}: DSAL_me or QVAP_me is absent")
        b = np.asarray(source["DSAL_me"], dtype=float)
        q = np.asarray(source["QVAP_me"], dtype=float)
    return b + gamma * q


def radial_shell_spectrum(field: np.ndarray, dx: float, dy: float, lx: float, ly: float):
    ny, nx = field.shape
    transform = np.fft.fft2(field)
    power = np.abs(transform) ** 2 / float(nx * ny) ** 2
    kx = 2.0 * np.pi * np.fft.fftfreq(nx, d=dx)
    ky = 2.0 * np.pi * np.fft.fftfreq(ny, d=dy)
    kh = np.hypot(ky[:, None], kx[None, :])
    dk = 2.0 * np.pi / max(lx, ly)
    shell = np.rint(kh / dk).astype(int)
    spectrum = np.bincount(shell.ravel(), weights=power.ravel())
    count = np.bincount(shell.ravel())
    k = np.arange(spectrum.size, dtype=float) * dk
    valid = count > 0
    return k[valid], spectrum[valid]


def radial_correlation(field: np.ndarray, dx: float, dy: float, lx: float, ly: float):
    ny, nx = field.shape
    transform = np.fft.fft2(field)
    corr2d = np.fft.fftshift(np.fft.ifft2(np.abs(transform) ** 2).real)
    center = corr2d[ny // 2, nx // 2]
    if center <= np.finfo(float).eps:
        return np.array([0.0]), np.array([1.0]), 0.0
    corr2d /= center

    xx = (np.arange(nx) - nx // 2) * dx
    yy = (np.arange(ny) - ny // 2) * dy
    radius = np.hypot(yy[:, None], xx[None, :])
    dr = min(dx, dy)
    index = np.floor(radius / dr + 0.5).astype(int)
    maximum = int(min(lx, ly) * 0.5 / dr)
    weight = np.bincount(index.ravel(), weights=corr2d.ravel(), minlength=maximum + 1)
    count = np.bincount(index.ravel(), minlength=maximum + 1)
    count = count[: maximum + 1]
    correlation = weight[: maximum + 1] / np.maximum(count, 1)
    r = np.arange(maximum + 1, dtype=float) * dr

    nonpositive = np.flatnonzero(correlation[1:] <= 0.0)
    stop = int(nonpositive[0] + 1) if nonpositive.size else correlation.size - 1
    if stop < 1:
        length = 0.0
    else:
        # Interpolate the first zero so the positive-lobe integral is not
        # biased by one radial bin.
        r_use = r[: stop + 1].copy()
        c_use = correlation[: stop + 1].copy()
        if c_use[-1] < 0.0 and c_use[-2] > 0.0:
            fraction = c_use[-2] / (c_use[-2] - c_use[-1])
            r_use[-1] = r_use[-2] + fraction * (r_use[-1] - r_use[-2])
            c_use[-1] = 0.0
        length = float(np.trapz(c_use, r_use))
    return r, correlation, length


def periodic_components(mask: np.ndarray):
    tiled = np.tile(mask, (3, 3))
    labelled, _ = ndimage.label(tiled, structure=np.ones((3, 3), dtype=int))
    ny, nx = mask.shape
    central = labelled[ny : 2 * ny, nx : 2 * nx]
    labels = np.unique(central)
    return central, labels[labels > 0]


def longest_false_run(occupied: np.ndarray) -> int:
    n = occupied.size
    if np.all(occupied):
        return 0
    if not np.any(occupied):
        return n
    doubled = np.r_[~occupied, ~occupied]
    best = run = 0
    for value in doubled:
        run = run + 1 if value else 0
        best = max(best, run)
    return min(best, n)


def cluster_records(mask: np.ndarray, dx: float, dy: float):
    central, labels = periodic_components(mask)
    records = []
    ny, nx = mask.shape
    for label in labels:
        component = central == label
        count = int(component.sum())
        area = count * dx * dy
        radius = math.sqrt(area / math.pi)
        rows, cols = np.nonzero(component)
        occupied_x = np.zeros(nx, dtype=bool)
        occupied_y = np.zeros(ny, dtype=bool)
        occupied_x[cols] = True
        occupied_y[rows] = True
        extent_x = (nx - longest_false_run(occupied_x)) * dx
        extent_y = (ny - longest_false_run(occupied_y)) * dy
        records.append((area, radius, max(extent_x, extent_y)))
    return records


def select_periodic_peaks(field: np.ndarray, sigma: float, min_distance: int):
    size = 2 * min_distance + 1
    local_max = field == ndimage.maximum_filter(field, size=size, mode="wrap")
    candidates = np.argwhere(local_max & (field > sigma))
    if candidates.size == 0:
        return np.empty((0, 2), dtype=int)
    order = np.argsort(field[candidates[:, 0], candidates[:, 1]])[::-1]
    candidates = candidates[order]
    selected: list[np.ndarray] = []
    ny, nx = field.shape
    for candidate in candidates:
        if selected:
            old = np.asarray(selected)
            dy = np.abs(old[:, 0] - candidate[0])
            dx = np.abs(old[:, 1] - candidate[1])
            dy = np.minimum(dy, ny - dy)
            dx = np.minimum(dx, nx - dx)
            if np.any(dx * dx + dy * dy < min_distance * min_distance):
                continue
        selected.append(candidate)
    return np.asarray(selected, dtype=int)


def polygon_area(vertices: np.ndarray) -> float:
    x = vertices[:, 0]
    y = vertices[:, 1]
    return 0.5 * abs(float(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1))))


def periodic_voronoi(points: np.ndarray, lx: float, ly: float):
    """Return central-cell Voronoi areas and polygons for periodic points."""
    if len(points) < 2:
        return np.empty(0), []
    tiled = []
    central_indices = []
    for sy in (-ly, 0.0, ly):
        for sx in (-lx, 0.0, lx):
            start = len(tiled)
            tiled.extend(points + np.array([sx, sy]))
            if sx == 0.0 and sy == 0.0:
                central_indices = list(range(start, start + len(points)))
    voronoi = Voronoi(np.asarray(tiled))
    areas = []
    polygons = []
    for point_index in central_indices:
        region_index = voronoi.point_region[point_index]
        region = voronoi.regions[region_index]
        if not region or -1 in region:
            areas.append(math.nan)
            polygons.append(np.empty((0, 2)))
            continue
        vertices = voronoi.vertices[region]
        areas.append(polygon_area(vertices))
        polygons.append(vertices)
    return np.asarray(areas, dtype=float), polygons


def style_axis(ax: mpl.axes.Axes, box_aspect: bool = True):
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


def save_figure(fig: mpl.figure.Figure, output: Path, stem: str):
    fig.savefig(output / f"{stem}.png", dpi=300)
    fig.savefig(output / f"{stem}.pdf")
    plt.close(fig)


def plot_outputs(
    output: Path,
    times: np.ndarray,
    z: np.ndarray,
    a_m: np.ndarray,
    var_z_time: np.ndarray,
    l_peak: np.ndarray,
    l_integral: np.ndarray,
    cluster_summary: dict[str, np.ndarray],
    last_m2d: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    last_points: np.ndarray,
    last_areas: np.ndarray,
    last_polygons: list[np.ndarray],
):
    mpl.rcParams.update(MPL_STYLE)
    time_label = times if np.all(np.isfinite(times)) else np.arange(len(times))
    xlabel = r"$t$" if np.all(np.isfinite(times)) else "snapshot index"

    fig, ax = single_axis()
    ax.plot(time_label, a_m, color=(0.0, 0.0, 1.0), lw=3.5)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"$A_m$")
    save_figure(fig, output, "Figure1_MSE_horizontal_variance")

    fig = plt.figure(figsize=(7.8, 5.84), facecolor="white")
    ax = fig.add_axes([0.15, 0.20, 0.68, 0.72])
    image = ax.pcolormesh(time_label, z, var_z_time.T, shading="auto", cmap="magma")
    style_axis(ax)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"$z/H$")
    cax = fig.add_axes([0.86, 0.20, 0.035, 0.72])
    colorbar = fig.colorbar(image, cax=cax)
    colorbar.set_label(r"${\rm Var}_{xy}(m)$", fontsize=20)
    colorbar.ax.tick_params(labelsize=13, width=1.2, length=7, direction="in")
    save_figure(fig, output, "Figure1b_MSE_variance_height_time")

    fig, ax = single_axis()
    ax.plot(time_label, l_peak, color=(0.74, 0.14, 0.18), lw=3.5, label=r"$L_{\rm peak}$")
    ax.plot(time_label, l_integral, color=(0.0, 0.0, 1.0), lw=3.5, ls="--", label=r"$L_{\rm int}$")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"aggregation scale")
    ax.legend(frameon=False, loc="best")
    save_figure(fig, output, "Figure2_MSE_aggregation_scales")

    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.84), facecolor="white")
    for ax, key, title in zip(
        axes,
        ("threshold0", "threshold_sigma"),
        (r"$m'>0$", r"$m'>\sigma_m$"),
    ):
        summary = cluster_summary[key]
        ax.plot(time_label, summary[:, 1], color=(0.0, 0.0, 1.0), lw=3.5, label="mean radius")
        ax.plot(time_label, summary[:, 2], color=(1.0, 0.0, 0.0), lw=3.5, ls="--", label="maximum radius")
        style_axis(ax)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(r"cluster radius")
        ax.set_title(title, fontsize=20)
        ax.legend(frameon=False, loc="best")
    fig.subplots_adjust(left=0.09, right=0.97, bottom=0.19, top=0.90, wspace=0.32)
    save_figure(fig, output, "Figure3_cluster_radius_evolution")

    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.84), facecolor="white")
    ax = axes[0]
    image = ax.imshow(
        last_m2d,
        origin="lower",
        extent=(0.0, x.size * np.median(np.diff(x)), 0.0, y.size * np.median(np.diff(y))),
        cmap="RdBu_r",
        aspect="equal",
    )
    segments = []
    for polygon in last_polygons:
        if polygon.size:
            segments.extend(np.stack([polygon, np.roll(polygon, -1, axis=0)], axis=1))
    if segments:
        ax.add_collection(LineCollection(segments, colors="black", linewidths=0.8))
    if len(last_points):
        ax.scatter(last_points[:, 0], last_points[:, 1], s=18, c="black")
    style_axis(ax)
    ax.set_xlabel(r"$x/H$")
    ax.set_ylabel(r"$y/H$")
    ax.set_xlim(0.0, x.size * np.median(np.diff(x)))
    ax.set_ylim(0.0, y.size * np.median(np.diff(y)))
    cbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(r"$m'_{2D}$", fontsize=20)

    ax = axes[1]
    finite = last_areas[np.isfinite(last_areas) & (last_areas > 0.0)]
    if finite.size:
        ax.hist(finite, bins=min(20, max(5, int(np.sqrt(finite.size)))), color=(0.0, 0.0, 1.0), alpha=0.85)
    style_axis(ax)
    ax.set_xlabel(r"Voronoi cell area")
    ax.set_ylabel(r"count")
    fig.subplots_adjust(left=0.08, right=0.97, bottom=0.19, top=0.94, wspace=0.35)
    save_figure(fig, output, "Figure4_periodic_Voronoi")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-root", type=Path, default=Path.cwd())
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--gamma", type=float, default=1.1)
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--last", type=int, default=None, help="Use only the latest N physical-time snapshots")
    parser.add_argument("--vertical-mode", choices=("mean", "height"), default="mean")
    parser.add_argument("--height", type=float, default=0.5)
    parser.add_argument("--min-peak-distance", type=int, default=3)
    parser.add_argument("--no-plots", action="store_true")
    args = parser.parse_args()

    case_root = args.case_root.expanduser().resolve()
    output = (args.output or (case_root / "aggregation_metrics")).expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    frames = discover_frames(case_root, args.stride, args.last)

    with h5py.File(frames[0].path, "r") as source:
        shape = tuple(source["DSAL_me"].shape)
    if len(shape) != 3:
        raise ValueError(f"Expected a 3-D field, found {shape}")
    grid_path = find_grid(frames[0].path, case_root)
    x, y, z, dx, dy, dz, lx, ly = read_grid(grid_path, shape)
    z_weights = dz / dz.sum()
    height_index = int(np.argmin(np.abs(z - args.height)))

    nt = len(frames)
    times = np.array([record.time for record in frames], dtype=float)
    if not np.all(np.isfinite(times)):
        times = np.arange(nt, dtype=float)
    var_z_time = np.empty((nt, shape[0]), dtype=float)
    a_m = np.empty(nt, dtype=float)
    l_peak = np.empty(nt, dtype=float)
    l_integral = np.empty(nt, dtype=float)
    spectra = []
    correlations = []
    cluster_all = {"threshold0": [], "threshold_sigma": []}
    cluster_summary = {
        "threshold0": np.full((nt, 3), np.nan),
        "threshold_sigma": np.full((nt, 3), np.nan),
    }
    voronoi_all = []
    last_m2d = np.empty((shape[1], shape[2]))
    last_points_xy = np.empty((0, 2))
    last_areas = np.empty(0)
    last_polygons: list[np.ndarray] = []

    for iframe, record in enumerate(frames):
        m = read_m(record.path, args.gamma)
        if m.shape != shape:
            raise ValueError(f"{record.path}: shape changed from {shape} to {m.shape}")
        mprime = m - m.mean(axis=(1, 2), keepdims=True)
        var_profile = np.mean(mprime * mprime, axis=(1, 2))
        var_z_time[iframe] = var_profile
        a_m[iframe] = float(np.dot(var_profile, z_weights))

        if args.vertical_mode == "mean":
            m2d = np.tensordot(z_weights, mprime, axes=(0, 0))
        else:
            m2d = mprime[height_index]
        m2d -= m2d.mean()
        sigma = float(np.sqrt(np.mean(m2d * m2d)))

        k, spectrum = radial_shell_spectrum(m2d, dx, dy, lx, ly)
        spectra.append(spectrum)
        if spectrum.size > 1 and np.any(spectrum[1:] > 0.0):
            peak_index = 1 + int(np.argmax(spectrum[1:]))
            l_peak[iframe] = 2.0 * np.pi / k[peak_index]
        else:
            l_peak[iframe] = np.nan
        r, corr, l_integral[iframe] = radial_correlation(m2d, dx, dy, lx, ly)
        correlations.append(corr)

        for key, threshold in (("threshold0", 0.0), ("threshold_sigma", sigma)):
            records_now = cluster_records(m2d > threshold, dx, dy)
            for icluster, (area, radius, extent) in enumerate(records_now):
                cluster_all[key].append(
                    (iframe, times[iframe], icluster, area, radius, extent)
                )
            cluster_summary[key][iframe, 0] = len(records_now)
            if records_now:
                radii = np.asarray([item[1] for item in records_now])
                cluster_summary[key][iframe, 1] = radii.mean()
                cluster_summary[key][iframe, 2] = radii.max()

        peak_indices = select_periodic_peaks(m2d, sigma, max(1, args.min_peak_distance))
        points_xy = np.c_[peak_indices[:, 1] * dx, peak_indices[:, 0] * dy]
        areas, polygons = periodic_voronoi(points_xy, lx, ly)
        for icore, (point, area) in enumerate(zip(points_xy, areas)):
            voronoi_all.append(
                (iframe, times[iframe], icore, point[0], point[1], area, math.sqrt(area) if area > 0 else math.nan)
            )

        if iframe == nt - 1:
            last_m2d = m2d
            last_points_xy = points_xy
            last_areas = areas
            last_polygons = polygons
        print(f"[{iframe + 1:4d}/{nt}] t={times[iframe]:.6g}  A_m={a_m[iframe]:.6e}  L_int={l_integral[iframe]:.6g}", flush=True)

    cluster_dtype = np.dtype(
        [
            ("frame", "i4"),
            ("time", "f8"),
            ("cluster", "i4"),
            ("area", "f8"),
            ("equivalent_radius", "f8"),
            ("maximum_periodic_extent", "f8"),
        ]
    )
    voronoi_dtype = np.dtype(
        [
            ("frame", "i4"),
            ("time", "f8"),
            ("core", "i4"),
            ("x", "f8"),
            ("y", "f8"),
            ("area", "f8"),
            ("sqrt_area", "f8"),
        ]
    )

    np.save(output / "time.npy", times)
    np.save(output / "z.npy", z)
    np.save(output / "MSE_variance.npy", a_m)
    np.save(output / "Var_m_z_time.npy", var_z_time)
    np.save(output / "L_peak.npy", l_peak)
    np.save(output / "L_integral.npy", l_integral)
    np.save(output / "MSE_spectrum.npy", np.asarray(spectra))
    np.save(output / "MSE_spectrum_k.npy", k)
    np.save(output / "MSE_correlation.npy", np.asarray(correlations))
    np.save(output / "MSE_correlation_r.npy", r)
    np.save(output / "cluster_radius_threshold0.npy", np.asarray(cluster_all["threshold0"], dtype=cluster_dtype))
    np.save(output / "cluster_radius_threshold_sigma.npy", np.asarray(cluster_all["threshold_sigma"], dtype=cluster_dtype))
    np.save(output / "cluster_summary_threshold0.npy", cluster_summary["threshold0"])
    np.save(output / "cluster_summary_threshold_sigma.npy", cluster_summary["threshold_sigma"])
    np.save(output / "voronoi_area.npy", np.asarray(voronoi_all, dtype=voronoi_dtype))

    final_count = min(10, nt)
    final = slice(nt - final_count, nt)
    summary = {
        "case_root": str(case_root),
        "output": str(output),
        "gamma": args.gamma,
        "definition": "m=b+gamma*q; mprime=m-mean_xy(m); frame-first diagnostics",
        "vertical_mode": args.vertical_mode,
        "requested_height": args.height if args.vertical_mode == "height" else None,
        "actual_height": float(z[height_index]) if args.vertical_mode == "height" else None,
        "n_frames": nt,
        "stride": args.stride,
        "shape": list(shape),
        "domain": {"Lx": lx, "Ly": ly, "z_min": float(z.min()), "z_max": float(z.max())},
        "final_average_frames": final_count,
        "A_m_final": float(np.nanmean(a_m[final])),
        "L_peak_final": float(np.nanmean(l_peak[final])),
        "L_integral_final": float(np.nanmean(l_integral[final])),
        "cluster_mean_radius_threshold0_final": float(np.nanmean(cluster_summary["threshold0"][final, 1])),
        "cluster_mean_radius_threshold_sigma_final": float(np.nanmean(cluster_summary["threshold_sigma"][final, 1])),
        "voronoi_mean_sqrt_area_last": float(np.nanmean(np.sqrt(last_areas))) if last_areas.size else math.nan,
        "voronoi_variance_sqrt_area_last": float(np.nanvar(np.sqrt(last_areas))) if last_areas.size else math.nan,
        "first_field": str(frames[0].path),
        "last_field": str(frames[-1].path),
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=True) + "\n")
    (output / "frames_used.txt").write_text(
        "\n".join(f"{time:.12g}\t{record.path}" for time, record in zip(times, frames)) + "\n"
    )

    if not args.no_plots:
        plot_outputs(
            output,
            times,
            z,
            a_m,
            var_z_time,
            l_peak,
            l_integral,
            cluster_summary,
            last_m2d,
            x,
            y,
            last_points_xy,
            last_areas,
            last_polygons,
        )
    print(json.dumps(summary, indent=2, allow_nan=True))


if __name__ == "__main__":
    main()
