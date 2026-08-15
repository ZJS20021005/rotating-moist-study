#!/usr/bin/env python3
"""Export latest-10-frame periodic MSE Voronoi diagnostics for mature AR10 cases."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
import re
import sys
from pathlib import Path

import h5py
import numpy as np


def load_helpers(path: Path):
    spec = importlib.util.spec_from_file_location("aggregation_helpers", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def canonical_case(case: Path) -> Path:
    if re.fullmatch(r"conti\d*", case.name, flags=re.IGNORECASE):
        return case.parent
    return case


def read_m2d(frame: Path, gamma: float, weights: np.ndarray) -> np.ndarray:
    with h5py.File(frame, "r") as source:
        b = source["DSAL_me"]
        q = source["QVAP_me"]
        nz, ny, nx = b.shape
        if nz != weights.size:
            raise ValueError(f"{frame}: nz={nz}, but weights={weights.size}")
        result = np.zeros((ny, nx), dtype=np.float64)
        for iz, weight in enumerate(weights):
            layer = np.asarray(b[iz], dtype=np.float64)
            layer += gamma * np.asarray(q[iz], dtype=np.float64)
            layer -= layer.mean()
            result += weight * layer
    result -= result.mean()
    return result


def polygon_segments(polygons):
    segments = []
    for polygon in polygons:
        if polygon.size:
            segments.extend(np.stack([polygon, np.roll(polygon, -1, axis=0)], axis=1))
    if not segments:
        return np.empty((0, 2, 2), dtype=float)
    return np.asarray(segments, dtype=float)


def parse_case(canonical: Path):
    text = str(canonical)
    ra_match = re.search(r"Ra([0-9.eE+-]+)", text)
    ek_match = re.search(r"Ek([0-9.eE+-]+)AR", text)
    ra = float(ra_match.group(1)) if ra_match else math.nan
    if "norotating" in str(canonical).lower():
        ek = math.inf
        label = "norotating"
    else:
        ek = float(ek_match.group(1)) if ek_match else math.nan
        label = f"Ek{ek_match.group(1)}" if ek_match else text
    return ra, ek, label


def scientific_tag(value: float) -> str:
    exponent = int(math.floor(math.log10(value)))
    coefficient = value / 10**exponent
    return f"{coefficient:g}e{exponent}"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--helpers", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--minimum-time", type=float, default=800.0)
    parser.add_argument("--latest", type=int, default=10)
    parser.add_argument("--gamma", type=float, default=1.1)
    parser.add_argument("--min-peak-distance", type=int, default=3)
    args = parser.parse_args()

    root = args.root.resolve()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    helpers = load_helpers(args.helpers.resolve())

    selected = {}
    for movie in sorted(root.rglob("run/movie")):
        case = movie.parent.parent
        if "AR10" not in str(case):
            continue
        try:
            frames = helpers.discover_frames(case, 1, None)
        except Exception:
            continue
        if len(frames) < args.latest or not np.isfinite(frames[-1].time):
            continue
        if frames[-1].time < args.minimum_time:
            continue
        key = canonical_case(case)
        old = selected.get(key)
        if old is None or frames[-1].time > old[1][-1].time:
            selected[key] = (case, frames)

    rows = []
    for canonical in sorted(selected, key=lambda p: str(p)):
        source_case, frames_all = selected[canonical]
        frames = frames_all[-args.latest :]
        first = frames[0].path
        with h5py.File(first, "r") as source:
            shape = tuple(source["DSAL_me"].shape)
        grid_path = helpers.find_grid(first, source_case)
        x, y, z, dx, dy, dz, lx, ly = helpers.read_grid(grid_path, shape)
        weights = dz / dz.sum()

        m2d_sum = np.zeros((shape[1], shape[2]), dtype=np.float64)
        frame_core_counts = []
        frame_mean_scales = []
        frame_sigmas = []
        for record in frames:
            m2d = read_m2d(record.path, args.gamma, weights)
            m2d_sum += m2d
            sigma = float(np.sqrt(np.mean(m2d * m2d)))
            peaks = helpers.select_periodic_peaks(
                m2d, sigma, max(1, args.min_peak_distance)
            )
            points = np.c_[peaks[:, 1] * dx, peaks[:, 0] * dy]
            areas, _ = helpers.periodic_voronoi(points, lx, ly)
            valid = areas[np.isfinite(areas) & (areas > 0.0)]
            frame_core_counts.append(len(points))
            frame_mean_scales.append(
                float(np.mean(np.sqrt(valid))) if valid.size else math.nan
            )
            frame_sigmas.append(sigma)

        mean_m2d = m2d_sum / len(frames)
        mean_m2d -= mean_m2d.mean()
        mean_sigma = float(np.sqrt(np.mean(mean_m2d * mean_m2d)))
        peak_indices = helpers.select_periodic_peaks(
            mean_m2d, mean_sigma, max(1, args.min_peak_distance)
        )
        points_xy = np.c_[peak_indices[:, 1] * dx, peak_indices[:, 0] * dy]
        areas, polygons = helpers.periodic_voronoi(points_xy, lx, ly)
        segments = polygon_segments(polygons)
        valid = areas[np.isfinite(areas) & (areas > 0.0)]

        ra, ek, ek_label = parse_case(canonical)
        case_id = f"Ra{scientific_tag(ra)}_{ek_label}_AR10"
        times = np.asarray([record.time for record in frames], dtype=float)
        np.savez_compressed(
            output / f"{case_id}_latest10_mean_voronoi.npz",
            mean_m2d=mean_m2d,
            normalized_m2d=mean_m2d / mean_sigma,
            points_xy=points_xy,
            areas=areas,
            segments=segments,
            x=x,
            y=y,
            z=z,
            times=times,
            lx=lx,
            ly=ly,
            sigma=mean_sigma,
            frame_core_counts=np.asarray(frame_core_counts),
            frame_mean_scales=np.asarray(frame_mean_scales),
            frame_sigmas=np.asarray(frame_sigmas),
        )
        rows.append(
            {
                "case_id": case_id,
                "Ra": ra,
                "Ek": "norotating" if math.isinf(ek) else ek,
                "canonical_case": str(canonical),
                "source_case": str(source_case),
                "first_time": float(times[0]),
                "last_time": float(times[-1]),
                "n_frames": len(times),
                "mean_field_sigma": mean_sigma,
                "mean_field_core_count": len(points_xy),
                "mean_field_mean_sqrt_area": float(np.mean(np.sqrt(valid))) if valid.size else math.nan,
                "mean_field_max_sqrt_area": float(np.max(np.sqrt(valid))) if valid.size else math.nan,
                "frame_mean_core_count": float(np.mean(frame_core_counts)),
                "frame_mean_sqrt_area": float(np.nanmean(frame_mean_scales)),
                "domain_Lx": lx,
                "domain_Ly": ly,
            }
        )
        print(
            f"{case_id}: t={times[0]:g}-{times[-1]:g}, "
            f"cores(mean field)={len(points_xy)}, frame mean cores={np.mean(frame_core_counts):.2f}",
            flush=True,
        )

    if not rows:
        raise RuntimeError("No AR10 cases reached the requested minimum time")
    with (output / "voronoi_latest10_metadata.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (output / "voronoi_latest10_metadata.json").write_text(
        json.dumps(rows, indent=2, allow_nan=True) + "\n"
    )


if __name__ == "__main__":
    main()
