#!/usr/bin/env python3
"""Compute reduced force and kinetic-energy budgets from remote movie fields.

This script is intended to run on the cluster.  It avoids h5py by using
h5dump's binary output and reads one latest common field/horizontal-velocity
snapshot per selected latest-program case.
"""

from __future__ import print_function

import argparse
import csv
import json
import math
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np


TIME_RE = re.compile(r'<Time\s+Value="([+\-0-9.eEdD]+)"')
FRAME_RE = re.compile(r"field(\d+)\.xmf$")


def dump_dataset(h5_path, dataset, output_path):
    subprocess.check_call(
        ["h5dump", "-d", "/" + dataset, "-b", "LE", "-o", str(output_path), str(h5_path)],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    return np.fromfile(str(output_path), dtype="<f8")


def read_time(xmf_path):
    try:
        text = xmf_path.read_text(errors="ignore")
    except Exception:
        return None
    match = TIME_RE.search(text)
    if not match:
        return None
    try:
        return float(match.group(1).replace("D", "e").replace("d", "e"))
    except Exception:
        return None


def latest_snapshot(run_paths):
    candidates = all_snapshots(run_paths)
    if not candidates:
        return None
    return max(candidates, key=lambda item: (item[0], item[1], item[2]))


def all_snapshots(run_paths):
    candidates = []
    for priority, run_text in enumerate(run_paths):
        movie = Path(run_text) / "movie"
        if not movie.is_dir():
            continue
        for xmf in movie.glob("field[0-9]*.xmf"):
            match = FRAME_RE.search(xmf.name)
            if not match:
                continue
            frame = match.group(1)
            field = movie / ("field" + frame + ".h5")
            hvel = movie / ("horizontal_velocity" + frame + ".h5")
            coord = movie / "cordin_info.h5"
            if not (field.exists() and hvel.exists() and coord.exists()):
                continue
            time = read_time(xmf)
            if time is None:
                # Production movies are written every 10 nondimensional time
                # units and use a continuous frame index across continuations.
                time = 10.0 * float(frame)
            # Some production XMF files omit the optional Time element.  The
            # frame index still provides a reliable ordering within a run.
            rank_time = time if time is not None and math.isfinite(time) else float(frame)
            candidates.append((rank_time, priority, int(frame), time, field, hvel, coord))
    # Deduplicate repeated frame/time pairs across continuation directories,
    # preferring the later run path.
    selected = {}
    for item in candidates:
        key = (item[0], item[2])
        if key not in selected or item[1] > selected[key][1]:
            selected[key] = item
    return sorted(selected.values(), key=lambda item: (item[0], item[2]))


def periodic_first(array, spacing, axis):
    return (np.roll(array, -1, axis=axis) - np.roll(array, 1, axis=axis)) / (2.0 * spacing)


def periodic_second(array, spacing, axis):
    return (np.roll(array, -1, axis=axis) - 2.0 * array + np.roll(array, 1, axis=axis)) / spacing**2


def vertical_first(array, z):
    return np.gradient(array, z, axis=0, edge_order=2)


def vertical_second(array, z):
    return np.gradient(np.gradient(array, z, axis=0, edge_order=2), z, axis=0, edge_order=2)


def plane_variance(component):
    anomaly = component - component.mean(axis=(1, 2), keepdims=True)
    return np.mean(anomaly * anomaly, axis=(1, 2))


def horizontal_plane_rms(components, index):
    total = 0.0
    for component in components[:2]:
        plane = component[index]
        anomaly = plane - np.mean(plane)
        total += float(np.mean(anomaly * anomaly))
    return math.sqrt(total)


def pressure_plane_diagnostics(field_path, frame, x, y, z, components, case_work):
    run = field_path.parent.parent
    pressure_path = run / "flowmov" / ("frame_{:05d}_zcut.h5".format(frame))
    if not pressure_path.exists():
        return {}
    pressure_raw = dump_dataset(pressure_path, "Pr", case_work / "pressure.bin")
    if pressure_raw.size != len(y) * len(x):
        return {}
    pressure = pressure_raw.reshape((len(y), len(x)))
    px = -periodic_first(pressure, float(np.mean(np.diff(x))), 1)
    py = -periodic_first(pressure, float(np.mean(np.diff(y))), 0)

    # PLANEMOVIE uses the one-based mov_zcut_k directly.  Current cases use
    # mov_zcut_k=5, corresponding to Python index 4 (z ~= 0.0483).
    bou = run / "bou.in"
    z_index = 4
    if bou.exists():
        lines = bou.read_text(errors="ignore").splitlines()
        for line_index, line in enumerate(lines):
            if "mov_zcut_k" in line:
                for values in lines[line_index + 1:]:
                    clean = values.split("!")[0].strip()
                    if clean:
                        try: z_index = max(0, min(len(z) - 1, int(clean.split()[0]) - 1))
                        except Exception: pass
                        break
                break

    pressure_force = math.sqrt(float(np.mean(px * px + py * py)))
    coriolis_x = components["coriolis"][0][z_index]
    coriolis_y = components["coriolis"][1][z_index]
    cp_residual = math.sqrt(float(np.mean((px + coriolis_x) ** 2 + (py + coriolis_y) ** 2)))
    return {
        "pressure_plane": pressure_force,
        "coriolis_plane": horizontal_plane_rms(components["coriolis"], z_index),
        "inertia_plane": horizontal_plane_rms(components["inertia"], z_index),
        "viscous_plane": horizontal_plane_rms(components["viscous"], z_index),
        "cp_residual_plane": cp_residual,
        "cp_residual_over_coriolis": cp_residual / max(horizontal_plane_rms(components["coriolis"], z_index), 1e-30),
        "pressure_plane_z": float(z[z_index]),
    }


def z_average(profile, z, mask):
    selected_z = z[mask]
    selected = profile[mask]
    if len(selected_z) < 2:
        return float(np.mean(selected))
    return float(np.trapz(selected, selected_z) / (selected_z[-1] - selected_z[0]))


def classify(forces):
    inertia, coriolis, viscous, buoyancy = (
        forces["inertia"], forces["coriolis"], forces["viscous"], forces["buoyancy"]
    )
    if coriolis <= 1.0e-15:
        return "nonrotating: inertia-buoyancy-viscous", math.inf
    local_ro = inertia / coriolis
    if local_ro >= 1.0:
        if coriolis / max(inertia, 1.0e-30) < 0.1:
            return "buoyancy/inertia dominated; nonrotating-like", local_ro
        return "rotation-affected; inertia-pressure leading expected", local_ro
    if inertia < 0.5 * max(viscous, buoyancy):
        return "geostrophic leading; VAC-like ageostrophic balance", local_ro
    if viscous < 0.5 * min(inertia, buoyancy) and 0.25 <= inertia / max(buoyancy, 1.0e-30) <= 4.0:
        return "geostrophic leading; CIA-like ageostrophic balance", local_ro
    return "geostrophic leading; CIA/VAC mixed ageostrophic balance", local_ro


def process_case(case, work_root, snapshot=None):
    snapshot = snapshot or latest_snapshot(case["run_paths"])
    if snapshot is None:
        raise RuntimeError("no common latest movie snapshot")
    _, _, frame, time, field_path, hvel_path, coord_path = snapshot
    case_work = Path(tempfile.mkdtemp(prefix="force_", dir=str(work_root)))
    try:
        x = dump_dataset(coord_path, "x", case_work / "x.bin")
        y = dump_dataset(coord_path, "y", case_work / "y.bin")
        z = dump_dataset(coord_path, "z", case_work / "z.bin")
        shape = (len(z), len(y), len(x))
        u = dump_dataset(hvel_path, "VX_me", case_work / "u.bin").reshape(shape)
        v = dump_dataset(hvel_path, "VY_me", case_work / "v.bin").reshape(shape)
        w = dump_dataset(field_path, "VZ_me", case_work / "w.bin").reshape(shape)
        b = dump_dataset(field_path, "DSAL_me", case_work / "b.bin").reshape(shape)

        dx = float(np.mean(np.diff(x)))
        dy = float(np.mean(np.diff(y)))
        nu = math.sqrt(float(case["Pr"]) / float(case["Ra"]))
        inv_ro = float(case.get("invRo") or 0.0)

        derivatives = {}
        for name, array in (("u", u), ("v", v), ("w", w)):
            derivatives[(name, "x")] = periodic_first(array, dx, 2)
            derivatives[(name, "y")] = periodic_first(array, dy, 1)
            derivatives[(name, "z")] = vertical_first(array, z)

        inertia_components = []
        viscous_components = []
        for name, array in (("u", u), ("v", v), ("w", w)):
            inertia_components.append(
                -(u * derivatives[(name, "x")] + v * derivatives[(name, "y")] + w * derivatives[(name, "z")])
            )
            viscous_components.append(
                nu * (
                    periodic_second(array, dx, 2)
                    + periodic_second(array, dy, 1)
                    + vertical_second(array, z)
                )
            )

        coriolis_components = (inv_ro * v, -inv_ro * u, np.zeros_like(w))
        b_anomaly = b - b.mean(axis=(1, 2), keepdims=True)
        buoyancy_components = (np.zeros_like(w), np.zeros_like(w), b_anomaly)
        component_groups = {
            "inertia": inertia_components, "coriolis": coriolis_components,
            "viscous": viscous_components, "buoyancy": buoyancy_components,
        }

        profiles = {"z": z}
        for label, components in component_groups.items():
            total_variance = sum(plane_variance(component) for component in components)
            profiles[label] = np.sqrt(total_variance)

        strain_sum = np.zeros_like(w)
        for derivative in derivatives.values():
            strain_sum += derivative * derivative
        profiles["buoyancy_power"] = np.mean(w * b_anomaly, axis=(1, 2))
        profiles["viscous_dissipation"] = nu * np.mean(strain_sum, axis=(1, 2))

        bulk_mask = (z >= 0.2) & (z <= 0.8)
        bulk = {name: z_average(profile, z, bulk_mask) for name, profile in profiles.items() if name != "z"}
        regime, local_ro = classify(bulk)
        summary = {
            "Ra": case["Ra"], "Pr": case["Pr"], "Ek": case.get("Ek"),
            "AR": case["AR"], "Nx": len(x), "Ny": len(y), "Nz": len(z),
            "time": "" if time is None else time, "frame": frame, "field_path": str(field_path),
            "inertia_bulk": bulk["inertia"], "coriolis_bulk": bulk["coriolis"],
            "viscous_bulk": bulk["viscous"], "buoyancy_bulk": bulk["buoyancy"],
            "local_Ro_force": local_ro,
            "buoyancy_power_bulk": bulk["buoyancy_power"],
            "viscous_dissipation_bulk": bulk["viscous_dissipation"],
            "regime": regime,
        }
        summary.update(pressure_plane_diagnostics(
            field_path, frame, x, y, z, component_groups, case_work
        ))
        profile_rows = []
        for index, height in enumerate(z):
            row = {"Ra": case["Ra"], "Ek": case.get("Ek"), "time": time, "z": float(height)}
            for name, profile in profiles.items():
                if name != "z":
                    row[name] = float(profile[index])
            profile_rows.append(row)
        return summary, profile_rows
    finally:
        shutil.rmtree(str(case_work), ignore_errors=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("cases_json")
    parser.add_argument("output_dir")
    parser.add_argument("--all-snapshots", action="store_true")
    parser.add_argument("--stride", type=int, default=1)
    args = parser.parse_args()
    cases = json.loads(Path(args.cases_json).read_text())
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    work = output / "temporary_binary"
    work.mkdir(exist_ok=True)
    summaries = []
    profile_rows = []
    errors = []
    for index, case in enumerate(cases, 1):
        label = case.get("Ek_label", "unknown")
        print("[{}/{}] {}".format(index, len(cases), label), flush=True)
        try:
            snapshots = all_snapshots(case["run_paths"]) if args.all_snapshots else [latest_snapshot(case["run_paths"])]
            if args.all_snapshots and args.stride > 1 and len(snapshots) > 1:
                sampled = snapshots[::args.stride]
                if sampled[-1] != snapshots[-1]:
                    sampled.append(snapshots[-1])
                snapshots = sampled
            for snapshot in snapshots:
                summary, rows = process_case(case, work, snapshot=snapshot)
                summaries.append(summary)
                if not args.all_snapshots:
                    profile_rows.extend(rows)
        except Exception as exc:
            errors.append({"Ek": case.get("Ek"), "anchor_path": case.get("anchor_path"), "error": str(exc)})
    shutil.rmtree(str(work), ignore_errors=True)

    if summaries:
        with (output / "force_balance_bulk_summary.csv").open("w", newline="") as handle:
            fieldnames = []
            for summary in summaries:
                for field in summary:
                    if field not in fieldnames:
                        fieldnames.append(field)
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader(); writer.writerows(summaries)
    if profile_rows:
        with (output / "force_balance_profiles.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(profile_rows[0].keys()))
            writer.writeheader(); writer.writerows(profile_rows)
    (output / "force_balance_errors.json").write_text(json.dumps(errors, indent=2))
    print("completed={} errors={}".format(len(summaries), len(errors)))


if __name__ == "__main__":
    main()
