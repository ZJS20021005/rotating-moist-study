from __future__ import annotations

import csv
import importlib.util
import json
import math
import shutil
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
PROJECT_ROOT = Path(r"E:\moist RB\rotating_case_inventory")
SOURCE_SCRIPT = HERE / "plot_highres_energy_mprime_scales.py"
OUTPUT_DIR = (
    PROJECT_ROOT
    / "04_outputs_and_figures"
    / f"target7_latest_program_timeseries_{datetime.now():%Y%m%d}"
)
ARCHIVE_LONG = (
    PROJECT_ROOT
    / "04_outputs_and_figures"
    / "high_resolution_timeseries_latest_program_20260801"
    / "high_resolution_timeseries_long.csv"
)

TARGET_EKS = (1.5e-4, 2.0e-4, 3.0e-3, 5.0e-3, 7.0e-3, 1.0e-2, 3.0e-2)


def load_module():
    spec = importlib.util.spec_from_file_location("latest_plot_source", SOURCE_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load {SOURCE_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def has_gap(rows: list[list[float]]) -> bool:
    if len(rows) < 3:
        return False
    time = np.asarray(rows, dtype=float)[:, 0]
    dt = float(np.median(np.diff(time)))
    return bool(np.any(np.diff(time) > 10.0 * dt)) if dt > 0.0 else False


def same_case(row: dict[str, str], record: dict) -> bool:
    try:
        if (int(row["Nx"]), int(row["Ny"]), int(row["Nz"])) != (
            int(record["Nx"]),
            int(record["Ny"]),
            int(record["Nz"]),
        ):
            return False
        for name in ("Ra", "Ek", "AR", "beta", "qbot"):
            if not math.isclose(
                float(row[name]),
                float(record[name]),
                rel_tol=1.0e-10,
                abs_tol=1.0e-12,
            ):
                return False
        return True
    except (KeyError, TypeError, ValueError):
        return False


def load_archive(record: dict, source) -> dict[str, list[list[float]]]:
    simple = {
        "kinetic",
        "mprime",
        "mse_peak",
        "mse_spectral",
        "mse_integral",
        "convective_mid",
        "convective_zavg",
        "velocity_z075",
    }
    values: dict[str, list[list[float]]] = {name: [] for name in source.SERIES_NAMES}
    l0_parts: dict[float, dict[int, float]] = defaultdict(dict)
    if not ARCHIVE_LONG.exists():
        return values
    with ARCHIVE_LONG.open("r", newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            metric = row.get("metric", "")
            if metric not in simple and metric not in {"l0_1", "l0_2"}:
                continue
            if not same_case(row, record):
                continue
            try:
                time = float(row["time"])
                value = float(row["value"])
            except (KeyError, TypeError, ValueError):
                continue
            if metric in simple:
                values[metric].append([time, value])
            else:
                component = int(metric.rsplit("_", 1)[1])
                l0_parts[time][component] = value
    values["l0"] = [
        [time, parts[1], parts[2]]
        for time, parts in sorted(l0_parts.items())
        if 1 in parts and 2 in parts
    ]
    return values


def backfill_record(record: dict, source) -> bool:
    if not any(has_gap(record["series"].get(name, [])) for name in source.SERIES_NAMES):
        return False
    archived = load_archive(record, source)
    used = False
    for name in source.SERIES_NAMES:
        old_rows = archived.get(name, [])
        if not old_rows:
            continue
        record["series"][name] = source.deduplicate_rows(
            old_rows + record["series"].get(name, [])
        )
        used = True
    return used


def main() -> None:
    source = load_module()
    source.configure_style()
    raw = source.extract_remote()
    selected, duplicates = source.choose_highest_resolution(raw)
    selected = [
        record
        for record in selected
        if record.get("Ek") is not None
        and math.isclose(float(record["Ra"]), 8.0e6, rel_tol=1.0e-10)
        and any(
            math.isclose(float(record["Ek"]), target, rel_tol=1.0e-8)
            for target in TARGET_EKS
        )
    ]
    selected.sort(key=lambda record: float(record["Ek"]))
    if len(selected) != len(TARGET_EKS):
        raise RuntimeError(
            f"Expected {len(TARGET_EKS)} target cases but found {len(selected)}."
        )

    archive_cases = []
    for record in selected:
        if backfill_record(record, source):
            archive_cases.append(float(record["Ek"]))
        all_times = [
            row[0]
            for name in source.SERIES_NAMES
            for row in record["series"].get(name, [])
        ]
        record["time_min"] = min(all_times, default=None)
        record["time_max"] = max(all_times, default=None)

    if OUTPUT_DIR.exists():
        shutil.rmtree(OUTPUT_DIR)
    OUTPUT_DIR.mkdir(parents=True)
    source.write_tables(selected, duplicates, OUTPUT_DIR)

    colors = source.case_colors(selected)
    figures = []
    figures.extend(source.plot_energy_mprime("Ra8e6", selected, colors, OUTPUT_DIR))
    figures.extend(source.plot_l0("Ra8e6", selected, colors, OUTPUT_DIR))
    figures.extend(source.plot_mse_scales("Ra8e6", selected, colors, OUTPUT_DIR))
    figures.extend(source.plot_velocity_scales("Ra8e6", selected, colors, OUTPUT_DIR))

    manifest = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "remote_host": source.REMOTE_HOST,
        "remote_root": source.REMOTE_ROOT,
        "target_Ek": list(TARGET_EKS),
        "latest_program_only": True,
        "definitions": {
            "kinetic": "K=0.5*avgvar one-based column 4",
            "mprime_variance": "<<[m-<m>xy]^2>xy>z with stretched-grid vertical weighting",
            "l0": "avgvar one-based columns 8 and 9 at z approximately 0.5 and 0.75",
            "mse_peak": "L_peak,m=2*pi/k_peak,m",
            "mse_spectral": "2*pi*sum(E_m)/sum(k_h E_m)",
            "convective_mid": "2*pi*sum(E_w)/sum(k_h E_w) at z approximately 0.5",
            "convective_zavg": "vertical weighted mean of per-height convective scale, multiplied by 2*pi",
        },
        "archive_gap_fill": {
            "source": str(ARCHIVE_LONG),
            "cases": archive_cases,
            "rule": "used only when the current remote merged series contains a physical-time gap; current remote rows override archived duplicate times",
        },
        "selected_cases": len(selected),
        "case_time_max": {
            f"{float(record['Ek']):.8g}": record.get("time_max") for record in selected
        },
        "figures": [str(path) for path in figures],
    }
    (OUTPUT_DIR / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps({"output_dir": str(OUTPUT_DIR), **manifest}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
