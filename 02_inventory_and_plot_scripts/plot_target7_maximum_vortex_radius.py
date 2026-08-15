from __future__ import annotations

import csv
import importlib.util
import json
import math
import subprocess
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
PROJECT_ROOT = Path(r"E:\moist RB\rotating_case_inventory")
SOURCE_SCRIPT = HERE / "plot_highres_energy_mprime_scales.py"
OUTPUT_DIR = (
    PROJECT_ROOT
    / "04_outputs_and_figures"
    / f"target7_latest_program_timeseries_{datetime.now():%Y%m%d}"
)
TARGET_EKS = (1.5e-4, 2.0e-4, 3.0e-3, 5.0e-3, 7.0e-3, 1.0e-2, 3.0e-2)

REMOTE_READER = r"""
import json, math
from pathlib import Path

payload = json.loads(r'''__PAYLOAD__''')
out = {}
for case in payload:
    merged = {}
    for run_path in case['source_run_paths']:
        path = Path(run_path) / 'diagnostics' / 'vortex' / 'maximum_radius.dat'
        if not path.exists() or path.stat().st_size == 0:
            continue
        with path.open('r', errors='ignore') as handle:
            for line in handle:
                parts = line.split()
                if len(parts) < 5:
                    continue
                try:
                    row = [float(parts[0]), int(parts[1]), float(parts[2]),
                           int(parts[3]), float(parts[4]), str(path)]
                except Exception:
                    continue
                if all(math.isfinite(value) for value in (row[0], row[2], row[4])):
                    merged[round(row[0], 8)] = row
    out[case['key']] = [merged[key] for key in sorted(merged)]
print(json.dumps(out, ensure_ascii=False))
"""


def load_module():
    spec = importlib.util.spec_from_file_location("latest_plot_source", SOURCE_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load {SOURCE_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def select_records(source) -> list[dict]:
    selected, _ = source.choose_highest_resolution(source.extract_remote())
    records = [
        record
        for record in selected
        if record.get("Ek") is not None
        and math.isclose(float(record["Ra"]), 8.0e6, rel_tol=1.0e-10)
        and any(
            math.isclose(float(record["Ek"]), target, rel_tol=1.0e-8)
            for target in TARGET_EKS
        )
    ]
    records.sort(key=lambda record: float(record["Ek"]))
    if len(records) != len(TARGET_EKS):
        raise RuntimeError(f"Expected 7 target cases, found {len(records)}")
    return records


def read_remote(source, records: list[dict]) -> None:
    payload = [
        {
            "key": f"{float(record['Ek']):.12g}",
            "source_run_paths": record.get("source_run_paths", [record["run_path"]]),
        }
        for record in records
    ]
    code = REMOTE_READER.replace("__PAYLOAD__", json.dumps(payload))
    process = subprocess.run(
        ["ssh", source.REMOTE_HOST, "python3", "-"],
        input=code,
        text=True,
        capture_output=True,
        encoding="utf-8",
        errors="replace",
        check=False,
    )
    if process.returncode != 0:
        raise RuntimeError(process.stderr)
    result = json.loads(process.stdout)
    for record in records:
        key = f"{float(record['Ek']):.12g}"
        record["maximum_radius"] = result.get(key, [])


def segments(time: np.ndarray) -> list[np.ndarray]:
    if len(time) < 2:
        return [np.arange(len(time), dtype=int)]
    dt = float(np.median(np.diff(time)))
    breaks = np.flatnonzero(np.diff(time) > 10.0 * dt) + 1 if dt > 0 else []
    return [chunk for chunk in np.split(np.arange(len(time)), breaks) if len(chunk)]


def plot_radius(source, records: list[dict], t_start: float | None, stem: str) -> list[Path]:
    figure = plt.figure(figsize=(16.5, 5.84), facecolor="white")
    axes = [
        figure.add_axes([0.055, 0.260, 0.312, 0.705]),
        figure.add_axes([0.425, 0.260, 0.312, 0.705]),
    ]
    for axis in axes:
        source.style_axis(axis)
    colors = source.case_colors(records)

    for record in records:
        rows = np.asarray(record["maximum_radius"], dtype=object)
        if not len(rows):
            continue
        time = rows[:, 0].astype(float)
        positive = rows[:, 2].astype(float)
        core = rows[:, 4].astype(float)
        mask = np.ones(len(time), dtype=bool) if t_start is None else time >= t_start
        time, positive, core = time[mask], positive[mask], core[mask]
        if not len(time):
            continue
        for number, segment in enumerate(segments(time)):
            axes[0].plot(
                time[segment], positive[segment], color=colors[id(record)], lw=source.LINE_WIDTH
            )
            axes[1].plot(
                time[segment], core[segment], color=colors[id(record)], lw=source.LINE_WIDTH
            )

    axes[0].set_xlabel(r"$t$")
    axes[1].set_xlabel(r"$t$")
    axes[0].set_ylabel(r"$R_{\max}^{\,m'_{2D}>0}$")
    axes[1].set_ylabel(r"$R_{\max}^{\,m'_{2D}>\sigma_m}$")
    axes[0].set_ylim(bottom=0.0)
    axes[1].set_ylim(bottom=0.0)
    if t_start is not None:
        end = max(float(np.asarray(record["maximum_radius"], dtype=object)[-1, 0]) for record in records)
        axes[0].set_xlim(t_start, end)
        axes[1].set_xlim(t_start, end)
    source.add_case_legend(figure, records, colors, (0.770, 0.52))
    return source.save_figure(figure, OUTPUT_DIR, stem)


def write_csv(records: list[dict]) -> Path:
    path = OUTPUT_DIR / "maximum_mse_vortex_radius_timeseries.csv"
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "Ek",
                "AR",
                "Nx",
                "Ny",
                "Nz",
                "time",
                "N_positive",
                "Rmax_positive",
                "N_core",
                "Rmax_core",
                "remote_source_file",
            ]
        )
        for record in records:
            for row in record["maximum_radius"]:
                writer.writerow(
                    [
                        record["Ek"],
                        record["AR"],
                        record["Nx"],
                        record["Ny"],
                        record["Nz"],
                        *row,
                    ]
                )
    return path


def main() -> None:
    source = load_module()
    source.configure_style()
    records = select_records(source)
    read_remote(source, records)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    csv_path = write_csv(records)
    figures = []
    figures.extend(
        plot_radius(
            source,
            records,
            None,
            "Ra8e6_maximum_mse_vortex_radius_timeseries_full",
        )
    )
    figures.extend(
        plot_radius(
            source,
            records,
            100.0,
            "Ra8e6_maximum_mse_vortex_radius_timeseries_tge100",
        )
    )
    manifest = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "remote_host": source.REMOTE_HOST,
        "target_Ek": list(TARGET_EKS),
        "definition": {
            "mprime_2D": "stretched-grid weighted vertical mean of m-<m>xy(z)",
            "positive_mask": "mprime_2D > 0",
            "core_mask": "mprime_2D > sigma(mprime_2D)",
            "equivalent_radius": "R_i=sqrt(A_i/pi)",
            "plotted_quantity": "Rmax=max_i R_i after the optional bou.in criterion R_i>2*vortex_lc",
        },
        "important_scope": "This is the largest MSE-anomaly connected-structure radius from mse_vortex_diagnostics, not a vertical-vorticity or Q-criterion dynamical vortex radius.",
        "known_gap": "Ek=7e-3 has no maximum_radius.dat history for 100<t<1000; the line is intentionally broken and no interpolation is used.",
        "csv": str(csv_path),
        "figures": [str(path) for path in figures],
    }
    (OUTPUT_DIR / "maximum_mse_vortex_radius_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps({"output_dir": str(OUTPUT_DIR), **manifest}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
