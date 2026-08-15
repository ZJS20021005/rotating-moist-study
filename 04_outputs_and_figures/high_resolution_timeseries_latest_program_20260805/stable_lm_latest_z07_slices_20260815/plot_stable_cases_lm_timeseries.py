from __future__ import annotations

import csv
import importlib.util
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("lc_plot", HERE / "plot_stable_cases_lc_timeseries.py")
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("Cannot load shared plotting definitions")
SHARED = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SHARED)

BASE = SHARED.BASE
LONG_CSV = SHARED.LONG_CSV
CASES = SHARED.CASES
COLORS = SHARED.COLORS
OUT = HERE / "lm_timeseries"
OUT.mkdir(parents=True, exist_ok=True)


def read_moist_scale(path: Path) -> list[tuple[float, float]]:
    rows: list[tuple[float, float]] = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            try:
                values = [float(item.replace("D", "E").replace("d", "e")) for item in line.split()]
            except ValueError:
                continue
            if len(values) >= 4 and np.all(np.isfinite(values[:4])) and values[1] > 0.0:
                rows.append((values[0], 2.0 * math.pi * values[1]))
    return rows


def load_merged() -> tuple[dict[str, np.ndarray], list[dict[str, object]]]:
    old = pd.read_csv(LONG_CSV)
    old = old[(old["Ra"] == 8.0e6) & (old["metric"] == "mse_spectral")]
    merged: dict[str, np.ndarray] = {}
    provenance: list[dict[str, object]] = []

    for label, (ek, root) in CASES.items():
        values: dict[float, tuple[float, str]] = {}
        historical = old[np.isclose(old["Ek"].to_numpy(dtype=float), ek, rtol=0.0, atol=1.0e-12)]
        for row in historical.itertuples(index=False):
            values[round(float(row.time), 8)] = (float(row.value), str(LONG_CSV))

        raw_files = [path for path in root.rglob("moist_integral_scale.dat") if SHARED.valid(path)]
        raw_files.sort(key=lambda path: (path.stat().st_mtime, len(path.parts), str(path)))
        for path in raw_files:
            for time, value in read_moist_scale(path):
                values[round(time, 8)] = (value, str(path))

        if not values:
            raise RuntimeError(f"No moist spectral scale found for {label}")
        rows = np.asarray([[time, values[time][0]] for time in sorted(values)], dtype=float)
        merged[label] = rows

        dt = float(np.median(np.diff(rows[:, 0]))) if len(rows) > 1 else math.nan
        gaps = np.diff(rows[:, 0]) if len(rows) > 1 else np.asarray([])
        gap_count = int(np.sum(gaps > 1.51 * dt)) if math.isfinite(dt) and dt > 0 else 0
        provenance.append(
            {
                "case": label,
                "Ek": ek,
                "n_rows": len(rows),
                "time_min": rows[0, 0],
                "time_max": rows[-1, 0],
                "latest_two_pi_Lm_H": rows[-1, 1],
                "median_dt": dt,
                "gap_count_above_1p5dt": gap_count,
                "historical_rows_before_override": len(historical),
                "local_raw_files": len(raw_files),
                "local_root": str(root),
            }
        )
    return merged, provenance


def new_figure() -> tuple[plt.Figure, plt.Axes]:
    fig = plt.figure(figsize=(6.5, 5.84), facecolor="white")
    ax = fig.add_axes([0.186, 0.260, 0.792, 0.705])
    SHARED.style_axis(ax)
    ax.set_ylabel(r"$2\pi L_m/H$")
    return fig, ax


def save_combined(series: dict[str, np.ndarray]) -> None:
    fig, ax = new_figure()
    for label, (ek, _) in CASES.items():
        data = series[label]
        ax.plot(data[:, 0], data[:, 1], color=COLORS[label], lw=3.5, label=SHARED.ek_label(ek))
    ax.set_xlim(left=0.0)
    ax.set_ylim(bottom=0.0)
    ax.legend(frameon=False, ncol=2, loc="lower right", handlelength=2.1, columnspacing=1.0)
    fig.savefig(OUT / "Ra8e6_stable_Lm_cases_two_pi_Lm_timeseries_combined.png")
    fig.savefig(OUT / "Ra8e6_stable_Lm_cases_two_pi_Lm_timeseries_combined.pdf")
    plt.close(fig)


def save_individual(series: dict[str, np.ndarray]) -> None:
    for label, (ek, _) in CASES.items():
        data = series[label]
        fig, ax = new_figure()
        ax.plot(data[:, 0], data[:, 1], color=(0.11, 0.44, 0.71), lw=3.5, label=SHARED.ek_label(ek))
        ax.set_xlim(left=0.0)
        ax.set_ylim(bottom=0.0)
        ax.legend(frameon=False, loc="lower right", handlelength=2.1)
        stem = OUT / f"Ra8e6_{label}_two_pi_Lm_timeseries"
        fig.savefig(stem.with_suffix(".png"))
        fig.savefig(stem.with_suffix(".pdf"))
        plt.close(fig)


def save_tables(series: dict[str, np.ndarray], provenance: list[dict[str, object]]) -> None:
    with (OUT / "Ra8e6_stable_Lm_cases_two_pi_Lm_timeseries.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(["case", "Ek", "time", "two_pi_Lm_H"])
        for label, data in series.items():
            ek = CASES[label][0]
            for time, value in data:
                writer.writerow([label, ek, time, value])
    with (OUT / "Ra8e6_stable_Lm_cases_two_pi_Lm_provenance.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(provenance[0]))
        writer.writeheader()
        writer.writerows(provenance)


def main() -> None:
    SHARED.configure_style()
    series, provenance = load_merged()
    save_combined(series)
    save_individual(series)
    save_tables(series, provenance)
    for row in provenance:
        print(row)
    print(OUT)


if __name__ == "__main__":
    main()
