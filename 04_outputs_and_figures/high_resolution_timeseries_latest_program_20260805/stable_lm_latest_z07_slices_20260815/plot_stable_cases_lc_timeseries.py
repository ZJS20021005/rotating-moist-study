from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


BASE = Path(r"E:\moist RB\rotating_case_inventory\04_outputs_and_figures\high_resolution_timeseries_latest_program_20260805")
OUT = BASE / "stable_lm_latest_z07_slices_20260815" / "lc_timeseries"
LONG_CSV = BASE / "high_resolution_timeseries_long.csv"
OUT.mkdir(parents=True, exist_ok=True)

CASES = {
    "Ek1e-3": (1.0e-3, Path(r"G:\moist convection\Ek1e-3")),
    "Ek3e-3": (3.0e-3, Path(r"G:\moist convection\Ek3e-3")),
    "Ek5e-3": (5.0e-3, Path(r"G:\moist convection\Ek5e-3")),
    "Ek7e-3": (7.0e-3, Path(r"G:\moist convection\Ek7e-3")),
    "Ek1e-2": (1.0e-2, Path(r"G:\moist convection\Ek1e-2")),
}

COLORS = {
    "Ek1e-3": (0.74, 0.14, 0.18),
    "Ek3e-3": (0.93, 0.32, 0.23),
    "Ek5e-3": (0.96, 0.58, 0.19),
    "Ek7e-3": (0.16, 0.56, 0.80),
    "Ek1e-2": (0.11, 0.44, 0.71),
}


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
            "legend.fontsize": 10,
            "figure.dpi": 180,
            "savefig.dpi": 300,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


def valid(path: Path) -> bool:
    text = str(path).lower()
    return "qbot0p5_qtop0p004978" in text and "\\run\\run\\" not in text


def read_convective(path: Path) -> list[tuple[float, float]]:
    rows: list[tuple[float, float]] = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            try:
                values = [float(item.replace("D", "E").replace("d", "e")) for item in line.split()]
            except ValueError:
                continue
            if len(values) >= 3 and np.all(np.isfinite(values[:3])) and values[1] > 0.0:
                rows.append((values[0], 2.0 * math.pi * values[1]))
    return rows


def load_merged() -> tuple[dict[str, np.ndarray], list[dict[str, object]]]:
    old = pd.read_csv(LONG_CSV)
    old = old[(old["Ra"] == 8.0e6) & (old["metric"] == "convective_mid")]
    merged: dict[str, np.ndarray] = {}
    provenance: list[dict[str, object]] = []

    for label, (ek, root) in CASES.items():
        values: dict[float, tuple[float, str]] = {}
        historical = old[np.isclose(old["Ek"].to_numpy(dtype=float), ek, rtol=0.0, atol=1.0e-12)]
        for row in historical.itertuples(index=False):
            values[round(float(row.time), 8)] = (float(row.value), str(LONG_CSV))

        raw_files = [path for path in root.rglob("convective_scale.dat") if valid(path)]
        raw_files.sort(key=lambda path: (path.stat().st_mtime, len(path.parts), str(path)))
        for path in raw_files:
            for time, value in read_convective(path):
                values[round(time, 8)] = (value, str(path))

        if not values:
            raise RuntimeError(f"No convective-scale history found for {label}")
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
                "latest_two_pi_lc_H": rows[-1, 1],
                "median_dt": dt,
                "gap_count_above_1p5dt": gap_count,
                "historical_rows_before_override": len(historical),
                "local_raw_files": len(raw_files),
                "local_root": str(root),
            }
        )
    return merged, provenance


def style_axis(ax: plt.Axes) -> None:
    ax.set_box_aspect(5.2 / 6.5)
    for spine in ax.spines.values():
        spine.set_linewidth(4.5)
    ax.tick_params(which="major", direction="in", length=12, width=1.2, top=True, right=True)
    ax.tick_params(which="minor", direction="in", length=6, width=1.0, top=True, right=True)
    ax.minorticks_on()
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$2\pi l_c(z\simeq0.5)/H$")


def new_figure() -> tuple[plt.Figure, plt.Axes]:
    fig = plt.figure(figsize=(6.5, 5.84), facecolor="white")
    ax = fig.add_axes([0.186, 0.260, 0.792, 0.705])
    style_axis(ax)
    return fig, ax


def ek_label(ek: float) -> str:
    exponent = int(math.floor(math.log10(ek)))
    coefficient = ek / 10.0**exponent
    if math.isclose(coefficient, 1.0, rel_tol=0.0, abs_tol=1.0e-10):
        return rf"$Ek=10^{{{exponent}}}$"
    return rf"$Ek={coefficient:g}\times10^{{{exponent}}}$"


def save_combined(series: dict[str, np.ndarray]) -> None:
    fig, ax = new_figure()
    for label in CASES:
        data = series[label]
        ek = CASES[label][0]
        ax.plot(data[:, 0], data[:, 1], color=COLORS[label], lw=3.5, label=ek_label(ek))
    ax.set_xlim(left=0.0)
    ax.set_ylim(bottom=0.0)
    ax.legend(frameon=False, ncol=2, loc="lower right", handlelength=2.1, columnspacing=1.0)
    fig.savefig(OUT / "Ra8e6_stable_Lm_cases_two_pi_lc_z05_timeseries_combined.png")
    fig.savefig(OUT / "Ra8e6_stable_Lm_cases_two_pi_lc_z05_timeseries_combined.pdf")
    plt.close(fig)


def save_individual(series: dict[str, np.ndarray]) -> None:
    for label, (ek, _) in CASES.items():
        data = series[label]
        fig, ax = new_figure()
        ax.plot(data[:, 0], data[:, 1], color=(0.11, 0.44, 0.71), lw=3.5, label=ek_label(ek))
        ax.set_xlim(left=0.0)
        ax.set_ylim(bottom=0.0)
        ax.legend(frameon=False, loc="lower right", handlelength=2.1)
        stem = OUT / f"Ra8e6_{label}_two_pi_lc_z05_timeseries"
        fig.savefig(stem.with_suffix(".png"))
        fig.savefig(stem.with_suffix(".pdf"))
        plt.close(fig)


def save_tables(series: dict[str, np.ndarray], provenance: list[dict[str, object]]) -> None:
    with (OUT / "Ra8e6_stable_Lm_cases_two_pi_lc_z05_timeseries.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(["case", "Ek", "time", "two_pi_lc_z05_H"])
        for label, data in series.items():
            ek = CASES[label][0]
            for time, value in data:
                writer.writerow([label, ek, time, value])
    with (OUT / "Ra8e6_stable_Lm_cases_two_pi_lc_z05_provenance.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(provenance[0]))
        writer.writeheader()
        writer.writerows(provenance)


def main() -> None:
    configure_style()
    series, provenance = load_merged()
    save_combined(series)
    save_individual(series)
    save_tables(series, provenance)
    for row in provenance:
        print(row)
    print(OUT)


if __name__ == "__main__":
    main()
