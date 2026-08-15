from __future__ import annotations

import csv
import importlib.util
import json
import math
import shutil
from datetime import datetime
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
PROJECT_ROOT = Path(r"E:\moist RB\rotating_case_inventory")
SOURCE_SCRIPT = HERE / "plot_highres_energy_mprime_scales.py"
FIT_SCRIPT = (
    PROJECT_ROOT
    / "04_outputs_and_figures"
    / "Ra8e6_Ek7e-3_qbot0p5_energy_convective_scale_20260728"
    / "plot_energy_convective_scale.py"
)
OUTPUT_ROOT = PROJECT_ROOT / "04_outputs_and_figures"
OUTPUT_DIR = OUTPUT_ROOT / "linear_growth_lc_steady_latest_program_20260801"
ARCHIVE_LONG = (
    OUTPUT_ROOT
    / "high_resolution_timeseries_latest_program_20260801"
    / "high_resolution_timeseries_long.csv"
)

FIGSIZE = (6.5, 5.84)
AX_RECT = [0.186, 0.260, 0.792, 0.705]
BOX_ASPECT = 5.2 / 6.5
FRAME_WIDTH = 4.5
LINE_WIDTH = 3.5

BLUE = (0.00, 0.32, 0.90)
LIGHT_BLUE = (0.18, 0.70, 0.95)
FIT_RED = (0.88, 0.08, 0.12)
SHADE = (0.72, 0.72, 0.72)

# The last quarter of each available record is used as a consistent terminal
# diagnostic window. It is never allowed to be shorter than 150 time units.
MIN_LATE_DURATION = 150.0
STRONG_DRIFT_THRESHOLD = 0.05
BORDERLINE_DRIFT_THRESHOLD = 0.03
TARGET_EKS = (1.5e-4, 2.0e-4, 3.0e-3, 5.0e-3, 7.0e-3, 1.0e-2, 3.0e-2)


def load_source_module():
    spec = importlib.util.spec_from_file_location("latest_source", SOURCE_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load source extractor: {SOURCE_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_fit_module():
    spec = importlib.util.spec_from_file_location("growth_fit", FIT_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load growth-fit implementation: {FIT_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Times New Roman",
            "mathtext.fontset": "stix",
            "axes.linewidth": FRAME_WIDTH,
            "axes.labelsize": 24,
            "xtick.labelsize": 13,
            "ytick.labelsize": 13,
            "legend.fontsize": 10,
            "figure.dpi": 180,
            "savefig.dpi": 300,
            "axes.unicode_minus": False,
        }
    )


def make_axes() -> tuple[plt.Figure, plt.Axes]:
    fig = plt.figure(figsize=FIGSIZE, facecolor="white")
    ax = fig.add_axes(AX_RECT)
    ax.set_box_aspect(BOX_ASPECT)
    for spine in ax.spines.values():
        spine.set_linewidth(FRAME_WIDTH)
    ax.tick_params(
        which="major", direction="in", length=12, width=1.2, top=True, right=True
    )
    ax.tick_params(
        which="minor", direction="in", length=6, width=1.0, top=True, right=True
    )
    ax.minorticks_on()
    return fig, ax


def save_figure(fig: plt.Figure, path: Path) -> None:
    fig.savefig(path.with_suffix(".png"), facecolor="white")
    fig.savefig(path.with_suffix(".pdf"), facecolor="white")
    plt.close(fig)


def as_array(rows: list[list[float]], columns: int = 2) -> np.ndarray:
    if not rows:
        return np.empty((0, columns), dtype=float)
    data = np.asarray(rows, dtype=float)
    data = data[np.all(np.isfinite(data), axis=1)]
    data = data[np.argsort(data[:, 0], kind="stable")]
    by_time = {round(float(row[0]), 8): row for row in data}
    return np.asarray([by_time[key] for key in sorted(by_time)], dtype=float)


def merge_arrays(archived: np.ndarray, current: np.ndarray) -> np.ndarray:
    """Merge same-definition histories, allowing current remote rows to win."""
    if not len(archived):
        return current
    if not len(current):
        return archived
    return as_array(np.vstack((archived, current)).tolist(), columns=current.shape[1])


def load_archived_gap_fill(record: dict) -> dict[str, np.ndarray]:
    metrics = {"kinetic", "convective_mid", "convective_zavg"}
    values = {metric: [] for metric in metrics}
    if not ARCHIVE_LONG.exists():
        return {metric: np.empty((0, 2), dtype=float) for metric in metrics}

    targets = {
        "Ra": float(record["Ra"]),
        "Ek": float(record["Ek"]),
        "AR": float(record["AR"]),
        "beta": float(record["beta"]),
        "qbot": float(record["qbot"]),
    }
    target_grid = (int(record["Nx"]), int(record["Ny"]), int(record["Nz"]))
    with ARCHIVE_LONG.open("r", newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            metric = row.get("metric", "")
            if metric not in metrics:
                continue
            try:
                grid = (int(row["Nx"]), int(row["Ny"]), int(row["Nz"]))
                same_case = grid == target_grid and all(
                    math.isclose(float(row[key]), target, rel_tol=1.0e-10, abs_tol=1.0e-12)
                    for key, target in targets.items()
                )
                if same_case:
                    values[metric].append([float(row["time"]), float(row["value"])])
            except (KeyError, TypeError, ValueError):
                continue
    return {metric: as_array(rows) for metric, rows in values.items()}


def contiguous_segments(time: np.ndarray) -> list[np.ndarray]:
    """Split a merged history at restart gaps before fitting or judging drift."""
    if len(time) < 2:
        return [np.arange(len(time), dtype=int)]
    dt = float(np.median(np.diff(time)))
    if not np.isfinite(dt) or dt <= 0.0:
        return [np.arange(len(time), dtype=int)]
    break_points = np.flatnonzero(np.diff(time) > 10.0 * dt) + 1
    return [
        chunk
        for chunk in np.split(np.arange(len(time), dtype=int), break_points)
        if len(chunk)
    ]


def terminal_window(time: np.ndarray, value: np.ndarray) -> dict[str, object]:
    segments = contiguous_segments(time)
    tail = segments[-1]
    tail_time = time[tail]
    tail_value = value[tail]
    end = float(tail_time[-1])
    span = max(float(tail_time[-1] - tail_time[0]), MIN_LATE_DURATION)
    duration = max(MIN_LATE_DURATION, 0.25 * span)
    start = max(float(tail_time[0]), end - duration)
    mask = tail_time >= start
    tw = tail_time[mask]
    vw = tail_value[mask]
    slope = float(np.polyfit(tw, vw, 1)[0]) if len(tw) >= 2 else float("nan")
    mean = float(np.mean(vw))
    normalized_slope = slope * (tw[-1] - tw[0]) / mean if mean != 0 else float("nan")
    split = start + 0.5 * (end - start)
    first = vw[tw < split]
    second = vw[tw >= split]
    half_drift = (
        float((np.mean(second) - np.mean(first)) / mean)
        if len(first) and len(second) and mean != 0
        else float("nan")
    )
    continuous_duration = float(tail_time[-1] - tail_time[0])
    has_sufficient_continuous_tail = continuous_duration >= MIN_LATE_DURATION
    diagnostic = max(abs(normalized_slope), abs(half_drift))
    if not has_sufficient_continuous_tail:
        status = "insufficient_continuous_late_record"
    elif diagnostic > STRONG_DRIFT_THRESHOLD:
        status = "clear_nonsteady"
    elif diagnostic > BORDERLINE_DRIFT_THRESHOLD:
        status = "borderline"
    else:
        status = "steady_by_terminal_drift"
    return {
        "start": start,
        "end": end,
        "duration": end - start,
        "continuous_tail_start": float(tail_time[0]),
        "continuous_tail_end": float(tail_time[-1]),
        "continuous_tail_duration": continuous_duration,
        "continuity_gap_detected": bool(len(segments) > 1),
        "sufficient_continuous_tail": has_sufficient_continuous_tail,
        "mean": mean,
        "std": float(np.std(vw)),
        "coefficient_of_variation": float(np.std(vw) / mean) if mean else float("nan"),
        "linear_slope": slope,
        "normalized_linear_drift": normalized_slope,
        "half_window_relative_drift": half_drift,
        "diagnostic_drift": diagnostic,
        "status": status,
    }


def scale_rows(record: dict, archived: dict[str, np.ndarray] | None = None) -> np.ndarray:
    rows = as_array(record["series"]["convective_mid"], columns=2)
    zavg = as_array(record["series"]["convective_zavg"], columns=2)
    if archived:
        rows = merge_arrays(archived.get("convective_mid", np.empty((0, 2))), rows)
        zavg = merge_arrays(archived.get("convective_zavg", np.empty((0, 2))), zavg)
    if not len(rows) or not len(zavg):
        return np.empty((0, 3), dtype=float)
    zmap = {round(float(row[0]), 8): float(row[1]) for row in zavg}
    merged = []
    for row in rows:
        key = round(float(row[0]), 8)
        if key in zmap:
            merged.append([float(row[0]), float(row[1]), zmap[key]])
    return np.asarray(merged, dtype=float)


def fit_growth(fit_module, time: np.ndarray, kinetic: np.ndarray) -> dict[str, float]:
    fit = fit_module.detect_exponential_interval(time, kinetic)
    return {key: float(value) for key, value in fit.items()}


def finite_range(values: np.ndarray) -> tuple[float, float]:
    values = values[np.isfinite(values)]
    if not len(values):
        return 0.0, 1.0
    low, high = float(np.min(values)), float(np.max(values))
    if high <= low:
        pad = max(abs(low) * 0.05, 1.0e-3)
    else:
        pad = 0.08 * (high - low)
    return max(0.0, low - pad), high + pad


def case_label(record: dict) -> str:
    ek = float(record["Ek"])
    ar = float(record["AR"])
    return f"Ek{ek:.2e}_AR{ar:g}"


def plot_kinetic(case_dir: Path, time: np.ndarray, kinetic: np.ndarray, fit: dict) -> None:
    fig, ax = make_axes()
    segments = contiguous_segments(time)
    for number, segment in enumerate(segments):
        ax.semilogy(
            time[segment],
            kinetic[segment],
            color=BLUE,
            lw=LINE_WIDTH,
            label=r"$K(t)$" if number == 0 else None,
        )
    fit_mask = (time >= fit["t_start"]) & (time <= fit["t_end"])
    fitted = np.exp(fit["intercept"] + fit["slope"] * time[fit_mask])
    ax.semilogy(
        time[fit_mask],
        fitted,
        color=FIT_RED,
        lw=LINE_WIDTH,
        ls="--",
        label=r"linear fit in $\ln K$",
    )
    ax.axvspan(
        fit["t_start"], fit["t_end"], color=SHADE, alpha=0.22, linewidth=0, zorder=0
    )
    ax.set_xlim(float(time[0]), float(time[-1]))
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$K=\frac{1}{2}\langle u^2+v^2+w^2\rangle_V$")
    gap_note = "\nrestart gap shown as a break" if len(segments) > 1 else ""
    ax.text(
        0.04,
        0.96,
        rf"${fit['t_start']:.1f}\leq t\leq {fit['t_end']:.1f}$"
        + "\n"
        + rf"$d\ln K/dt={fit['slope']:.3f}$"
        + gap_note,
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=12,
    )
    ax.legend(frameon=False, loc="lower right")
    save_figure(fig, case_dir / "kinetic_energy_linear_growth")


def plot_lc_window(
    case_dir: Path,
    scale: np.ndarray,
    start: float,
    end: float,
    stem: str,
    label_suffix: str,
) -> dict[str, float]:
    mask = (scale[:, 0] >= start) & (scale[:, 0] <= end)
    selected = scale[mask]
    if not len(selected):
        raise RuntimeError(f"No convective-scale rows in {label_suffix} window.")
    fig, ax = make_axes()
    ax.plot(selected[:, 0], selected[:, 1], color=BLUE, lw=LINE_WIDTH, label=r"$z\simeq0.5$")
    ax.plot(
        selected[:, 0],
        selected[:, 2],
        color=LIGHT_BLUE,
        lw=LINE_WIDTH,
        ls="--",
        label=r"$z$-weighted mean",
    )
    ax.set_xlim(float(selected[0, 0]), float(selected[-1, 0]))
    low, high = finite_range(selected[:, 1:])
    ax.set_ylim(low, high)
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$2\pi l_c$")
    ax.legend(frameon=False, loc="best")
    save_figure(fig, case_dir / stem)
    return {
        "start": float(selected[0, 0]),
        "end": float(selected[-1, 0]),
        "midheight_mean": float(np.mean(selected[:, 1])),
        "midheight_std": float(np.std(selected[:, 1])),
        "vertical_mean_mean": float(np.mean(selected[:, 2])),
        "vertical_mean_std": float(np.std(selected[:, 2])),
        "midheight_min": float(np.min(selected[:, 1])),
        "midheight_max": float(np.max(selected[:, 1])),
        "vertical_mean_min": float(np.min(selected[:, 2])),
        "vertical_mean_max": float(np.max(selected[:, 2])),
    }


def write_case_data(case_dir: Path, time: np.ndarray, kinetic: np.ndarray, scale: np.ndarray) -> None:
    with (case_dir / "reduced_timeseries.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["time", "kinetic_energy", "convective_wavelength_z05_2pi", "convective_wavelength_z_weighted_2pi"])
        scale_map = {round(float(row[0]), 8): row for row in scale}
        for t, k in zip(time, kinetic):
            row = scale_map.get(round(float(t), 8))
            writer.writerow([t, k, "" if row is None else row[1], "" if row is None else row[2]])


def main() -> None:
    configure_style()
    source = load_source_module()
    fit_module = load_fit_module()
    raw_records = source.extract_remote()
    records, duplicates = source.choose_highest_resolution(raw_records)
    records = [
        record
        for record in records
        if record.get("Ek") is not None
        and any(
            math.isclose(float(record["Ek"]), ek, rel_tol=1.0e-8)
            for ek in TARGET_EKS
        )
    ]
    records = sorted(records, key=lambda record: float(record["Ek"]))

    if OUTPUT_DIR.exists():
        shutil.rmtree(OUTPUT_DIR)
    OUTPUT_DIR.mkdir(parents=True)

    summary_rows = []
    summary_json = []
    for record in records:
        kinetic_rows = as_array(record["series"]["kinetic"])
        kinetic_rows = kinetic_rows[(kinetic_rows[:, 1] >= 0.0) & (kinetic_rows[:, 1] < 0.1)]
        if len(kinetic_rows) < 3:
            continue
        remote_gap_detected = len(contiguous_segments(kinetic_rows[:, 0])) > 1
        archived = load_archived_gap_fill(record) if remote_gap_detected else {}
        archive_backfill_used = bool(
            remote_gap_detected and len(archived.get("kinetic", np.empty((0, 2))))
        )
        if archive_backfill_used:
            kinetic_rows = merge_arrays(archived["kinetic"], kinetic_rows)
        time, kinetic = kinetic_rows[:, 0], kinetic_rows[:, 1]
        fit_segment = contiguous_segments(time)[0]
        fit = fit_growth(fit_module, time[fit_segment], kinetic[fit_segment])
        terminal = terminal_window(time, kinetic)
        scale = scale_rows(record, archived)
        if len(scale) < 3:
            continue
        case_dir = OUTPUT_DIR / case_label(record)
        case_dir.mkdir()
        write_case_data(case_dir, time, kinetic, scale)
        early_stats = plot_lc_window(
            case_dir,
            scale,
            fit["t_start"],
            fit["t_end"],
            "convective_wavelength_linear_growth",
            "linear-growth",
        )
        late_stats = plot_lc_window(
            case_dir,
            scale,
            float(terminal["start"]),
            float(terminal["end"]),
            "convective_wavelength_late_window",
            "late",
        )
        plot_kinetic(case_dir, time, kinetic, fit)

        row = {
            "Ek": float(record["Ek"]),
            "Ra": float(record["Ra"]),
            "Pr": float(record["Pr"]),
            "AR": float(record["AR"]),
            "Nx": int(record["Nx"]),
            "Ny": int(record["Ny"]),
            "Nz": int(record["Nz"]),
            "time_end": float(time[-1]),
            "remote_gap_detected_before_backfill": remote_gap_detected,
            "archive_backfill_used": archive_backfill_used,
            "archive_backfill_source": str(ARCHIVE_LONG) if archive_backfill_used else "",
            "source_run_paths": " | ".join(record.get("source_run_paths", [record["run_path"]])),
            "linear_t_start": fit["t_start"],
            "linear_t_end": fit["t_end"],
            "linear_duration": fit["duration"],
            "dlnK_dt": fit["slope"],
            "linear_R2": fit["r_squared"],
            "linear_max_abs_log_residual": fit["max_abs_log_residual"],
            "late_t_start": terminal["start"],
            "late_t_end": terminal["end"],
            "continuous_tail_start": terminal["continuous_tail_start"],
            "continuous_tail_end": terminal["continuous_tail_end"],
            "continuous_tail_duration": terminal["continuous_tail_duration"],
            "continuity_gap_detected": terminal["continuity_gap_detected"],
            "sufficient_continuous_tail": terminal["sufficient_continuous_tail"],
            "late_K_mean": terminal["mean"],
            "late_K_std": terminal["std"],
            "late_K_normalized_drift": terminal["normalized_linear_drift"],
            "late_K_half_window_relative_drift": terminal["half_window_relative_drift"],
            "late_K_diagnostic_drift": terminal["diagnostic_drift"],
            "steady_status": terminal["status"],
            **{f"linear_lc_{key}": value for key, value in early_stats.items() if key != "start" and key != "end"},
            "linear_lc_t_start": early_stats["start"],
            "linear_lc_t_end": early_stats["end"],
            **{f"late_lc_{key}": value for key, value in late_stats.items() if key != "start" and key != "end"},
            "late_lc_t_start": late_stats["start"],
            "late_lc_t_end": late_stats["end"],
        }
        summary_rows.append(row)
        summary_json.append({**row, "case_label": case_label(record)})

    fields = list(summary_rows[0]) if summary_rows else []
    with (OUTPUT_DIR / "case_summary.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(summary_rows)

    manifest = {
        "generated_from": str(SOURCE_SCRIPT),
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "remote_host": "c01n0006",
        "remote_root": "/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case",
        "latest_program_only": True,
        "target_Ek": list(TARGET_EKS),
        "growth_fit_source": str(FIT_SCRIPT),
        "selection": "Only cases with mprime_square.dat, convective_scale.dat, moist_integral_scale.dat and peak_scale.dat; actual Nx,Ny >=257; highest resolution kept for identical parameters; base and continuation segments merged by physical time.",
        "gap_fill": f"If a current remote restart history has a physical-time gap, fill only that gap from the same latest-program reduced archive: {ARCHIVE_LONG}; current remote rows override duplicate archived times.",
        "kinetic_definition": "K(t)=0.5*avgvar one-based column 4",
        "convective_scale_definition": "plotted wavelength = 2*pi*l_c from diagnostics/scale/convective_scale.dat; solid is z approximately 0.5 and dashed is vertically averaged l_c",
        "linear_growth_fit": "straight-line least-squares fit to ln(K) over the longest accepted positive-slope interval before the first dominant K peak, using the established R2 >= 0.9999 and max log residual <= 0.06 criteria",
        "terminal_window": "last quarter of the final contiguous K segment, with minimum duration 150 time units; restart gaps are not bridged",
        "steady_classification": {
            "insufficient_continuous_late_record": "the final contiguous segment is shorter than 150 time units, so a late steady-state claim is not possible",
            "clear_nonsteady": "max(abs(normalized terminal linear drift), abs(first-half to second-half relative drift)) > 0.05",
            "borderline": "same diagnostic > 0.03 and <= 0.05",
            "steady_by_terminal_drift": "same diagnostic <= 0.03",
            "interpretation": "borderline is not claimed as fully statistically steady; it indicates that the current record is close but still has measurable systematic drift",
        },
        "selected_cases": len(summary_rows),
        "excluded_lower_resolution_duplicates": len(duplicates),
        "cases": summary_json,
    }
    (OUTPUT_DIR / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")

    print(json.dumps({"output_dir": str(OUTPUT_DIR), "cases": summary_json}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
