from __future__ import annotations

import csv
import json
import math
import re
import subprocess
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D


REMOTE_HOST = "c01n0011"
REMOTE_ROOT = "/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case"
PROJECT_ROOT = Path(r"E:\moist RB\rotating_case_inventory")
OUTPUT_ROOT = PROJECT_ROOT / "04_outputs_and_figures"
HISTORICAL_PATCH_ROOT = (
    PROJECT_ROOT / "03_inventory_tables" / "historical_time_series_patches"
)

# The Ek=7e-3 run lost its remote t=10--1000 histories after an accidental
# file operation.  This immutable reduced-data patch was exported from the
# complete 2026-08-01 figures.  Patch rows are merged first, so current remote
# rows always win at duplicate physical times.
HISTORICAL_PATCH_RULES = (
    {
        "name": "Ra8e6_Ek7e-3_t10_1000",
        "path": HISTORICAL_PATCH_ROOT
        / "Ra8e6_Ek7e-3_t10_1000_latest_program.csv",
        "Ra": 8.0e6,
        "Ek": 7.0e-3,
        "AR": 16.0,
        "beta": 1.02,
        "qbot": 0.5,
        "Nx": 385,
        "Ny": 385,
        "Nz": 65,
        "time_min": 10.0,
        "time_max": 1000.0,
    },
)

MIN_HORIZONTAL_RESOLUTION = 257
LATEST_PROGRAM_ONLY = True
FRAME_WIDTH = 4.5
LINE_WIDTH = 3.5
AXIS_LABEL_SIZE = 24
TICK_LABEL_SIZE = 13
LEGEND_SIZE = 10
BOX_ASPECT = 5.2 / 6.5
KINETIC_RED = (0.90, 0.05, 0.08)
MSE_BLUE = (0.00, 0.25, 0.90)


# The cluster login node uses an old Python version. Keep this embedded
# extractor compatible with Python 3.6 and return only small text histories.
REMOTE_CODE = r"""
import json, math, re
from pathlib import Path

ROOT = Path(r'''__REMOTE_ROOT__''')

def as_float(token):
    try:
        return float(str(token).replace('D', 'e').replace('d', 'e'))
    except Exception:
        return None

def label_float(token, prefix):
    if not str(token).startswith(prefix):
        return None
    return as_float(str(token)[len(prefix):].replace('p', '.'))

def row_after(lines, marker):
    for i, line in enumerate(lines):
        if marker in line:
            for j in range(i + 1, len(lines)):
                clean = lines[j].split('!')[0].strip()
                if clean:
                    return clean.split()
    return []

def read_bou(run):
    path = run / 'bou.in'
    with path.open('r', errors='ignore') as handle:
        lines = handle.readlines()
    grid = row_after(lines, 'N1      N2')
    geometry = row_after(lines, 'ALX3D')
    control = row_after(lines, 'Ra       Pr')
    boundary = row_after(lines, 'A_stopmod')
    if len(grid) < 3 or len(geometry) < 3 or len(control) < 9 or len(boundary) < 8:
        return None
    values = {
        'Nx': int(as_float(grid[0])),
        'Ny': int(as_float(grid[1])),
        'Nz': int(as_float(grid[2])),
        'AR': as_float(geometry[1]),
        'Ra': as_float(control[0]),
        'Pr': as_float(control[1]),
        'invRo': as_float(control[2]),
        'control_alpha': as_float(control[3]),
        'gamma': as_float(control[4]),
        'Sm': as_float(control[5]),
        'alphaqs': as_float(control[6]),
        'beta': as_float(control[7]),
        'tau': as_float(control[8]),
        'qtop': as_float(boundary[6]),
        'qbot': as_float(boundary[7]),
        'dsaltop': as_float(boundary[4]),
        'dsalbot': as_float(boundary[5]),
    }
    if values['invRo'] is not None and abs(values['invRo']) > 1.0e-15:
        values['Ek'] = math.sqrt(values['Pr'] / values['Ra']) / values['invRo']
    else:
        values['Ek'] = None
    return values

def read_rows(path, minimum_columns):
    if not path.exists() or path.stat().st_size == 0:
        return []
    rows = []
    with path.open('r', errors='ignore') as handle:
        for line in handle:
            parts = line.split()
            if len(parts) < minimum_columns:
                continue
            values = [as_float(item) for item in parts]
            if any(item is None or not math.isfinite(item) for item in values):
                continue
            rows.append(values)
    return rows

def convert_series(run, params):
    data = run / 'data'
    diagnostics = run / 'diagnostics'

    avg = read_rows(data / 'avgvar.out', 4)
    kinetic = []
    l0 = []
    for row in avg:
        kinetic_value = 0.5 * row[3]
        if 0.0 <= kinetic_value < 0.1:
            kinetic.append([row[0], kinetic_value])
        if len(row) >= 9 and row[7] > 0.0 and row[8] > 0.0:
            l0.append([row[0], row[7], row[8]])

    # nu_profiles.out stores one row per face and time:
    # time, z_face, NuT(z), Num(z), Nuq(z).  Integrate the instantaneous
    # Num profile over the nonuniform vertical coordinate before merging
    # continuation segments.
    num_by_time = {}
    for row in read_rows(data / 'nu_profiles.out', 5):
        num_by_time.setdefault(round(row[0], 8), {})[round(row[1], 12)] = row[3]
    delta_m = (
        float(params.get('dsalbot', 0.0) - params.get('dsaltop', 0.0))
        + float(params.get('gamma', 0.0))
        * float(params.get('qbot', 0.0) - params.get('qtop', 0.0))
    )
    num = []
    for time_key in sorted(num_by_time):
        profile = sorted(num_by_time[time_key].items())
        if len(profile) < 2:
            continue
        integral = 0.0
        for left, right in zip(profile[:-1], profile[1:]):
            integral += 0.5 * (left[1] + right[1]) * (right[0] - left[0])
        height = profile[-1][0] - profile[0][0]
        if height > 0.0:
            raw_flux = integral / height
            if delta_m > 0.0:
                num.append([time_key, raw_flux / delta_m])

    mprime = []
    new_mprime = read_rows(diagnostics / 'thermo' / 'mprime_square.dat', 5)
    old_mprime = read_rows(data / 'mse_aggregation.out', 3)
    if new_mprime:
        mprime = [[row[0], row[4]] for row in new_mprime if row[4] >= 0.0]
        mprime_source = 'diagnostics/thermo/mprime_square.dat'
    else:
        mprime = [[row[0], row[1]] for row in old_mprime if row[1] >= 0.0]
        mprime_source = 'data/mse_aggregation.out' if old_mprime else ''

    mse_peak = []
    mse_spectral = []
    mse_integral = []
    new_moist = read_rows(diagnostics / 'scale' / 'moist_integral_scale.dat', 4)
    new_peak = read_rows(diagnostics / 'scale' / 'peak_scale.dat', 4)
    old_scales = read_rows(data / 'mse_aggregation_scales.out', 11)
    if new_moist:
        mse_spectral = [[row[0], 2.0 * math.pi * row[1]] for row in new_moist if row[1] > 0.0]
    elif old_scales:
        mse_spectral = [[row[0], row[4]] for row in old_scales if row[4] > 0.0]
    if new_peak:
        mse_peak = [[row[0], row[2]] for row in new_peak if row[2] > 0.0]
    elif old_scales:
        mse_peak = [[row[0], row[1]] for row in old_scales if row[1] > 0.0]
    if old_scales:
        mse_integral = [[row[0], row[2]] for row in old_scales if row[2] > 0.0]

    convective_mid = []
    convective_zavg = []
    velocity_z075 = []
    new_convective = read_rows(diagnostics / 'scale' / 'convective_scale.dat', 3)
    old_velocity = read_rows(data / 'w_z075_spectral_length.out', 3)
    if new_convective:
        convective_mid = [[row[0], 2.0 * math.pi * row[1]] for row in new_convective if row[1] > 0.0]
        convective_zavg = [[row[0], 2.0 * math.pi * row[2]] for row in new_convective if row[2] > 0.0]
    if old_velocity:
        velocity_z075 = [[row[0], row[2]] for row in old_velocity if row[2] > 0.0]

    return {
        'kinetic': kinetic,
        'num': num,
        'mprime': mprime,
        'l0': l0,
        'mse_peak': mse_peak,
        'mse_spectral': mse_spectral,
        'mse_integral': mse_integral,
        'convective_mid': convective_mid,
        'convective_zavg': convective_zavg,
        'velocity_z075': velocity_z075,
        'mprime_source': mprime_source,
    }

records = []
for run in sorted(ROOT.rglob('run')):
    if not (run / 'bou.in').exists():
        continue
    rel = run.relative_to(ROOT).parts
    if len(rel) < 8 or not rel[0].startswith('Pr') or not rel[1].startswith('Ra'):
        continue
    params = read_bou(run)
    if params is None:
        continue
    latest_program_files = [
        run / 'diagnostics' / 'thermo' / 'mprime_square.dat',
        run / 'diagnostics' / 'scale' / 'convective_scale.dat',
        run / 'diagnostics' / 'scale' / 'moist_integral_scale.dat',
        run / 'diagnostics' / 'scale' / 'peak_scale.dat',
    ]
    if not all(
        path.exists() and path.stat().st_size > 0 for path in latest_program_files
    ):
        continue
    series = convert_series(run, params)
    if not series['mprime']:
        continue
    params.update({
        'run_path': str(run),
        'anchor_path': str(ROOT.joinpath(*rel[:7])),
        'Pr_label': rel[0],
        'Ra_label': rel[1],
        'Ek_label': rel[2],
        'AR_label': rel[3],
        'Beta_label': rel[4],
        'q_label': rel[5],
        'grid_directory_label': rel[6],
        'segment': '/'.join(rel[7:-1]),
        'series': series,
    })
    records.append(params)

print(json.dumps(records, ensure_ascii=False))
"""


SERIES_NAMES = (
    "kinetic",
    "num",
    "mprime",
    "l0",
    "mse_peak",
    "mse_spectral",
    "mse_integral",
    "convective_mid",
    "convective_zavg",
    "velocity_z075",
)


def extract_remote() -> list[dict]:
    code = REMOTE_CODE.replace("__REMOTE_ROOT__", REMOTE_ROOT)
    process = subprocess.run(
        ["ssh", REMOTE_HOST, "python3", "-"],
        input=code,
        text=True,
        capture_output=True,
        encoding="utf-8",
        errors="replace",
        check=False,
    )
    if process.returncode != 0:
        raise RuntimeError(
            "Remote extraction failed.\n"
            f"STDOUT:\n{process.stdout}\nSTDERR:\n{process.stderr}"
        )
    return json.loads(process.stdout)


def rounded(value: float | None, digits: int = 12):
    return None if value is None else round(float(value), digits)


def segment_priority(segment: str) -> tuple[int, int, str]:
    if not segment:
        return (0, 0, "")
    match = re.search(r"conti(?:nuation)?(\d+)", segment, re.IGNORECASE)
    if match:
        return (1, int(match.group(1)), segment)
    return (2, 0, segment)


def deduplicate_rows(rows: list[list[float]]) -> list[list[float]]:
    by_time: dict[float, list[float]] = {}
    for row in rows:
        if not row:
            continue
        values = [float(value) for value in row]
        if all(math.isfinite(value) for value in values):
            by_time[round(values[0], 8)] = values
    return [by_time[key] for key in sorted(by_time)]


def close_value(value: float | None, target: float, tolerance: float = 1.0e-8) -> bool:
    if value is None:
        return False
    scale = max(1.0, abs(float(value)), abs(float(target)))
    return abs(float(value) - float(target)) <= tolerance * scale


def patch_matches(record: dict, rule: dict) -> bool:
    for name in ("Ra", "Ek", "AR", "beta", "qbot"):
        if not close_value(record.get(name), float(rule[name])):
            return False
    return all(int(record.get(name, -1)) == int(rule[name]) for name in ("Nx", "Ny", "Nz"))


def read_historical_patch(path: Path, time_min: float, time_max: float) -> dict[str, list[list[float]]]:
    components: dict[str, dict[float, dict[int, float]]] = defaultdict(
        lambda: defaultdict(dict)
    )
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            try:
                time = float(row["time"])
                value = float(row["value"])
            except (KeyError, TypeError, ValueError):
                continue
            if not (math.isfinite(time) and math.isfinite(value)):
                continue
            if time < time_min - 1.0e-8 or time > time_max + 1.0e-8:
                continue
            metric = str(row.get("metric", ""))
            if metric in SERIES_NAMES:
                base, component = metric, 1
            else:
                match = re.fullmatch(r"(.+)_(\d+)", metric)
                if not match or match.group(1) not in SERIES_NAMES:
                    continue
                base, component = match.group(1), int(match.group(2))
            components[base][round(time, 8)][component] = value

    series = {name: [] for name in SERIES_NAMES}
    for name, by_time in components.items():
        component_count = max(
            (max(values) for values in by_time.values() if values), default=0
        )
        for time in sorted(by_time):
            values = by_time[time]
            if component_count and all(index in values for index in range(1, component_count + 1)):
                series[name].append(
                    [float(time)] + [values[index] for index in range(1, component_count + 1)]
                )
    return series


def apply_historical_patches(records: list[dict]) -> list[dict]:
    reports = []
    for rule in HISTORICAL_PATCH_RULES:
        path = Path(rule["path"])
        if not path.exists():
            raise FileNotFoundError(f"Required historical patch is missing: {path}")
        matches = [record for record in records if patch_matches(record, rule)]
        if len(matches) != 1:
            raise RuntimeError(
                f"Historical patch {rule['name']} matched {len(matches)} selected cases; expected one."
            )
        record = matches[0]
        patch_series = read_historical_patch(
            path, float(rule["time_min"]), float(rule["time_max"])
        )
        added_by_metric = {}
        for name in SERIES_NAMES:
            current = record["series"].get(name, [])
            patched = deduplicate_rows(patch_series.get(name, []) + current)
            added_by_metric[name] = len(patched) - len(current)
            record["series"][name] = patched

        report = {
            "name": rule["name"],
            "path": str(path),
            "time_min": float(rule["time_min"]),
            "time_max": float(rule["time_max"]),
            "priority": "current remote rows override duplicate historical rows",
            "added_rows_by_metric": added_by_metric,
        }
        record.setdefault("historical_patches", []).append(report)
        all_times = [
            row[0]
            for name in SERIES_NAMES
            for row in record["series"].get(name, [])
        ]
        record["time_min"] = min(all_times, default=None)
        record["time_max"] = max(all_times, default=None)
        reports.append(report)
    return reports


def branch_identity(record: dict) -> tuple:
    return (
        record["anchor_path"],
        int(record["Nx"]),
        int(record["Ny"]),
        int(record["Nz"]),
    )


def physical_identity(record: dict) -> tuple:
    return (
        rounded(record.get("Pr")),
        rounded(record.get("Ra")),
        rounded(record.get("Ek")),
        rounded(record.get("AR")),
        rounded(record.get("beta")),
        rounded(record.get("gamma")),
        rounded(record.get("alphaqs")),
        rounded(record.get("tau")),
        rounded(record.get("qbot")),
        rounded(record.get("qtop")),
    )


def merge_branch(segments: list[dict]) -> dict:
    segments = sorted(segments, key=lambda record: segment_priority(record["segment"]))
    merged = {key: value for key, value in segments[0].items() if key != "series"}
    merged["series"] = {name: [] for name in SERIES_NAMES}
    merged["source_run_paths"] = []
    merged["mprime_sources"] = []
    for segment in segments:
        merged["source_run_paths"].append(segment["run_path"])
        source = segment["series"].get("mprime_source")
        if source and source not in merged["mprime_sources"]:
            merged["mprime_sources"].append(source)
        for name in SERIES_NAMES:
            merged["series"][name].extend(segment["series"].get(name, []))
    for name in SERIES_NAMES:
        merged["series"][name] = deduplicate_rows(merged["series"][name])
    all_times = [
        row[0]
        for name in SERIES_NAMES
        for row in merged["series"][name]
    ]
    merged["time_max"] = max(all_times, default=None)
    merged["time_min"] = min(all_times, default=None)
    return merged


def choose_highest_resolution(raw_records: list[dict]) -> tuple[list[dict], list[dict]]:
    branches: dict[tuple, list[dict]] = defaultdict(list)
    for record in raw_records:
        branches[branch_identity(record)].append(record)
    merged_branches = [merge_branch(segments) for segments in branches.values()]

    eligible = [
        record
        for record in merged_branches
        if int(record["Nx"]) >= MIN_HORIZONTAL_RESOLUTION
        and int(record["Ny"]) >= MIN_HORIZONTAL_RESOLUTION
    ]
    grouped: dict[tuple, list[dict]] = defaultdict(list)
    for record in eligible:
        grouped[physical_identity(record)].append(record)

    selected = []
    excluded_duplicates = []
    for candidates in grouped.values():
        ordered = sorted(
            candidates,
            key=lambda record: (
                min(int(record["Nx"]), int(record["Ny"])),
                int(record["Nx"]) * int(record["Ny"]) * int(record["Nz"]),
                float(record.get("time_max") or -1.0),
                len(record["series"]["mprime"]),
            ),
            reverse=True,
        )
        selected.append(ordered[0])
        excluded_duplicates.extend(ordered[1:])
    return selected, excluded_duplicates


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Times New Roman",
            "mathtext.fontset": "stix",
            "axes.linewidth": FRAME_WIDTH,
            "axes.labelsize": AXIS_LABEL_SIZE,
            "xtick.labelsize": TICK_LABEL_SIZE,
            "ytick.labelsize": TICK_LABEL_SIZE,
            "legend.fontsize": LEGEND_SIZE,
            "axes.formatter.use_mathtext": True,
            "figure.dpi": 180,
            "savefig.dpi": 300,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def style_axis(axis: plt.Axes) -> None:
    for spine in axis.spines.values():
        spine.set_linewidth(FRAME_WIDTH)
    axis.tick_params(
        which="major", direction="in", length=12, width=1.2, top=True, right=True
    )
    axis.minorticks_on()
    axis.tick_params(
        which="minor", direction="in", length=6, width=1.0, top=True, right=True
    )
    axis.ticklabel_format(
        axis="y", style="sci", scilimits=(-2, 2), useMathText=True
    )
    axis.yaxis.get_offset_text().set_fontsize(TICK_LABEL_SIZE)
    axis.set_box_aspect(BOX_ASPECT)


def sci_label(value: float | None) -> str:
    if value is None:
        return "NR"
    # Ek is reconstructed from Ra, Pr, and invRo and can lie a few ulps below
    # an exact decade (for example 9.99999999999935e-3).  The small offset
    # keeps clean decades from being rendered as 10 x 10^(n-1).
    exponent = int(math.floor(math.log10(abs(value)) + 1.0e-10))
    mantissa = value / 10.0**exponent
    if abs(mantissa - round(mantissa)) < 1.0e-8:
        mantissa_text = str(int(round(mantissa)))
    else:
        mantissa_text = f"{mantissa:.2g}"
    if mantissa_text == "1":
        return rf"10^{{{exponent}}}"
    return rf"{mantissa_text}\times10^{{{exponent}}}"


def plain_label(value: float | None, digits: int = 3) -> str:
    if value is None:
        return ""
    return f"{float(value):.{digits}g}"


def scientific_slug(value: float) -> str:
    mantissa, exponent = f"{float(value):.8e}".split("e")
    mantissa = mantissa.rstrip("0").rstrip(".").replace(".", "p")
    return f"{mantissa}e{int(exponent)}"


def individual_case_stem(record: dict) -> str:
    ek = "NR" if record.get("Ek") is None else scientific_slug(record["Ek"])
    return "_".join(
        (
            f"Ra{scientific_slug(record['Ra'])}",
            f"Ek{ek}",
            f"AR{plain_label(record.get('AR')).replace('.', 'p')}",
            f"Beta{plain_label(record.get('beta')).replace('.', 'p')}",
            f"qbot{plain_label(record.get('qbot')).replace('.', 'p')}",
            f"N{record['Nx']}x{record['Ny']}x{record['Nz']}",
        )
    )


def individual_case_title(record: dict) -> str:
    ek = r"$\mathrm{NR}$" if record.get("Ek") is None else rf"$Ek={sci_label(record['Ek'])}$"
    return ", ".join(
        (
            rf"$Ra={sci_label(record['Ra'])}$",
            ek,
            rf"$\Gamma={plain_label(record.get('AR'))}$",
            rf"$\beta={plain_label(record.get('beta'))}$",
            rf"$q_{{\mathrm{{bot}}}}={plain_label(record.get('qbot'))}$",
            rf"${record['Nx']}\!\times\!{record['Ny']}\!\times\!{record['Nz']}$",
        )
    )


def case_label(record: dict, records: list[dict]) -> str:
    parts = ["NR" if record.get("Ek") is None else rf"$Ek={sci_label(record['Ek'])}$"]
    ar_values = {rounded(item.get("AR")) for item in records}
    beta_values = {rounded(item.get("beta")) for item in records}
    qbot_values = {rounded(item.get("qbot")) for item in records}
    if len(ar_values) > 1:
        parts.append(rf"$\Gamma={plain_label(record.get('AR'))}$")
    if len(beta_values) > 1:
        parts.append(rf"$\beta={plain_label(record.get('beta'))}$")
    if len(qbot_values) > 1 or abs(float(record.get("qbot") or 0.0) - 1.0) > 1.0e-12:
        parts.append(rf"$q_{{\mathrm{{bot}}}}={plain_label(record.get('qbot'))}$")
    parts.append(rf"${record['Nx']}\!\times\!{record['Ny']}\!\times\!{record['Nz']}$")
    return ", ".join(parts)


def record_sort_key(record: dict) -> tuple:
    ek = -1.0 if record.get("Ek") is None else float(record["Ek"])
    return (
        float(record.get("qbot") or 0.0),
        float(record.get("AR") or 0.0),
        float(record.get("beta") or 0.0),
        ek,
    )


def case_colors(records: list[dict]) -> dict[int, tuple]:
    colors = plt.cm.turbo(np.linspace(0.05, 0.95, max(len(records), 2)))
    return {id(record): tuple(colors[index]) for index, record in enumerate(records)}


def add_case_legend(
    figure: plt.Figure,
    records: list[dict],
    colors: dict[int, tuple],
    anchor: tuple[float, float],
) -> None:
    handles = [
        Line2D([0], [0], color=colors[id(record)], lw=LINE_WIDTH)
        for record in records
    ]
    labels = [case_label(record, records) for record in records]
    figure.legend(
        handles,
        labels,
        frameon=False,
        loc="center left",
        bbox_to_anchor=anchor,
        handlelength=2.4,
        labelspacing=0.55,
        borderaxespad=0.0,
    )


def plot_series(axis: plt.Axes, records: list[dict], colors: dict[int, tuple], name: str) -> int:
    plotted = 0
    for record in records:
        rows = record["series"].get(name, [])
        if not rows:
            continue
        values = np.asarray(rows, dtype=float)
        axis.plot(values[:, 0], values[:, 1], color=colors[id(record)], lw=LINE_WIDTH)
        plotted += 1
    return plotted


def save_figure(figure: plt.Figure, output_dir: Path, stem: str) -> list[Path]:
    png = output_dir / f"{stem}.png"
    pdf = output_dir / f"{stem}.pdf"
    figure.savefig(png, facecolor="white")
    figure.savefig(pdf, facecolor="white")
    plt.close(figure)
    return [png, pdf]


PLUME_BURST_WINDOWS = {
    2.0e-4: [(353.3, 432.6), (622.9, 732.1), (939.4, 1029.8)],
    5.0e-4: [(390.5, 429.0), (451.8, 488.5)],
    7.0e-4: [(447.9, 469.6), (491.1, 523.1), (580.6, 684.6)],
}


def add_plume_burst_windows(axis: plt.Axes, record: dict) -> None:
    ek = record.get("Ek")
    if ek is None:
        return
    for target, windows in PLUME_BURST_WINDOWS.items():
        if math.isclose(float(ek), target, rel_tol=1.0e-8, abs_tol=1.0e-12):
            for start, stop in windows:
                axis.axvspan(start, stop, color=(1.0, 0.82, 0.12), alpha=0.20, lw=0, zorder=0)
            return


def plot_energy_mprime(
    ra_name: str,
    records: list[dict],
    colors: dict[int, tuple],
    output_dir: Path,
) -> list[Path]:
    figure = plt.figure(figsize=(16.5, 5.84), facecolor="white")
    axes = [
        figure.add_axes([0.055, 0.260, 0.312, 0.705]),
        figure.add_axes([0.425, 0.260, 0.312, 0.705]),
    ]
    for axis in axes:
        style_axis(axis)
    plot_series(axes[0], records, colors, "kinetic")
    plot_series(axes[1], records, colors, "mprime")
    axes[0].set_xlabel(r"$t$")
    axes[0].set_ylabel(r"$K$")
    axes[1].set_xlabel(r"$t$")
    axes[1].set_ylabel(r"$\left\langle\left\langle m'^2\right\rangle_{xy}\right\rangle_z$")
    axes[0].set_ylim(bottom=0.0)
    axes[1].set_ylim(bottom=0.0)
    add_case_legend(figure, records, colors, (0.770, 0.52))
    return save_figure(figure, output_dir, f"{ra_name}_kinetic_mprime_timeseries")


def plot_l0(
    ra_name: str,
    records: list[dict],
    colors: dict[int, tuple],
    output_dir: Path,
) -> list[Path]:
    figure = plt.figure(figsize=(16.5, 5.84), facecolor="white")
    axes = [
        figure.add_axes([0.055, 0.260, 0.312, 0.705]),
        figure.add_axes([0.425, 0.260, 0.312, 0.705]),
    ]
    for axis in axes:
        style_axis(axis)
    for record in records:
        rows = record["series"].get("l0", [])
        if not rows:
            continue
        values = np.asarray(rows, dtype=float)
        axes[0].plot(values[:, 0], values[:, 1], color=colors[id(record)], lw=LINE_WIDTH)
        axes[1].plot(values[:, 0], values[:, 2], color=colors[id(record)], lw=LINE_WIDTH)
    for axis in axes:
        axis.set_xlabel(r"$t$")
        axis.set_ylim(bottom=0.0)
    axes[0].set_ylabel(r"$l_0(z\simeq0.5)$")
    axes[1].set_ylabel(r"$l_0(z\simeq0.75)$")
    add_case_legend(figure, records, colors, (0.770, 0.52))
    return save_figure(figure, output_dir, f"{ra_name}_horizontal_velocity_l0_timeseries")


def plot_mse_scales(
    ra_name: str,
    records: list[dict],
    colors: dict[int, tuple],
    output_dir: Path,
) -> list[Path]:
    figure = plt.figure(figsize=(16.5, 5.84), facecolor="white")
    axes = [
        figure.add_axes([0.055, 0.260, 0.312, 0.705]),
        figure.add_axes([0.425, 0.260, 0.312, 0.705]),
    ]
    names = ("mse_peak", "mse_spectral")
    labels = (
        r"$L_{\mathrm{peak},m}$",
        r"$2\pi L_m$",
    )
    for axis, name, label in zip(axes, names, labels):
        style_axis(axis)
        plot_series(axis, records, colors, name)
        axis.set_xlabel(r"$t$")
        axis.set_ylabel(label)
        axis.set_ylim(bottom=0.0)
    add_case_legend(figure, records, colors, (0.770, 0.52))
    return save_figure(figure, output_dir, f"{ra_name}_mse_length_scales_timeseries")


def plot_velocity_scales(
    ra_name: str,
    records: list[dict],
    colors: dict[int, tuple],
    output_dir: Path,
) -> list[Path]:
    figure = plt.figure(figsize=(16.5, 5.84), facecolor="white")
    axes = [
        figure.add_axes([0.055, 0.260, 0.312, 0.705]),
        figure.add_axes([0.425, 0.260, 0.312, 0.705]),
    ]
    names = ("convective_mid", "convective_zavg")
    labels = (
        r"$2\pi l_c(z\simeq0.5)$",
        r"$2\pi\left\langle l_c\right\rangle_z$",
    )
    for axis, name, label in zip(axes, names, labels):
        style_axis(axis)
        plot_series(axis, records, colors, name)
        axis.set_xlabel(r"$t$")
        axis.set_ylabel(label)
        axis.set_ylim(bottom=0.0)
    add_case_legend(figure, records, colors, (0.770, 0.52))
    return save_figure(figure, output_dir, f"{ra_name}_velocity_length_scales_timeseries")


def make_individual_figure(record: dict) -> tuple[plt.Figure, plt.Axes]:
    # The wider outer canvas leaves enough room for a second y-axis without
    # changing the approved physical size or aspect ratio of the black frame.
    figure = plt.figure(figsize=(7.80, 5.84), facecolor="white")
    axis = figure.add_axes([0.160, 0.260, 0.660, 0.705])
    style_axis(axis)
    return figure, axis


def plot_individual_energy_mprime(
    record: dict,
    output_dir: Path,
) -> list[Path]:
    figure, axis = make_individual_figure(record)
    right_axis = axis.twinx()
    style_axis(right_axis)
    axis.yaxis.set_offset_position("left")
    right_axis.yaxis.set_offset_position("right")

    kinetic_rows = record["series"].get("kinetic", [])
    handles = []
    labels = []
    if kinetic_rows:
        kinetic = np.asarray(kinetic_rows, dtype=float)
        handles.append(
            axis.plot(
                kinetic[:, 0], kinetic[:, 1], color=KINETIC_RED, lw=LINE_WIDTH
            )[0]
        )
        labels.append(r"$K$")
    mprime_rows = record["series"].get("mprime", [])
    if mprime_rows:
        mprime = np.asarray(mprime_rows, dtype=float)
        handles.append(
            right_axis.plot(
                mprime[:, 0], mprime[:, 1], color=MSE_BLUE, lw=LINE_WIDTH
            )[0]
        )
        labels.append(
            r"$\left\langle\left\langle m'^2\right\rangle_{xy}\right\rangle_z$"
        )

    axis.set_xlabel(r"$t$")
    axis.set_ylabel(r"$K$")
    right_axis.set_ylabel(
        r"$\left\langle\left\langle m'^2\right\rangle_{xy}\right\rangle_z$"
    )
    axis.set_ylim(bottom=0.0)
    right_axis.set_ylim(bottom=0.0)
    add_plume_burst_windows(axis, record)
    if handles:
        axis.legend(
            handles,
            labels,
            frameon=False,
            loc="lower right",
            handlelength=2.4,
        )
    return save_figure(
        figure,
        output_dir,
        f"{individual_case_stem(record)}_kinetic_mprime_dual_axis_timeseries",
    )


def plot_individual_lm(
    record: dict,
    output_dir: Path,
) -> list[Path]:
    figure, axis = make_individual_figure(record)
    rows = record["series"].get("mse_spectral", [])
    if rows:
        values = np.asarray(rows, dtype=float)
        axis.plot(values[:, 0], values[:, 1], color=MSE_BLUE, lw=LINE_WIDTH)
    axis.set_xlabel(r"$t$")
    axis.set_ylabel(r"$2\pi L_m$")
    axis.set_ylim(bottom=0.0)
    add_plume_burst_windows(axis, record)
    return save_figure(
        figure,
        output_dir,
        f"{individual_case_stem(record)}_lm_timeseries",
    )


def write_tables(records: list[dict], duplicates: list[dict], output_dir: Path) -> None:
    metadata_fields = [
        "Pr",
        "Ra",
        "Ek",
        "AR",
        "beta",
        "gamma",
        "alphaqs",
        "tau",
        "qbot",
        "qtop",
        "Nx",
        "Ny",
        "Nz",
        "time_min",
        "time_max",
        "anchor_path",
        "source_run_paths",
        "mprime_sources",
        "historical_patches",
    ]
    with (output_dir / "selected_case_metadata.csv").open(
        "w", newline="", encoding="utf-8-sig"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=metadata_fields)
        writer.writeheader()
        for record in records:
            writer.writerow({field: record.get(field) for field in metadata_fields})

    with (output_dir / "excluded_lower_resolution_duplicates.csv").open(
        "w", newline="", encoding="utf-8-sig"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=metadata_fields)
        writer.writeheader()
        for record in duplicates:
            writer.writerow({field: record.get(field) for field in metadata_fields})

    long_fields = ["Ra", "Ek", "AR", "beta", "qbot", "Nx", "Ny", "Nz", "metric", "time", "value"]
    with (output_dir / "high_resolution_timeseries_long.csv").open(
        "w", newline="", encoding="utf-8-sig"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=long_fields)
        writer.writeheader()
        for record in records:
            base = {field: record.get(field) for field in long_fields[:8]}
            for metric in SERIES_NAMES:
                for row in record["series"].get(metric, []):
                    for component, value in enumerate(row[1:], start=1):
                        metric_name = metric if len(row) == 2 else f"{metric}_{component}"
                        writer.writerow({**base, "metric": metric_name, "time": row[0], "value": value})


def main() -> None:
    configure_style()
    output_dir = OUTPUT_ROOT / f"high_resolution_timeseries_latest_program_{datetime.now():%Y%m%d}"
    output_dir.mkdir(parents=True, exist_ok=True)
    individual_dir = output_dir / "individual_cases"
    individual_dir.mkdir(parents=True, exist_ok=True)
    for pattern in ("*.png", "*.pdf"):
        for old_figure in individual_dir.glob(pattern):
            old_figure.unlink()

    raw_records = extract_remote()
    selected, duplicates = choose_highest_resolution(raw_records)
    historical_patch_reports = apply_historical_patches(selected)
    selected = sorted(selected, key=lambda record: (float(record["Ra"]), record_sort_key(record)))
    write_tables(selected, duplicates, output_dir)

    groups: dict[float, list[dict]] = defaultdict(list)
    for record in selected:
        groups[float(record["Ra"])].append(record)

    figures = []
    individual_index = []
    for ra, records in sorted(groups.items()):
        records.sort(key=record_sort_key)
        colors = case_colors(records)
        ra_name = "Ra" + f"{ra:.0e}".replace("+", "").replace("e0", "e")
        figures.extend(plot_energy_mprime(ra_name, records, colors, output_dir))
        figures.extend(plot_l0(ra_name, records, colors, output_dir))
        figures.extend(plot_mse_scales(ra_name, records, colors, output_dir))
        figures.extend(plot_velocity_scales(ra_name, records, colors, output_dir))
        for record in records:
            case_figures = []
            case_figures.extend(
                plot_individual_energy_mprime(record, individual_dir)
            )
            case_figures.extend(
                plot_individual_lm(record, individual_dir)
            )
            figures.extend(case_figures)
            individual_index.append(
                {
                    "Ra": record.get("Ra"),
                    "Ek": record.get("Ek"),
                    "AR": record.get("AR"),
                    "beta": record.get("beta"),
                    "qbot": record.get("qbot"),
                    "Nx": record.get("Nx"),
                    "Ny": record.get("Ny"),
                    "Nz": record.get("Nz"),
                    "time_min": record.get("time_min"),
                    "time_max": record.get("time_max"),
                    "kinetic_mprime_dual_axis_png": str(case_figures[0]),
                    "lm_png": str(case_figures[2]),
                }
            )

    individual_fields = [
        "Ra",
        "Ek",
        "AR",
        "beta",
        "qbot",
        "Nx",
        "Ny",
        "Nz",
        "time_min",
        "time_max",
        "kinetic_mprime_dual_axis_png",
        "lm_png",
    ]
    with (individual_dir / "individual_case_figure_index.csv").open(
        "w", newline="", encoding="utf-8-sig"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=individual_fields)
        writer.writeheader()
        writer.writerows(individual_index)

    manifest = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "remote_host": REMOTE_HOST,
        "remote_root": REMOTE_ROOT,
        "selection": (
            "Actual bou.in Nx and Ny are both >=257. For identical physical "
            "parameters, retain the highest actual grid; break equal-grid ties "
            "by longest available physical time. Only cases with the latest "
            "diagnostics/thermo/mprime_square.dat and all three latest scale "
            "files are included."
        ),
        "definitions": {
            "kinetic": "K=0.5*avgvar one-based column 4",
            "num": "vertical trapezoidal average of nu_profiles.out column 4 divided by Delta_m=(b_bot-b_top)+gamma(q_bot-q_top)",
            "mprime": "< <(m-<m>xy)^2>xy >z with stretched-grid vertical weighting",
            "l0": "avgvar one-based columns 8 and 9 at z approximately 0.5 and 0.75",
            "mse_peak": "L_peak,m=2*pi/k_peak,m",
            "mse_spectral": "2*pi*sum(E_m)/sum(k_h E_m)",
            "mse_integral": "positive-lobe radial autocorrelation integral; legacy online output only",
            "convective_mid": "2*pi*sum(E_w)/sum(k_h E_w) at z approximately 0.5",
            "convective_zavg": "stretched-grid average of the per-height convective scale, multiplied by 2*pi",
            "velocity_z075": "2*pi*sum(E_w)/sum(k_h E_w) at z=0.75; legacy online output only",
        },
        "latest_program_only": LATEST_PROGRAM_ONLY,
        "historical_patches": historical_patch_reports,
        "selected_cases": len(selected),
        "excluded_lower_resolution_duplicates": len(duplicates),
        "groups": {str(ra): len(records) for ra, records in groups.items()},
        "individual_case_figure_directory": str(individual_dir),
        "individual_case_count": len(individual_index),
        "figures": [str(path) for path in figures],
    }
    (output_dir / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps({"output_dir": str(output_dir), **manifest}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
