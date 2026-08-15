#!/usr/bin/env python3
"""Reduce current Ra=8e6 continuation runs to Re and Num time series.

Run on c01n0011.  Only text diagnostics are read remotely:
  data/avgvar.out
  data/nu_profiles.out
  bou.in

Re = sqrt(<u^2+v^2+w^2>) * sqrt(Ra/Pr)
Num(t) = integral_z [nu_profiles.out column 4] dz / Delta_m

The last 200 physical-time units are summarized with 10-unit block averages.
No scheduler command is called.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path


DEFAULT_ROOT = Path(
    "/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case/Pr0p7/Ra8e6"
)
DEFAULT_OUTPUT = Path(
    "/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case/postprocess/re_num_20260805"
)


def ff(text: str) -> float:
    return float(text.replace("D", "e").replace("d", "e"))


def numeric_rows(path: Path, min_columns: int) -> list[list[float]]:
    rows: list[list[float]] = []
    if not path.exists():
        return rows
    for line in path.read_text(errors="ignore").splitlines():
        parts = line.split()
        if len(parts) < min_columns:
            continue
        try:
            values = [ff(item) for item in parts]
        except ValueError:
            continue
        if all(math.isfinite(value) for value in values):
            rows.append(values)
    return rows


def parse_bou(path: Path) -> dict[str, float]:
    lines = path.read_text(errors="ignore").splitlines()
    out: dict[str, float] = {}
    for i, line in enumerate(lines):
        if "Ra" in line and "Pr" in line and "invRo" in line and i + 1 < len(lines):
            vals = lines[i + 1].split("!")[0].split()
            out["Ra"], out["Pr"], out["invRo"] = map(ff, vals[:3])
        if "A_stopmod" in line and "qvapbot" in line and i + 1 < len(lines):
            vals = lines[i + 1].split("!")[0].split()
            out["dsaltop"], out["dsalbot"], out["qtop"], out["qbot"] = (
                ff(vals[4]),
                ff(vals[5]),
                ff(vals[6]),
                ff(vals[7]),
            )
    out["Delta_m"] = (
        out["dsalbot"]
        - out["dsaltop"]
        + out["gamma"] * (out["qbot"] - out["qtop"])
        if "gamma" in out
        else 0.0
    )
    # gamma is fixed by the maintained Ra8e6 setup; read it from the control row.
    for i, line in enumerate(lines):
        if "Ra" in line and "Pr" in line and "invRo" in line and i + 1 < len(lines):
            vals = lines[i + 1].split("!")[0].split()
            out["gamma"] = ff(vals[4])
            break
    out["Delta_m"] = (
        out["dsalbot"] - out["dsaltop"]
        + out["gamma"] * (out["qbot"] - out["qtop"])
    )
    return out


def ek_from_path(path: Path) -> float:
    match = re.search(r"/Ek([0-9epP+\-\.]+)", str(path).replace("\\", "/"))
    if not match:
        return math.nan
    return ff(match.group(1).replace("p", "."))


def block_stats(times: list[float], values: list[float], window: float = 200.0) -> dict[str, float | str]:
    if not times:
        return {"status": "no_data", "mean": math.nan}
    pairs = sorted(zip(times, values))
    end = pairs[-1][0]
    start = max(pairs[0][0], end - window)
    tail = [(t, v) for t, v in pairs if t >= start]
    blocks: dict[int, list[float]] = {}
    for t, v in tail:
        blocks.setdefault(int(math.floor((t - start) / 10.0)), []).append(v)
    bx = []
    by = []
    for key in sorted(blocks):
        bx.append(start + (key + 0.5) * 10.0)
        by.append(sum(blocks[key]) / len(blocks[key]))
    mean = sum(by) / len(by) if by else math.nan
    if len(by) >= 2 and mean != 0.0:
        slope = sum((x - sum(bx) / len(bx)) * (y - mean) for x, y in zip(bx, by)) / sum(
            (x - sum(bx) / len(bx)) ** 2 for x in bx
        )
        drift = slope * (max(bx) - min(bx)) / mean
        mid = start + 0.5 * (end - start)
        left = [v for t, v in tail if t < mid]
        right = [v for t, v in tail if t >= mid]
        half_shift = (sum(right) / len(right) - sum(left) / len(left)) / mean if left and right else math.nan
        cv = (sum((v - mean) ** 2 for v in by) / (len(by) - 1)) ** 0.5 / abs(mean)
    else:
        drift = half_shift = cv = math.nan
    stable = all(math.isfinite(x) and abs(x) <= 0.10 for x in (drift, half_shift))
    return {
        "status": "stable" if stable else "not_stable",
        "mean": mean,
        "drift_200": drift,
        "half_shift_200": half_shift,
        "block_cv_200": cv,
        "window_start": start,
        "window_end": end,
        "n_blocks": len(by),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    ap.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = ap.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, object]] = []
    summaries: list[dict[str, object]] = []
    runs = sorted(args.root.glob("**/conti_strict_force_500/run"))
    # Use the highest-resolution qbot=0.5 nonrotating baseline.  Its current
    # branch is base + conti1 and has no Ek directory token.
    norotating = sorted(
        args.root.glob(
            "norotating/AR16/Beta1p02/qbot0p5_qtop0p004978/N*/conti1/run"
        )
    )
    if norotating:
        runs.append(norotating[-1])
    for run in runs:
        bou = run / "bou.in"
        if not bou.exists():
            continue
        params = parse_bou(bou)
        if not params or params.get("Delta_m", 0.0) <= 0.0:
            continue
        ek = ek_from_path(run)
        # Merge the base run and all continuation segments belonging to the
        # same actual grid/case root.  Later/deeper continuation paths win at
        # duplicate physical times.
        case_root = run.parent.parent
        branch_runs = sorted(
            (
                item
                for item in case_root.glob("**/run")
                if (item / "bou.in").exists()
                and "conti_strict_force_500" not in item.parts
            ),
            key=lambda item: (len(item.parts), str(item)),
        )
        branch_runs.append(run)

        avg_by_time: dict[float, list[float]] = {}
        profile_by_time: dict[float, dict[float, float]] = {}
        for branch in branch_runs:
            for item in numeric_rows(branch / "data" / "avgvar.out", 4):
                avg_by_time[round(item[0], 8)] = item
            for item in numeric_rows(branch / "data" / "nu_profiles.out", 5):
                profile_by_time.setdefault(round(item[0], 8), {})[item[1]] = item[3]

        avg = [avg_by_time[key] for key in sorted(avg_by_time)]
        re_times: list[float] = []
        re_values: list[float] = []
        for item in avg:
            # avgvar column 4 (one-based) is <u^2+v^2+w^2>_V.
            re_times.append(item[0])
            re_values.append(math.sqrt(max(item[3], 0.0)) * math.sqrt(params["Ra"] / params["Pr"]))

        num_times: list[float] = []
        num_values: list[float] = []
        for time in sorted(profile_by_time):
            profile = sorted(profile_by_time[time].items())
            if len(profile) < 2:
                continue
            integral = sum(
                0.5 * (left[1] + right[1]) * (right[0] - left[0])
                for left, right in zip(profile[:-1], profile[1:])
            )
            height = profile[-1][0] - profile[0][0]
            if height > 0.0:
                num_times.append(time)
                num_values.append((integral / height) / params["Delta_m"])

        for metric, times, values in (
            ("Re", re_times, re_values),
            ("Num", num_times, num_values),
        ):
            for time, value in zip(times, values):
                rows.append(
                    {
                        "Ek": "" if not math.isfinite(ek) else ek,
                        "Ra": params["Ra"],
                        "Pr": params["Pr"],
                        "run": str(run),
                        "metric": metric,
                        "time": time,
                        "value": value,
                    }
                )

        re_stats = block_stats(re_times, re_values)
        num_stats = block_stats(num_times, num_values)
        summaries.append(
            {
                "Ek": "" if not math.isfinite(ek) else ek,
                "Ra": params["Ra"],
                "Pr": params["Pr"],
                "run": str(run),
                "branch_runs": ";".join(str(item) for item in branch_runs),
                "time_Re_max": max(re_times) if re_times else math.nan,
                "time_Num_max": max(num_times) if num_times else math.nan,
                "Re_status": re_stats["status"],
                "Re_mean_200": re_stats.get("mean", math.nan),
                "Re_drift_200": re_stats.get("drift_200", math.nan),
                "Re_half_shift_200": re_stats.get("half_shift_200", math.nan),
                "Num_status": num_stats["status"],
                "Num_mean_200": num_stats.get("mean", math.nan),
                "Num_drift_200": num_stats.get("drift_200", math.nan),
                "Num_half_shift_200": num_stats.get("half_shift_200", math.nan),
                "stable_for_Ek_scaling": re_stats["status"] == "stable" and num_stats["status"] == "stable",
            }
        )

    with (args.output / "re_num_timeseries.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["Ek", "Ra", "Pr", "run", "metric", "time", "value"])
        writer.writeheader()
        writer.writerows(rows)
    fields = list(summaries[0].keys()) if summaries else ["Ek"]
    with (args.output / "re_num_stability_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(summaries)
    print(f"runs={len(summaries)} rows={len(rows)} output={args.output}")


if __name__ == "__main__":
    main()
