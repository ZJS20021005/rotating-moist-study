from pathlib import Path
import math

import numpy as np
import pandas as pd

ROOT = Path(r"E:\moist RB\rotating_case_inventory\04_outputs_and_figures\high_resolution_timeseries_latest_program_20260805")
INPUT = ROOT / "high_resolution_timeseries_long.csv"
OUTPUT = ROOT / "latest_program_case_stability_summary.csv"
METRICS = ["num", "kinetic", "mprime"]


def label_ek(value):
    return "Nonrotating" if not np.isfinite(value) else f"{value:.3g}"


def analyse(frame):
    frame = frame.sort_values("time").drop_duplicates("time", keep="last")
    tmax = float(frame.time.max())
    tail_start = max(float(frame.time.min()), tmax - max(300.0, 0.30 * (tmax - float(frame.time.min()))))
    tail = frame[frame.time >= tail_start].copy()
    tail["block"] = np.floor((tail.time - tail_start) / 5.0).astype(int)
    tail = tail.groupby("block", as_index=False).agg(time=("time", "mean"), value=("value", "mean"))
    if len(tail) < 8:
        return dict(time_max=tmax, tail_start=tail_start, mean=np.nan, cv=np.nan,
                    trend_change=np.nan, half_change=np.nan, status="数据不足")
    time = tail.time.to_numpy(float); value = tail.value.to_numpy(float)
    mean = float(np.mean(value)); scale = max(abs(mean), 1e-30)
    slope = float(np.polyfit(time, value, 1)[0])
    trend_change = slope * (time[-1] - time[0]) / scale
    middle = len(value) // 2
    half_change = (float(np.mean(value[middle:])) - float(np.mean(value[:middle]))) / scale
    cv = float(np.std(value) / scale)
    drift = max(abs(trend_change), abs(half_change))
    if drift <= 0.08:
        status = "统计稳定（间歇波动）" if cv > 0.15 else "统计稳定"
    elif drift <= 0.15:
        status = "接近稳定"
    else:
        status = "仍明显演化"
    return dict(time_max=tmax, tail_start=tail_start, mean=mean, cv=cv,
                trend_change=trend_change, half_change=half_change, status=status)


def main():
    data = pd.read_csv(INPUT)
    keys = ["Ra", "Ek", "AR", "beta", "qbot", "Nx", "Ny", "Nz"]
    rows = []
    for key, case in data.groupby(keys, dropna=False):
        row = dict(zip(keys, key)); row["Ek_label"] = label_ek(float(row["Ek"]))
        row["time_max_all"] = float(case.time.max())
        for metric in METRICS:
            result = analyse(case[case.metric == metric])
            for name, value in result.items(): row[f"{metric}_{name}"] = value
        statuses = [row[f"{metric}_status"] for metric in METRICS]
        if all(s.startswith("统计稳定") for s in statuses): overall = "三项均统计稳定"
        elif any(s == "仍明显演化" for s in statuses): overall = "至少一项仍明显演化"
        else: overall = "接近联合稳定"
        row["overall_status"] = overall; rows.append(row)
    output = pd.DataFrame(rows).sort_values("Ek", na_position="last")
    output.to_csv(OUTPUT, index=False, encoding="utf-8-sig")
    print(OUTPUT)
    print(output[["Ek_label", "time_max_all", "num_status", "kinetic_status", "mprime_status", "overall_status"]].to_string(index=False))


if __name__ == "__main__": main()
