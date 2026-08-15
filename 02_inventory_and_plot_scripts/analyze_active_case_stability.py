from __future__ import annotations

import json
import math
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(r"E:\moist RB\rotating_case_inventory")
ACTIVE_JSON = ROOT / "03_inventory_tables" / "active_running_cases_latest.json"
SERIES_PARENT = ROOT / "04_outputs_and_figures"
SERIES_DIR = max(
    (
        path
        for path in SERIES_PARENT.glob("high_resolution_timeseries_latest_program_*")
        if (path / "high_resolution_timeseries_long.csv").exists()
    ),
    key=lambda path: path.name,
)
LONG_CSV = SERIES_DIR / "high_resolution_timeseries_long.csv"
OUTPUT_CSV = SERIES_DIR / "active_case_stability_latest.csv"
OUTPUT_JSON = SERIES_DIR / "active_case_stability_latest.json"

TAIL_LONG = 200.0
TAIL_SHORT = 100.0
BLOCK_WIDTH = 10.0


def close(left: float, right: float, rtol: float = 2.0e-6) -> bool:
    if not (math.isfinite(left) and math.isfinite(right)):
        return False
    return abs(left - right) <= rtol * max(1.0, abs(left), abs(right))


def metric_stats(frame: pd.DataFrame, duration: float) -> dict:
    if frame.empty:
        return {
            "mean": math.nan,
            "normalized_drift": math.nan,
            "half_shift": math.nan,
            "block_cv": math.nan,
            "window_start": math.nan,
            "window_end": math.nan,
            "n_blocks": 0,
        }
    frame = frame.sort_values("time")
    end = float(frame["time"].max())
    start = max(float(frame["time"].min()), end - duration)
    tail = frame.loc[frame["time"] >= start].copy()
    tail["block"] = np.floor((tail["time"] - start) / BLOCK_WIDTH).astype(int)
    blocked = tail.groupby("block", as_index=False).agg(time=("time", "mean"), value=("value", "mean"))
    blocked = blocked[np.isfinite(blocked["value"])].copy()
    if len(blocked) < 2:
        return {
            "mean": float(tail["value"].mean()),
            "normalized_drift": math.nan,
            "half_shift": math.nan,
            "block_cv": math.nan,
            "window_start": start,
            "window_end": end,
            "n_blocks": len(blocked),
        }
    x = blocked["time"].to_numpy(dtype=float)
    y = blocked["value"].to_numpy(dtype=float)
    mean = float(np.mean(y))
    slope = float(np.polyfit(x, y, 1)[0])
    span = float(x.max() - x.min())
    drift = slope * span / mean if mean != 0.0 else math.nan
    midpoint = start + 0.5 * (end - start)
    first = tail.loc[tail["time"] < midpoint, "value"]
    second = tail.loc[tail["time"] >= midpoint, "value"]
    half_shift = (
        (float(second.mean()) - float(first.mean())) / mean
        if len(first) and len(second) and mean != 0.0
        else math.nan
    )
    cv = float(np.std(y, ddof=1) / abs(mean)) if len(y) > 1 and mean != 0.0 else math.nan
    return {
        "mean": mean,
        "normalized_drift": drift,
        "half_shift": half_shift,
        "block_cv": cv,
        "window_start": start,
        "window_end": end,
        "n_blocks": len(blocked),
    }


def absmax(*values: float) -> float:
    finite = [abs(float(value)) for value in values if value is not None and math.isfinite(float(value))]
    return max(finite) if finite else math.nan


def classify(k_long: dict, k_short: dict, m_long: dict, lm_long: dict) -> tuple[str, str]:
    k_amp = k_short["mean"]
    m_amp = m_long["mean"]
    if math.isfinite(k_amp) and k_amp < 1.0e-6:
        return "未进入非线性稳态", "动能仍处于极小振幅，不能把平坦曲线当作成熟稳态"

    k200 = absmax(k_long["normalized_drift"], k_long["half_shift"])
    k100 = absmax(k_short["normalized_drift"], k_short["half_shift"])
    m200 = absmax(m_long["normalized_drift"], m_long["half_shift"])
    lm200 = absmax(lm_long["normalized_drift"], lm_long["half_shift"])

    if all(math.isfinite(value) and value <= 0.10 for value in (k200, m200, lm200)):
        return "已基本稳定", "最近200时间内K、m'^2和2πLm的系统漂移均不超过约10%"
    if (
        math.isfinite(k200)
        and k200 <= 0.10
        and math.isfinite(m200)
        and m200 <= 0.20
        and math.isfinite(lm200)
        and lm200 <= 0.20
    ):
        return "快稳定", "动能已成平台，湿静能方差和尺度仅有弱慢漂移"
    if (
        math.isfinite(k100)
        and k100 <= 0.10
        and math.isfinite(m200)
        and m200 <= 0.25
        and math.isfinite(lm200)
        and lm200 <= 0.25
    ):
        return "可能接近稳定", "最近100时间的动能已变平，但结构量仍需更长窗口确认"
    return "仍明显演化", "最近窗口仍有显著系统漂移或重组"


def main() -> None:
    active_payload = json.loads(ACTIVE_JSON.read_text(encoding="utf-8"))
    active = [row for row in active_payload["records"] if row.get("active_during_observation")]
    data = pd.read_csv(LONG_CSV)
    rows = []
    for case in active:
        if not case.get("latest_program"):
            continue
        mask = (
            np.isclose(data["Ra"], float(case["Ra"]), rtol=2e-6)
            & np.isclose(data["AR"], float(case["AR"]), rtol=2e-6)
            & np.isclose(data["beta"], float(case["beta"]), rtol=2e-6)
            & np.isclose(data["qbot"], float(case["qbot"]), rtol=2e-6)
            & (data["Nx"] == int(case["Nx"]))
            & (data["Ny"] == int(case["Ny"]))
            & (data["Nz"] == int(case["Nz"]))
        )
        if case.get("Ek") is None:
            mask &= data["Ek"].isna()
        else:
            mask &= np.isclose(data["Ek"], float(case["Ek"]), rtol=2e-6)
        subset = data.loc[mask].copy()

        stats = {}
        for metric in ("kinetic", "mprime", "mse_spectral"):
            metric_frame = subset.loc[subset["metric"] == metric, ["time", "value"]]
            stats[(metric, "200")] = metric_stats(metric_frame, TAIL_LONG)
            stats[(metric, "100")] = metric_stats(metric_frame, TAIL_SHORT)

        status, reason = classify(
            stats[("kinetic", "200")],
            stats[("kinetic", "100")],
            stats[("mprime", "200")],
            stats[("mse_spectral", "200")],
        )
        row = {
            "Ra": case["Ra"],
            "Ek": case.get("Ek"),
            "Ek_label": "NR" if case.get("Ek") is None else f"{case['Ek']:.6g}",
            "AR": case["AR"],
            "beta": case["beta"],
            "qbot": case["qbot"],
            "Nx": case["Nx"],
            "Ny": case["Ny"],
            "Nz": case["Nz"],
            "time_at_active_scan": case.get("latest_diagnostic_time"),
            "time_in_reduced_series": float(subset["time"].max()) if not subset.empty else math.nan,
            "seconds_since_write_at_scan": case.get("seconds_since_latest_write"),
            "status": status,
            "reason": reason,
            "run_path": case["run_path"],
        }
        for metric, short_name in (
            ("kinetic", "K"),
            ("mprime", "mprime2"),
            ("mse_spectral", "two_pi_Lm"),
        ):
            for window in ("100", "200"):
                values = stats[(metric, window)]
                for name in ("mean", "normalized_drift", "half_shift", "block_cv", "window_start", "window_end", "n_blocks"):
                    row[f"{short_name}_{name}_{window}"] = values[name]
        rows.append(row)

    result = pd.DataFrame(rows)
    order = {"已基本稳定": 0, "快稳定": 1, "可能接近稳定": 2, "仍明显演化": 3, "未进入非线性稳态": 4}
    result["_order"] = result["status"].map(order).fillna(99)
    result["_ek_order"] = result["Ek"].fillna(1.0)
    result = result.sort_values(["_order", "_ek_order"]).drop(columns=["_order", "_ek_order"])
    result.to_csv(OUTPUT_CSV, index=False, encoding="utf-8-sig")
    payload = {
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "active_scan": str(ACTIVE_JSON),
        "series_source": str(LONG_CSV),
        "criterion": {
            "block_width": BLOCK_WIDTH,
            "long_window": TAIL_LONG,
            "short_window": TAIL_SHORT,
            "normalized_drift": "linear-fit change across the window divided by the window mean",
            "half_shift": "second-half mean minus first-half mean, divided by the window mean",
            "stable": "max absolute drift/half-shift <=10% for K, m'^2, and 2*pi*Lm",
            "near": "K <=10%, m'^2 and 2*pi*Lm <=20% over 200 time units",
        },
        "records": result.replace({np.nan: None}).to_dict(orient="records"),
    }
    OUTPUT_JSON.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    print(result[[
        "Ek_label", "time_in_reduced_series", "status",
        "K_normalized_drift_200", "mprime2_normalized_drift_200",
        "two_pi_Lm_normalized_drift_200",
    ]].to_string(index=False))
    print(f"saved {OUTPUT_CSV}")
    print(f"saved {OUTPUT_JSON}")


if __name__ == "__main__":
    main()
