from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


ROOT = Path(r"E:\moist RB\rotating_case_inventory")
SOURCE = ROOT / "04_outputs_and_figures" / "high_resolution_timeseries_latest_program_20260805"
NUM_SOURCE = SOURCE / "num_timeseries" / "Ra8e6_latest_program_Num_timeseries_normalized.csv"
OUT = SOURCE / "force_num_vs_Ek_20260812"

FORCE_EK2E4 = Path(
    r"G:\moist convection\Ek2e-4\AR4\Beta1p02\qbot0p5_qtop0p004978"
    r"\N257x257x65\conti1\conti_strict_force_500\run\diagnostics"
    r"\force_balance_strict_fp3d_20260806\strict_force_balance_bulk_timeseries.csv"
)
FORCE_EK2E3 = Path(
    r"E:\moist RB\rotating_case_inventory\04_outputs_and_figures"
    r"\high_resolution_timeseries_latest_program_20260805\force_num_vs_Ek_20260812"
    r"\downloaded_force_balance\Ek2e-3\force_balance.out"
)

BURST_WINDOWS = {
    2.0e-4: [(353.3, 432.6), (622.9, 732.1), (939.4, 1029.8)],
    5.0e-4: [(390.5, 429.0), (451.8, 488.5)],
    7.0e-4: [(447.9, 469.6), (491.1, 523.1), (580.6, 684.6)],
}

FRAME = 4.5
LINE = 3.5
BLUE = (0.00, 0.28, 0.82)
YELLOW = (1.00, 0.72, 0.00)
GRAY = (0.35, 0.35, 0.35)
AX_RECT = [0.186, 0.260, 0.792, 0.705]

FORCE_STYLE = {
    "F_I": (r"$F_I$", (0.74, 0.14, 0.18), "o"),
    "F_C": (r"$F_C$", (0.11, 0.44, 0.71), "s"),
    "F_P": (r"$F_P$", (0.58, 0.05, 0.90), "D"),
    "F_V": (r"$F_V$", (0.96, 0.58, 0.19), "^"),
    "F_B": (r"$F_B$", (0.00, 0.65, 0.12), "v"),
    "F_T": (r"$F_T$", (0.20, 0.20, 0.20), "P"),
}


REMOTE_RUNS = {
    1.0e-1: "/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case/Pr0p7/Ra8e6/Ek1e-1/AR16/Beta1p02/qbot0p5_qtop0p004978/N385x385x65/conti1/conti_strict_force_500/run",
    5.0e-2: "/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case/Pr0p7/Ra8e6/Ek5e-2/AR16/Beta1p02/qbot0p5_qtop0p004978/N385x385x65/conti_strict_force_500/run",
    3.0e-2: "/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case/Pr0p7/Ra8e6/Ek3e-2/AR16/Beta1p02/qbot0p5_qtop0p004978/N385x385x65/conti2/conti1/conti_strict_force_500/run",
    1.0e-2: "/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case/Pr0p7/Ra8e6/Ek1e-2/AR16/Beta1p02/qbot0p5_qtop0p004978/N257x257x65/conti2/conti1/conti_strict_force_500/run",
    7.0e-3: "/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case/Pr0p7/Ra8e6/Ek7e-3/AR16/Beta1p02/qbot0p5_qtop0p004978/N257x257x65/conti2/conti1/conti_strict_force_500/run",
    5.0e-3: "/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case/Pr0p7/Ra8e6/Ek5e-3/AR16/Beta1p02/qbot0p5_qtop0p004978/N257x257x65/conti2/conti1/conti_strict_force_500/run",
    3.0e-3: "/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case/Pr0p7/Ra8e6/Ek3e-3/AR16/Beta1p02/qbot0p5_qtop0p004978/N257x257x65/conti3/conti_strict_force_500/run",
    2.0e-3: "/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case/Pr0p7/Ra8e6/Ek2e-3/AR16/Beta1p02/qbot0p5_qtop0p004978/N385x385x65/conti_strict_force_500/run",
    1.0e-3: "/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case/Pr0p7/Ra8e6/Ek1e-3/AR16/Beta1p02/qbot0p5_qtop0p004978/N257x257x65/conti1/conti_strict_force_500/run",
    7.0e-4: "/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case/Pr0p7/Ra8e6/Ek7e-4/AR16/Beta1p02/qbot0p5_qtop0p004978/N257x257x65/conti_strict_force_500/run",
    5.0e-4: "/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case/Pr0p7/Ra8e6/Ek5e-4/AR16/Beta1p02/qbot0p5_qtop0p004978/N257x257x65/conti_strict_force_500/run",
    2.0e-4: "/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case/Pr0p7/Ra8e6/Ek2e-4/AR4/Beta1p02/qbot0p5_qtop0p004978/N257x257x65/conti1/conti_strict_force_500/run",
    1.5e-4: "/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case/Pr0p7/Ra8e6/Ek1p5e-4/AR4/Beta1p02/qbot0p5_qtop0p004978/N257x257x65/conti1/conti_strict_force_500/run",
    np.nan: "/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case/Pr0p7/Ra8e6/norotating/AR16/Beta1p02/qbot0p5_qtop0p004978/N257x257x65/conti1/conti_strict_force_500/run",
}


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Times New Roman",
            "mathtext.fontset": "stix",
            "axes.linewidth": FRAME,
            "axes.labelsize": 24,
            "xtick.labelsize": 13,
            "ytick.labelsize": 13,
            "legend.fontsize": 10,
            "figure.dpi": 180,
            "savefig.dpi": 300,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def make_axis() -> tuple[plt.Figure, plt.Axes]:
    fig = plt.figure(figsize=(6.5, 5.84), facecolor="white")
    ax = fig.add_axes(AX_RECT)
    for spine in ax.spines.values():
        spine.set_linewidth(FRAME)
    ax.tick_params(which="major", direction="in", length=12, width=1.2, top=True, right=True)
    ax.minorticks_on()
    ax.tick_params(which="minor", direction="in", length=6, width=1.0, top=True, right=True)
    ax.set_box_aspect(5.2 / 6.5)
    return fig, ax


def save(fig: plt.Figure, stem: str) -> None:
    fig.savefig(OUT / f"{stem}.png", facecolor="white")
    fig.savefig(OUT / f"{stem}.pdf", facecolor="white")
    plt.close(fig)


def load_num_summary() -> tuple[pd.DataFrame, float]:
    data = pd.read_csv(NUM_SOURCE)
    data["Ek"] = pd.to_numeric(data["Ek"], errors="coerce")
    data["time"] = pd.to_numeric(data["time"], errors="coerce")
    data["value"] = pd.to_numeric(data["value"], errors="coerce")
    rows = []
    baseline = np.nan
    for ek, group in data.groupby("Ek", dropna=False):
        group = group.sort_values("time").drop_duplicates("time", keep="last")
        if not np.isfinite(ek):
            tmax = float(group["time"].max())
            chosen = group[group["time"] >= tmax - 200.0]
            baseline = float(chosen["value"].mean())
            rows.append(
                {"Ek": np.nan, "Num": baseline, "selection": "nonrotating late 200", "time_windows": f"{tmax-200.0:g}-{tmax:g}", "samples": len(chosen)}
            )
            continue
        matching = next((windows for target, windows in BURST_WINDOWS.items() if np.isclose(ek, target)), None)
        if matching:
            masks = [(group["time"] >= start) & (group["time"] <= stop) for start, stop in matching]
            selected = group[np.logical_or.reduce(masks)]
            selection = "plume-burst windows"
            windows_text = "; ".join(f"{start:g}-{stop:g}" for start, stop in matching)
        else:
            tmax = float(group["time"].max())
            selected = group[group["time"] >= tmax - 200.0]
            selection = "late 200"
            windows_text = f"{max(float(group['time'].min()), tmax-200.0):g}-{tmax:g}"
        rows.append(
            {"Ek": float(ek), "Num": float(selected["value"].mean()), "selection": selection, "time_windows": windows_text, "samples": len(selected)}
        )
    return pd.DataFrame(rows), baseline


def plot_num(summary: pd.DataFrame, baseline: float) -> None:
    finite = summary[np.isfinite(summary["Ek"])].sort_values("Ek")
    normal = finite[finite["selection"].ne("plume-burst windows")]
    burst = finite[finite["selection"].eq("plume-burst windows")]
    fig, ax = make_axis()
    ax.plot(finite["Ek"], finite["Num"], color=BLUE, lw=2.2, zorder=1)
    ax.scatter(normal["Ek"], normal["Num"], s=100, marker="o", facecolor=BLUE, edgecolor="black", linewidth=1.5, zorder=3)
    ax.scatter(burst["Ek"], burst["Num"], s=110, marker="o", facecolor=YELLOW, edgecolor="black", linewidth=1.5, zorder=4)
    ax.axhline(baseline, color=GRAY, ls="--", lw=2.5, zorder=0)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1.2e-4, 1.3e-1)
    ax.set_ylim(0.9 * min(finite["Num"].min(), baseline), 1.18 * max(finite["Num"].max(), baseline))
    ax.set_xlabel(r"$Ek$")
    ax.set_ylabel(r"$Nu_m$")
    handles = [
        Line2D([], [], ls="None", marker="o", ms=9, mfc=BLUE, mec="black", mew=1.2, label="Late-time mean"),
        Line2D([], [], ls="None", marker="o", ms=9, mfc=YELLOW, mec="black", mew=1.2, label="Plume-burst phase"),
        Line2D([], [], color=GRAY, ls="--", lw=2.5, label="Nonrotating"),
    ]
    ax.legend(handles=handles, frameon=False, loc="upper left")
    save(fig, "Ra8e6_Num_vs_Ek_plume_burst_windows")


def load_available_force_points() -> pd.DataFrame:
    rows = []
    if FORCE_EK2E4.exists():
        frame = pd.read_csv(FORCE_EK2E4)
        row = {"Ek": 2.0e-4, "time_min": frame["time"].min(), "time_max": frame["time"].max(), "samples": len(frame), "status": "complete strict offline"}
        for field in FORCE_STYLE:
            row[field] = pd.to_numeric(frame[field], errors="coerce").mean()
        rows.append(row)
    if FORCE_EK2E3.exists():
        names = ["time", "F_I", "F_C", "F_P", "F_V", "F_B", "F_T"] + [f"unused_{i}" for i in range(27)]
        frame = pd.read_csv(FORCE_EK2E3, sep=r"\s+", header=None, names=names)
        row = {"Ek": 2.0e-3, "time_min": frame["time"].min(), "time_max": frame["time"].max(), "samples": len(frame), "status": "partial online strict"}
        for field in FORCE_STYLE:
            row[field] = pd.to_numeric(frame[field], errors="coerce").mean()
        rows.append(row)
    return pd.DataFrame(rows).sort_values("Ek")


def plot_force(points: pd.DataFrame) -> None:
    fig, ax = make_axis()
    ax.axvspan(2.0e-4, 7.0e-4, color=YELLOW, alpha=0.10, lw=0, zorder=0)
    for field, (label, color, marker) in FORCE_STYLE.items():
        valid = np.isfinite(points[field]) & (points[field] > 0.0)
        ax.plot(
            points.loc[valid, "Ek"], points.loc[valid, field], ls="None", marker=marker,
            ms=9, mfc=color, mec="black", mew=1.1, color=color, label=label, zorder=3,
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1.2e-4, 3.0e-3)
    ax.set_xlabel(r"$Ek$")
    ax.set_ylabel(r"Force magnitude")
    ax.legend(frameon=False, ncol=2, loc="best")
    save(fig, "Ra8e6_strict_forces_vs_Ek_currently_available")


def write_download_checklist(points: pd.DataFrame) -> pd.DataFrame:
    available = {float(value) for value in points.loc[points["status"].eq("complete strict offline"), "Ek"]}
    rows = []
    for ek, run in REMOTE_RUNS.items():
        finite = np.isfinite(ek)
        if finite and float(ek) in available:
            status = "complete locally"
            action = "none"
        elif finite and np.isclose(float(ek), 2.0e-3):
            status = "partial locally: t=660.1-684.2 only"
            action = "download updated complete files"
        else:
            status = "strict bulk time series missing locally"
            action = "download"
        rows.append(
            {
                "case": "norotating" if not finite else f"Ek={float(ek):.6g}",
                "status": status,
                "remote_run_directory": run,
                "required_file": f"{run}/data/force_balance.out",
                "recommended_profile_file": f"{run}/data/force_balance_z.out",
                "action": action,
                "note": "pressure_force.out is legacy horizontal pressure only and cannot replace strict F_P=-grad(p)",
            }
        )
    table = pd.DataFrame(rows)
    table.to_csv(OUT / "strict_force_balance_download_checklist.csv", index=False, encoding="utf-8-sig")
    missing = table[table["action"].ne("none")]
    lines = [
        "# 严格力平衡下载清单（2026-08-12）",
        "",
        "最终跨 case 图需要同一续算段中的 `data/force_balance.out`。建议同时下载 `data/force_balance_z.out`，以后可直接画高度-时间图。",
        "",
        "`pressure_force.out` 只有旧定义的水平压力梯度，不能代替包含垂直压力梯度的严格三维 `F_P=-grad(p)`。",
        "",
        f"当前完整：Ek=2e-4；当前部分可用：Ek=2e-3（仅 t=660.1-684.2）；需要下载/刷新：{len(missing)} 个 case。",
        "",
    ]
    for row in missing.itertuples(index=False):
        lines += [
            f"## {row.case}",
            "",
            f"远端 run：`{row.remote_run_directory}`",
            "",
            f"必须下载：`{row.required_file}`",
            "",
            f"建议同时下载：`{row.recommended_profile_file}`",
            "",
        ]
    (OUT / "strict_force_balance_download_checklist.md").write_text("\n".join(lines), encoding="utf-8-sig")
    return table


def main() -> None:
    configure_style()
    OUT.mkdir(parents=True, exist_ok=True)
    num_summary, baseline = load_num_summary()
    num_summary.to_csv(OUT / "Ra8e6_Num_vs_Ek_values_and_windows.csv", index=False, encoding="utf-8-sig")
    plot_num(num_summary, baseline)
    force_points = load_available_force_points()
    force_points.to_csv(OUT / "Ra8e6_strict_force_points_currently_available.csv", index=False, encoding="utf-8-sig")
    plot_force(force_points)
    checklist = write_download_checklist(force_points)
    print(OUT)
    print(num_summary.to_string(index=False))
    print(force_points.to_string(index=False))
    print(checklist[["case", "status", "action"]].to_string(index=False))


if __name__ == "__main__":
    main()
