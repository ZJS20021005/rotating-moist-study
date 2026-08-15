from __future__ import annotations

import json
import math
import subprocess
from datetime import datetime
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


HOST = "c01n0006"
PROJECT_ROOT = Path(r"E:\moist RB\rotating_case_inventory")
OUT_ROOT = PROJECT_ROOT / "04_outputs_and_figures" / "case_timeseries"

CASES = [
    {
        "key": "Ra8e6_Ek3e-3_AR16_Beta1p02_qbot0p5",
        "title": r"$Ra=8\times10^6,\ Ek=3\times10^{-3},\ \Gamma=16,\ \beta=1.02,\ q_{\rm bot}=0.5$",
        "short": r"$Ek=3\times10^{-3},\ q_{\rm bot}=0.5$",
        "color": (0.00, 0.45, 1.00),
        "remote_run": "/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case/Pr0p7/Ra8e6/Ek3e-3/AR16/Beta1p02/qbot0p5_qtop0p004978/N257x257x65/run",
        "plot_lengths": True,
    },
    {
        "key": "Ra8e6_AR10_norotating_Beta1p02_botsat",
        "title": r"$Ra=8\times10^6,\ {\rm nonrotating},\ \Gamma=10,\ \beta=1.02$",
        "short": r"${\rm nonrotating},\ \Gamma=10$",
        "color": (0.22, 0.22, 0.22),
        "remote_run": "/share/org/SHUTUANL/shu_zhangjs/rainy model/ns/transition_study/beta_not1/Ra8e6/AR10/norotating/Beta1p02/run",
        "plot_lengths": False,
    },
]


def apply_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Times New Roman",
            "mathtext.fontset": "stix",
            "axes.linewidth": 4.5,
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


def style_axis(ax) -> None:
    for spine in ax.spines.values():
        spine.set_linewidth(4.5)
    ax.tick_params(direction="in", length=12, width=1.2, top=True, right=True)
    ax.minorticks_on()
    ax.tick_params(which="minor", direction="in", length=6, width=1.0, top=True, right=True)
    ax.set_box_aspect(5.2 / 6.5)


def fetch_remote_file(remote_path: str, local_path: Path) -> bool:
    local_path.parent.mkdir(parents=True, exist_ok=True)
    code = (
        "from pathlib import Path\n"
        "import sys\n"
        f"p=Path({remote_path!r})\n"
        "sys.exit(2) if not p.exists() or p.stat().st_size==0 else None\n"
        "sys.stdout.buffer.write(p.read_bytes())\n"
    )
    proc = subprocess.run(
        ["ssh", HOST, "python3", "-"],
        input=code.encode("utf-8"),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if proc.returncode != 0:
        return False
    local_path.write_bytes(proc.stdout)
    return True


def read_table(path: Path) -> np.ndarray:
    rows = []
    with path.open("r", errors="ignore") as f:
        for line in f:
            parts = line.split()
            if not parts:
                continue
            vals = []
            ok = True
            for p in parts:
                try:
                    vals.append(float(p.replace("D", "E").replace("d", "E")))
                except ValueError:
                    ok = False
                    break
            if ok and all(math.isfinite(v) for v in vals):
                rows.append(vals)
    if not rows:
        return np.empty((0, 0), dtype=float)
    ncol = max(len(r) for r in rows)
    out = np.full((len(rows), ncol), np.nan, dtype=float)
    for i, r in enumerate(rows):
        out[i, : len(r)] = r
    return out


def deduplicate_time(arr: np.ndarray) -> np.ndarray:
    if arr.size == 0:
        return arr
    by_t: dict[float, np.ndarray] = {}
    for row in arr:
        if not np.isfinite(row[0]):
            continue
        by_t[round(float(row[0]), 10)] = row
    if not by_t:
        return np.empty((0, arr.shape[1]), dtype=float)
    return np.vstack([by_t[k] for k in sorted(by_t)])


def plot_energy_mprime(case: dict, data: dict[str, np.ndarray], out_dir: Path) -> list[Path]:
    avg = data.get("avgvar", np.empty((0, 0)))
    mse = data.get("mse_aggregation", np.empty((0, 0)))
    avg = deduplicate_time(avg)
    mse = deduplicate_time(mse)

    fig = plt.figure(figsize=(13.0, 5.84), facecolor="white")
    ax1 = fig.add_axes([0.085, 0.260, 0.365, 0.705])
    ax2 = fig.add_axes([0.575, 0.260, 0.365, 0.705])
    for ax in (ax1, ax2):
        style_axis(ax)

    color = case["color"]
    if avg.shape[1] >= 4:
        t = avg[:, 0]
        K = 0.5 * avg[:, 3]
        mask = np.isfinite(t) & np.isfinite(K) & (K >= 0)
        ax1.plot(t[mask], K[mask], color=color, lw=3.5)
    if mse.shape[1] >= 3:
        t = mse[:, 0]
        sqrt_am = mse[:, 2]
        mask = np.isfinite(t) & np.isfinite(sqrt_am) & (sqrt_am >= 0)
        ax2.plot(t[mask], sqrt_am[mask], color=color, lw=3.5)

    ax1.set_xlabel(r"$t$")
    ax1.set_ylabel(r"$K$")
    ax2.set_xlabel(r"$t$")
    ax2.set_ylabel(r"$\left\langle\left\langle m'^2\right\rangle_{xy}\right\rangle_z^{1/2}$")
    ax1.set_title(case["short"], fontsize=13, pad=8)
    ax2.set_title(case["short"], fontsize=13, pad=8)

    png = out_dir / f"{case['key']}_K_mprime_timeseries.png"
    pdf = out_dir / f"{case['key']}_K_mprime_timeseries.pdf"
    fig.savefig(png)
    fig.savefig(pdf)
    plt.close(fig)

    # Save reduced CSVs next to the figure.
    if avg.shape[1] >= 4:
        pd.DataFrame({"time": avg[:, 0], "K": 0.5 * avg[:, 3]}).to_csv(
            out_dir / f"{case['key']}_kinetic_energy.csv", index=False, encoding="utf-8-sig"
        )
    if mse.shape[1] >= 3:
        pd.DataFrame({"time": mse[:, 0], "A_m": mse[:, 1], "sqrt_A_m": mse[:, 2]}).to_csv(
            out_dir / f"{case['key']}_mprime_mse_aggregation.csv", index=False, encoding="utf-8-sig"
        )
    return [png, pdf]


def plot_lengths(case: dict, data: dict[str, np.ndarray], out_dir: Path) -> list[Path]:
    scales = deduplicate_time(data.get("mse_aggregation_scales", np.empty((0, 0))))
    wz = deduplicate_time(data.get("w_z075_spectral_length", np.empty((0, 0))))
    figures: list[Path] = []
    if scales.shape[1] < 11:
        return figures

    colors = {
        "L_peak": (0.00, 0.20, 0.95),
        "L_integral": (0.85, 0.10, 0.10),
        "lambda_m": (0.95, 0.45, 0.00),
        "lambda_w": (0.00, 0.55, 0.18),
    }
    fig = plt.figure(figsize=(6.5, 5.84), facecolor="white")
    ax = fig.add_axes([0.186, 0.260, 0.792, 0.705])
    style_axis(ax)
    ax.plot(scales[:, 0], scales[:, 1], lw=3.5, color=colors["L_peak"], label=r"$L_{\rm peak}$")
    ax.plot(scales[:, 0], scales[:, 2], lw=3.5, color=colors["L_integral"], label=r"$L_{\rm int}$")
    ax.plot(scales[:, 0], scales[:, 4], lw=3.5, color=colors["lambda_m"], label=r"$2\pi\ell_m$")
    if wz.shape[1] >= 3:
        ax.plot(wz[:, 0], wz[:, 2], lw=3.5, color=colors["lambda_w"], label=r"$2\pi\ell_w(z=0.75)$")
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$L/H$")
    ax.legend(frameon=False, loc="best", handlelength=2.5)
    png = out_dir / f"{case['key']}_aggregation_lengths_timeseries.png"
    pdf = out_dir / f"{case['key']}_aggregation_lengths_timeseries.pdf"
    fig.savefig(png)
    fig.savefig(pdf)
    plt.close(fig)
    figures.extend([png, pdf])

    fig = plt.figure(figsize=(6.5, 5.84), facecolor="white")
    ax = fig.add_axes([0.186, 0.260, 0.792, 0.705])
    style_axis(ax)
    ax.plot(scales[:, 0], scales[:, 6], lw=3.5, color=(0.00, 0.45, 1.00), label=r"$\overline{R}_{m'>0}$")
    ax.plot(scales[:, 0], scales[:, 7], lw=3.5, color=(0.00, 0.45, 1.00), ls="--", label=r"$R_{\max,m'>0}$")
    ax.plot(scales[:, 0], scales[:, 9], lw=3.5, color=(0.95, 0.25, 0.05), label=r"$\overline{R}_{m'>\sigma}$")
    ax.plot(scales[:, 0], scales[:, 10], lw=3.5, color=(0.95, 0.25, 0.05), ls="--", label=r"$R_{\max,m'>\sigma}$")
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$R/H$")
    ax.legend(frameon=False, loc="best", handlelength=2.5)
    png2 = out_dir / f"{case['key']}_cluster_radii_timeseries.png"
    pdf2 = out_dir / f"{case['key']}_cluster_radii_timeseries.pdf"
    fig.savefig(png2)
    fig.savefig(pdf2)
    plt.close(fig)
    figures.extend([png2, pdf2])

    pd.DataFrame(
        {
            "time": scales[:, 0],
            "L_peak": scales[:, 1],
            "L_integral": scales[:, 2],
            "ell_m": scales[:, 3],
            "lambda_m_2pi_ell_m": scales[:, 4],
            "cluster_count_mprime_gt0": scales[:, 5],
            "mean_radius_mprime_gt0": scales[:, 6],
            "max_radius_mprime_gt0": scales[:, 7],
            "cluster_count_mprime_gtsigma": scales[:, 8],
            "mean_radius_mprime_gtsigma": scales[:, 9],
            "max_radius_mprime_gtsigma": scales[:, 10],
        }
    ).to_csv(out_dir / f"{case['key']}_mse_aggregation_scales.csv", index=False, encoding="utf-8-sig")
    if wz.shape[1] >= 3:
        pd.DataFrame({"time": wz[:, 0], "ell_w_z075": wz[:, 1], "lambda_w_z075": wz[:, 2]}).to_csv(
            out_dir / f"{case['key']}_w_z075_spectral_length.csv", index=False, encoding="utf-8-sig"
        )
    return figures


def main() -> None:
    apply_style()
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = OUT_ROOT / f"selected_qbot05_and_beta_not1_{stamp}"
    raw_dir = out_dir / "remote_raw_out_files"
    out_dir.mkdir(parents=True, exist_ok=True)
    raw_dir.mkdir(parents=True, exist_ok=True)

    outputs: list[Path] = []
    metadata = []
    file_map = {
        "avgvar": "avgvar.out",
        "mse_aggregation": "mse_aggregation.out",
        "mse_aggregation_scales": "mse_aggregation_scales.out",
        "w_z075_spectral_length": "w_z075_spectral_length.out",
    }

    for case in CASES:
        data: dict[str, np.ndarray] = {}
        case_meta = {k: v for k, v in case.items() if k not in {"color"}}
        case_meta["files"] = {}
        for key, fname in file_map.items():
            remote_file = case["remote_run"] + "/data/" + fname
            local_file = raw_dir / case["key"] / fname
            ok = fetch_remote_file(remote_file, local_file)
            case_meta["files"][key] = {"remote": remote_file, "local": str(local_file), "available": ok}
            if ok:
                arr = read_table(local_file)
                data[key] = arr
                arr = deduplicate_time(arr)
                if arr.size:
                    case_meta["files"][key]["n_rows"] = int(arr.shape[0])
                    case_meta["files"][key]["last_time"] = float(np.nanmax(arr[:, 0]))
        outputs.extend(plot_energy_mprime(case, data, out_dir))
        if case.get("plot_lengths"):
            outputs.extend(plot_lengths(case, data, out_dir))
        metadata.append(case_meta)

    (out_dir / "metadata.json").write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps({"output_dir": str(out_dir), "figures": [str(p) for p in outputs]}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
