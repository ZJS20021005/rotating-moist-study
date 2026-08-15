from __future__ import annotations

import json
import math
import subprocess
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


REMOTE_HOST = "c01n0006"
REMOTE_ROOT = "/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case"
PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUT_ROOT = PROJECT_ROOT / "04_outputs_and_figures" / "case_timeseries"


REMOTE_CODE = r"""
import json, math, re
from pathlib import Path

ROOT = Path(r'''__REMOTE_ROOT__''')

def label_float(s, prefix):
    if not s.startswith(prefix): return None
    return float(s[len(prefix):].replace('p','.'))

def sci_float(s, prefix):
    if not s.startswith(prefix): return None
    return float(s[len(prefix):].replace('p','.'))

def parse_run(run):
    rel = run.relative_to(ROOT).parts
    info = {
        "Pr_label": rel[0],
        "Ra_label": rel[1],
        "Ek_label": rel[2],
        "AR_label": rel[3],
        "Beta_label": rel[4],
        "q_label": rel[5],
        "grid_label": rel[6],
        "run_path": str(run),
        "case_path": str(run.parent),
    }
    info["Pr"] = label_float(rel[0], "Pr")
    info["Ra"] = sci_float(rel[1], "Ra")
    info["Ek"] = None if rel[2] == "norotating" else sci_float(rel[2], "Ek")
    info["rotation"] = "norotating" if info["Ek"] is None else "rotating"
    info["AR"] = label_float(rel[3], "AR")
    info["beta"] = label_float(rel[4], "Beta")
    gm = re.match(r"N(\d+)[x×](\d+)[x×](\d+)", rel[6])
    if not gm:
        gm = re.match(r"N(\d+)[xX×](\d+)[xX×](\d+)", rel[6])
    if gm:
        info["Nx"], info["Ny"], info["Nz"] = map(int, gm.groups())
    return info

def read_timeseries(path, cols):
    if not path.exists():
        return []
    out = []
    with path.open("r", errors="ignore") as f:
        for line in f:
            parts = line.split()
            if len(parts) <= max(cols):
                continue
            try:
                vals = [float(parts[i].replace("D","e").replace("d","e")) for i in cols]
            except Exception:
                continue
            if all(math.isfinite(v) for v in vals):
                out.append(vals)
    # remove duplicated physical time, keep latest row
    by_t = {}
    for vals in out:
        by_t[round(vals[0], 10)] = vals
    return [by_t[k] for k in sorted(by_t)]

raw_records = []
for run in sorted(ROOT.rglob("run")):
    if not (run / "bou.in").exists():
        continue
    info = parse_run(run)
    data = run / "data"
    avg = read_timeseries(data / "avgvar.out", [0,3])
    # avgvar one-based column 4 is <u^2+v^2+w^2>_V; K=0.5*column4.
    avg_clean = []
    for t, col4 in avg:
        K = 0.5 * col4
        if 0.0 <= K < 0.1:
            avg_clean.append([t, K])
    mse = read_timeseries(data / "mse_aggregation.out", [0,1,2])
    if not avg_clean and not mse:
        continue
    rec = dict(info)
    rec["kinetic"] = avg_clean
    rec["mprime"] = mse
    raw_records.append(rec)

merged = {}
for rec in raw_records:
    key = (
        rec.get("Pr_label"), rec.get("Ra_label"), rec.get("Ek_label"),
        rec.get("AR_label"), rec.get("Beta_label"), rec.get("q_label"),
        rec.get("grid_label")
    )
    if key not in merged:
        out = dict(rec)
        out["kinetic"] = []
        out["mprime"] = []
        out["source_run_paths"] = []
        merged[key] = out
    merged[key]["kinetic"].extend(rec.get("kinetic", []))
    merged[key]["mprime"].extend(rec.get("mprime", []))
    merged[key]["source_run_paths"].append(rec.get("run_path"))

records = []
for rec in merged.values():
    by_t = {}
    for vals in rec.get("kinetic", []):
        by_t[round(vals[0], 10)] = vals
    rec["kinetic"] = [by_t[k] for k in sorted(by_t)]
    by_t = {}
    for vals in rec.get("mprime", []):
        by_t[round(vals[0], 10)] = vals
    rec["mprime"] = [by_t[k] for k in sorted(by_t)]
    rec["max_T_kinetic"] = max([x[0] for x in rec["kinetic"]], default=None)
    rec["max_T_mprime"] = max([x[0] for x in rec["mprime"]], default=None)
    rec["n_kinetic"] = len(rec["kinetic"])
    rec["n_mprime"] = len(rec["mprime"])
    records.append(rec)

print(json.dumps(records, ensure_ascii=False))
"""


def run_remote() -> list[dict]:
    code = REMOTE_CODE.replace("__REMOTE_ROOT__", REMOTE_ROOT)
    proc = subprocess.run(
        ["ssh", REMOTE_HOST, "python3", "-"],
        input=code,
        text=True,
        capture_output=True,
        encoding="utf-8",
        errors="replace",
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"remote extraction failed\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}")
    return json.loads(proc.stdout)


def apply_style():
    mpl.rcParams.update({
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
    })


def style_axis(ax):
    for spine in ax.spines.values():
        spine.set_linewidth(4.5)
    ax.tick_params(direction="in", length=12, width=1.2, top=True, right=True)
    ax.minorticks_on()
    ax.tick_params(which="minor", direction="in", length=6, width=1.0, top=True, right=True)
    ax.set_box_aspect(5.2 / 6.5)


def fmt_num(x):
    if x is None:
        return "NR"
    if x == 0:
        return "0"
    exp = int(np.floor(np.log10(abs(x))))
    mant = x / 10**exp
    if abs(mant - round(mant)) < 1e-8:
        mant_s = str(int(round(mant)))
    else:
        mant_s = f"{mant:.2g}".rstrip("0").rstrip(".")
    if mant_s == "1":
        return rf"10^{{{exp}}}"
    return rf"{mant_s}\times10^{{{exp}}}"


def fmt_plain(x, nd=3):
    if x is None:
        return ""
    return f"{float(x):.{nd}g}".rstrip("0").rstrip(".")


def parse_q_label(q_label):
    if not q_label or not isinstance(q_label, str):
        return None, None
    import re
    m = re.match(r"qbot(.+)_qtop(.+)", q_label)
    if not m:
        return None, None
    return float(m.group(1).replace("p", ".")), float(m.group(2).replace("p", "."))


def case_label(rec, group):
    parts = []
    if rec.get("Ek") is None:
        parts.append("NR")
    else:
        parts.append(rf"$Ek={fmt_num(rec['Ek'])}$")
    beta_values = {round(float(r.get("beta", -999)), 10) for r in group if r.get("beta") is not None}
    ar_values = {round(float(r.get("AR", -999)), 10) for r in group if r.get("AR") is not None}
    grid_values = {(r.get("Nx"), r.get("Ny"), r.get("Nz")) for r in group}
    q_values = {r.get("q_label") for r in group if r.get("q_label")}
    if len(beta_values) > 1 or abs(float(rec.get("beta", 0)) - 1.02) > 1e-10:
        parts.append(rf"$\beta={fmt_plain(rec.get('beta'), 3)}$")
    if len(ar_values) > 1:
        parts.append(rf"$\Gamma={rec.get('AR'):.0f}$")
    if len(grid_values) > 1:
        parts.append(rf"${rec.get('Nx')}\!\times\!{rec.get('Ny')}\!\times\!{rec.get('Nz')}$")
    qbot, qtop = parse_q_label(rec.get("q_label"))
    if qbot is not None and (len(q_values) > 1 or abs(qbot - 1.0) > 1e-12):
        parts.append(rf"$q_{{bot}}={fmt_plain(qbot, 3)}$")
    return ", ".join(parts)


def case_sort_key(rec):
    ek = rec.get("Ek")
    ek_sort = -1.0 if ek is None else float(ek)
    return (ek_sort, float(rec.get("beta") or 0), float(rec.get("AR") or 0), rec.get("Nx") or 0, rec.get("run_path", ""))


def plot_group(ra_label: str, group: list[dict], out_dir: Path):
    apply_style()
    group = sorted(group, key=case_sort_key)
    colors = plt.cm.turbo(np.linspace(0.06, 0.94, max(len(group), 2)))[:len(group)]

    # two-panel: kinetic energy and m' standard deviation
    fig = plt.figure(figsize=(13.0, 5.84), facecolor="white")
    ax1 = fig.add_axes([0.085, 0.260, 0.365, 0.705])
    ax2 = fig.add_axes([0.575, 0.260, 0.365, 0.705])
    for ax in (ax1, ax2):
        style_axis(ax)
    plotted_energy = 0
    plotted_mprime = 0
    for rec, color in zip(group, colors):
        label = case_label(rec, group)
        if rec.get("kinetic"):
            arr = np.asarray(rec["kinetic"], float)
            ax1.plot(arr[:, 0], arr[:, 1], lw=3.5, color=color, label=label)
            plotted_energy += 1
        if rec.get("mprime"):
            arr = np.asarray(rec["mprime"], float)
            ax2.plot(arr[:, 0], arr[:, 2], lw=3.5, color=color, label=label)
            plotted_mprime += 1
    ax1.set_xlabel(r"$t$")
    ax1.set_ylabel(r"$K$")
    ax2.set_xlabel(r"$t$")
    ax2.set_ylabel(r"$\left\langle\left\langle m'^2\right\rangle_{xy}\right\rangle_z^{1/2}$")
    if plotted_energy:
        ax1.legend(frameon=False, loc="best", handlelength=2.4)
    if plotted_mprime:
        ax2.legend(frameon=False, loc="best", handlelength=2.4)
    safe_ra = ra_label.replace(".", "p")
    png = out_dir / f"{safe_ra}_kinetic_mprime_std_timeseries.png"
    pdf = out_dir / f"{safe_ra}_kinetic_mprime_std_timeseries.pdf"
    fig.savefig(png)
    fig.savefig(pdf)
    plt.close(fig)

    # single-panel MSE horizontal variance for aggregation amplitude
    fig = plt.figure(figsize=(6.5, 5.84), facecolor="white")
    ax = fig.add_axes([0.186, 0.260, 0.792, 0.705])
    style_axis(ax)
    plotted = 0
    for rec, color in zip(group, colors):
        if rec.get("mprime"):
            arr = np.asarray(rec["mprime"], float)
            ax.plot(arr[:, 0], arr[:, 1], lw=3.5, color=color, label=case_label(rec, group))
            plotted += 1
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$\left\langle\left\langle m'^2\right\rangle_{xy}\right\rangle_z$")
    if plotted:
        ax.legend(frameon=False, loc="best", handlelength=2.4)
    png2 = out_dir / f"{safe_ra}_mse_variance_Am_timeseries.png"
    pdf2 = out_dir / f"{safe_ra}_mse_variance_Am_timeseries.pdf"
    fig.savefig(png2)
    fig.savefig(pdf2)
    plt.close(fig)
    return [png, pdf, png2, pdf2]


def write_csvs(records: list[dict], out_dir: Path):
    meta_rows = []
    energy_rows = []
    m_rows = []
    for rec in records:
        base = {k: rec.get(k) for k in [
            "Ra_label", "Ra", "Pr", "Ek_label", "Ek", "rotation", "AR", "beta",
            "q_label", "Nx", "Ny", "Nz", "run_path", "case_path", "max_T_kinetic",
            "max_T_mprime", "n_kinetic", "n_mprime", "source_run_paths"
        ]}
        meta_rows.append(base)
        for t, K in rec.get("kinetic", []):
            r = dict(base)
            r.update({"time": t, "K": K})
            energy_rows.append(r)
        for t, A_m, sqrt_A_m in rec.get("mprime", []):
            r = dict(base)
            r.update({"time": t, "A_m": A_m, "sqrt_A_m": sqrt_A_m})
            m_rows.append(r)
    pd.DataFrame(meta_rows).to_csv(out_dir / "timeseries_case_metadata.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(energy_rows).to_csv(out_dir / "kinetic_energy_timeseries.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(m_rows).to_csv(out_dir / "mprime_mse_aggregation_timeseries.csv", index=False, encoding="utf-8-sig")


def main():
    stamp = datetime.now().strftime("%Y%m%d")
    out_dir = OUT_ROOT / f"already_run_by_Ra_{stamp}"
    out_dir.mkdir(parents=True, exist_ok=True)
    records = run_remote()
    (out_dir / "raw_remote_timeseries.json").write_text(json.dumps(records, ensure_ascii=False, indent=2), encoding="utf-8")
    write_csvs(records, out_dir)
    groups = defaultdict(list)
    for rec in records:
        groups[rec["Ra_label"]].append(rec)
    figs = []
    for ra_label, group in sorted(groups.items(), key=lambda kv: float(kv[1][0].get("Ra") or 0)):
        figs.extend(plot_group(ra_label, group, out_dir))
    print(json.dumps({
        "output_dir": str(out_dir),
        "n_cases_with_any_timeseries": len(records),
        "groups": {k: len(v) for k, v in groups.items()},
        "figures": [str(p) for p in figs],
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
