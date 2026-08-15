from __future__ import annotations

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
import pandas as pd


REMOTE_HOST = "c01n0006"
PROJECT_ROOT = Path(r"E:\moist RB\rotating_case_inventory")
OUT_ROOT = PROJECT_ROOT / "04_outputs_and_figures" / "case_timeseries"

REMOTE_ROOTS = {
    "rotating_case": "/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case",
    "ar10_old": "/share/org/SHUTUANL/shu_zhangjs/rainy model/ns/transition_study/beta1/lowrestest/aspect_ratio_study",
}


REMOTE_CODE = r"""
import json, math, re
from pathlib import Path

ROOTS = __REMOTE_ROOTS__

def parse_float_token(token):
    if token is None:
        return None
    token = str(token)
    token = token.replace("p", ".")
    try:
        return float(token)
    except Exception:
        return None

def value_after_prefix(token, prefix):
    if token is None or not str(token).startswith(prefix):
        return None
    return parse_float_token(str(token)[len(prefix):])

def parse_sci_label(token, prefix):
    return value_after_prefix(token, prefix)

def parse_grid(token):
    if not token:
        return None, None, None
    s = str(token)
    m = re.match(r"N?(\d+)[xX×脳](\d+)[xX×脳](\d+)$", s)
    if not m:
        m = re.match(r"(\d{3})(\d{3})(\d{2,3})$", s)
    if not m:
        return None, None, None
    return tuple(int(x) for x in m.groups())

def read_timeseries(path, cols):
    if not path.exists() or path.stat().st_size == 0:
        return []
    rows = []
    with path.open("r", errors="ignore") as f:
        for line in f:
            parts = line.split()
            if len(parts) <= max(cols):
                continue
            try:
                vals = [float(parts[i].replace("D", "e").replace("d", "e")) for i in cols]
            except Exception:
                continue
            if all(math.isfinite(v) for v in vals):
                rows.append(vals)
    by_t = {}
    for vals in rows:
        by_t[round(vals[0], 10)] = vals
    return [by_t[k] for k in sorted(by_t)]

def read_case_info_json(run):
    for p in [run.parent / "case_info.json", run / "case_info.json"]:
        if p.exists():
            try:
                return json.loads(p.read_text(errors="ignore"))
            except Exception:
                return {}
    return {}

def parse_rotating_case(run, root):
    rel = run.relative_to(root).parts
    info = {
        "source_family": "rotating_case",
        "source_root": str(root),
        "run_path": str(run),
        "case_path": str(run.parent),
        "variant": "",
    }
    if len(rel) < 8:
        return None
    info["Pr_label"] = rel[0]
    info["Ra_label"] = rel[1]
    info["Ek_label"] = rel[2]
    info["AR_label"] = rel[3]
    info["Beta_label"] = rel[4]
    info["q_label"] = rel[5]
    info["grid_label"] = rel[6]
    if len(rel) > 8:
        info["variant"] = "/".join(rel[7:-1])
    info["Pr"] = parse_sci_label(rel[0], "Pr")
    info["Ra"] = parse_sci_label(rel[1], "Ra")
    info["Ek"] = None if rel[2] == "norotating" else parse_sci_label(rel[2], "Ek")
    info["rotation"] = "norotating" if info["Ek"] is None else "rotating"
    info["AR"] = parse_sci_label(rel[3], "AR")
    info["beta"] = parse_sci_label(rel[4], "Beta")
    nx, ny, nz = parse_grid(rel[6])
    info["Nx"], info["Ny"], info["Nz"] = nx, ny, nz
    qm = re.match(r"qbot(.+)_qtop(.+)", rel[5])
    if qm:
        info["qbot"] = parse_float_token(qm.group(1))
        info["qtop"] = parse_float_token(qm.group(2))
    return info

def parse_beta_not1(run, root):
    rel = run.relative_to(root).parts
    info = {
        "source_family": "beta_not1",
        "source_root": str(root),
        "run_path": str(run),
        "case_path": str(run.parent),
        "Pr_label": "Pr0p7",
        "Pr": 0.7,
        "Ek": None,
        "Ek_label": "norotating",
        "rotation": "norotating",
        "AR": None,
        "AR_label": "",
        "beta": None,
        "Beta_label": "",
        "qbot": 1.0,
        "qtop": 0.004978,
        "q_label": "qbot1_qtop0p004978",
        "grid_label": "",
        "variant": "",
    }
    if not rel:
        return None
    # Common shape:
    # Ra1e8/AR10/norotating/Beta1p20/run
    # Ra1e8/AR10/norotating/Beta1p20/384384192/run
    # Ra1e8/AR10/norotating/Beta1p20/bot05/run
    info["Ra_label"] = rel[0] if rel[0].startswith("Ra") else "Ra?"
    info["Ra"] = parse_sci_label(info["Ra_label"], "Ra")
    if len(rel) >= 2 and rel[1].startswith("AR"):
        info["AR_label"] = rel[1]
        info["AR"] = parse_sci_label(rel[1], "AR")
    if len(rel) >= 3:
        if rel[2] == "norotating":
            info["Ek_label"] = "norotating"
            info["Ek"] = None
            info["rotation"] = "norotating"
        elif rel[2].startswith("Ek"):
            info["Ek_label"] = rel[2]
            info["Ek"] = parse_sci_label(rel[2], "Ek")
            info["rotation"] = "rotating"
    if len(rel) >= 4 and rel[3].startswith("Beta"):
        info["Beta_label"] = rel[3]
        info["beta"] = parse_sci_label(rel[3], "Beta")
    variant_parts = []
    if len(rel) > 5:
        variant_parts = list(rel[4:-1])
    elif len(rel) > 4:
        variant_parts = list(rel[4:-1])
    info["variant"] = "/".join(variant_parts)
    if "bot05" in variant_parts:
        info["qbot"] = 0.5
        info["q_label"] = "qbot0p5_qtop0p004978"
    if "botsat" in variant_parts:
        info["qbot"] = 1.0
        info["q_label"] = "qbot1_qtop0p004978"
        info["variant"] = "botsat"
    for v in variant_parts:
        nx, ny, nz = parse_grid(v)
        if nx:
            info["Nx"], info["Ny"], info["Nz"] = nx, ny, nz
            info["grid_label"] = "N{}x{}x{}".format(nx, ny, nz)
    if "Nx" not in info:
        info["Nx"], info["Ny"], info["Nz"] = None, None, None
    if len(rel) == 2 and rel[-1] == "run":
        # Fallback for a short test path like Ra8e6/run.
        info["variant"] = "root_run"
    return info

def parse_ar10_old(run, root):
    rel = run.relative_to(root).parts
    info = {
        "source_family": "ar10_old",
        "source_root": str(root),
        "run_path": str(run),
        "case_path": str(run.parent),
        "Pr_label": "Pr0p7",
        "Pr": 0.7,
        "beta": 1.0,
        "Beta_label": "Beta1",
        "qbot": 1.0,
        "qtop": 0.004978,
        "q_label": "qbot1_qtop0p004978",
        "variant": "",
    }
    if len(rel) < 4:
        return None
    info["Ra_label"] = rel[0]
    info["Ra"] = parse_sci_label(rel[0], "Ra")
    info["AR_label"] = rel[1]
    info["AR"] = parse_sci_label(rel[1], "AR")
    folder = rel[2]
    if folder.startswith("Ra"):
        m = re.search(r"Ek(.+?)AR", folder)
        info["Ek_label"] = "Ek" + m.group(1) if m else ""
        info["Ek"] = parse_sci_label(info["Ek_label"], "Ek") if m else None
        info["rotation"] = "rotating"
    elif folder == "norotating":
        info["Ek_label"] = "norotating"
        info["Ek"] = None
        info["rotation"] = "norotating"
        if len(rel) >= 4:
            info["variant"] = "" if rel[3] == "run" else rel[3]
            nx, ny, nz = parse_grid(rel[3])
            if nx:
                info["Nx"], info["Ny"], info["Nz"] = nx, ny, nz
                info["grid_label"] = "N{}x{}x{}".format(nx, ny, nz)
    else:
        return None
    if "Nx" not in info:
        info["Nx"], info["Ny"], info["Nz"] = None, None, None
        info["grid_label"] = ""
    return info

def parse_run(run, family, root):
    if family == "rotating_case":
        return parse_rotating_case(run, root)
    if family == "beta_not1":
        return parse_beta_not1(run, root)
    if family == "ar10_old":
        return parse_ar10_old(run, root)
    return None

raw_records = []
for family, root_s in ROOTS.items():
    root = Path(root_s)
    if not root.exists():
        continue
    # Use mse_aggregation.out as the inclusion criterion so each plotted case has m'.
    for mse_file in sorted(root.rglob("mse_aggregation.out")):
        if mse_file.stat().st_size == 0:
            continue
        run = mse_file.parent.parent
        data = run / "data"
        info = parse_run(run, family, root)
        if info is None:
            continue
        avg = read_timeseries(data / "avgvar.out", [0, 3])
        avg_clean = []
        for t, col4 in avg:
            K = 0.5 * col4
            if 0.0 <= K < 0.1:
                avg_clean.append([t, K])
        mse = read_timeseries(data / "mse_aggregation.out", [0, 1, 2])
        if not mse and not avg_clean:
            continue
        info["kinetic"] = avg_clean
        info["mprime"] = mse
        raw_records.append(info)

merged = {}
for rec in raw_records:
    # Merge continuations under the structured rotating_case root only; keep
    # explicitly different legacy variants (384 grid, bot05, botsat) separate.
    key = (
        rec.get("source_family"), rec.get("Pr_label"), rec.get("Ra_label"),
        rec.get("Ek_label"), rec.get("AR_label"), rec.get("Beta_label"),
        rec.get("q_label"), rec.get("grid_label"), rec.get("variant")
    )
    if rec.get("source_family") == "rotating_case":
        key = (
            rec.get("source_family"), rec.get("Pr_label"), rec.get("Ra_label"),
            rec.get("Ek_label"), rec.get("AR_label"), rec.get("Beta_label"),
            rec.get("q_label"), rec.get("grid_label")
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
    code = REMOTE_CODE.replace("__REMOTE_ROOTS__", repr(REMOTE_ROOTS))
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


def fmt_num(x) -> str:
    if x is None:
        return "NR"
    x = float(x)
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


def fmt_plain(x, nd=3) -> str:
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return ""
    return f"{float(x):.{nd}g}".rstrip("0").rstrip(".")


def source_tag(source_family: str) -> str:
    return {
        "rotating_case": "new",
        "beta_not1": "β-sweep",
        "ar10_old": "AR10-old",
    }.get(source_family, source_family)


def case_label(rec: dict, group: list[dict]) -> str:
    parts = []
    if rec.get("Ek") is None:
        parts.append("NR")
    else:
        parts.append(rf"$Ek={fmt_num(rec['Ek'])}$")
    if rec.get("AR") is not None:
        parts.append(rf"$\Gamma={fmt_plain(rec.get('AR'), 3)}$")
    if rec.get("beta") is not None:
        parts.append(rf"$\beta={fmt_plain(rec.get('beta'), 3)}$")
    qbot = rec.get("qbot")
    if qbot is not None and abs(float(qbot) - 1.0) > 1e-12:
        parts.append(rf"$q_{{bot}}={fmt_plain(qbot, 3)}$")
    grid_values = {(r.get("Nx"), r.get("Ny"), r.get("Nz")) for r in group}
    if len(grid_values) > 1 and rec.get("Nx"):
        parts.append(rf"${rec.get('Nx')}\!\times\!{rec.get('Ny')}\!\times\!{rec.get('Nz')}$")
    variant = rec.get("variant") or ""
    if variant and variant not in {"run"}:
        parts.append(variant.replace("_", "-"))
    parts.append(source_tag(rec.get("source_family", "")))
    return ", ".join(parts)


def case_sort_key(rec: dict):
    family_order = {"rotating_case": 0, "beta_not1": 1, "ar10_old": 2}
    ek = rec.get("Ek")
    ek_sort = -1.0 if ek is None else float(ek)
    beta = -1 if rec.get("beta") is None else float(rec.get("beta"))
    ar = -1 if rec.get("AR") is None else float(rec.get("AR"))
    qbot = -1 if rec.get("qbot") is None else float(rec.get("qbot"))
    grid = rec.get("Nx") or 0
    return (family_order.get(rec.get("source_family"), 99), ek_sort, beta, ar, qbot, grid, rec.get("variant") or "")


def safe_name(label: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", label).strip("_")


def plot_group(ra_label: str, group: list[dict], out_dir: Path) -> list[Path]:
    apply_style()
    group = sorted(group, key=case_sort_key)
    colors = plt.cm.turbo(np.linspace(0.05, 0.95, max(len(group), 2)))[: len(group)]

    fig = plt.figure(figsize=(16.0, 5.84), facecolor="white")
    ax1 = fig.add_axes([0.065, 0.260, 0.315, 0.705])
    ax2 = fig.add_axes([0.455, 0.260, 0.315, 0.705])
    for ax in (ax1, ax2):
        style_axis(ax)

    handles = []
    labels = []
    for rec, color in zip(group, colors):
        label = case_label(rec, group)
        h = None
        if rec.get("kinetic"):
            arr = np.asarray(rec["kinetic"], float)
            h, = ax1.plot(arr[:, 0], arr[:, 1], lw=3.5, color=color, label=label)
        if rec.get("mprime"):
            arr = np.asarray(rec["mprime"], float)
            if h is None:
                h, = ax2.plot(arr[:, 0], arr[:, 2], lw=3.5, color=color, label=label)
            else:
                ax2.plot(arr[:, 0], arr[:, 2], lw=3.5, color=color, label=label)
        if h is not None:
            handles.append(h)
            labels.append(label)

    ax1.set_xlabel(r"$t$")
    ax1.set_ylabel(r"$K$")
    ax2.set_xlabel(r"$t$")
    ax2.set_ylabel(r"$\left\langle\left\langle m'^2\right\rangle_{xy}\right\rangle_z^{1/2}$")
    if handles:
        fig.legend(
            handles,
            labels,
            frameon=False,
            loc="center left",
            bbox_to_anchor=(0.800, 0.52),
            handlelength=2.2,
            labelspacing=0.50,
            borderaxespad=0.0,
        )

    out_paths: list[Path] = []
    base = f"{safe_name(ra_label)}_all_current_kinetic_mprime_timeseries"
    png = out_dir / f"{base}.png"
    pdf = out_dir / f"{base}.pdf"
    fig.savefig(png)
    fig.savefig(pdf)
    plt.close(fig)
    out_paths.extend([png, pdf])

    # Variance-only companion figure, using the same source cases.
    fig = plt.figure(figsize=(9.2, 5.84), facecolor="white")
    ax = fig.add_axes([0.150, 0.260, 0.560, 0.705])
    style_axis(ax)
    handles = []
    labels = []
    for rec, color in zip(group, colors):
        if not rec.get("mprime"):
            continue
        arr = np.asarray(rec["mprime"], float)
        h, = ax.plot(arr[:, 0], arr[:, 1], lw=3.5, color=color)
        handles.append(h)
        labels.append(case_label(rec, group))
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$\left\langle\left\langle m'^2\right\rangle_{xy}\right\rangle_z$")
    if handles:
        fig.legend(
            handles,
            labels,
            frameon=False,
            loc="center left",
            bbox_to_anchor=(0.740, 0.52),
            handlelength=2.2,
            labelspacing=0.50,
            borderaxespad=0.0,
        )
    png2 = out_dir / f"{safe_name(ra_label)}_all_current_mse_variance_timeseries.png"
    pdf2 = out_dir / f"{safe_name(ra_label)}_all_current_mse_variance_timeseries.pdf"
    fig.savefig(png2)
    fig.savefig(pdf2)
    plt.close(fig)
    out_paths.extend([png2, pdf2])
    return out_paths


def write_csvs(records: list[dict], out_dir: Path) -> None:
    meta_rows = []
    energy_rows = []
    m_rows = []
    for rec in records:
        base = {
            k: rec.get(k)
            for k in [
                "source_family",
                "Ra_label",
                "Ra",
                "Pr",
                "Ek_label",
                "Ek",
                "rotation",
                "AR",
                "beta",
                "qbot",
                "qtop",
                "q_label",
                "variant",
                "Nx",
                "Ny",
                "Nz",
                "run_path",
                "case_path",
                "max_T_kinetic",
                "max_T_mprime",
                "n_kinetic",
                "n_mprime",
                "source_run_paths",
            ]
        }
        meta_rows.append(base)
        for t, k_val in rec.get("kinetic", []):
            row = dict(base)
            row.update({"time": t, "K": k_val})
            energy_rows.append(row)
        for t, a_m, sqrt_a_m in rec.get("mprime", []):
            row = dict(base)
            row.update({"time": t, "A_m": a_m, "sqrt_A_m": sqrt_a_m})
            m_rows.append(row)

    pd.DataFrame(meta_rows).to_csv(out_dir / "all_current_timeseries_case_metadata.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(energy_rows).to_csv(out_dir / "all_current_kinetic_energy_timeseries.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(m_rows).to_csv(out_dir / "all_current_mprime_mse_aggregation_timeseries.csv", index=False, encoding="utf-8-sig")


def main() -> None:
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = OUT_ROOT / f"all_current_energy_mprime_no_beta_not1_{stamp}"
    out_dir.mkdir(parents=True, exist_ok=True)
    records = run_remote()
    records = sorted(records, key=lambda r: (float(r.get("Ra") or -1), case_sort_key(r)))
    (out_dir / "all_current_raw_remote_timeseries.json").write_text(
        json.dumps(records, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    write_csvs(records, out_dir)
    groups = defaultdict(list)
    for rec in records:
        groups[rec.get("Ra_label", "Ra_unknown")].append(rec)
    figures = []
    for ra_label, group in sorted(groups.items(), key=lambda kv: float(kv[1][0].get("Ra") or -1)):
        figures.extend(plot_group(ra_label, group, out_dir))
    print(
        json.dumps(
            {
                "output_dir": str(out_dir),
                "n_cases": len(records),
                "groups": {k: len(v) for k, v in groups.items()},
                "figures": [str(p) for p in figures],
                "csv": [
                    str(out_dir / "all_current_timeseries_case_metadata.csv"),
                    str(out_dir / "all_current_kinetic_energy_timeseries.csv"),
                    str(out_dir / "all_current_mprime_mse_aggregation_timeseries.csv"),
                ],
            },
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
