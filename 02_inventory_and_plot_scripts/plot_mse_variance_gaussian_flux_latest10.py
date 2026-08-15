from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(
    r"E:\moist RB\rotating_case_inventory\04_outputs_and_figures"
    r"\Ra8e6_mse_variance_gaussian_flux_20260802\final_latest10"
)
FIGDIR = ROOT.parent / "figures"

CASES = [
    ("Ek1p5e-4", 1.5e-4),
    ("Ek2e-4", 2.0e-4),
    ("Ek3e-3", 3.0e-3),
    ("Ek5e-3", 5.0e-3),
    ("Ek7e-3", 7.0e-3),
    ("Ek1e-2", 1.0e-2),
    ("Ek3e-2", 3.0e-2),
]

FRAME_WIDTH = 4.5
LINE_WIDTH = 3.5
AXIS_LABEL_SIZE = 24
TICK_LABEL_SIZE = 13
LEGEND_SIZE = 10
FIGSIZE = (6.5, 5.84)
AX_RECT = [0.186, 0.260, 0.792, 0.705]


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
            "figure.dpi": 180,
            "savefig.dpi": 300,
        }
    )


def axis() -> tuple[plt.Figure, plt.Axes]:
    figure = plt.figure(figsize=FIGSIZE, facecolor="white")
    ax = figure.add_axes(AX_RECT)
    for spine in ax.spines.values():
        spine.set_linewidth(FRAME_WIDTH)
    ax.tick_params(direction="in", length=12, width=1.2, top=True, right=True)
    ax.minorticks_on()
    ax.tick_params(
        which="minor", direction="in", length=6, width=1.0, top=True, right=True
    )
    ax.set_box_aspect(5.2 / 6.5)
    return figure, ax


def colors() -> dict[str, tuple[float, float, float, float]]:
    values = plt.cm.turbo(np.linspace(0.06, 0.94, len(CASES)))
    return {label: color for (label, _), color in zip(CASES, values)}


def ek_label(ek: float) -> str:
    exponent = int(math.floor(math.log10(ek)))
    coefficient = ek / 10.0**exponent
    if abs(coefficient - 1.0) < 1.0e-10:
        return rf"$10^{{{exponent}}}$"
    if abs(coefficient - round(coefficient)) < 1.0e-10:
        coefficient_text = str(int(round(coefficient)))
    else:
        coefficient_text = f"{coefficient:g}"
    return rf"${coefficient_text}\times10^{{{exponent}}}$"


def load() -> dict[str, dict[str, object]]:
    output: dict[str, dict[str, object]] = {}
    for label, ek in CASES:
        directory = ROOT / label
        data = np.load(directory / "mse_variance_spectrum_budget_gaussian_flux.npz")
        manifest = json.loads((directory / "manifest.json").read_text(encoding="utf-8"))
        output[label] = {"ek": ek, "data": data, "manifest": manifest}
    return output


def save(figure: plt.Figure, stem: str) -> list[Path]:
    FIGDIR.mkdir(parents=True, exist_ok=True)
    paths = [FIGDIR / f"{stem}.png", FIGDIR / f"{stem}.pdf"]
    for path in paths:
        figure.savefig(path, facecolor="white")
    plt.close(figure)
    return paths


def add_case_legend(ax: plt.Axes, location: str = "best", ncol: int = 1) -> None:
    ax.legend(
        title=r"$Ek$",
        frameon=False,
        loc=location,
        ncol=ncol,
        handlelength=2.8,
        columnspacing=1.0,
    )


def plot_spectrum(records, palette) -> list[Path]:
    figure, ax = axis()
    for label, ek in CASES:
        data = records[label]["data"]
        k = data["k_shell"]
        energy = data["spectrum_frames"].mean(axis=0)
        mask = (k > 0.0) & (energy > 0.0)
        ax.loglog(
            k[mask],
            energy[mask],
            color=palette[label],
            lw=LINE_WIDTH,
            label=ek_label(ek),
        )
    ax.set_xlabel(r"$k_h$")
    ax.set_ylabel(r"$E_m(k_h)$")
    add_case_legend(ax, "lower left", ncol=2)
    return save(figure, "Ra8e6_MSE_variance_spectrum_latest10")


def plot_flux(records, palette, strong_only: bool = False) -> list[Path]:
    figure, ax = axis()
    selected = CASES[:2] if strong_only else CASES
    for label, ek in selected:
        data = records[label]["data"]
        k = data["k_c"]
        flux = data["pi_gaussian_frames"].mean(axis=0)
        ax.semilogx(
            k,
            flux,
            color=palette[label],
            lw=LINE_WIDTH,
            label=ek_label(ek),
        )
    ax.axhline(0.0, color=(0.25, 0.25, 0.25), lw=1.6, ls="--", zorder=0)
    ax.set_xlabel(r"$k_c=1/\ell$")
    ax.set_ylabel(r"$\Pi_m(k_c)$")
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0), useMathText=True)
    add_case_legend(ax, "best", ncol=2 if not strong_only else 1)
    stem = (
        "Ra8e6_MSE_variance_Gaussian_flux_strong_rotation_latest10"
        if strong_only
        else "Ra8e6_MSE_variance_Gaussian_flux_all_cases_latest10"
    )
    return save(figure, stem)


def plot_flux_normalized(records, palette) -> list[Path]:
    figure, ax = axis()
    for label, ek in CASES:
        data = records[label]["data"]
        k = data["k_c"]
        flux = data["pi_gaussian_frames"].mean(axis=0)
        amplitude = np.max(np.abs(flux))
        if amplitude <= 0.0:
            continue
        ax.semilogx(
            k,
            flux / amplitude,
            color=palette[label],
            lw=LINE_WIDTH,
            label=ek_label(ek),
        )
    ax.axhline(0.0, color=(0.25, 0.25, 0.25), lw=1.6, ls="--", zorder=0)
    ax.set_xlabel(r"$k_c=1/\ell$")
    ax.set_ylabel(r"$\Pi_m/\max|\Pi_m|$")
    ax.set_ylim(-1.12, 1.12)
    add_case_legend(ax, "best", ncol=2)
    return save(figure, "Ra8e6_MSE_variance_Gaussian_flux_normalized_latest10")


def plot_case_budgets(label: str, ek: float, record, palette) -> list[Path]:
    output: list[Path] = []
    data = record["data"]

    figure, ax = axis()
    k = data["k_shell"]
    mask = k > 0.0
    terms = [
        (data["transfer_frames"].mean(axis=0), r"$T_m$", "-"),
        (data["production_frames"].mean(axis=0), r"$P_m$", "--"),
        (data["dissipation_frames"].mean(axis=0), r"$D_m$", "-."),
        (data["spectral_dedt"], r"$\partial_tE_m$", ":"),
    ]
    term_colors = [
        (0.00, 0.00, 1.00),
        (0.90, 0.05, 0.05),
        (0.00, 0.62, 0.12),
        (0.72, 0.10, 0.85),
    ]
    for (values, legend, linestyle), color in zip(terms, term_colors):
        ax.semilogx(k[mask], values[mask], color=color, lw=LINE_WIDTH, ls=linestyle, label=legend)
    ax.axhline(0.0, color=(0.25, 0.25, 0.25), lw=1.4, ls="--", zorder=0)
    ax.set_xlabel(r"$k_h$")
    ax.set_ylabel(r"$\partial_t E_m(k_h)$ terms")
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0), useMathText=True)
    ax.legend(frameon=False, loc="best", ncol=2, handlelength=2.8)
    output.extend(save(figure, f"{label}_MSE_variance_spectral_budget_latest10"))

    figure, ax = axis()
    kc = data["k_c"]
    gaussian_terms = [
        (data["production_gaussian_frames"].mean(axis=0), r"$P_m^{<}$", "--"),
        (-data["dissipation_gaussian_frames"].mean(axis=0), r"$-\chi_m^{<}$", "-."),
        (-data["pi_gaussian_frames"].mean(axis=0), r"$-\Pi_m$", "-"),
        (data["dvariance_dt"], r"$\partial_tV_m^{<}$", ":"),
    ]
    for (values, legend, linestyle), color in zip(gaussian_terms, term_colors[1:] + term_colors[:1]):
        ax.semilogx(kc, values, color=color, lw=LINE_WIDTH, ls=linestyle, label=legend)
    ax.axhline(0.0, color=(0.25, 0.25, 0.25), lw=1.4, ls="--", zorder=0)
    ax.set_xlabel(r"$k_c=1/\ell$")
    ax.set_ylabel(r"$\partial_t V_m^{<}$ terms")
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0), useMathText=True)
    ax.legend(frameon=False, loc="best", ncol=2, handlelength=2.8)
    output.extend(save(figure, f"{label}_MSE_variance_Gaussian_budget_latest10"))
    return output


def write_summary(records) -> Path:
    rows = []
    for label, ek in CASES:
        data = records[label]["data"]
        manifest = records[label]["manifest"]
        k = data["k_c"]
        flux = data["pi_gaussian_frames"].mean(axis=0)
        negative = flux < 0.0
        if np.any(negative):
            negative_indices = np.where(negative)[0]
            inverse_k_min = float(k[negative_indices[0]])
            inverse_k_max = float(k[negative_indices[-1]])
        else:
            inverse_k_min = math.nan
            inverse_k_max = math.nan
        p = data["production_gaussian_frames"].mean(axis=0)
        chi = data["dissipation_gaussian_frames"].mean(axis=0)
        residual = data["gaussian_closure"]
        scale = max(np.max(np.abs(p)), np.max(np.abs(chi)), np.max(np.abs(flux)), 1.0e-30)
        energy = data["spectrum_frames"].mean(axis=0)
        kh = data["k_shell"]
        peak_index = 1 + int(np.argmax(energy[1:]))
        rows.append(
            {
                "case": label,
                "Ek": ek,
                "AR": manifest["parameters"]["Lx"],
                "Nx_movie": manifest["movie_shape_zyx"][2],
                "t_first": float(data["times"][0]),
                "t_last": float(data["times"][-1]),
                "nframes": len(data["times"]),
                "Pi_min": float(np.min(flux)),
                "Pi_max": float(np.max(flux)),
                "inverse_k_min": inverse_k_min,
                "inverse_k_max": inverse_k_max,
                "inverse_length_max_2pi_over_kmin": (
                    2.0 * np.pi / inverse_k_min if math.isfinite(inverse_k_min) else math.nan
                ),
                "inverse_length_min_2pi_over_kmax": (
                    2.0 * np.pi / inverse_k_max if math.isfinite(inverse_k_max) else math.nan
                ),
                "spectrum_peak_k": float(kh[peak_index]),
                "spectrum_peak_length_2pi_over_k": float(2.0 * np.pi / kh[peak_index]),
                "max_relative_gaussian_closure_residual": float(np.max(np.abs(residual)) / scale),
                "mprime_reconstruction_max_abs_error": manifest[
                    "mprime_reconstruction_max_abs_error_first_frame"
                ],
            }
        )
    path = FIGDIR / "Ra8e6_MSE_variance_flux_summary_latest10.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def main() -> None:
    configure_style()
    records = load()
    palette = colors()
    generated = []
    generated.extend(plot_spectrum(records, palette))
    generated.extend(plot_flux(records, palette, strong_only=False))
    generated.extend(plot_flux(records, palette, strong_only=True))
    generated.extend(plot_flux_normalized(records, palette))
    for label, ek in CASES:
        generated.extend(plot_case_budgets(label, ek, records[label], palette))
    summary = write_summary(records)
    print(json.dumps({"figures": [str(path) for path in generated], "summary": str(summary)}, indent=2))


if __name__ == "__main__":
    main()
