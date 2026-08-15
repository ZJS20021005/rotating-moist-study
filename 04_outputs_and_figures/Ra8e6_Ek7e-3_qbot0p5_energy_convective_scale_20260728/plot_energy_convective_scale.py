from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import AutoMinorLocator, ScalarFormatter


HERE = Path(__file__).resolve().parent
AVGVAR_FILE = HERE / "avgvar.out"
CONVECTIVE_SCALE_FILE = HERE / "convective_scale.dat"

FIGSIZE = (6.5, 5.84)
AX_RECT = [0.186, 0.260, 0.792, 0.705]
BOX_ASPECT = 5.2 / 6.5
FRAME_WIDTH = 4.5
LINE_WIDTH = 3.5

BLUE = (0.00, 0.32, 0.90)
LIGHT_BLUE = (0.18, 0.70, 0.95)
FIT_RED = (0.88, 0.08, 0.12)
SHADE = (0.72, 0.72, 0.72)


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


def clean_rows(array: np.ndarray) -> np.ndarray:
    array = np.atleast_2d(array)
    array = array[np.all(np.isfinite(array), axis=1)]
    order = np.argsort(array[:, 0], kind="stable")
    array = array[order]
    # Restart files can contain duplicate physical times. Keep the latest row.
    latest_by_time = {float(row[0]): row for row in array}
    return np.asarray([latest_by_time[t] for t in sorted(latest_by_time)])


def linear_fit(x: np.ndarray, y: np.ndarray) -> dict[str, float]:
    slope, intercept = np.polyfit(x, y, 1)
    fitted = slope * x + intercept
    residual = y - fitted
    ss_res = float(np.sum(residual**2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0.0 else 1.0
    return {
        "slope": float(slope),
        "intercept": float(intercept),
        "r_squared": float(r_squared),
        "max_abs_log_residual": float(np.max(np.abs(residual))),
    }


def detect_exponential_interval(time: np.ndarray, kinetic: np.ndarray) -> dict[str, float]:
    """
    Select the longest contiguous interval that is nearly linear in log(K).

    The first five nondimensional time units are excluded as initialization
    transients. A useful interval must last at least ten time units, have
    R^2 >= 0.9999, and deviate from its exponential fit by less than 0.06 in
    log(K), approximately six percent in K.
    """
    finite = np.isfinite(time) & np.isfinite(kinetic) & (kinetic > 0.0)
    finite_indices = np.flatnonzero(finite)
    if finite_indices.size < 3:
        raise RuntimeError("Not enough positive kinetic-energy samples.")

    # Initial linear instability must precede the first dominant nonlinear
    # peak. If the record is still monotonic, retain the full available span.
    finite_kinetic = kinetic[finite_indices]
    peak_position = int(np.argmax(finite_kinetic))
    if peak_position < finite_indices.size - 1:
        peak_time = float(time[finite_indices[peak_position]])
    else:
        peak_time = float(np.max(time[finite]))

    valid = finite & (time >= 5.0) & (time <= peak_time)
    t = time[valid]
    log_k = np.log(kinetic[valid])
    if t.size < 3:
        raise RuntimeError("Not enough samples before the first energy peak.")

    # Search candidate boundaries every approximately 0.5 time units. The
    # selected interval is then refined at the original output cadence. This
    # retains the original strict criteria without an O(N^3) exhaustive scan.
    dt = float(np.median(np.diff(t)))
    coarse_step = max(1, int(round(0.5 / dt)))
    coarse = np.arange(0, t.size, coarse_step, dtype=int)
    if coarse[-1] != t.size - 1:
        coarse = np.append(coarse, t.size - 1)

    sx = np.concatenate(([0.0], np.cumsum(t)))
    sy = np.concatenate(([0.0], np.cumsum(log_k)))
    sxx = np.concatenate(([0.0], np.cumsum(t * t)))
    syy = np.concatenate(([0.0], np.cumsum(log_k * log_k)))
    sxy = np.concatenate(([0.0], np.cumsum(t * log_k)))

    def prefix_fit(i: int, j: int) -> tuple[float, float, float]:
        n = float(j - i + 1)
        sum_x = sx[j + 1] - sx[i]
        sum_y = sy[j + 1] - sy[i]
        centered_xx = sxx[j + 1] - sxx[i] - sum_x * sum_x / n
        centered_yy = syy[j + 1] - syy[i] - sum_y * sum_y / n
        centered_xy = sxy[j + 1] - sxy[i] - sum_x * sum_y / n
        if centered_xx <= 0.0 or centered_yy <= 0.0:
            return 0.0, 0.0, -np.inf
        slope = centered_xy / centered_xx
        intercept = (sum_y - slope * sum_x) / n
        r_squared = slope * centered_xy / centered_yy
        return float(slope), float(intercept), float(r_squared)

    coarse_best: tuple[int, int] | None = None
    for gap in range(coarse.size - 1, 0, -1):
        accepted: list[tuple[float, int, int]] = []
        for start_pos in range(coarse.size - gap):
            i = int(coarse[start_pos])
            j = int(coarse[start_pos + gap])
            if t[j] - t[i] < 10.0:
                continue
            slope, intercept, r_squared = prefix_fit(i, j)
            if slope <= 0.0 or r_squared < 0.9999:
                continue
            residual = log_k[i : j + 1] - (slope * t[i : j + 1] + intercept)
            max_residual = float(np.max(np.abs(residual)))
            if max_residual <= 0.06:
                accepted.append((r_squared, i, j))
        if accepted:
            _, i_best, j_best = max(accepted, key=lambda item: item[0])
            coarse_best = (i_best, j_best)
            break

    if coarse_best is None:
        raise RuntimeError("No reliable exponential-growth interval was found.")

    i_coarse, j_coarse = coarse_best
    refined: list[dict[str, float]] = []
    for i in range(max(0, i_coarse - coarse_step), min(t.size, i_coarse + coarse_step + 1)):
        for j in range(
            max(i + 2, j_coarse - coarse_step),
            min(t.size, j_coarse + coarse_step + 1),
        ):
            duration = float(t[j] - t[i])
            if duration < 10.0:
                continue
            fit = linear_fit(t[i : j + 1], log_k[i : j + 1])
            if (
                fit["slope"] > 0.0
                and fit["r_squared"] >= 0.9999
                and fit["max_abs_log_residual"] <= 0.06
            ):
                refined.append(
                    {
                        **fit,
                        "t_start": float(t[i]),
                        "t_end": float(t[j]),
                        "duration": duration,
                        "n_points": int(j - i + 1),
                    }
                )

    if not refined:
        raise RuntimeError("The coarse exponential interval failed refinement.")
    best = max(
        refined,
        key=lambda item: (
            item["duration"],
            item["r_squared"],
            -item["max_abs_log_residual"],
        ),
    )
    best["search_peak_time"] = peak_time
    best["amplitude_growth_rate"] = 0.5 * best["slope"]
    best["energy_e_folding_time"] = 1.0 / best["slope"]
    best["amplitude_e_folding_time"] = 1.0 / best["amplitude_growth_rate"]
    return best


def make_axes() -> tuple[plt.Figure, plt.Axes]:
    fig = plt.figure(figsize=FIGSIZE, facecolor="white")
    ax = fig.add_axes(AX_RECT)
    ax.set_box_aspect(BOX_ASPECT)
    for spine in ax.spines.values():
        spine.set_linewidth(FRAME_WIDTH)
    ax.tick_params(
        which="major",
        direction="in",
        length=12,
        width=1.2,
        top=True,
        right=True,
    )
    ax.tick_params(
        which="minor",
        direction="in",
        length=6,
        width=1.0,
        top=True,
        right=True,
    )
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    return fig, ax


def save_figure(fig: plt.Figure, stem: str) -> None:
    fig.savefig(HERE / f"{stem}.png", facecolor="white")
    fig.savefig(HERE / f"{stem}.pdf", facecolor="white")
    plt.close(fig)


def plot_kinetic_linear(
    time: np.ndarray, kinetic: np.ndarray, fit: dict[str, float]
) -> None:
    fig, ax = make_axes()
    ax.plot(time, kinetic, color=BLUE, lw=LINE_WIDTH)
    ax.axvspan(
        fit["t_start"],
        fit["t_end"],
        color=SHADE,
        alpha=0.22,
        linewidth=0,
        zorder=0,
    )
    ax.set_xlim(float(time[0]), float(time[-1]))
    ax.set_ylim(bottom=0.0)
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$K=\frac{1}{2}\langle u^2+v^2+w^2\rangle_V$")
    formatter = ScalarFormatter(useMathText=True)
    formatter.set_powerlimits((0, 0))
    ax.yaxis.set_major_formatter(formatter)
    ax.yaxis.get_offset_text().set_fontsize(13)
    save_figure(fig, "kinetic_energy_timeseries")


def plot_kinetic_semilog(
    time: np.ndarray, kinetic: np.ndarray, fit: dict[str, float]
) -> None:
    fig, ax = make_axes()
    ax.semilogy(time, kinetic, color=BLUE, lw=LINE_WIDTH, label=r"$K(t)$")
    mask = (time >= fit["t_start"]) & (time <= fit["t_end"])
    fit_values = np.exp(fit["intercept"] + fit["slope"] * time[mask])
    ax.semilogy(
        time[mask],
        fit_values,
        color=FIT_RED,
        lw=LINE_WIDTH,
        ls="--",
        label="Exponential fit",
    )
    ax.axvspan(
        fit["t_start"],
        fit["t_end"],
        color=SHADE,
        alpha=0.22,
        linewidth=0,
        zorder=0,
    )
    ax.set_xlim(float(time[0]), float(time[-1]))
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$K=\frac{1}{2}\langle u^2+v^2+w^2\rangle_V$")
    ax.yaxis.set_minor_locator(mpl.ticker.LogLocator(subs=np.arange(2, 10) * 0.1))
    ax.text(
        0.055,
        0.94,
        (
            rf"${fit['t_start']:.1f}\leq t\leq {fit['t_end']:.1f}$"
            "\n"
            rf"$d\ln K/dt={fit['slope']:.3f}$"
            "\n"
            rf"$R^2={fit['r_squared']:.6f}$"
        ),
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=12,
    )
    ax.legend(frameon=False, loc="lower right")
    save_figure(fig, "kinetic_energy_linear_instability")


def plot_kinetic_log(time: np.ndarray, kinetic: np.ndarray) -> None:
    fig, ax = make_axes()
    ax.semilogy(time, kinetic, color=BLUE, lw=LINE_WIDTH)
    ax.set_xlim(float(time[0]), float(time[-1]))
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$K=\frac{1}{2}\langle u^2+v^2+w^2\rangle_V$")
    ax.yaxis.set_minor_locator(mpl.ticker.LogLocator(subs=np.arange(2, 10) * 0.1))
    save_figure(fig, "kinetic_energy_timeseries_logy")


def plot_convective_scale(
    time: np.ndarray,
    wavelength_mid: np.ndarray,
    wavelength_vertical: np.ndarray,
    fit: dict[str, float],
) -> None:
    fig, ax = make_axes()
    ax.plot(
        time,
        wavelength_mid,
        color=BLUE,
        lw=LINE_WIDTH,
        label=r"$z\simeq0.5$",
    )
    ax.plot(
        time,
        wavelength_vertical,
        color=LIGHT_BLUE,
        lw=LINE_WIDTH,
        ls="--",
        label=r"$z$-weighted mean",
    )
    ax.axvspan(
        fit["t_start"],
        fit["t_end"],
        color=SHADE,
        alpha=0.22,
        linewidth=0,
        zorder=0,
    )
    ax.set_xlim(float(time[0]), float(time[-1]))
    ax.set_ylim(bottom=0.0)
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"Convective wavelength, $2\pi l_c$")
    ax.legend(frameon=False, loc="lower right")
    save_figure(fig, "convective_wavelength_timeseries")


def write_outputs(
    avg: np.ndarray,
    scale: np.ndarray,
    kinetic: np.ndarray,
    fit: dict[str, float],
    case_info: dict[str, object],
) -> None:
    scale_by_time = {float(row[0]): row for row in scale}
    with (HERE / "reduced_timeseries.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "time",
                "kinetic_energy",
                "lc_w_z05",
                "lc_w_z_weighted",
                "convective_wavelength_z05_2pi_lc",
                "convective_wavelength_z_weighted_2pi_lc",
            ]
        )
        for time, k_value in zip(avg[:, 0], kinetic):
            row = scale_by_time.get(float(time))
            if row is None:
                writer.writerow([time, k_value, "", "", "", ""])
            else:
                writer.writerow(
                    [
                        time,
                        k_value,
                        row[1],
                        row[2],
                        2.0 * np.pi * row[1],
                        2.0 * np.pi * row[2],
                    ]
                )

    mask = (scale[:, 0] >= fit["t_start"]) & (scale[:, 0] <= fit["t_end"])
    summary = {
        "case": case_info,
        "definitions": {
            "kinetic_energy": "K=0.5*avgvar one-based column 4",
            "lc_w": "sum(E_w)/sum(k_h*E_w), excluding k_h=0",
            "plotted_convective_wavelength": "2*pi*lc_w",
        },
        "available_time": {
            "avgvar_start": float(avg[0, 0]),
            "avgvar_end": float(avg[-1, 0]),
            "convective_scale_start": float(scale[0, 0]),
            "convective_scale_end": float(scale[-1, 0]),
        },
        "linear_instability_fit": fit,
        "linear_interval_convective_wavelength": {
            "z05_mean": float(np.mean(2.0 * np.pi * scale[mask, 1])),
            "z05_std": float(np.std(2.0 * np.pi * scale[mask, 1])),
            "z_weighted_mean": float(np.mean(2.0 * np.pi * scale[mask, 2])),
            "z_weighted_std": float(np.std(2.0 * np.pi * scale[mask, 2])),
        },
    }
    with (HERE / "linear_instability_summary.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(summary, handle, indent=2, ensure_ascii=False)


def main(case_info: dict[str, object]) -> None:
    configure_style()
    avg = clean_rows(np.loadtxt(AVGVAR_FILE))
    scale = clean_rows(np.loadtxt(CONVECTIVE_SCALE_FILE))
    if avg.shape[1] < 4:
        raise ValueError("avgvar.out must contain at least four columns.")
    if scale.shape[1] < 3:
        raise ValueError("convective_scale.dat must contain three columns.")

    kinetic = 0.5 * avg[:, 3]
    fit = detect_exponential_interval(avg[:, 0], kinetic)
    wavelength_mid = 2.0 * np.pi * scale[:, 1]
    wavelength_vertical = 2.0 * np.pi * scale[:, 2]

    plot_kinetic_linear(avg[:, 0], kinetic, fit)
    plot_kinetic_semilog(avg[:, 0], kinetic, fit)
    plot_kinetic_log(avg[:, 0], kinetic)
    plot_convective_scale(scale[:, 0], wavelength_mid, wavelength_vertical, fit)
    write_outputs(avg, scale, kinetic, fit, case_info)

    print(json.dumps(fit, indent=2))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot kinetic energy, convective scale, and linear growth."
    )
    parser.add_argument(
        "--directory",
        type=Path,
        default=HERE,
        help="Directory containing avgvar.out and convective_scale.dat.",
    )
    parser.add_argument("--ra", type=float, default=8.0e6)
    parser.add_argument("--pr", type=float, default=0.7)
    parser.add_argument("--ek", type=float, default=7.0e-3)
    parser.add_argument("--ar", type=float, default=16.0)
    parser.add_argument("--beta", type=float, default=1.02)
    parser.add_argument("--qbot", type=float, default=0.5)
    parser.add_argument("--qtop", type=float, default=0.004978)
    parser.add_argument("--resolution", nargs=3, type=int, default=[257, 257, 65])
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    HERE = args.directory.resolve()
    AVGVAR_FILE = HERE / "avgvar.out"
    CONVECTIVE_SCALE_FILE = HERE / "convective_scale.dat"
    main(
        {
            "Ra": args.ra,
            "Pr": args.pr,
            "Ek": args.ek,
            "AR": args.ar,
            "beta": args.beta,
            "qbot": args.qbot,
            "qtop": args.qtop,
            "resolution": args.resolution,
        }
    )
