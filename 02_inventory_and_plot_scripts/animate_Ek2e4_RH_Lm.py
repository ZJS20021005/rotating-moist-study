from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.animation import FuncAnimation, PillowWriter


RUN = Path(
    r"G:\moist convection\Ek2e-4\AR4\Beta1p02\qbot0p5_qtop0p004978"
    r"\N257x257x65\conti1\conti_strict_force_500\run"
)
OUTPUT = Path(
    r"E:\moist RB\rotating_case_inventory\04_outputs_and_figures"
    r"\high_resolution_timeseries_latest_program_20260805"
    r"\force_num_vs_Ek_20260812\Ek2e-4_q_RH_profile_animation_20260815"
)
PROFILE_CSV = OUTPUT / "q_RH_profiles_time_z.csv"
PHASE_MEAN_CSV = OUTPUT / "q_RH_phase_mean_profiles.csv"
LM_FILE = RUN / "diagnostics" / "scale" / "moist_integral_scale.dat"
AVGVAR_FILE = RUN / "data" / "avgvar.out"

BURST_WINDOWS = [(1219.4, 1321.8), (1523.1, 1597.9)]
TIME_MIN = 1210.0
TIME_MAX = 1700.0

FRAME_WIDTH = 4.5
LINE_WIDTH = 3.5
ACCUMULATION_BLUE = (0.11, 0.44, 0.71)
BURST_ORANGE = (0.96, 0.58, 0.19)
BURST_GRAY = (0.55, 0.55, 0.55)
LM_BLUE = (0.00, 0.25, 0.90)
KINETIC_RED = (0.74, 0.14, 0.18)


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
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def style_axis(axis: plt.Axes) -> None:
    for spine in axis.spines.values():
        spine.set_linewidth(FRAME_WIDTH)
        spine.set_color("black")
    axis.tick_params(
        which="major",
        direction="in",
        length=12,
        width=1.2,
        top=True,
        right=True,
        colors="black",
    )
    axis.minorticks_on()
    axis.tick_params(
        which="minor",
        direction="in",
        length=6,
        width=1.0,
        top=True,
        right=True,
        colors="black",
    )
    axis.set_box_aspect(5.2 / 6.5)


def is_burst(time: float) -> bool:
    return any(start <= time <= stop for start, stop in BURST_WINDOWS)


def load_data() -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    profiles = pd.read_csv(PROFILE_CSV)
    times = np.sort(profiles["time"].unique().astype(float))
    z = np.sort(profiles["z"].unique().astype(float))
    rh_profiles = (
        profiles.pivot(index="time", columns="z", values="RH_xy_mean")
        .reindex(index=times, columns=z)
        .to_numpy(float)
    )

    phase_mean = pd.read_csv(PHASE_MEAN_CSV).sort_values("z")
    rh_accumulation = phase_mean["RH_accumulation_mean"].to_numpy(float)
    rh_burst = phase_mean["RH_plume_burst_mean"].to_numpy(float)

    moist = np.loadtxt(LM_FILE)
    lm_mask = (moist[:, 0] >= TIME_MIN) & (moist[:, 0] <= TIME_MAX)
    lm_time = moist[lm_mask, 0]
    two_pi_lm = 2.0 * np.pi * moist[lm_mask, 1]
    avgvar = np.loadtxt(AVGVAR_FILE)
    kinetic_mask = (avgvar[:, 0] >= TIME_MIN) & (avgvar[:, 0] <= TIME_MAX)
    kinetic_time = avgvar[kinetic_mask, 0]
    kinetic_energy = 0.5 * avgvar[kinetic_mask, 3]
    return (
        z,
        times,
        rh_profiles,
        rh_accumulation,
        rh_burst,
        lm_time,
        two_pi_lm,
        kinetic_time,
        kinetic_energy,
    )


def make_animation(series_kind: str) -> tuple[Path, Path]:
    (
        z,
        times,
        rh_profiles,
        rh_accumulation,
        rh_burst,
        lm_time,
        two_pi_lm,
        kinetic_time,
        kinetic_energy,
    ) = load_data()

    if series_kind == "lm":
        series_time = lm_time
        series_values = two_pi_lm
        series_color = LM_BLUE
        series_ylabel = r"$2\pi L_m/H$"
        stem = "Ek2e-4_RH_profile_and_Lm_timeseries"
    elif series_kind == "kinetic":
        series_time = kinetic_time
        series_values = kinetic_energy
        series_color = KINETIC_RED
        series_ylabel = r"$K$"
        stem = "Ek2e-4_RH_profile_and_kinetic_energy_timeseries"
    else:
        raise ValueError(f"Unknown series kind: {series_kind}")

    figure, axes = plt.subplots(1, 2, figsize=(13.0, 5.84), facecolor="white")
    profile_axis, lm_axis = axes
    for axis in axes:
        style_axis(axis)

    rh_min = float(np.nanmin(rh_profiles))
    rh_max = float(np.nanmax(rh_profiles))
    rh_span = max(rh_max - rh_min, 1.0e-3)
    profile_axis.set_xlim(max(0.0, rh_min - 0.05 * rh_span), rh_max + 0.06 * rh_span)
    profile_axis.set_ylim(0.0, 1.0)
    profile_axis.set_xlabel(r"$\langle RH\rangle_{xy}$")
    profile_axis.set_ylabel(r"$z$")

    profile_axis.plot(
        rh_accumulation,
        z,
        color=ACCUMULATION_BLUE,
        lw=2.5,
        ls="--",
        alpha=0.72,
        label="accumulation mean",
    )
    profile_axis.plot(
        rh_burst,
        z,
        color=BURST_ORANGE,
        lw=2.5,
        ls="--",
        alpha=0.72,
        label="plume-burst mean",
    )
    current_profile, = profile_axis.plot(
        rh_profiles[0],
        z,
        color=ACCUMULATION_BLUE,
        lw=LINE_WIDTH + 0.5,
        solid_capstyle="round",
        zorder=5,
    )
    profile_axis.legend(frameon=False, loc="lower left", handlelength=2.8)

    lm_axis.plot(
        series_time,
        series_values,
        color=series_color,
        lw=LINE_WIDTH,
        solid_capstyle="round",
        zorder=3,
    )
    for start, stop in BURST_WINDOWS:
        lm_axis.axvspan(start, stop, color=BURST_GRAY, alpha=0.18, lw=0, zorder=0)
    lm_axis.set_xlim(TIME_MIN, TIME_MAX)
    series_pad = 0.06 * (
        float(np.nanmax(series_values)) - float(np.nanmin(series_values))
    )
    series_bottom = (
        0.0
        if series_kind == "kinetic"
        else float(np.nanmin(series_values)) - series_pad
    )
    lm_axis.set_ylim(series_bottom, float(np.nanmax(series_values)) + series_pad)
    lm_axis.set_xlabel(r"$t$")
    lm_axis.set_ylabel(series_ylabel)
    if series_kind == "kinetic":
        lm_axis.ticklabel_format(
            axis="y", style="sci", scilimits=(-2, 2), useMathText=True
        )
        lm_axis.yaxis.get_offset_text().set_fontsize(13)

    initial_lm = float(np.interp(times[0], series_time, series_values))
    moving_point, = lm_axis.plot(
        [times[0]],
        [initial_lm],
        marker="o",
        ms=11,
        mec="black",
        mew=1.5,
        mfc=ACCUMULATION_BLUE,
        ls="none",
        zorder=6,
    )

    phase_text = figure.text(
        0.5,
        0.985,
        "",
        ha="center",
        va="top",
        fontsize=20,
        color="black",
    )
    figure.subplots_adjust(left=0.10, right=0.98, bottom=0.20, top=0.90, wspace=0.22)

    def update(index: int):
        time = float(times[index])
        burst = is_burst(time)
        phase_color = BURST_ORANGE if burst else ACCUMULATION_BLUE
        phase_label = "plume burst" if burst else "accumulation"
        current_profile.set_data(rh_profiles[index], z)
        current_profile.set_color(phase_color)
        current_lm = float(np.interp(time, series_time, series_values))
        moving_point.set_data([time], [current_lm])
        moving_point.set_markerfacecolor(phase_color)
        phase_text.set_text(rf"$t={time:.0f}$    {phase_label}")
        phase_text.set_color(phase_color)
        return current_profile, moving_point, phase_text

    animation = FuncAnimation(
        figure,
        update,
        frames=len(times),
        interval=180,
        blit=False,
        repeat=True,
    )
    gif_path = OUTPUT / f"{stem}.gif"
    mp4_path = OUTPUT / f"{stem}.mp4"
    animation.save(gif_path, writer=PillowWriter(fps=6), dpi=120)

    import imageio_ffmpeg

    mpl.rcParams["animation.ffmpeg_path"] = imageio_ffmpeg.get_ffmpeg_exe()
    animation.save(mp4_path, writer="ffmpeg", fps=6, dpi=150, bitrate=2600)
    plt.close(figure)
    return gif_path, mp4_path


def main() -> None:
    configure_style()
    for series_kind in ("lm", "kinetic"):
        gif_path, mp4_path = make_animation(series_kind)
        print(gif_path)
        print(mp4_path)


if __name__ == "__main__":
    main()
