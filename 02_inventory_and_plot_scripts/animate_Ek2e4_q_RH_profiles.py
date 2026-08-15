from __future__ import annotations

import re
from pathlib import Path

import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.animation import FuncAnimation, PillowWriter


RUN = Path(
    r"G:\moist convection\Ek2e-4\AR4\Beta1p02\qbot0p5_qtop0p004978"
    r"\N257x257x65\conti1\conti_strict_force_500\run"
)
MOVIE = RUN / "movie"
OUTPUT = Path(
    r"E:\moist RB\rotating_case_inventory\04_outputs_and_figures"
    r"\high_resolution_timeseries_latest_program_20260805"
    r"\force_num_vs_Ek_20260812\Ek2e-4_q_RH_profile_animation_20260815"
)

ALPHA_QS = 3.0
BETA_QS = 1.02
BURST_WINDOWS = [(1219.4, 1321.8), (1523.1, 1597.9)]

FRAME_WIDTH = 4.5
LINE_WIDTH = 3.5
ACCUMULATION_BLUE = (0.11, 0.44, 0.71)
BURST_ORANGE = (0.96, 0.58, 0.19)
REFERENCE_GRAY = (0.30, 0.30, 0.30)


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


def frame_time(xmf_path: Path) -> float:
    text = xmf_path.read_text(encoding="utf-8", errors="replace")
    match = re.search(r'<Time\s+Value="\s*([^\"]+)"', text)
    if not match:
        raise ValueError(f"No XDMF time in {xmf_path}")
    return float(match.group(1))


def in_burst(time: float) -> bool:
    return any(start <= time <= stop for start, stop in BURST_WINDOWS)


def load_profiles() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, list[str]]:
    with h5py.File(MOVIE / "cordin_info.h5", "r") as source:
        z = np.asarray(source["z"], dtype=float)

    records: list[tuple[float, np.ndarray, np.ndarray, str]] = []
    skipped: list[str] = []
    for field_path in sorted(MOVIE.glob("field*.h5")):
        xmf_path = field_path.with_suffix(".xmf")
        if field_path.stat().st_size == 0 or not xmf_path.exists():
            skipped.append(field_path.name)
            continue
        try:
            time = frame_time(xmf_path)
            with h5py.File(field_path, "r") as source:
                b = np.asarray(source["DSAL_me"], dtype=float)
                rh = np.asarray(source["RH_me"], dtype=float)
        except (OSError, KeyError, ValueError) as exc:
            skipped.append(f"{field_path.name}: {exc}")
            continue

        temperature = b - BETA_QS * z[:, None, None]
        qsat = np.exp(ALPHA_QS * temperature)
        q = rh * qsat
        q_profile = np.mean(q, axis=(1, 2))
        rh_profile = np.mean(rh, axis=(1, 2))
        phase = "plume burst" if in_burst(time) else "accumulation"
        records.append((time, q_profile, rh_profile, phase))

    records.sort(key=lambda item: item[0])
    times = np.asarray([item[0] for item in records], dtype=float)
    q_profiles = np.stack([item[1] for item in records])
    rh_profiles = np.stack([item[2] for item in records])
    phases = np.asarray([item[3] for item in records], dtype=object)
    return z, times, q_profiles, rh_profiles, phases.tolist(), skipped


def padded_limits(values: np.ndarray, lower_zero: bool = False) -> tuple[float, float]:
    lo = float(np.nanmin(values))
    hi = float(np.nanmax(values))
    span = max(hi - lo, abs(hi) * 0.08, 1.0e-8)
    left = 0.0 if lower_zero and lo >= 0.0 else lo - 0.06 * span
    return left, hi + 0.08 * span


def save_reduced_data(
    z: np.ndarray,
    times: np.ndarray,
    q_profiles: np.ndarray,
    rh_profiles: np.ndarray,
    phases: list[str],
    skipped: list[str],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    burst_mask = np.asarray([phase == "plume burst" for phase in phases])
    accumulation_mask = ~burst_mask
    q_acc = np.mean(q_profiles[accumulation_mask], axis=0)
    q_burst = np.mean(q_profiles[burst_mask], axis=0)
    rh_acc = np.mean(rh_profiles[accumulation_mask], axis=0)
    rh_burst = np.mean(rh_profiles[burst_mask], axis=0)

    long_rows = []
    for index, time in enumerate(times):
        for level, height in enumerate(z):
            long_rows.append(
                {
                    "time": time,
                    "phase": phases[index],
                    "z": height,
                    "q_xy_mean": q_profiles[index, level],
                    "RH_xy_mean": rh_profiles[index, level],
                }
            )
    pd.DataFrame(long_rows).to_csv(OUTPUT / "q_RH_profiles_time_z.csv", index=False)
    pd.DataFrame(
        {
            "z": z,
            "q_accumulation_mean": q_acc,
            "q_plume_burst_mean": q_burst,
            "RH_accumulation_mean": rh_acc,
            "RH_plume_burst_mean": rh_burst,
        }
    ).to_csv(OUTPUT / "q_RH_phase_mean_profiles.csv", index=False)
    pd.DataFrame(BURST_WINDOWS, columns=["burst_start", "burst_end"]).to_csv(
        OUTPUT / "burst_windows.csv", index=False
    )
    (OUTPUT / "processing_note.txt").write_text(
        "Horizontal means are computed snapshot-first.\n"
        "q = RH * exp(alpha_qs * (b - beta_qs*z)), alpha_qs=3, beta_qs=1.02.\n"
        "Burst windows are inherited from the low-flat 2*pi*Lm plateau diagnostic.\n"
        f"Skipped frames: {skipped}\n",
        encoding="utf-8",
    )
    return q_acc, q_burst, rh_acc, rh_burst


def make_phase_mean_figure(
    z: np.ndarray,
    q_acc: np.ndarray,
    q_burst: np.ndarray,
    rh_acc: np.ndarray,
    rh_burst: np.ndarray,
) -> None:
    figure, axes = plt.subplots(1, 2, figsize=(13.0, 5.84), sharey=True)
    figure.patch.set_facecolor("white")
    for axis in axes:
        style_axis(axis)
        axis.set_ylim(0.0, 1.0)
    axes[0].plot(q_acc, z, color=ACCUMULATION_BLUE, lw=LINE_WIDTH, label="accumulation")
    axes[0].plot(q_burst, z, color=BURST_ORANGE, lw=LINE_WIDTH, label="plume burst")
    axes[1].plot(rh_acc, z, color=ACCUMULATION_BLUE, lw=LINE_WIDTH)
    axes[1].plot(rh_burst, z, color=BURST_ORANGE, lw=LINE_WIDTH)
    axes[0].set_xlabel(r"$\langle q\rangle_{xy,t}$")
    axes[1].set_xlabel(r"$\langle RH\rangle_{xy,t}$")
    axes[0].set_ylabel(r"$z$")
    axes[0].legend(frameon=False, loc="best")
    figure.subplots_adjust(left=0.10, right=0.98, bottom=0.20, top=0.97, wspace=0.18)
    figure.savefig(OUTPUT / "Ek2e-4_q_RH_phase_mean_profiles.png", facecolor="white")
    figure.savefig(OUTPUT / "Ek2e-4_q_RH_phase_mean_profiles.pdf", facecolor="white")
    plt.close(figure)


def make_animation(
    z: np.ndarray,
    times: np.ndarray,
    q_profiles: np.ndarray,
    rh_profiles: np.ndarray,
    phases: list[str],
    q_acc: np.ndarray,
    q_burst: np.ndarray,
    rh_acc: np.ndarray,
    rh_burst: np.ndarray,
) -> Path:
    figure, axes = plt.subplots(1, 2, figsize=(13.0, 6.5), sharey=True)
    figure.patch.set_facecolor("white")
    for axis in axes:
        style_axis(axis)
        axis.set_ylim(0.0, 1.0)

    axes[0].set_xlim(*padded_limits(q_profiles, lower_zero=True))
    axes[1].set_xlim(*padded_limits(rh_profiles, lower_zero=True))
    axes[0].set_xlabel(r"$\langle q\rangle_{xy}$")
    axes[1].set_xlabel(r"$\langle RH\rangle_{xy}$")
    axes[0].set_ylabel(r"$z$")

    axes[0].plot(q_acc, z, color=ACCUMULATION_BLUE, lw=2.4, ls="--", alpha=0.75)
    axes[0].plot(q_burst, z, color=BURST_ORANGE, lw=2.4, ls="--", alpha=0.75)
    axes[1].plot(rh_acc, z, color=ACCUMULATION_BLUE, lw=2.4, ls="--", alpha=0.75)
    axes[1].plot(rh_burst, z, color=BURST_ORANGE, lw=2.4, ls="--", alpha=0.75)
    current_q, = axes[0].plot([], [], color=ACCUMULATION_BLUE, lw=LINE_WIDTH + 0.5)
    current_rh, = axes[1].plot([], [], color=ACCUMULATION_BLUE, lw=LINE_WIDTH + 0.5)

    phase_text = figure.text(
        0.5,
        0.965,
        "",
        ha="center",
        va="top",
        fontsize=20,
        color="black",
    )
    figure.text(
        0.50,
        0.035,
        "dashed: phase mean",
        ha="center",
        va="bottom",
        fontsize=12,
        color=REFERENCE_GRAY,
    )

    timeline = figure.add_axes([0.18, 0.09, 0.64, 0.045])
    timeline.set_xlim(times.min(), times.max())
    timeline.set_ylim(0.0, 1.0)
    timeline.set_yticks([])
    timeline.tick_params(axis="x", direction="in", length=5, width=1.0, labelsize=10)
    for spine in timeline.spines.values():
        spine.set_linewidth(1.2)
    timeline.axvspan(times.min(), times.max(), color=ACCUMULATION_BLUE, alpha=0.12)
    for start, stop in BURST_WINDOWS:
        timeline.axvspan(start, stop, color=BURST_ORANGE, alpha=0.30)
    time_marker = timeline.axvline(times[0], color="black", lw=2.0)
    timeline.set_xlabel(r"$t$", fontsize=13, labelpad=-1)

    figure.subplots_adjust(left=0.10, right=0.98, bottom=0.20, top=0.90, wspace=0.18)

    def update(index: int):
        phase = phases[index]
        color = BURST_ORANGE if phase == "plume burst" else ACCUMULATION_BLUE
        current_q.set_data(q_profiles[index], z)
        current_rh.set_data(rh_profiles[index], z)
        current_q.set_color(color)
        current_rh.set_color(color)
        time_marker.set_xdata([times[index], times[index]])
        phase_text.set_text(rf"$t={times[index]:.0f}$    {phase}")
        phase_text.set_color(color)
        return current_q, current_rh, time_marker, phase_text

    animation = FuncAnimation(
        figure,
        update,
        frames=len(times),
        interval=180,
        blit=False,
        repeat=True,
    )
    gif_path = OUTPUT / "Ek2e-4_q_RH_profiles_accumulation_burst.gif"
    animation.save(gif_path, writer=PillowWriter(fps=6), dpi=120)

    try:
        import imageio_ffmpeg

        mpl.rcParams["animation.ffmpeg_path"] = imageio_ffmpeg.get_ffmpeg_exe()
        mp4_path = OUTPUT / "Ek2e-4_q_RH_profiles_accumulation_burst.mp4"
        animation.save(mp4_path, writer="ffmpeg", fps=6, dpi=150, bitrate=2600)
    except (ImportError, RuntimeError, FileNotFoundError):
        pass
    plt.close(figure)
    return gif_path


def main() -> None:
    configure_style()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    z, times, q_profiles, rh_profiles, phases, skipped = load_profiles()
    q_acc, q_burst, rh_acc, rh_burst = save_reduced_data(
        z, times, q_profiles, rh_profiles, phases, skipped
    )
    make_phase_mean_figure(z, q_acc, q_burst, rh_acc, rh_burst)
    animation_path = make_animation(
        z,
        times,
        q_profiles,
        rh_profiles,
        phases,
        q_acc,
        q_burst,
        rh_acc,
        rh_burst,
    )
    print(f"frames={len(times)} time={times.min():.1f}-{times.max():.1f}")
    print(f"burst_frames={sum(phase == 'plume burst' for phase in phases)}")
    print(f"skipped={skipped}")
    print(animation_path)


if __name__ == "__main__":
    main()
