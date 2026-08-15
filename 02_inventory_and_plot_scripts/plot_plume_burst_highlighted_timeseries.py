from pathlib import Path
import math

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(r"E:\moist RB\rotating_case_inventory\04_outputs_and_figures\high_resolution_timeseries_latest_program_20260805")
DATA = ROOT / "high_resolution_timeseries_long.csv"
META = ROOT / "selected_case_metadata.csv"
WINDOWS = {
    2e-4: [(353.3, 432.6), (622.9, 732.1), (939.4, 1029.8)],
    5e-4: [(390.5, 429.0), (451.8, 488.5)],
    7e-4: [(447.9, 469.6), (491.1, 523.1), (580.6, 684.6)],
}
FRAME, LINE = 4.5, 3.5
RED, BLUE, YELLOW = (0.74, 0.14, 0.18), (0.00, 0.25, 0.90), (1.0, 0.82, 0.12)


def configure():
    mpl.rcParams.update({"font.family":"Times New Roman", "mathtext.fontset":"stix",
        "axes.linewidth":FRAME, "axes.labelsize":24, "xtick.labelsize":13,
        "ytick.labelsize":13, "legend.fontsize":10, "figure.dpi":180,
        "savefig.dpi":300, "pdf.fonttype":42})


def style(ax):
    for spine in ax.spines.values(): spine.set_linewidth(FRAME)
    ax.tick_params(which="major", direction="in", length=12, width=1.2, top=True, right=True)
    ax.minorticks_on(); ax.tick_params(which="minor", direction="in", length=6, width=1.0, top=True, right=True)
    ax.set_box_aspect(5.2/6.5)


def shade(ax, windows):
    for start, stop in windows: ax.axvspan(start, stop, color=YELLOW, alpha=0.20, lw=0, zorder=0)


def axis(wide=False):
    fig = plt.figure(figsize=((7.8 if wide else 6.5), 5.84), facecolor="white")
    ax = fig.add_axes(([0.160,0.260,0.660,0.705] if wide else [0.186,0.260,0.792,0.705])); style(ax)
    return fig, ax


def save(fig, stem):
    fig.savefig(stem.with_suffix('.png'), facecolor='white'); fig.savefig(stem.with_suffix('.pdf'), facecolor='white'); plt.close(fig)


def main():
    configure(); data=pd.read_csv(DATA); meta=pd.read_csv(META)
    for target, windows in WINDOWS.items():
        mrow=meta[np.isclose(meta.Ek, target)].iloc[0]
        stem=(f"Ra8e6_Ek{target:.0e}_AR{int(mrow.AR)}_Beta1p02_qbot0p5_"
              f"N{int(mrow.Nx)}x{int(mrow.Ny)}x{int(mrow.Nz)}")
        case=data[np.isclose(data.Ek, target)]

        fig, left=axis(wide=True); right=left.twinx(); style(right)
        k=case[case.metric.eq('kinetic')].sort_values('time'); m=case[case.metric.eq('mprime')].sort_values('time')
        left.plot(k.time,k.value,color=RED,lw=LINE,label=r'$K$'); right.plot(m.time,m.value,color=BLUE,lw=LINE,label=r'$\langle\langle m\prime^2\rangle_{xy}\rangle_z$')
        shade(left,windows); left.set_xlabel(r'$t$'); left.set_ylabel(r'$K$'); right.set_ylabel(r'$\langle\langle m\prime^2\rangle_{xy}\rangle_z$')
        left.set_ylim(bottom=0); right.set_ylim(bottom=0)
        left.legend(left.lines[:1]+right.lines[:1],[r'$K$',r'$\langle\langle m\prime^2\rangle_{xy}\rangle_z$'],frameon=False,loc='lower right')
        save(fig, ROOT/'individual_cases'/(stem+'_kinetic_mprime_dual_axis_timeseries'))

        fig, ax=axis(); lm=case[case.metric.eq('mse_spectral')].sort_values('time')
        ax.plot(lm.time,lm.value,color=BLUE,lw=LINE); shade(ax,windows); ax.set_xlabel(r'$t$'); ax.set_ylabel(r'$2\pi L_m$'); ax.set_ylim(bottom=0)
        save(fig, ROOT/'individual_cases'/(stem+'_lm_timeseries'))

        fig, ax=axis(); nu=case[case.metric.eq('num')].sort_values('time')
        ax.plot(nu.time,nu.value,color=BLUE,lw=LINE); shade(ax,windows); ax.set_xlabel(r'$t$'); ax.set_ylabel(r'$Nu_m(t)$')
        save(fig, ROOT/'num_timeseries'/'individual_cases'/(f'Ra8e6_Ek{str(target).replace(".","p")}_Num_timeseries'))


if __name__=='__main__': main()
