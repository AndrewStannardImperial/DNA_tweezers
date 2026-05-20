from pathlib import Path

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


mpl.rcParams.update(
    {
        "text.usetex": False,
        "mathtext.fontset": "dejavusans",
        "font.family": "helvetica",
        "font.size": 15,
        "axes.titlesize": 22,
        "axes.labelsize": 15,
        "xtick.labelsize": 16,
        "ytick.labelsize": 16,
        "axes.linewidth": 1.2,
        "xtick.major.width": 1.5,
        "ytick.major.width": 1.5,
        "xtick.major.size": 6,
        "xtick.minor.size": 3,
        "ytick.major.size": 6,
    }
)


conc_exp = np.array([0, 1, 2, 5, 10, 20, 50, 100, 200])

rel_FR_Mg = np.array([1, 0.995, 0.994, 0.991, 0.986, 0.983, 0.976, 0.969, 0.963])
rel_FR_Ca = np.array([1, 0.993, 0.990, 0.984, 0.979, 0.972, 0.966, 0.956, 0.950])
rel_FR_Mg_uncert = np.array([0, 0.002, 0.002, 0.003, 0.002, 0.002, 0.002, 0.002, 0.002])
rel_FR_Ca_uncert = np.array([0, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.003])

FR_Mg = np.array([5.71, 5.68, 5.67, 5.66, 5.63, 5.61, 5.57, 5.53, 5.50]) * 10
FR_Ca = np.array([5.71, 5.67, 5.65, 5.62, 5.59, 5.55, 5.52, 5.46, 5.42]) * 10
FR_Mg_uncert = np.array([0.4, 0.4, 0.4, 0.4, 0.39, 0.39, 0.39, 0.39, 0.39]) * 10
FR_Ca_uncert = np.array([0.4, 0.4, 0.4, 0.39, 0.39, 0.39, 0.39, 0.38, 0.38]) * 10


def R_c(c, a, k):
    return 1 - a * np.log(1 + c / k)


def main():
    params_Mg = curve_fit(R_c, conc_exp, rel_FR_Mg)[0]
    params_Ca = curve_fit(R_c, conc_exp, rel_FR_Ca)[0]

    fig, ax = plt.subplots(1, 2, figsize=(8, 4))
    conc_interp = np.linspace(0, 250, 10000)

    ax[0].plot(
        conc_interp,
        R_c(conc_interp, params_Mg[0], params_Mg[1]),
        ls="--",
        c="#07baec",
        label=r"$\mathrm{Mg}^{2+}$ fit",
    )
    ax[0].plot(
        conc_interp,
        R_c(conc_interp, params_Ca[0], params_Ca[1]),
        ls="--",
        c="#f59621",
        label=r"$\mathrm{Ca}^{2+}$ fit",
    )
    ax[1].plot(
        conc_interp,
        FR_Mg[0] * R_c(conc_interp, params_Mg[0], params_Mg[1]),
        ls="--",
        c="#07baec",
    )
    ax[1].plot(
        conc_interp,
        FR_Ca[0] * R_c(conc_interp, params_Ca[0], params_Ca[1]),
        ls="--",
        c="#f59621",
    )

    ax[0].errorbar(
        conc_exp,
        rel_FR_Mg,
        yerr=rel_FR_Mg_uncert,
        marker="^",
        markersize=8,
        c="#07baec",
        ls="",
        capsize=3,
        elinewidth=1.5,
        capthick=2.5,
        markeredgewidth=1,
        markeredgecolor="black",
        label=r"$\mathrm{Mg}^{2+}$ expt.",
    )
    ax[0].errorbar(
        conc_exp,
        rel_FR_Ca,
        yerr=rel_FR_Ca_uncert,
        marker="^",
        markersize=8,
        c="#f59621",
        ls="",
        capsize=3,
        elinewidth=1.5,
        capthick=2.5,
        markeredgewidth=1,
        markeredgecolor="black",
        label=r"$\mathrm{Ca}^{2+}$ expt.",
    )
    ax[1].errorbar(
        conc_exp,
        FR_Mg,
        yerr=FR_Mg_uncert,
        marker="^",
        markersize=8,
        c="#07baec",
        ls="",
        capsize=3,
        elinewidth=1.5,
        capthick=2.5,
        markeredgewidth=1,
        markeredgecolor="black",
    )
    ax[1].errorbar(
        conc_exp,
        FR_Ca,
        yerr=FR_Ca_uncert,
        marker="^",
        markersize=8,
        c="#f59621",
        ls="",
        capsize=3,
        elinewidth=1.5,
        capthick=2.5,
        markeredgewidth=1,
        markeredgecolor="black",
    )

    ax[0].spines[["top", "right"]].set_visible(False)
    ax[1].spines[["top", "right"]].set_visible(False)

    ax[0].set_xlabel(r"[M$^{2+}$] (mM)")
    ax[1].set_xlabel(r"[M$^{2+}$] (mM)")
    ax[0].set_ylabel("Relative Forster Radius")
    ax[1].set_ylabel(r"Forster Radius ($\AA$)")

    fig.subplots_adjust(wspace=0.4)
    ax[0].set_xscale("symlog")
    ax[1].set_xscale("symlog")
    ax[0].get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    ax[1].get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    ax[1].set_xticks([0, 1, 10, 100])
    ax[0].set_xticks([0, 1, 10, 100])
    ax[0].tick_params(which="both", direction="in", top=False, right=False)
    ax[1].tick_params(which="both", direction="in", top=False, right=False)
    ax[0].set_xlim(-0.5, 250)
    ax[1].set_xlim(-0.5, 250)
    fig.legend(fontsize="x-small", fancybox=False, edgecolor="white", bbox_to_anchor=(1.1, 0.55))

    out_path = Path(__file__).with_name("FR_fit.pdf")
    plt.savefig(out_path, format="pdf", dpi=300, bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    main()
