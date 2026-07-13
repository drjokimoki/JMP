"""
k* vs T -- exogenous-only figure.

Replaces plot_kstar_vs_T_v2.py (which plotted endogenous Toeplitz/Factor
RMSFE curves). The endogenous case has been removed from the paper; this
figure uses only the exogenous population-grid and Monte-Carlo data already
reported in Table~\\ref{tab:exo_scaling_h6} (B=3, H~6.08, m=50, 500 reps),
alongside the T^{1/3} and T^{1/2} reference rates that bracket the rigorous
rate window of Theorem~\\ref{thm:ratewindow}.
"""

import numpy as np
import matplotlib.pyplot as plt

# -- Data from Table tab:exo_scaling_h6 (B=3, H~6.08, m=50, 500 reps) --------
T_exo         = np.array([80, 240, 640, 1000])
kstar_grid    = np.array([2, 3, 4, 5])
kstar_formula = np.array([2.0, 2.9, 4.0, 4.6])

# Fitted power law from the paper: k*(grid) ~ C * T^0.352 (R^2=0.994)
slope_fit = 0.352
C_fit = kstar_grid[1] / T_exo[1] ** slope_fit   # anchor at T=240

# -- Reference curves ---------------------------------------------------------
T_ref = np.linspace(60, 1100, 500)
sqrt_ref = np.sqrt(T_ref)
cube_ref = C_fit * T_ref ** (1 / 3)
fit_ref = C_fit * T_ref ** slope_fit

COL_DATA = "#0072B2"
COL_FORMULA = "#D55E00"
COL_SQRT = "#009E73"
COL_CUBE = "#CC79A7"

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 10.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 180,
})

fig, axes = plt.subplots(1, 2, figsize=(11, 5))

for ax, log_scale in zip(axes, [False, True]):
    if log_scale:
        ax.set_xscale("log")
        ax.set_yscale("log")

    ax.plot(T_ref, sqrt_ref, "--", color=COL_SQRT, lw=1.5, alpha=0.85,
            label=r"$\sqrt{T}$ reference", zorder=1)
    ax.plot(T_ref, cube_ref, ":", color=COL_CUBE, lw=1.5, alpha=0.85,
            label=r"$T^{1/3}$ reference (fitted constant)", zorder=1)
    ax.plot(T_ref, fit_ref, "-", color=COL_DATA, lw=1.0, alpha=0.35, zorder=2)

    ax.plot(T_exo, kstar_formula, "D--", color=COL_FORMULA, lw=1.3, ms=6,
            alpha=0.85, label=r"Exogenous $k^*$ (power-law formula)", zorder=4)
    ax.plot(T_exo, kstar_grid, "o-", color=COL_DATA, lw=2.0, ms=7,
            label=r"Exogenous $k^*$ (population grid minimiser)", zorder=5)

    if log_scale:
        ax.annotate(r"fitted slope $0.352\approx 1/3$",
                    xy=(400, C_fit * 400 ** slope_fit),
                    xytext=(150, 8), fontsize=9, color=COL_DATA,
                    arrowprops=dict(arrowstyle="->", color=COL_DATA, lw=0.8))
        ax.annotate(r"slope = 1/2",
                    xy=(700, np.sqrt(700)),
                    xytext=(500, 40), fontsize=9, color=COL_SQRT,
                    arrowprops=dict(arrowstyle="->", color=COL_SQRT, lw=0.8))
        ax.annotate(r"slope = 1/3",
                    xy=(900, C_fit * 900 ** (1 / 3)),
                    xytext=(700, 3.2), fontsize=9, color=COL_CUBE,
                    arrowprops=dict(arrowstyle="->", color=COL_CUBE, lw=0.8))

    ax.set_xlabel(r"Sample size $T$", fontsize=12)
    ax.set_ylabel(r"Optimal $k^*$", fontsize=12)
    ax.set_title("Linear scale" if not log_scale else "Log--log scale",
                 fontsize=11, pad=6)
    ax.set_xlim(55, 1150)
    if not log_scale:
        ax.set_ylim(0, 34)
        ax.fill_between(T_ref, 0.8 * sqrt_ref, 1.2 * sqrt_ref,
                         color=COL_SQRT, alpha=0.06, label=r"$\pm20\%$ of $\sqrt{T}$")
    else:
        ax.set_ylim(1.5, 34)
    ax.grid(True, which="both" if log_scale else "major",
            ls=":", lw=0.4, alpha=0.45)

handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="lower center", ncol=2,
           fontsize=9, frameon=True, framealpha=0.92,
           bbox_to_anchor=(0.5, -0.14))

fig.suptitle(
    r"Optimal subset size $k^*$ vs sample size $T$ (exogenous, $H\approx 6.08$, $m=50$)",
    fontsize=11, y=1.02)

fig.tight_layout()

for ext in ("pdf", "png"):
    path = f"/Users/boris/Desktop/RSShapley/Combination_Monte/Final/codess/JMP/kstar_vs_T_v2.{ext}"
    plt.savefig(path, bbox_inches="tight", dpi=180)
    print("Saved:", path)

plt.close()
