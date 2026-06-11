"""
Plot optimal k* vs T for exogenous (H=1.05, H=6 placeholder) and
endogenous (Toeplitz, Factor) cases.
Generates two panels: log-log and linear scale.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# ── Data ─────────────────────────────────────────────────────────────────────

# Exogenous H=1.05 (population grid k*, always = 1; empirical k* shown)
T_exo_105   = [80,  240,  640,  1000]
kstar_exo_105_pop = [1, 1, 1, 1]         # population grid minimiser
kstar_exo_105_emp = [2, 2, 4,  3]        # empirical (finite-sample MSPE)

# Exogenous H=6 — only T=80 available; rest pending (H=6 scaling still running)
T_exo_6_avail   = [80]
kstar_exo_6_pop  = [3]   # from joint_B_H at T=240: k*=3; T=80 from H6 CSV
# We compute what the cube-root formula predicts across T for H=6
# From joint_B_H: C_delta=0.1006, A=s+Q_opt=2.06, B=3, alpha_hat=1.72
C_delta_H6 = 0.1006
A_H6       = 2.06
alpha_H6   = 1.72  # fitted from MC log-log

# Endogenous Toeplitz (B4, empirical RMSFE minimiser)
T_endo      = [100, 150, 200, 300, 500, 800]
kstar_toep  = [10,  10,  14,  14,  18,  22]
kstar_fact  = [14,  14,  14,  18,  26,  26]

# Fitted log-log slopes and intercepts (from OLS, computed previously)
# Toeplitz: log(k*) = 0.444 + 0.394*log(T)
# Factor:   log(k*) = 0.836 + 0.367*log(T)
slope_toep, intercept_toep = 0.394, np.exp(0.444)
slope_fact, intercept_fact = 0.367, np.exp(0.836)

# ── Reference curves ─────────────────────────────────────────────────────────
T_ref = np.linspace(60, 1100, 400)

# sqrt(T) reference
sqrt_T = np.sqrt(T_ref)

# T^{1/3} reference — scaled to match endogenous Toeplitz at T=200
C_13 = 14 / 200**(1/3)
cube_T = C_13 * T_ref**(1/3)

# Exogenous H=6 cube-root formula: k* = (2*C_delta*T/A)^{1/3}
exo_H6_formula = (2 * C_delta_H6 * T_ref / A_H6)**(1/3)

# Exogenous H=6 power-law formula: k* = (C*alpha*T/A)^{1/(alpha+1)}
exo_H6_powerlaw = (C_delta_H6 * alpha_H6 * T_ref / A_H6)**(1/(alpha_H6 + 1))

# Endogenous fitted lines
endo_toep_fit = intercept_toep * T_ref**slope_toep
endo_fact_fit = intercept_fact * T_ref**slope_fact

# ── Colours & markers ────────────────────────────────────────────────────────
COL_EXO105  = "#888888"
COL_EXO6    = "#E69F00"
COL_TOEP    = "#0072B2"
COL_FACT    = "#D55E00"
COL_SQRT    = "#009E73"
COL_CUBE    = "#CC79A7"

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

fig = plt.figure(figsize=(13, 5.5))
gs  = gridspec.GridSpec(1, 2, wspace=0.38)

for ax_idx, (ax, log_scale) in enumerate(
        [(fig.add_subplot(gs[0]), False),
         (fig.add_subplot(gs[1]), True)]):

    if log_scale:
        ax.set_xscale("log"); ax.set_yscale("log")

    # ── Reference lines ──────────────────────────────────────────────────
    ax.plot(T_ref, sqrt_T,   "--", color=COL_SQRT,  lw=1.2,
            label=r"$\sqrt{T}$", zorder=1)
    ax.plot(T_ref, cube_T,   ":",  color=COL_CUBE,  lw=1.2,
            label=r"$2.19\cdot T^{1/3}$ (endo. scaling)", zorder=1)

    # Exo H=6 power-law prediction
    ax.plot(T_ref, exo_H6_powerlaw, "-.", color=COL_EXO6, lw=1.0, alpha=0.7,
            label=r"Exo. $H=6$ power-law pred. ($\hat\alpha=1.72$)", zorder=1)

    # ── Endogenous MC data ────────────────────────────────────────────────
    ax.plot(T_endo, kstar_toep, "o-", color=COL_TOEP, lw=1.8, ms=6,
            label=r"Endo. Toeplitz $k^*$ (RMSFE)", zorder=3)
    ax.plot(T_endo, kstar_fact, "s-", color=COL_FACT, lw=1.8, ms=6,
            label=r"Endo. Factor $k^*$ (RMSFE)", zorder=3)

    # Fitted regression lines for endogenous
    ax.plot(T_ref, endo_toep_fit, "-", color=COL_TOEP, lw=0.8, alpha=0.4, zorder=2)
    ax.plot(T_ref, endo_fact_fit, "-", color=COL_FACT, lw=0.8, alpha=0.4, zorder=2)

    # ── Exogenous H=1.05 population k* (always 1, for reference) ─────────
    ax.axhline(1, color=COL_EXO105, lw=1.0, ls="--", alpha=0.5,
               label=r"Exo. $H=1.05$ pop. $k^*=1$ (flat)")

    # Exo H=1.05 empirical k*
    ax.plot(T_exo_105, kstar_exo_105_emp, "^", color=COL_EXO105, ms=6,
            label=r"Exo. $H=1.05$ emp. $k^*$", zorder=3)

    # ── Exogenous H=6: only T=80 available; mark as pending ──────────────
    ax.plot(T_exo_6_avail, kstar_exo_6_pop, "D", color=COL_EXO6, ms=7,
            zorder=4, label=r"Exo. $H=6$ pop. $k^*$ (T=80 only; rest pending)")

    # Annotation: pending
    ax.annotate("H=6 scaling\npending",
                xy=(80, 3), xytext=(140 if not log_scale else 110, 6),
                fontsize=7.5, color=COL_EXO6,
                arrowprops=dict(arrowstyle="->", color=COL_EXO6, lw=0.8))

    # ── Slope annotations (log-log panel only) ────────────────────────────
    if log_scale:
        # Annotate slope = 0.394 near endogenous Toeplitz
        ax.annotate(r"slope$\approx$0.39",
                    xy=(300, intercept_toep*300**slope_toep),
                    xytext=(180, 20), fontsize=8, color=COL_TOEP,
                    arrowprops=dict(arrowstyle="->", color=COL_TOEP, lw=0.7))
        ax.annotate(r"slope$\approx$0.37",
                    xy=(300, intercept_fact*300**slope_fact),
                    xytext=(400, 22), fontsize=8, color=COL_FACT,
                    arrowprops=dict(arrowstyle="->", color=COL_FACT, lw=0.7))
        ax.annotate(r"slope$=\frac{1}{2}$",
                    xy=(600, np.sqrt(600)),
                    xytext=(650, 32), fontsize=8, color=COL_SQRT,
                    arrowprops=dict(arrowstyle="->", color=COL_SQRT, lw=0.7))
        ax.annotate(r"slope$=\frac{1}{3}$",
                    xy=(600, C_13*600**(1/3)),
                    xytext=(400, 10), fontsize=8, color=COL_CUBE,
                    arrowprops=dict(arrowstyle="->", color=COL_CUBE, lw=0.7))

    # ── Axes ──────────────────────────────────────────────────────────────
    ax.set_xlabel(r"Sample size $T$", fontsize=11)
    ax.set_ylabel(r"Optimal $k^*$", fontsize=11)
    title = ("Linear scale" if not log_scale else "Log–log scale")
    ax.set_title(title, fontsize=11, pad=6)
    ax.set_xlim(55, 1100)
    if not log_scale:
        ax.set_ylim(-1, 35)
        ax.axhspan(0, 1.5, alpha=0.04, color=COL_EXO105)  # floor region
    else:
        ax.set_ylim(0.7, 50)
    ax.grid(True, which="both" if log_scale else "major",
            ls=":", lw=0.4, alpha=0.5)

# ── Shared legend below ───────────────────────────────────────────────────────
handles, labels = fig.axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="lower center", ncol=3,
           fontsize=8.5, frameon=True, framealpha=0.9,
           bbox_to_anchor=(0.5, -0.18))

fig.suptitle(
    r"Optimal subset size $k^*$ vs sample size $T$: exogenous and endogenous cases"
    "\n" r"($m=50$; endo. curves = RMSFE minimiser over 500 MC reps)",
    fontsize=11, y=1.01)

outpath = "/Users/boris/Desktop/RSShapley/Combination_Monte/Final/codess/JMP/kstar_vs_T.pdf"
plt.savefig(outpath, bbox_inches="tight", dpi=180)
outpath_png = outpath.replace(".pdf", ".png")
plt.savefig(outpath_png, bbox_inches="tight", dpi=180)
print("Saved:", outpath)
plt.close()
