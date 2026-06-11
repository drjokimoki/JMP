"""
k* vs T — clean two-panel figure.
Left panel: linear scale, endogenous cases only with reference curves.
Right panel: log-log scale, same data + fitted slope annotations.
Exogenous H=6 shown only as formula prediction (no data pending).
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# ── Data ─────────────────────────────────────────────────────────────────────
T_endo     = np.array([100, 150, 200, 300, 500, 800])
kstar_toep = np.array([10,  10,  14,  14,  18,  22])
kstar_fact = np.array([14,  14,  14,  18,  26,  26])

# Fitted free-OLS slopes and intercepts for endogenous
slope_toep, C_toep = 0.394, np.exp(0.444)   # k* ~ 1.559 * T^{0.394}
slope_fact, C_fact = 0.367, np.exp(0.836)   # k* ~ 2.307 * T^{0.367}

# Exogenous H=6: only power-law formula prediction (no MC data yet for T>80)
# From joint_B_H: alpha_hat=1.72, C_delta=0.1006, A=2.06
alpha_H6 = 1.72; C_delta_H6 = 0.1006; A_H6 = 2.06
# Exo H=6 single data point at T=80 (population grid k*)
T_exo6_data = np.array([80])
k_exo6_data = np.array([3])

# ── Reference curves ─────────────────────────────────────────────────────────
T_ref = np.linspace(60, 900, 500)

sqrt_ref  = np.sqrt(T_ref)
cube_ref  = C_toep * T_ref**(1/3)           # C≈2.19 at slope=1/3 constrained
fit_toep  = C_toep * T_ref**slope_toep
fit_fact  = C_fact * T_ref**slope_fact
exo_pred  = (C_delta_H6 * alpha_H6 * T_ref / A_H6)**(1/(alpha_H6+1))

# ── Palette ──────────────────────────────────────────────────────────────────
COL_TOEP  = "#0072B2"
COL_FACT  = "#D55E00"
COL_SQRT  = "#009E73"
COL_CUBE  = "#CC79A7"
COL_EXO6  = "#E69F00"

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 10.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 180,
})

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

for ax, log_scale in zip(axes, [False, True]):
    if log_scale:
        ax.set_xscale("log"); ax.set_yscale("log")

    # ── Reference lines ──────────────────────────────────────────────────
    ax.plot(T_ref, sqrt_ref, "--", color=COL_SQRT, lw=1.5, alpha=0.85,
            label=r"$\sqrt{T}$ reference", zorder=1)
    ax.plot(T_ref, cube_ref,  ":", color=COL_CUBE, lw=1.5, alpha=0.85,
            label=r"$2.19\cdot T^{1/3}$ reference", zorder=1)

    # ── Exogenous H=6 formula (dashed, no data) ──────────────────────────
    ax.plot(T_ref, exo_pred, "-.", color=COL_EXO6, lw=1.2, alpha=0.75,
            label=r"Exo.\ $H=6$ power-law formula ($\hat\alpha=1.72$)", zorder=2)
    ax.plot(T_exo6_data, k_exo6_data, "D", color=COL_EXO6, ms=7, zorder=4,
            label=r"Exo.\ $H=6$, $T=80$ (only point available)")

    # ── Endogenous fitted lines (thin, behind) ────────────────────────────
    ax.plot(T_ref, fit_toep, "-",  color=COL_TOEP, lw=0.9, alpha=0.35, zorder=2)
    ax.plot(T_ref, fit_fact, "-",  color=COL_FACT, lw=0.9, alpha=0.35, zorder=2)

    # ── Endogenous MC data ────────────────────────────────────────────────
    ax.plot(T_endo, kstar_toep, "o-", color=COL_TOEP, lw=2.0, ms=6.5,
            label=r"Endo.\ Toeplitz $k^*$ (RMSFE, 500 reps)", zorder=5)
    ax.plot(T_endo, kstar_fact, "s-", color=COL_FACT, lw=2.0, ms=6.5,
            label=r"Endo.\ Factor $k^*$ (RMSFE, 500 reps)", zorder=5)

    # ── Slope labels (log-log panel only) ────────────────────────────────
    if log_scale:
        ax.annotate(r"slope $0.39$",
                    xy=(400, C_toep*400**slope_toep),
                    xytext=(200, 20), fontsize=9, color=COL_TOEP,
                    arrowprops=dict(arrowstyle="->", color=COL_TOEP, lw=0.8))
        ax.annotate(r"slope $0.37$",
                    xy=(400, C_fact*400**slope_fact),
                    xytext=(500, 25), fontsize=9, color=COL_FACT,
                    arrowprops=dict(arrowstyle="->", color=COL_FACT, lw=0.8))
        ax.annotate(r"slope = 1/2",
                    xy=(700, np.sqrt(700)),
                    xytext=(500, 32), fontsize=9, color=COL_SQRT,
                    arrowprops=dict(arrowstyle="->", color=COL_SQRT, lw=0.8))
        ax.annotate(r"slope = 1/3",
                    xy=(500, C_toep*500**(1/3)),
                    xytext=(200, 12), fontsize=9, color=COL_CUBE,
                    arrowprops=dict(arrowstyle="->", color=COL_CUBE, lw=0.8))

    # ── Axes ──────────────────────────────────────────────────────────────
    ax.set_xlabel(r"Sample size $T$", fontsize=12)
    ax.set_ylabel(r"Optimal $k^*$", fontsize=12)
    ax.set_title("Linear scale" if not log_scale else "Log--log scale",
                 fontsize=11, pad=6)
    ax.set_xlim(55, 920)
    if not log_scale:
        ax.set_ylim(0, 34)
        # Shaded region: sqrt(T) ± 20%
        ax.fill_between(T_ref, 0.8*sqrt_ref, 1.2*sqrt_ref,
                        color=COL_SQRT, alpha=0.06, label=r"$\pm20\%$ of $\sqrt{T}$")
    else:
        ax.set_ylim(1.5, 45)
    ax.grid(True, which="both" if log_scale else "major",
            ls=":", lw=0.4, alpha=0.45)

# ── Legend below both panels ──────────────────────────────────────────────────
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="lower center", ncol=3,
           fontsize=9, frameon=True, framealpha=0.92,
           bbox_to_anchor=(0.5, -0.22))

fig.suptitle(
    r"Optimal subset size $k^*$ vs sample size $T$" "\n"
    r"Endogenous cases: RMSFE minimiser ($m=50$, 500 MC reps); "
    r"exogenous $H=6$: formula only (data pending)",
    fontsize=11, y=1.02)

fig.tight_layout()

for ext in ("pdf", "png"):
    path = f"/Users/boris/Desktop/RSShapley/Combination_Monte/Final/codess/JMP/kstar_vs_T_v2.{ext}"
    plt.savefig(path, bbox_inches="tight", dpi=180)
    print("Saved:", path)

plt.close()
