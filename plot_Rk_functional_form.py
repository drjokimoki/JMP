import csv, math, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

OUTDIR = "/Users/boris/Desktop/RSShapley/Combination_Monte/Final/codess/JMP/"
m = 50

# ── helpers ────────────────────────────────────────────────────────────────

def read_csv(path, kcol="k", rcol="R_bag_inf_mc"):
    rows = []
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            k = int(row[kcol])
            R = float(row[rcol])
            if R > 0:
                rows.append((k, R))
    return rows

def ols_alpha(data, m=50, kmin=2):
    """Global log-log OLS: log R ~ alpha * log z_k."""
    pts = [(k, R) for k, R in data if k >= kmin and R > 0 and (1/k - 1/m) > 0]
    if len(pts) < 3:
        return float("nan"), float("nan")
    lz = [math.log(1/k - 1/m) for k, R in pts]
    lr = [math.log(R)          for k, R in pts]
    n = len(lz); xb = sum(lz)/n; yb = sum(lr)/n
    sxx = sum((x-xb)**2 for x in lz)
    sxy = sum((x-xb)*(y-yb) for x,y in zip(lz,lr))
    slope = sxy/sxx
    res   = [lr[i] - (yb + slope*(lz[i]-xb)) for i in range(n)]
    ss_res = sum(r**2 for r in res)
    ss_tot = sum((y-yb)**2 for y in lr)
    r2 = 1 - ss_res/ss_tot if ss_tot > 0 else float("nan")
    return slope, r2

def normalize(data):
    """Normalize so R_k / R at smallest k = 1."""
    R0 = data[0][1]
    return [(k, R/R0) for k, R in data]

# ── load data ──────────────────────────────────────────────────────────────

base = "/Users/boris/Desktop/RSShapley/Combination_Monte/Final/codess/JMP/"

exo_files = {
    "Exo B3 H=1.05": base+"RSM_exogenous_results/joint_BH_m50_T240_B3_H1p05_population_weight_distance_by_k.csv",
    "Exo B3 H=6":    base+"RSM_exogenous_results/joint_BH_m50_T240_B3_H6_population_weight_distance_by_k.csv",
    "Exo B8 H=1.05": base+"RSM_exogenous_results/joint_BH_m50_T240_B8_H1p05_population_weight_distance_by_k.csv",
    "Exo B8 H=6":    base+"RSM_exogenous_results/joint_BH_m50_T240_B8_H6_population_weight_distance_by_k.csv",
}
endo_files = {
    "Toeplitz B2": base+"RSM_endogenous_results/endo_toeplitz_m_50_T_200_R2_0.3_rho_0.8_beta_B2_phif_0.3_phie_0.7_pop_mc_general.csv",
    "Toeplitz B4": base+"RSM_endogenous_results/endo_toeplitz_m_50_T_200_R2_0.3_rho_0.8_beta_B4_phif_0.3_phie_0.7_pop_mc_general.csv",
    "Factor B2":   base+"RSM_endogenous_results/endo_factor_m_50_T_200_R2_0.3_r_2_load_random_fs_0.5_id_0.25_beta_B2_phif_0.3_phie_0.7_pop_mc_general.csv",
    "Factor B4":   base+"RSM_endogenous_results/endo_factor_m_50_T_200_R2_0.3_r_2_load_random_fs_0.5_id_0.25_beta_B4_phif_0.3_phie_0.7_pop_mc_general.csv",
}

exo_data  = {lbl: read_csv(f) for lbl, f in exo_files.items()}
endo_data = {lbl: read_csv(f) for lbl, f in endo_files.items()}

# ── colour / style maps ────────────────────────────────────────────────────
exo_styles = {
    "Exo B3 H=1.05": dict(color="#1f77b4", ls="--",  marker="o", ms=4, label="B=3, H=1.05"),
    "Exo B3 H=6":    dict(color="#1f77b4", ls="-",   marker="o", ms=4, label="B=3, H=6"),
    "Exo B8 H=1.05": dict(color="#ff7f0e", ls="--",  marker="s", ms=4, label="B=8, H=1.05"),
    "Exo B8 H=6":    dict(color="#ff7f0e", ls="-",   marker="s", ms=4, label="B=8, H=6"),
}
endo_styles = {
    "Toeplitz B2": dict(color="#2ca02c", ls="--", marker="o", ms=5, label="Toeplitz β=B2"),
    "Toeplitz B4": dict(color="#2ca02c", ls="-",  marker="o", ms=5, label="Toeplitz β=B4"),
    "Factor B2":   dict(color="#d62728", ls="--", marker="s", ms=5, label="Factor β=B2"),
    "Factor B4":   dict(color="#d62728", ls="-",  marker="s", ms=5, label="Factor β=B4"),
}

# ══════════════════════════════════════════════════════════════════════════
# FIGURE 1 — Level plots: R_k vs k  (linear and log-y)
# ══════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(2, 2, figsize=(13, 9))
fig.suptitle("Excess risk $R_k = Q_k^{\\mathrm{bag},\\infty} - Q^{\\mathrm{opt}}$ vs subset size $k$",
             fontsize=14, fontweight="bold")

titles = ["Exogenous — level", "Exogenous — log scale",
          "Endogenous — level", "Endogenous — log scale"]

for col, (data_dict, styles) in enumerate([(exo_data, exo_styles),
                                            (endo_data, endo_styles)]):
    for row, logy in enumerate([False, True]):
        ax = axes[row][col]
        for lbl, d in data_dict.items():
            ks = [k for k, R in d]
            Rs = [R for k, R in d]
            s  = styles[lbl]
            ax.plot(ks, Rs, color=s["color"], ls=s["ls"],
                    marker=s["marker"], ms=s["ms"],
                    markevery=max(1, len(ks)//12),
                    label=s["label"], lw=1.6, alpha=0.85)
        if logy:
            ax.set_yscale("log")
            ax.set_ylabel("$R_k$ (log scale)", fontsize=10)
        else:
            ax.set_ylabel("$R_k$", fontsize=10)
        ax.set_xlabel("$k$", fontsize=10)
        ax.set_title(titles[row*2 + col], fontsize=11)
        ax.legend(fontsize=8, framealpha=0.7)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(left=1)

plt.tight_layout()
fig.savefig(OUTDIR + "Rk_level_plots.pdf", bbox_inches="tight")
fig.savefig(OUTDIR + "Rk_level_plots.png", dpi=150, bbox_inches="tight")
plt.close()
print("Saved: Rk_level_plots.pdf/.png")


# ══════════════════════════════════════════════════════════════════════════
# FIGURE 2 — Log-log: log R_k vs log z_k  (power-law check)
# ══════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
fig.suptitle("Log-log plot: $\\log R_k$ vs $\\log z_k$,  $z_k = 1/k - 1/m$\n"
             "(power-law $R_k \\propto z_k^{\\alpha}$ appears as straight line)",
             fontsize=13, fontweight="bold")

for col, (data_dict, styles, title) in enumerate([
        (exo_data,  exo_styles,  "Exogenous (m=50, T=240)"),
        (endo_data, endo_styles, "Endogenous (m=50, T=200)"),
]):
    ax = axes[col]
    for lbl, d in data_dict.items():
        pts = [(k, R) for k, R in d if k >= 2 and R > 0 and (1/k - 1/m) > 0]
        lz = [math.log(1/k - 1/m) for k, R in pts]
        lr = [math.log(R)          for k, R in pts]
        alpha, r2 = ols_alpha(d, m=m)
        s = styles[lbl]
        leg = f"{s['label']}  (α̂={alpha:.2f}, R²={r2:.2f})"
        ax.plot(lz, lr, color=s["color"], ls=s["ls"],
                marker=s["marker"], ms=3.5,
                markevery=max(1, len(lz)//10),
                label=leg, lw=1.5, alpha=0.85)

    # Reference slopes
    lz_ref = np.linspace(-7, -0.2, 80)
    for alpha_ref, col_ref, lab_ref in [(1.0, "grey", "slope 1"), (2.0, "black", "slope 2")]:
        lr_ref = alpha_ref * (lz_ref - lz_ref[-1]) + np.log(d[-1][1] if d else 1e-4)
        # anchor at z_k = 0.8 (k≈1.25), choose intercept for visibility
        intercept = -2.5 - alpha_ref * (-0.2)
        ax.plot(lz_ref, alpha_ref * lz_ref + intercept,
                color=col_ref, ls=":", lw=1.2, alpha=0.6, label=lab_ref)

    ax.set_xlabel("$\\log z_k = \\log(1/k - 1/m)$", fontsize=10)
    ax.set_ylabel("$\\log R_k$", fontsize=10)
    ax.set_title(title, fontsize=11)
    ax.legend(fontsize=7.5, framealpha=0.75, loc="upper left")
    ax.grid(True, alpha=0.3)

plt.tight_layout()
fig.savefig(OUTDIR + "Rk_loglog_zk.pdf", bbox_inches="tight")
fig.savefig(OUTDIR + "Rk_loglog_zk.png", dpi=150, bbox_inches="tight")
plt.close()
print("Saved: Rk_loglog_zk.pdf/.png")


# ══════════════════════════════════════════════════════════════════════════
# FIGURE 3 — Normalised R_k / R_{k0}  (shape comparison, common scale)
# ══════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
fig.suptitle("Normalised excess risk $R_k / R_{k_0}$ vs $k/m$\n"
             "(all curves anchored at 1 at their first observed $k$)",
             fontsize=13, fontweight="bold")

for col, (data_dict, styles, title) in enumerate([
        (exo_data,  exo_styles,  "Exogenous"),
        (endo_data, endo_styles, "Endogenous"),
]):
    ax = axes[col]
    for lbl, d in data_dict.items():
        nd = normalize(d)
        ks_norm = [k/m for k, R in nd]
        Rs_norm = [R   for k, R in nd]
        s = styles[lbl]
        ax.plot(ks_norm, Rs_norm, color=s["color"], ls=s["ls"],
                marker=s["marker"], ms=4,
                markevery=max(1, len(ks_norm)//10),
                label=s["label"], lw=1.8, alpha=0.85)
    ax.set_xlabel("$k/m$", fontsize=10)
    ax.set_ylabel("$R_k / R_{k_0}$", fontsize=10)
    ax.set_title(title, fontsize=11)
    ax.legend(fontsize=9, framealpha=0.75)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 1)
    ax.set_ylim(bottom=0)

plt.tight_layout()
fig.savefig(OUTDIR + "Rk_normalised.pdf", bbox_inches="tight")
fig.savefig(OUTDIR + "Rk_normalised.png", dpi=150, bbox_inches="tight")
plt.close()
print("Saved: Rk_normalised.pdf/.png")


# ══════════════════════════════════════════════════════════════════════════
# FIGURE 4 — Local slope alpha(k) = -d log R / d log z_k  vs k
# ══════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
fig.suptitle("Local log-log slope $\\hat{\\alpha}(k)$ as a function of $k$\n"
             "$\\hat\\alpha(k) = [\\log R_k - \\log R_{k+1}] / [\\log z_k - \\log z_{k+1}]$",
             fontsize=13, fontweight="bold")

for col, (data_dict, styles, title) in enumerate([
        (exo_data,  exo_styles,  "Exogenous"),
        (endo_data, endo_styles, "Endogenous"),
]):
    ax = axes[col]
    for lbl, d in data_dict.items():
        pts = [(k, R) for k, R in d if R > 0 and (1/k - 1/m) > 0]
        alphas = []
        ks_mid = []
        for i in range(len(pts)-1):
            k0, R0 = pts[i];  k1, R1 = pts[i+1]
            z0 = 1/k0 - 1/m; z1 = 1/k1 - 1/m
            if R0 > 0 and R1 > 0 and z0 > 0 and z1 > 0:
                a = (math.log(R0) - math.log(R1)) / (math.log(z0) - math.log(z1))
                alphas.append(a)
                ks_mid.append((k0+k1)/2)
        s = styles[lbl]
        ax.plot(ks_mid, alphas, color=s["color"], ls=s["ls"],
                marker=s["marker"], ms=3.5,
                label=s["label"], lw=1.4, alpha=0.8)

    ax.axhline(1.0, color="grey",  ls=":", lw=1.2, label="α = 1")
    ax.axhline(2.0, color="black", ls=":", lw=1.2, label="α = 2")
    ax.set_ylim(-4, 8)
    ax.set_xlabel("$k$", fontsize=10)
    ax.set_ylabel("Local $\\hat{\\alpha}(k)$", fontsize=10)
    ax.set_title(title, fontsize=11)
    ax.legend(fontsize=8, framealpha=0.75)
    ax.grid(True, alpha=0.3)

plt.tight_layout()
fig.savefig(OUTDIR + "Rk_local_alpha.pdf", bbox_inches="tight")
fig.savefig(OUTDIR + "Rk_local_alpha.png", dpi=150, bbox_inches="tight")
plt.close()
print("Saved: Rk_local_alpha.pdf/.png")


# ══════════════════════════════════════════════════════════════════════════
# FIGURE 5 — Marginal gain (R_k - R_{k+1}) / R_k  vs k
# ══════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
fig.suptitle("Relative marginal gain $(R_k - R_{k+1})/R_k$ vs $k$\n"
             "(the FOC for optimal $k^*$ equates this to $1/(T-k)$)",
             fontsize=13, fontweight="bold")

for col, (data_dict, styles, title, T) in enumerate([
        (exo_data,  exo_styles,  "Exogenous (T=240)",  240),
        (endo_data, endo_styles, "Endogenous (T=200)", 200),
]):
    ax = axes[col]
    for lbl, d in data_dict.items():
        pts = [(k, R) for k, R in d if R > 0]
        drops = []
        ks    = []
        for i in range(len(pts)-1):
            k0, R0 = pts[i]; k1, R1 = pts[i+1]
            step = k1 - k0
            if step > 0 and R0 > 0:
                # normalise per unit k for comparability
                drops.append((R0 - R1) / (R0 * step))
                ks.append(k0)
        s = styles[lbl]
        ax.plot(ks, drops, color=s["color"], ls=s["ls"],
                marker=s["marker"], ms=3.5,
                label=s["label"], lw=1.4, alpha=0.8)

    # 1/(T-k) reference curve
    ks_ref = np.arange(1, m)
    ax.plot(ks_ref, 1/(T - ks_ref), color="black", ls="-", lw=1.8,
            label=f"$1/(T-k)$  [T={T}]")
    ax.set_yscale("log")
    ax.set_xlabel("$k$", fontsize=10)
    ax.set_ylabel("$(R_k - R_{k+1})/R_k$ per unit $k$  (log scale)", fontsize=10)
    ax.set_title(title, fontsize=11)
    ax.legend(fontsize=8, framealpha=0.75)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(1, m-1)

plt.tight_layout()
fig.savefig(OUTDIR + "Rk_marginal_gain.pdf", bbox_inches="tight")
fig.savefig(OUTDIR + "Rk_marginal_gain.png", dpi=150, bbox_inches="tight")
plt.close()
print("Saved: Rk_marginal_gain.pdf/.png")


# ══════════════════════════════════════════════════════════════════════════
# FIGURE 6 — All scenarios, single log-log panel (R_k vs k, not z_k)
# ══════════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(9, 6))
fig.suptitle("$\\log R_k$ vs $\\log k$  — all scenarios\n"
             "(straight line ↔ $R_k \\propto k^{-\\gamma}$)",
             fontsize=13, fontweight="bold")

all_styles = {**exo_styles, **endo_styles}
# merge endo style labels for clarity
endo_labels = {
    "Toeplitz B2": dict(color="#2ca02c", ls="--", marker="o", ms=5, label="Endo: Toeplitz β=B2"),
    "Toeplitz B4": dict(color="#2ca02c", ls="-",  marker="o", ms=5, label="Endo: Toeplitz β=B4"),
    "Factor B2":   dict(color="#d62728", ls="--", marker="s", ms=5, label="Endo: Factor β=B2"),
    "Factor B4":   dict(color="#d62728", ls="-",  marker="s", ms=5, label="Endo: Factor β=B4"),
}
exo_labels = {
    "Exo B3 H=1.05": dict(color="#aec7e8", ls="--", marker="o", ms=4, label="Exo: B3 H=1.05"),
    "Exo B3 H=6":    dict(color="#1f77b4", ls="-",  marker="o", ms=4, label="Exo: B3 H=6"),
    "Exo B8 H=1.05": dict(color="#ffbb78", ls="--", marker="s", ms=4, label="Exo: B8 H=1.05"),
    "Exo B8 H=6":    dict(color="#ff7f0e", ls="-",  marker="s", ms=4, label="Exo: B8 H=6"),
}

for data_dict, styles_map in [(exo_data, exo_labels), (endo_data, endo_labels)]:
    for lbl, d in data_dict.items():
        pts = [(k, R) for k, R in d if k >= 2 and R > 0]
        lk = [math.log(k) for k, R in pts]
        lr = [math.log(R) for k, R in pts]
        s = styles_map[lbl]
        ax.plot(lk, lr, color=s["color"], ls=s["ls"],
                marker=s["marker"], ms=3, markevery=max(1, len(lk)//8),
                label=s["label"], lw=1.5, alpha=0.85)

# Reference slopes anchored at log k = 1
lk_ref = np.linspace(0.5, 3.8, 80)
for gamma, col_ref, lab_ref in [(2.0, "grey", "slope −2"), (3.0, "black", "slope −3")]:
    intercept = -2.5
    ax.plot(lk_ref, -gamma * lk_ref + intercept,
            color=col_ref, ls=":", lw=1.3, label=lab_ref)

ax.set_xlabel("$\\log k$", fontsize=11)
ax.set_ylabel("$\\log R_k$", fontsize=11)
ax.legend(fontsize=7.5, framealpha=0.8, ncol=2, loc="lower left")
ax.grid(True, alpha=0.3)

plt.tight_layout()
fig.savefig(OUTDIR + "Rk_loglog_k.pdf", bbox_inches="tight")
fig.savefig(OUTDIR + "Rk_loglog_k.png", dpi=150, bbox_inches="tight")
plt.close()
print("Saved: Rk_loglog_k.pdf/.png")

print("\nAll figures done.")
