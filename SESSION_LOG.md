# Session Log — RSM Forecast Combination Project
**Date:** 3 June 2026  
**GitHub repo:** https://github.com/drjokimoki/JMP  
**Working directory:** `/Users/boris/Desktop/RSShapley/Combination_Monte/Final/codess/JMP/`

---

## 1. What this project is

A JMP paper on **Random Subset Methods (RSM) for forecast combination** in the common-loading heteroskedastic model. The paper derives:

- Exact and approximate expressions for the **infinite-bagged RSM** population risk $Q_k^{\text{bag},\infty}$ and the **average-subset RSM** risk $Q_k^{\text{sub}}$
- A **small-heterogeneity expansion** giving the endpoint-corrected approximation $Q_k^{\text{bag},\infty} - Q^{\text{opt}} \approx C_\delta(1/k - 1/m)^2$
- **Cube-root optimal subset size** $k^* \propto (T/m)^{1/3}$ for the infinite-bagged case vs **square-root** $k^* \propto (T/A_0)^{1/2}$ for the average-subset case
- **Dominance regions** showing when RSM beats full estimated OLS or equal weights

---

## 2. Files in the repo

| File | Purpose |
|---|---|
| `rsm_common_loading_small_heterogeneity_proposal.pdf` | Main paper + appended 6-page exact-dominance appendix (22 pp total) |
| `rsm_common_loading_grid_weights_general_distance_reviewfix_final.R` | **Exogenous** simulation code (verifies theory) |
| `endogenous_rsm_toeplitz_factor.R` | **Endogenous** simulation code (empirical study) |
| `RSM_endogenous_results/` | Output folder — **currently being written by a live run** |
| `SESSION_LOG.md` | This file |

Other folders (`GDP_constrained/`, `IP_Weighted/`, `Simulation/`, etc.) contain empirical application code for GDP and IP forecasting — not the focus of this session.

---

## 3. Key model notation

**Common-loading model:**  
$\Sigma_e = q\iota_m\iota_m' + D$, $D = \text{diag}(\sigma_1^2,\ldots,\sigma_m^2)$, $\tau_i = \sigma_i^{-2}$

**Two RSM population objects:**
- $Q_k^{\text{sub}} = E_S[v_S'\Sigma_e v_S] = q + E_S[1/A_S]$ — average-subset risk
- $Q_k^{\text{bag},\infty} = \bar v_k'\Sigma_e\bar v_k$ — infinite-bagged risk (Jensen: $\leq Q_k^{\text{sub}}$)

**Key identities (exact, no small-het assumption):**
- $Q(S) = q + 1/A_S$, $A_S = \sum_{i\in S}\tau_i$
- $w_{i,S} = \tau_i/A_S$, $w_i^{\text{opt}} = \tau_i/A_m$, $Q^{\text{opt}} = q + 1/A_m$
- $R_k^{\text{bag},\infty} := Q_k^{\text{bag},\infty} - Q^{\text{opt}} = \sum_i(\bar v_{i,k} - w_i^{\text{opt}})^2/\tau_i \geq 0$ ← **exact for all H**

**Small-heterogeneity ($\tau_i = \bar\tau(1+\delta_i)$, $\bar\delta=0$, small $V_\delta$):**
- $\bar v_{i,k} - w_i^{\text{opt}} \approx -\frac{\delta_i}{m}\cdot\frac{m-k}{k(m-1)}$
- $R_k^{\text{bag},\infty} \approx C_\delta(1/k - 1/m)^2$, $C_\delta = \frac{BV_\delta}{m}\left(\frac{m}{m-1}\right)^2$, $B = 1/\bar\tau$

**Finite-bagging bridge (Proposition 1):**  
$Q_k^{\text{bag},G} = Q_k^{\text{bag},\infty} + \frac{1}{G}(Q_k^{\text{sub}} - Q_k^{\text{bag},\infty})$  
→ $Q_k^{\text{bag},\infty} \leq Q_k^{\text{bag},G} \leq Q_k^{\text{sub}}$

**Inflation factor:** $I(T,k) = (T-1)/(T-k)$ (conservative for bagging, using $d^{\text{bag}}(k) = k$)

---

## 4. Optimal k results

| Object | Optimal $k^*$ | Rule |
|---|---|---|
| Average-subset | $(-B + \sqrt{B^2 + A_0 BT})/A_0 \sim (BT/A_0)^{1/2}$ | **Square-root** |
| Infinite-bagged | $(2C_\delta T/A)^{1/3}$ for $k \ll m$ | **Cube-root** |
| Finite-bagged | Interpolates; approaches square-root as $G \to 1$ | Transition |

**Scaling:** $k^*_{\text{bag},\infty} \propto (T/m)^{1/3}$ when $V_\delta, B, s+q$ stable.

---

## 5. Exact dominance regions (appended to PDF, Section 11)

The main paper (Section 10) derives dominance under small-het. The appendix derives the **exact** conditions for the common-loading model with no further assumption:

| Comparison | Exact condition |
|---|---|
| Avg-sub vs $L^{\text{opt}}$ | $(Q_k^{\text{sub}} - Q^{\text{opt}})(T-m) < A(m-k)$ |
| Avg-sub vs $L^{\text{ave}}$ | $(s + Q_k^{\text{sub}})(T-1) < M_{\text{ave}}(T-k)$ |
| $\infty$-bag vs $L^{\text{opt}}$ | $R_k^{\text{bag},\infty}(T-m) < A(m-k)$ |
| Finite-bag vs $L^{\text{opt}}$ | $(R_k^{\text{bag},\infty} + \Delta_k/G)(T-m) < A(m-k)$ |
| $\infty$-bag vs $L^{\text{ave}}$ | $(A + R_k^{\text{bag},\infty})(T-1) < M_{\text{ave}}(T-k)$ |

where $A = s + Q^{\text{opt}}$, $M_{\text{ave}} = s + Q^{\text{ave}}$, $\Delta_k = Q_k^{\text{sub}} - Q_k^{\text{bag},\infty}$.

Small-het formulas from Section 10 are recovered as corollaries.  
**Exact threshold** $\kappa^{\text{exact}}$: smallest $k$ where $R_k^{\text{bag},\infty}(T-m) < A(m-k)$ — computable directly from MC estimates of $R_k^{\text{bag},\infty}$ (stored in `*_pop_mc_general.csv` as `R_bag_inf_mc`).

**For general factor structure** (endogenous code): identity becomes  
$R_k^{\text{bag},\infty} = (\bar v_k - w^{\text{opt}})'\Sigma_e(\bar v_k - w^{\text{opt}})$ — still exact, computed by `q_population_mc_general()`.

---

## 6. Exogenous simulation code — key functions

**File:** `rsm_common_loading_grid_weights_general_distance_reviewfix_final.R`

```
make_tau_lognormal(m, B_scale, H_target)   # generate precision design
precision_summary(tau, q)                  # B, H, V_delta, C_delta, Q_opt, Q_ave
simulate_common_loading(...)               # DGP: mu_t AR(1), common q, idio D
constrained_ols_weights(y, Fmat)           # sum-to-one OLS (reparameterised)
oracle_weights_common(tau)                 # tau_i / A_m

# Population theory
q_subset_delta_approx(tau, q, k_grid)     # eq (69): B/k + B*(m-k)V_delta/k^2(m-1)
q_bag_inf_smallhet(tau, q, k_grid)        # Q_opt + C_delta*(1/k-1/m)^2
q_bag_G_smallhet(tau, q, k_grid, G_bag)   # adds finite-bagging (B/G)*(1/k-1/m) term
q_population_mc_common(tau, q, k_grid)    # exact MC: Q_sub, Q_bag_inf, weight errors
q_population_mc_general(Sigma_e, k_grid)  # general Sigma_e version

# k* calculations
kstar_subset_leading(...)                 # continuous positive root of quadratic
kstar_bag_cube_root_local(...)            # (2 C_delta T / A)^(1/3)

loss_from_Q(s, Q, T_train, k)             # (s+Q)*(T-1)/(T-k)
```

**Experiments run by default (`RUN_EXPERIMENTS`):**
1. `"baseline"` — MSPE curves for grid of (B_scale, H_target) values
2. `"scaling_T_m"` — verifies cube-root scaling $k^* \propto (T/m)^{1/3}$ via log-log regression
3. `"joint_B_H_grid"` — heatmaps over joint (B, H) grid

**Key parameters (paper profile):**
- `m=50`, `T_train=240`, `s=1`, `q=1`
- `B_scale ∈ {3, 8}`, `H_target ∈ {1.05, 6}`
- `N_REP_BASE=500`, `G_BAG_MC=800`, `N_SUBSETS_POP=3000`

**Fixes applied this session:**
- Added `"baseline"` and `"scaling_T_m"` to `RUN_EXPERIMENTS` (were missing)

---

## 7. Endogenous simulation code — key structure

**File:** `endogenous_rsm_toeplitz_factor.R`

**DGP:**
- $x_t = F_t^{\text{common}}\cdot\iota_m + \eta_t$, $F^{\text{common}} \sim \text{AR}(1)(\phi_f)$
- $\eta_t \sim \text{VAR}(1)(\phi_\eta, \Sigma_e)$ — either Toeplitz or factor structure
- $y_{t+1} = \beta'x_t + u_{t+1}$, $\beta$ scaled to target $R^2$
- Individual forecasters: $\hat y_{t+1|t} = \hat\gamma_i x_{it}$ (no-intercept OLS)
- RSM uses estimated $\hat\Sigma_u$ (forecast error covariance from training residuals)

**Two Sigma_e types:**
- `"toeplitz"`: $\Sigma_e[i,j] = \rho^{|i-j|}$, default `rho=0.8`
- `"factor"`: $\Sigma_e = \Lambda\Omega\Lambda' + D_e$, scaled to unit diagonal

**Beta designs:**
- `"B2"`: $\beta_i = m+1-i$ (decreasing)
- `"B4"` / `"sparse10"`: $\beta_i = 1$ for $i \leq 10$, else 0

**Fixes applied this session:**
1. **Inflation formula**: `1 + (k-1)/T` → `(T-1)/(T-k)` for RSM subsets (matches paper)
2. **Full-OLS inflation**: `Q*(1+(m-1)/T)` → `Q*(T-1)/(T-m)` (matches paper)
3. **Population MC added**: `q_population_mc_general(Sigma_e, z_grid)` called per scenario — computes true $Q_k^{\text{bag},\infty}$ on true $\Sigma_e$, exact dominance threshold $\kappa^{\text{exact}}$, saves `*_pop_mc_general.csv`
4. **T-scaling experiment added**: `run_scaling_T_endo()` — fixes $\Sigma_e$, varies $T$

**Profiles:**
- `"fast"`: 10 reps, quick logic check
- `"run"`: 200 reps, 200 bags, OOS=30 — current live run
- `"paper"`: 500 reps, 500 bags, OOS=50

---

## 8. Live simulation — currently running

**PID:** 66691 (as of session end — may have changed)  
**Started:** ~21:25 on 3 June 2026  
**Profile:** `"run"` (N_REP=200, K_RS=200, OOS_LEN=30)  
**Output:** `RSM_endogenous_results/`  
**Status at session end:** Scenario 1 of 4 (toeplitz, B2), replication 40/200 at elapsed 5:01

**Expected completion:** ~22:40–22:50 (≈75–85 min total from start)

**To check status:**
```bash
cd /Users/boris/Desktop/RSShapley/Combination_Monte/Final/codess/JMP
tail -20 RSM_endogenous_results/run_log.txt
ls -lh RSM_endogenous_results/
ps aux | grep "[e]ndogenous_rsm" | awk '{print "CPU:",$3"% | elapsed:",$10}'
```

**4 main scenarios:**
1. toeplitz, rho=0.8, beta=B2 — **running now**
2. toeplitz, rho=0.8, beta=B4
3. factor (r=2, random loadings, strength=0.5, idio=0.25), beta=B2
4. factor, beta=B4

**T-scaling experiment** (runs after scenario 4):
- Fixes $\Sigma_e$ with seed, varies $T \in \{100, 150, 200, 300, 500, 800\}$
- For both toeplitz (rho=0.8) and factor $\Sigma_e$
- 120 reps, 150 bags per T value
- Uses fast bag loop (`rsm_fast_one_origin`, no `eigen()` per bag)

---

## 9. Output files expected in RSM_endogenous_results/

**Per scenario (4 scenarios):**
| File pattern | Content |
|---|---|
| `*_rmsfe_by_z.png` | RMSFE vs subset size k — core empirical result |
| `*_Q_by_z.png` | Q diagnostics: bagged / subset-avg / full OLS |
| `*_L_by_z.png` | Inflation-adjusted L curves |
| `*_pop_vs_empirical_Q.png` | True population $Q_k^{\text{bag},\infty}$ vs estimated |
| `*_by_z.csv` | Full per-k table (RMSFE, Q, L, bias, variance) |
| `*_pop_mc_general.csv` | Population MC: $Q_k^{\text{sub}}$, $Q_k^{\text{bag},\infty}$, $R_k^{\text{bag},\infty}$, $\kappa^{\text{exact}}$ |

**Summary:**
| File | Content |
|---|---|
| `endogenous_rsm_summary.csv` | One row per scenario, all diagnostics |

**T-scaling:**
| File | Content |
|---|---|
| `scaling_T_rmsfe_curves_toeplitz.png` | All RMSFE curves overlaid by T |
| `scaling_T_rmsfe_curves_factor.png` | Same for factor $\Sigma_e$ |
| `scaling_T_kstar_toeplitz.png` | $k^*$ vs $T^{1/3}$ + log-log panel |
| `scaling_T_kstar_factor.png` | Same for factor |
| `scaling_T_summary.csv` | $k^*$ and relative RMSFE for each $(T, \Sigma_e)$ |
| `scaling_T_loglog_regression.txt` | OLS slope (theory predicts ≈ 1/3) |
| `scaling_*.csv` | Per-$(T, \Sigma_e)$ RMSFE curves |

---

## 10. Key finding already visible (from population MC, scenario 1)

**Toeplitz, rho=0.8, B2 beta, m=50:**

| k | $Q_k^{\text{bag},\infty}$ | $R_k^{\text{bag},\infty}$ | $\|\bar v_k - w^*\|_2$ |
|---|---|---|---|
| 2 | 0.1645 | 0.0093 | 0.098 |
| 11 | 0.1592 | 0.0041 | 0.081 |
| 20 | 0.1572 | 0.0020 | 0.069 |
| 29 | 0.1560 | 0.0008 | 0.047 |
| 50 | 0.1552 | 0.0000 | 0.000 |

$Q^{\text{opt}} = 0.1552$, $Q^{\text{ave}} = 0.164$  
$\kappa^{\text{exact}} = 2$ — even a 2-predictor bagged subset beats full estimated OLS!

This makes intuitive sense for rho=0.8: high cross-correlation means subsets already capture most of the benefit; the inflation penalty on the full 50-predictor OLS is severe.

---

## 11. Things to do in the next session

1. **Read and interpret results** from `RSM_endogenous_results/` once simulation completes
2. **Compare kappa_exact** across sigma types and beta designs
3. **Check T-scaling slope** — does $k^* \propto T^{1/3}$ hold empirically in the endogenous case?
4. **Consider running the exogenous code** to verify the theoretical derivations (add `"baseline"` and `"scaling_T_m"` already in `RUN_EXPERIMENTS`)
5. **Paper writing**: the exact dominance section is in the PDF appendix — decide whether to integrate it as a main section
6. **Dominance without small-het**: for the general factor structure, dominance conditions are in the appendix (Section 11, eqs 168–175) but only computable numerically via `R_bag_inf_mc` from `*_pop_mc_general.csv`

---

## 12. Git history (this session)

```
a0506d0  Add run profile, T-scaling experiment, q_population_mc_general
c4e4bba  Add exact dominance-region appendix and three simulation fixes
ea5eee9  Initial commit: JMP project files
```

All pushed to: https://github.com/drjokimoki/JMP

---

## 13. R package requirements

The scripts use **base R only** — no external packages. Tested on R ≥ 4.2.

LaTeX tools used for PDF: TeX Live 2024 (`/usr/local/texlive/2024/bin/universal-darwin/pdflatex`), poppler (`pdftotext`, `pdfunite`).
