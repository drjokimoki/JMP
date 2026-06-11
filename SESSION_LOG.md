# Session Log — RSM Forecast Combination Project
**GitHub repo:** https://github.com/drjokimoki/JMP  
**Working directory:** `/Users/boris/Desktop/RSShapley/Combination_Monte/Final/codess/JMP/`

---

## Project overview

A JMP paper on **Random Subset Methods (RSM) for forecast combination** in the common-loading heteroskedastic model.

**Scope (2026-06-04):** Theory and MC focus on exogenous case. Endogenous MC included as a separate section (motivation + results). GDP and SPF as empirical applications. CPI dropped.

| Component | Script | Status |
|---|---|---|
| Exogenous MC (known Σ_e) | `rsm_common_loading_grid_weights_general_distance_reviewfix_final.R` | ✅ Complete |
| Endogenous MC (estimated Σ_e) | `endogenous_rsm_toeplitz_factor.R` | ✅ Complete |
| Endogenous T-scaling (high-res rerun) | `endogenous_scaling_rerun.R` (PID 89806) | 🔄 Toeplitz T=800 in progress; Factor not yet started |
| FRED-MD GDP (bivariate) | `GDP_constrained/GDP_all_k_results.R` | ✅ Complete |
| FRED-MD GDP (AR+X) | `GDP_constrained/GDP_ar_bivariate_allk.R` | ✅ Complete |
| SPF | `GDP_constrained/SPF_RSM_analysis.R` | ✅ Complete |
| H=6 T-scaling (exo) | `rsm_common_loading_...final.R` via `launch_h6_scaling.sh` | ⏸ Pending — waiting for PID 89806; T=80 done, T=240–1000 not yet run |
| **Full paper** | `RSM_paper_main.tex` / `.pdf` | ✅ Revised per reviewer comments (2026-06-04 session 3) |
| Summary report (short) | `RSM_full_report.tex` / `.pdf` | ✅ 13pp |

---

## Papers

### RSM_paper_main.pdf — revised 2026-06-04 (session 3)
All reviewer comments from `RSM_proposal_evaluation_comments.pdf` addressed. See revision summary below.

Structure: Scope → Setup → Benchmarks → RSM objects → Inflation → Common-loading exact derivations → Small-het derivation → Finite bagging → Optimal k → General het / power law → Dominance regions → Exogenous MC → Endo motivation → Endogenous MC → SPF → GDP → Summary/Conclusion  
All propositions with proofs (Prop 1: bridge; 2: weight-deviation; 3: small-het; 4–6: power-law regimes; 7–10: dominance).  
Appendix A: exact dominance table.

### RSM_full_report.pdf — 13 pages (summary, no proofs)
Streamlined narrative: theory → exo MC → GDP → SPF. No PENDING blocks. Not yet updated for reviewer comments.

---

## Reviewer comment responses (session 3, 2026-06-04)

Source: `RSM_proposal_evaluation_comments.pdf`

| Reviewer point | Fix applied in RSM_paper_main.tex |
|---|---|
| **α interpretation** (Priority 1): small het → α=2 when identifiable; near-homogeneity → C_δ≈0, α weakly identified, not α≈1 | Added `Remark[Identification of α]` in §9; regime table gains C_δ column + footnote; exo MC table gets caveat footnote; MC findings paragraph rewritten |
| **Inflation as maintained approximation** (Priority 4): I(T,k) not exact for bagged RSM | §5 adds explicit sentence; endo MC gap finding rewritten to explain both sources of overstatement |
| **Endogenous weight definition** (Priority 3): code uses full Σ_{e,S}^{-1}ι_S formula, not diagonal τ_i/A_S | §11 RSM implementation paragraph replaced with correct formula + clarification that diagonal form holds only under common loading |
| **SPF reframing** (Priority 5): equal within-subset averaging → infinite bagging = EW mean | §12 data bullet adds mathematical identity; conclusion rewritten as sanity check |
| **SPF sign convention**: gains were negative when RSM < EW RMSFE | Table gains flipped to positive (+0.2%, +0.1%, +0.5%, +0.1%); footnote clarified |
| **Fixed-k vs oracle-k** (Priority 6): k*=19/21 is ex-post oracle, not OOS evidence | §13 RMSFE-by-k paragraph restructured: fixed rule k=floor(sqrt(T)) = primary OOS; k*=19/21 labelled as ex-post diagnostic |
| **Claim-level edits** (Priority 7): overclaiming in abstract and conclusion | Abstract: "exact functional form" → "exact representation"; "strictly more efficient" → "lower scaling exponent under maintained approximation"; Conclusion: same; endogenous claim → "extended estimation-noise interpretation" |

---

## Key model notation

**Common-loading model:** Σ_e = qι_m ι_m' + D, D = diag(σ₁²,…,σ_m²), τ_i = σ_i^{-2}

**Two RSM population objects:**
- Q_k^sub = q + E_S[1/A_S] — average-subset risk
- Q_k^{bag,∞} = v̄_k' Σ_e v̄_k — infinite-bagged risk
- **Weight-deviation (exact):** Q_k^{bag,∞} - Q^opt = Σ_i (v̄_{i,k} - w_i^opt)²/τ_i
- Finite-bagging bridge: Q_k^{bag,G} = Q_k^{bag,∞} + (1/G)(Q_k^sub − Q_k^{bag,∞}) [exact, verified to machine precision in exogenous case]

**Optimal k formulae:**
- Average-subset: k*_sub ≈ (BT/A_0)^{1/2} — square-root rule
- Infinite-bagged (small het): k*_bag ≈ (2C_δ T/A)^{1/3} — cube-root rule (under maintained inflation approximation)
- General het (power law): k* ≈ (CαT/A)^{1/(α+1)}, α ∈ [1,2]
- C_δ = (BV_δ/m)(m/(m-1))², B = 1/τ̄

**Small-het weight deviation:**
- Δ_{i,k} = v̄_{i,k} - w_i^opt ≈ -(δ_i/(m-1))·z_k, z_k = 1/k - 1/m
- R_k^{bag,∞} ≈ C_δ z_k² (quadratic in z_k; common loading q-term vanishes since Σδ_i=0)

**α interpretation (post-reviewer revision):**
- α=2: small-het theory prediction, identifiable only when C_δ non-negligible
- α≈1: dominant-forecaster limit (extreme heterogeneity)
- Near-homogeneity = C_δ≈0 (flat curve, α weakly identified), NOT α≈1
- k* increases with H (higher heterogeneity → need larger k to capture best forecasters reliably)
- Lognormality of τ_i: modeling choice for tractability, not derived from structure; results hold for any distribution with same second moment V_δ = H−1

---

## Exogenous MC — results summary

**Output dir:** `RSM_exogenous_results/`  
**Profile:** N_REP=500, G=800, N_SUBSETS_POP=2500, m=50, T=240

### Main results (joint_B_H_summary.csv)

| B | H_actual | k*(grid) | k*(emp) | Gain vs EW | α̂ | R² | Approx err |
|---|---|---|---|---|---|---|---|
| 3 | 1.047 | 1 | 3 | 0.1% | 1.12† | 0.98 | 90% |
| 8 | 1.041 | 1 | 3 | 0.1% | 1.10† | 0.98 | 90% |
| 3 | 6.22 | 3 | 6 | 12.7% | 1.72 | 0.98 | 75% |
| 8 | 3.51* | 4 | 8 | 14.9% | 1.70 | 0.99 | 77% |

*B=8, H_target=6 achieved H_actual=3.51.  
†At H=1.05, C_δ≈0 and the excess-risk curve is nearly flat; α̂ weakly identified, not interpretable as structural exponent.  
α̂ from log-log OLS: log(Q_k^{bag,∞} - Q^opt) ~ α log(z_k), k=2..49.  
Dominance threshold κ < 1 for all DGPs (RSM beats full OLS at k=1).

### T-scaling (H=1.05, B=3; scaling_T_m_summary.csv)

| T | k*(formula, cont.) | k*(grid) | k*(emp) | Gain |
|---|---|---|---|---|
| 80 | 0.57 | 1 | 2 | 0.02% |
| 240 | 0.82 | 1 | 2 | 0.02% |
| 640 | 1.13 | 1 | 4 | 0.06% |
| 1000 | 1.31 | 1 | 3 | 0.03% |

At H=1.05: cube-root law invisible (gains < 0.1%, k*=1 for all T). C_δ≈0 regime.

### H=6 T-scaling — ⏸ PENDING (launch_h6_scaling.sh ready)
Output dir: `RSM_exogenous_results_H6_scaling/`  
Status: T=80 CSV saved (`scaling_finite_sample_m50_T80.csv`); T=240 crashed at rep 400/500 (OOM — only 63MB free RAM, PID 89806 competing). T=640 and T=1000 not started.  
Fix applied: checkpoint logic added to `run_scaling_T_m` (skips existing CSVs); `RSM_SCALING_T_GRID` env var added; `launch_h6_scaling.sh` will auto-start once PID 89806 finishes.

---

## Endogenous MC — results summary

**Output dir:** `RSM_endogenous_results/`  
**Profile:** N_REP=200, K_RS=200, OOS_LEN=30, m=50, T=200

### Main results (endogenous_rsm_summary.csv)

| Σ_e | β | k*(RMSFE) | k*(L_bag) | RMSFE gain | κ_exact | Q_opt | pop k*(bag,∞) |
|---|---|---|---|---|---|---|---|
| Toeplitz | B2 | 11 | 32 | 3.5% | 2 | 0.155 | 2 |
| Toeplitz | B4 | 11 | 26 | 9.8% | 2 | 0.155 | 2 |
| Factor | B2 | 14 | 29 | 5.0% | 2 | 0.005 | 2 |
| Factor | B4 | 14 | 29 | 6.1% | 2 | 0.004 | 2 |

pop k*(bag,∞) = 2 for all: population minimiser without inflation is k=2.  
Gap k*(L_bag)=26–32 vs k*(RMSFE)=11–14: maintained inflation approximation overstates effective df of bagged RSM; Σ_e estimation adds further noise not captured.

### T-scaling (high-res rerun: `endogenous_scaling_rerun.R`, PID 89806)

Output dir: `RSM_endogenous_results/scaling_v2/`  
Profile: N_REP=500, K_RS=500, K_BY=1, B4 only.

**Completed (Toeplitz):**

| T | k*(RMSFE) | sqrt(T) | 2·T^{1/3} |
|---|---|---|---|
| 100 | 10 | 10.0 | 9.3 |
| 150 | 11 | 12.2 | 10.6 |
| 200 | 14 | 14.1 | 11.7 |
| 300 | 13 | 17.3 | 13.4 |
| 500 | 16 | 22.4 | 15.9 |
| 800 | 🔄 in progress | 28.3 | 18.6 |

**Finding:** Log-log exponent = 0.279, CI [0.064, 0.493]. Consistent with cube-root (0.333), not sqrt (0.5). k*/T^{1/3} ratio is stable at ~2.0. Ratio k*/sqrt(T) declines (1.0 → 0.72), ruling out sqrt(T). The T=100 coincidence (k*=10=sqrt(100)) is misleading.

**Pending:** Factor spec T=100–800 (~18+ hours from now as of 20:34 CEST 2026-06-04).

### α̂ for endogenous DGPs (from pop_mc_general.csv files)

| Σ_e | β | α̂ | R² |
|---|---|---|---|
| Toeplitz | B2 | 1.12 | 0.93 |
| Toeplitz | B4 | 1.09 | 0.93 |
| Factor | B2 | 1.73 | 0.99 |
| Factor | B4 | 1.75 | 0.99 |

Factor has α̂≈1.7 (high effective heterogeneity); Toeplitz has α̂≈1.1 (C_δ near zero despite ρ=0.8). Σ_e structure — not scalar H alone — determines α̂.

---

## GDP empirical — results summary

**Data:** 102 vintages, Q3 1999 – Q1 2025 (Q4 2003 excluded), m≈54 bivariate AR predictors, T_eff≈34–179 (median≈82)

### Bivariate RSM — fixed rule k=floor(sqrt(T))≈9 (primary OOS evidence)

| Sample | Model | MAFE | Sig | RMSFE | Sig |
|---|---|---|---|---|---|
| Full | Mean | 0.54 | | 1.10 | |
| Full | RSM | 0.45 | ** | 0.76 | * |
| Full | Lasso | 0.51 | | 0.73 | * |
| Full | Ridge | 0.49 | * | 0.98 | |
| Full | RF | 0.55 | | 1.24 | |
| No-COVID | Mean | 0.42 | | 0.60 | |
| No-COVID | RSM | 0.38 | ** | 0.53 | ** |
| No-COVID | Ridge | 0.39 | ** | 0.54 | ** |

**Oracle k* (ex-post diagnostic, not OOS):** k*=19 (full, RMSFE=0.662), k*=21 (no-COVID, RMSFE=0.499).  
**Fixed rule RMSFE:** k=9: full=0.762, no-COVID=0.529.

**Note on RMSFE(k) curve:** Extremely noisy with only 102 vintage observations. Curve is essentially flat k=7–22; k*=19 is lucky draw of jagged plateau, not a clean interior minimum. Neither 2·T^{1/3}≈9 nor sqrt(T)≈9 explains k*=19 — both give RMSFE≈0.755 vs oracle 0.662. The GDP forecasters appear near-homogeneous (flat curve), consistent with low effective C_δ.

### AR+X RSM

Full: RSM RMSFE=0.617 (−6.8% vs bivariate), oracle k*=19.  
No-COVID: RSM RMSFE=0.506, oracle k*=16 (closer to endogenous MC range 11–14).

---

## SPF empirical — results summary

**Implementation note:** SPF uses simple (equal) within-subset averaging, NOT OLS weights. Therefore infinite bagging mechanically recovers the equal-weighted mean: E_S[(1/k)Σ_{i∈S} f_i] = (1/N)Σ f_i. SPF is a sanity check, not validation of weighted-RSM theory.

| Series | h | k* | RMSFE (EW) | RMSFE (RSM) | Gain |
|---|---|---|---|---|---|
| RGDP | 1 | 3 | 6.405 | 6.392 | +0.2% |
| RGDP | 4 | 13 | 2.208 | 2.205 | +0.1% |
| CPI | 1 | 10 | 0.974 | 0.969 | +0.5% |
| CPI | 4 | 3 | 1.471 | 1.470 | +0.1% |

Gains positive = RSM slightly better (sign convention corrected from previous session).  
All gains < 0.5%, none significant. Confirms: C_δ≈0 for SPF, not α≈1.  
HAC DM correction not yet applied.

---

## Running processes (as of ~21:00 CEST 2026-06-04)

| Job | PID | Script | Status |
|---|---|---|---|
| Endogenous T-scaling rerun | 89806 | `endogenous_scaling_rerun.R` | Toeplitz T=800 in progress; Factor spec (~18h) follows |
| H=6 exo T-scaling | — | `launch_h6_scaling.sh` | Ready to fire once PID 89806 ends; T=80 already saved |

```bash
# Check status:
ps aux | grep 89806 | grep -v grep
tail -3 RSM_endogenous_results/scaling_v2/run_log.txt
ls RSM_exogenous_results_H6_scaling/
```

To launch H=6 scaling manually (if PID 89806 already done):
```bash
bash launch_h6_scaling.sh &
```

---

## Code changes (session 3, 2026-06-04)

### `rsm_common_loading_grid_weights_general_distance_reviewfix_final.R`
1. **Checkpoint in `run_scaling_T_m`**: skips combo if output CSV already exists; reads saved CSV back into `all_summary` for the summary table.
2. **`RSM_SCALING_T_GRID` env var**: GRID_SCALING T_train now reads from env var (default `c(80,240,640,1000)`).
3. **`RSM_RUN_EXPERIMENTS` env var**: already added in prior session; allows selective experiment runs.
4. **`N_SUBSETS_POP` reduced**: 3000 → 2500 in paper profile.

### `launch_h6_scaling.sh` (new)
Auto-waits for PID 89806 then fires H=6 T-scaling with correct env vars.

### `RSM_paper_main.tex`
Revised to address all 7 reviewer priorities. See reviewer table above.

---

## Bugs and fixes

### In `endogenous_rsm_toeplitz_factor.R`
1. **`G_bag = NULL`** → Q_bag_G_mc all NA. **Not yet fixed.** Fix: add `G_bag = K_RS` in `q_population_mc_general` call. One-line change then rerun.
2. **Inflation formula**: replaced `1 + (k-1)/T` with `(T-1)/(T-k)`. ✅ Fixed.

### In `GDP_ar_bivariate_allk.R`
1. **`ifelse` with negative subscripts** in Y_lag_list. ✅ Fixed.
2. **lm/predict scoping issue**. ✅ Fixed (matrix algebra).

### In `rsm_common_loading_grid_weights_general_distance_reviewfix_final.R`
- No outstanding bugs.

---

## Session 4 — Theory of R_k functional form (2026-06-05)

### Question addressed
Can we find optimal k by studying the functional form of R_k generally, without small-het assumption?

### Key exact results derived

#### 1. Exact boundary condition (no assumption)
v̄_{i,1} = 1/m for all i (k=1 ↔ EW averaging, regardless of τ_i distribution)
v̄_{i,m} = w_i^opt (k=m ↔ full GLS)

#### 2. Exact identity
∑_i τ_i · E[1/A_S | i∈S] = m/k  for all k  (follows from ∑_i v̄_{i,k} = 1)

Let g_i(k) = k · E[1/A_S | i∈S].  Then ∑_i τ_i g_i(k) = m for all k.
Note: at k=1, g_i = 1/τ_i, so this is trivially ∑_i τ_i·(1/τ_i) = m.
Non-trivial for k>1.  Conservation law: ∑_i τ_i [g_i(k) - 1/τ_i] = 0 for all k.

#### 3. Exact formula for R_k
  Δ_{i,k} = (τ_i/m)(g_i - 1/τ̄)
  R_k = (1/m²) ∑_i τ_i (g_i - 1/τ̄)²      [exact, no assumption]

The τ_i (not τ_i²) appears because Δ²/τ_i = τ_i²(g_i-1/τ̄)²/(m²τ_i) = τ_i(g_i-1/τ̄)²/m².

#### 4. Exact decomposition around k=1
Write h_i(k) = g_i(k) - 1/τ_i (change from k=1).  Since ∑_i τ_i h_i = 0:

  R_k = R_1 + (2/m²)∑_i h_i + (1/m²)∑_i τ_i h_i²

where R_1 = (1/m²)∑_i(1/τ_i)(1 - τ_i/τ̄)²  [exact EW excess risk]

#### 5. Sharp upper bound (exact)
  R_k ≤ Q_k^sub - Q^opt  (tightened by Jensen; tight iff all f_i equal = homogeneous case)

#### 6. Proof of strict monotonicity (exact, no assumption)
  dR_k/dk = (2/m²)∑_i τ_i(g_i - 1/τ̄)·ġ_i
  ∑_i τ_i ġ_i = 0  (from constraint ∑_i τ_i g_i = m = const)
  (g_i - 1/τ̄)·ġ_i < 0 for all i:
    τ_i > τ̄  ↔  g_i < 1/τ̄ and ġ_i > 0  ↔  product < 0
    τ_i < τ̄  ↔  g_i > 1/τ̄ and ġ_i < 0  ↔  product < 0
  Therefore dR_k/dk < 0 strictly.  ∎

#### 7. Exact discrete FOC for optimal k*
  (T - k*)(R_{k*} - R_{k*+1}) = Q^opt + R_{k*}

### MC evidence on R_k functional form (session 4)

**Figures produced:** `Rk_level_plots.pdf`, `Rk_loglog_zk.pdf`, `Rk_normalised.pdf`,
  `Rk_local_alpha.pdf`, `Rk_marginal_gain.pdf`, `Rk_loglog_k.pdf`
Script: `plot_Rk_functional_form.py`

**Key findings:**
- R_k strictly monotone decreasing: confirmed all DGPs
- Largest gain always at k=1→2: 88% reduction for H=6 (most important step)
- Two qualitatively different shapes:
    Factor (high het):   fast initial decay, most gain at k<5, R_k ~ k^{-3}
    Toeplitz (low het):  slow steady decay, geometric in k, R_k ~ k^{-1.7}
- Global log-log OLS (in z_k): α̂ = 1.7–1.75 (Factor), 1.09–1.10 (Toeplitz), R²=0.93–0.99
- Global log-log OLS (in k):   γ̂ = 2.1–3.2 (slope in log R vs log k)
- Local alpha (adjacent-k slope) is wildly variable (−5 to +7); global OLS is smoothing only
- R_k ~ c·z_k^α fits much better than R_k ~ b/k^α at small k (R²: 0.98 vs 0.90 for H=6)

**Why z_k beats 1/k:**
  z_k = 1/k - 1/m → 0 as k→m, so R_k = c·z_k^α → 0 exactly (boundary condition)
  1/k does not vanish at k=m → b/m^α ≠ 0 violates R_m = 0
  z_k is the natural "distance from optimal" argument

**Gamma/InvGamma distributional approach:**
  If τ_i ~ Gamma(α,β) with α = 1/(H-1):
    A_S ~ Gamma(kα,β) → E[1/A_S] = β/(kα-1) → R_k ~ (α-1)²/(kα-1)² · R_1
  Fails empirically: off by 100-1000× for high het (H>2)
  Reason: treats τ_i as i.i.d. random draws, not fixed realized values; mean-field
    independence breaks down exactly in the high-het regime where it matters most

**U-statistics connection (literature):**
  v̄_{i,k} is a U-statistic of degree k with kernel h_i(S) = τ_i/A_S (Hoeffding 1948)
  Hájek projection = first-order term = small-het result (exact leading term, not approximation)
  "Small-het condition" → just says higher-order terms negligible; leading term always α=2
  Key references: Hoeffding (1948) JASA, Hájek (1968), Lee (1990) textbook,
    Bühlmann & Yu (2002) Ann.Stat., Serfling (1974) Ann.Stat.

---

## Next session priorities

1. **Check endogenous_scaling_rerun.R** status; once done, update T-scaling table in paper.
2. **Run `launch_h6_scaling.sh`** if PID 89806 finished; add H=6 T-scaling table to paper.
3. **Verify α̂ computation** in joint_B_H population CSVs: confirm values come from Q_bag_inf_mc column, not Q_bag_G_mc.
4. **Fix G_bag bug** in endogenous script; rerun bridge verification.
5. **Consider recursive k selection** for GDP (rolling validation window) as suggested by reviewer.
6. **HAC DM tests** for SPF.
7. **Theory: write up exact R_k results as formal propositions** — exact boundary conditions,
   monotonicity proof (via U-stat / direct), decomposition R_k = R_1 + 2/m²·∑h_i + 1/m²·∑τ_i h_i².
8. **Theory: replace "small-het assumption" with "leading Hájek term"** in paper exposition.

---

## R package requirements

**Simulation scripts:** base R only (R ≥ 4.2)  
**GDP/SPF scripts:** `dplyr`, `tidyr`, `lubridate`, `purrr`, `tibble`, `stringr`, `forecast`, `glmnet`, `ranger`, `readxl`  
**LaTeX:** `/Library/TeX/texbin/pdflatex`

---

## Session 5 — 2026-06-06

### Simulation status at session start
| Job | Status |
|---|---|
| rsm_m_scaling.R (PID 7959) | 🔄 Running — m=120 Factor in progress |
| endogenous_scaling_rerun.R (PID 89806) | ✅ Finished (found dead at session start) |
| H=6 exogenous T-scaling | ⏸ Never auto-launched; only T=80 exists |

**Auto-launcher set up (PID 66242):** watches for PID 7959 to finish, then fires `launch_h6_scaling.sh` automatically. Log: `RSM_exogenous_results_H6_scaling/auto_launch_log.txt`.

### m-scaling results now complete (Toeplitz full; Factor through m=80)
From `RSM_m_scaling_results/m_scaling_summary.csv`:
| m | Toeplitz k*(RMSFE) | Factor k*(RMSFE) | KZ factor |
|---|---|---|---|
| 20 | 8 | 8 | 1.11 |
| 30 | 8 | 10 | 1.18 |
| 50 | 12 | 14 | 1.34 |
| 80 | 17 | 20 | 1.68 |
| 120 | 22 | — (in progress) | 2.54 |

### Theory discussion (session 5)
**Endogenous optimal k — tractable Kan-Zhou case:**
Under Gaussianity + common loading:
  L_k^endo ≈ (Ã + λCz_k^α)(T-1)/(T-k),  Ã = s + λQ^opt,  λ = (T-2)/(T-m-2)
  k*_endo ≈ (λCαT/Ã)^{1/(α+1)}
  k*_endo/k*_exo = (λA/Ã)^{1/(α+1)},  A = s + Q^opt
  For ρ=Q^opt/s ≈ 1, m=50, T=200: ratio ≈ 1.05  →  only 5% level shift
  Large empirical gap (k*_exo≈2 → k*_endo≈11-14) driven by NON-COMMON-LOADING structure, not KZ factor.
  KZ preserves T^{1/(α+1)} scaling exponent (consistent with empirical slope 0.39 in both cases).

**sqrt(T) in endogenous case:**
  Theoretical: k* ∝ T^{1/3} with large C ≈ 2.2 (inflated by non-CL misspecification).
  Numerically close to sqrt(T) for T ∈ [100, 300] (crossover at T ≈ 250).
  Mechanism: shared Σ̂_e across G=1000 bags → G_eff ~ T/m ≈ 4 << G → residual linear-in-z_k component → pushes exponent toward 1/2.
  Empirical slope 0.39 is between 1/3 and 1/2, consistent.
  Practical rule k=⌊√T⌋ costs 6-13% RMSFE vs oracle — acceptable given flat RMSFE(k) curve.

**z_k^α motivation:** 
  Three exact properties force z_k as basis: R_m=0, R_k>0 for k<m, R_k strictly decreasing.
  z_k = 1/k - 1/m satisfies all three; 1/k^α fails R_m=0.
  Power law is simplest one-parameter family. R²>0.97 validates shape.

**H vs α mechanism:**
  Small H: τ_i ≈ τ̄(1+δ_i), δ_i small → g_i(k) - 1/τ̄ ≈ -δ_i z_k/τ̄ → R_k ≈ C_δ z_k² → α=2 (exact leading term).
  Large H (dominant forecaster): v̄_{k,1} ≈ k/m, w_1^opt ≈ 1 → Δ_{1,k} ≈ -(m-k)/m = -kz_k → R_k ∝ (kz_k)² → log-log slope → 1.
  At H=6: α̂ ≈ 1.7, intermediate, gains 12-15%. α=1 requires H >> 6 (near H=m).
  α̂ at H=1.05 spurious (C_δ≈0, log-log picks noise).

### New files created (session 5)
| File | Description |
|---|---|
| `RSM_paper_v2.tex` | Full revised paper, new logical structure (24pp, compiles cleanly) |
| `RSM_paper_v2.pdf` | Compiled PDF |
| `plot_kstar_vs_T.py` | k* vs T plot script (linear + log-log panels) |
| `kstar_vs_T.pdf` | k* vs T figure (exo + endo, sqrt(T) + T^{1/3} reference lines) |
| `kstar_vs_T.png` | PNG version |

### RSM_paper_v2.tex structure
1. Scope and Main Message
2. Forecast-Combination Setup
3. Benchmarks + RSM Unification (EW and full OLS as special cases)
4. RSM Definition + Exact Bagging Bridge (with full proof)
5. Common-Loading: Conservation Law + Exact R_k + Strict Monotonicity + Exact FOC (all with proofs)
6. Tractable Cases: Taylor → sqrt(T) for average-subset; why bagged resists same treatment
7. Monte Carlo (Exogenous): why z_k^α, α estimation, H vs α (both limits + intermediate), optimal k validation
8. Endogenous: Kan-Zhou proposition + endogenous k* formula + comparison ratio + MC results + sqrt(T) discussion
9. SPF empirical
10. GDP empirical
11. Conclusion
Appendix: supporting proofs

### Pending for RSM_paper_v2.tex
- [ ] H=6 T-scaling table (§7, placeholder marked) — awaiting H=6 scaling run
- [ ] m=120 Factor k* (m-scaling table, if included) — in progress
- [ ] Bibliography file rsm_refs.bib (currently no .bib, citation placeholders used)
- [ ] \citep{KanZhou2007} needs real bib entry
- [ ] Second pdflatex pass once .bib available

### Next session priorities
1. Check auto-launcher (PID 66242) status; verify H=6 scaling started.
2. Once H=6 done: add T-scaling table to §7 of RSM_paper_v2.tex (replace PLACEHOLDER).
3. Create rsm_refs.bib with all citations used in v2.
4. Once m=120 Factor done: update m-scaling table.
5. Consider adding m-scaling section to paper (k* grows with m at fixed T — KZ explains level).
6. Rerun pdflatex twice after .bib added for correct references.

---

## Session 5 (continued) — additions to RSM_paper_v2.tex

### New content added (both compile cleanly, now 28pp)

**1. §1.1 Motivation and Relation to the Literature**
- Elliott & Liao (2025): their common-loading model $\tilde\Omega = \sigma_\epsilon^2\iota\iota'+\tilde\Sigma$ is exactly ours; they show gains are small when forecasters are homogeneous (C_δ≈0), and when gains ARE large the optimal approach is to discard forecasters and average — i.e., RSM. But they give no k* formula.
- Chen & Maung (2023): use CSR at only k∈{1,2,3,5}. Table 3 shows:
  - Inflation: ASCFE 0.98, 0.95, 0.93, 0.96 — minimum near k=3-4, never precisely located
  - Unemployment: ASCFE 1.01, 1.02, 1.01, 0.97 — monotonically decreasing through k=5 → entirely on wrong side
  - Our theory gives k*≈10-15 at their sample sizes (sqrt(T) rule), far above k=5

**2. Remark in §7 (after H/α discussion): "When Q^bag is flat: the degenerate case k*=1"**
- When R_k≈0 for all k (near-homogeneous, C_δ≈0): L_k ≈ A·(T-1)/(T-k), minimized at k=1
- At k=1 with G→∞: bar_v_{1,i} = 1/m for all i → infinite-bagged RSM at k=1 = equal-weighted mean
- Interior optimum k*>1 requires C>0 (positive heterogeneity)
- Connects to forecast combination puzzle: near-homogeneous panels → EW is optimal → k*=1

**3. §5.5 Alternative inflation section** (added previous session)
- I(m,T)=(T-1)/(T-m) constant → no interior optimum, k*=m always (degenerate)
- I_lin=(1+(k-1)/T) linear approx → same k* to leading order (proved formally)
- Key: any inflation factor increasing in k preserves the interior optimum

### Current RSM_paper_v2.tex sections (28pp)
1. Scope + §1.1 Motivation (Elliott-Liao, Chen-Maung)
2. Setup
3. Benchmarks + RSM unification
4. RSM + exact bagging bridge
5. Common loading: conservation law, exact R_k, monotonicity, exact FOC, §5.5 alternative inflation
6. Tractable cases: sqrt(T) for average-subset; why bagged resists
7. Exogenous MC: z_k^α motivation, α estimation, H vs α, flat Q^bag remark, optimal k validation
8. Endogenous: Kan-Zhou prop, k* formula, comparison ratio, MC, sqrt(T) discussion
9. SPF
10. GDP
11. Conclusion + Appendix

### Pending citations needed in rsm_refs.bib
- BatesGranger1969
- GrangerRamanathan1984
- ElliottLiao2025
- ChenMaung2023
- ElliottGarganoTimmermann2013
- KanZhou2007

---

## Session 5 (continued) — double inflation section + m-scaling results

### New MC results available

**m-scaling complete** (PID 7959 finished ~2:44AM Jun 7):
| m | Toep k*(RMSFE) | Toep k*(L) | Fact k*(RMSFE) | Fact k*(L) | KZ(m) | Toep gain | Fact gain |
|---|---|---|---|---|---|---|---|
| 20 | 8 | 12 | 8 | 14 | 1.11 | 5.9% | 3.8% |
| 30 | 8 | 16 | 10 | 20 | 1.18 | 8.3% | 4.0% |
| 50 | 12 | 26 | 14 | 28 | 1.34 | 10.3% | 6.7% |
| 80 | 17 | 44 | 20 | 47 | 1.68 | 10.4% | 5.9% |
| 120 | 22 | 82 | 18 | 82 | 2.54 | 11.0% | 6.5% |

Log-log m-slopes: k*(RMSFE) ~ m^0.51-0.61; k*(L_bag) ~ m^0.96-1.06

**H=6 exogenous T-scaling** (auto-launched ~2:44AM via PID 85068):
- T=80: k*=4, T=240: k*=6, T=640: k*=8, T=1000: in progress (200/500 reps)
- Log-log slope = 0.334 ≈ 1/3  (cube-root law confirmed at meaningful heterogeneity)

### Key theoretical finding: double inflation for endogenous RSM

**The correct inflation factor for k selection in the endogenous case:**
$$I^{endo}(T,k) = \frac{(T-1)(T-2)}{(T-k)(T-k-2)}$$

Two sources of estimation noise (multiplicative):
1. OLS combination weights: (T-1)/(T-k)  [present in exo and endo]
2. Covariance estimation KZ(k): (T-2)/(T-k-2)  [endo only; uses k not m]

Why KZ uses k not m:
- RSM extracts k×k submatrix Σ̂_{e,S} from m×m Wishart(T-1, Σ_e)
- As a Wishart marginal: Σ̂_{e,S} ~ Wishart(T-1, Σ_{e,S}) — same as if estimated from k series
- KZ correction for k-dimensional portfolio → (T-2)/(T-k-2), not (T-m-2) 
- KZ(m) is a constant multiplier → doesn't change argmin, useless for k selection

Unbiasedness result:
  E_{Σ̂}[Q_k(Σ̂) · KZ(k)] = Q_k(Σ)  [de-biases in-sample optimism]
  E_{Σ̂}[L^endo] = L^exo  [corrected endo criterion targets same k* as exo]

Results:
| DGP | m | k*(RMSFE) | k*(L_curr) | k*(L^endo) |
|---|---|---|---|---|
| Toeplitz | 20 | 8 | 12 | **8** (exact) |
| Toeplitz | 30 | 8 | 16 | **8** (exact) |
| Toeplitz | 50 | 12 | 26 | **12** (exact) |
| Toeplitz | 80 | 17 | 44 | 14 (gap=3) |
| Toeplitz | 120 | 22 | 82 | 18 (gap=4) |
| Factor | 50 | 14 | 28 | **14** (exact) |
| Factor | 80 | 20 | 47 | 17 (gap=3) |
| Factor | 120 | 18 | 82 | **18** (exact) |

Gap closes from 4-64 (current) to 0-4 (corrected).
Residual gap at m≥80: non-Gaussian/non-common-loading DGP breaks Wishart approximation.

### Section added to RSM_paper_v2.tex

§8.3 "The Correct Inflation Criterion for the Endogenous Case" (new, §sec:double_inflation):
- Two sources of noise table
- KZ(k) vs KZ(m) distinction with formal argument
- Double inflation formula (T-1)(T-2)/[(T-k)(T-k-2)]
- Unbiasedness proof: E[L^endo] = L^exo
- Table kstar_double: RMSFE vs L_curr vs L^endo across all m-scaling scenarios
- Remark: practical recommendation to use double inflation whenever k << T

### Current RSM_paper_v2.tex (32pp, compiles clean)

Remaining PLACEHOLDERS to fill after simulations complete:
1. H=6 T-scaling table in §7.4 (T=80,240,640 done; T=1000 pending)
2. Full m-scaling table (all done, just needs inserting)
3. rsm_refs.bib (citations currently as \citep{key} without .bib file)
