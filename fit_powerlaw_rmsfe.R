#!/usr/bin/env Rscript
###############################################################################
# fit_powerlaw_rmsfe.R
#
# Fits R_k = C * (1/k - 1/m)^alpha directly to the empirical RMSFE curve
# from the endogenous Monte Carlo simulations, using nonlinear least squares.
#
# Two targets are fitted:
#   (1) Excess MSPE:  R_k = MSPE(k) - MSPE(k=m)
#   (2) Excess RMSFE: R_k = RMSFE(k) - RMSFE(k=m)   [cruder, more intuitive]
#
# For each file, also runs the log-log OLS regression on log(R_k) ~ log(z_k)
# for direct comparison.
#
# Outputs a summary table and fit-quality plots.
###############################################################################

library(minpack.lm)  # nlsLM — more robust than base nls()

INDIR  <- file.path(getwd(), "RSM_endogenous_results")
OUTDIR <- file.path(getwd(), "powerlaw_fit_results")
dir.create(OUTDIR, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Note on sign conventions for the endogenous case
# ---------------------------------------------------------------------------
# In the exogenous case, MSPE(k) is monotonically decreasing in k, so
# R_k = MSPE(k) - MSPE(k=m) > 0 and vanishes at k=m.
#
# In the endogenous case (estimated Sigma_e), the full-pool k=m estimate
# suffers from the curse of dimensionality: MSPE(k=m) is the HIGHEST
# value, not the lowest. The RMSFE curve is U-shaped with minimum at k*.
#
# We therefore fit two objects:
#
#  (A) In-sample Q_rs_bag(k): this IS monotonically decreasing because
#      more forecasters always lower the in-sample fit. The floor is
#      Q_rs_bag(m). This directly parallels the exogenous excess risk.
#      R_k^Q = Q_rs_bag(k) - Q_rs_bag(k=m)  > 0.
#
#  (B) RMSFE(k) on the DECREASING branch only (k <= k*(RMSFE)):
#      R_k^rmsfe = RMSFE(k) - RMSFE(k*)  > 0.
#      This captures how much the RMSFE exceeds the optimum as k shrinks.
#      The power-law fit here tests whether the approach-to-optimal has
#      the same functional form C*z_k^alpha.

# ---------------------------------------------------------------------------
# 1. Helper: fit C*(1/k - 1/m)^alpha via NLS and via log-log OLS
# ---------------------------------------------------------------------------

fit_powerlaw <- function(k_grid, R_k, m, label = "") {
  # z_k = 1/k - 1/m  (the natural scaling variable)
  z <- 1 / k_grid - 1 / m

  # Keep only interior k (exclude k=m where z=0 and R_k should be 0)
  keep <- z > 0 & R_k > 0 & is.finite(R_k) & is.finite(z)
  # Also exclude k=1 if present (boundary effects)
  keep <- keep & k_grid > 1

  if (sum(keep) < 3) {
    warning(label, ": fewer than 3 usable points, skipping.")
    return(NULL)
  }

  z_fit <- z[keep]
  R_fit <- R_k[keep]

  # --- NLS fit: R = C * z^alpha -------------------------------------------
  # Starting values from log-log regression
  lfit   <- lm(log(R_fit) ~ log(z_fit))
  # unname() is critical: coef() returns named scalars and those names
  # would contaminate the nlsLM start list, producing garbled coef names.
  C0     <- unname(exp(coef(lfit)[1]))
  alpha0 <- unname(coef(lfit)[2])

  # Note: avoid "C" as parameter name — it shadows stats::C() and causes
  # nlsLM to misname the returned coefficients.
  df_fit  <- data.frame(z_fit = z_fit, R_fit = R_fit)
  nls_fit <- tryCatch(
    nlsLM(
      R_fit ~ Cpow * z_fit^alph,
      data  = df_fit,
      start = list(Cpow = C0, alph = alpha0),
      lower = c(Cpow = 1e-10, alph = 0.1),
      upper = c(Cpow = 1e6,   alph = 5),
      control = nls.lm.control(maxiter = 200, ftol = 1e-8, ptol = 1e-8)
    ),
    error = function(e) {
      message("NLS failed [", label, "]: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(nls_fit)) return(NULL)

  C_hat     <- coef(nls_fit)["Cpow"]
  alpha_hat <- coef(nls_fit)["alph"]

  # Fitted values and residuals (NLS)
  R_fitted_nls <- C_hat * z_fit^alpha_hat
  ss_res_nls   <- sum((R_fit - R_fitted_nls)^2)
  ss_tot       <- sum((R_fit - mean(R_fit))^2)
  R2_nls       <- 1 - ss_res_nls / ss_tot
  rmse_nls     <- sqrt(mean((R_fit - R_fitted_nls)^2))

  # --- Log-log OLS ----------------------------------------------------------
  alpha_ols <- coef(lfit)[2]
  C_ols     <- exp(coef(lfit)[1])
  R_fitted_ols <- C_ols * z_fit^alpha_ols
  ss_res_ols   <- sum((R_fit - R_fitted_ols)^2)
  R2_ols_linear <- 1 - ss_res_ols / ss_tot  # R2 in original (not log) space
  R2_ols_log    <- summary(lfit)$r.squared  # R2 in log-log space

  # --- Implied k* from NLS fit -------------------------------------------
  # k* minimises L_k = (A + C*z^alpha) * (T-1)/(T-k)   [approx for k << m]
  # FOC: A/T = C*alpha / k^(alpha+1)  =>  k* = (C*alpha*T/A)^{1/(alpha+1)}
  # We don't have A (= s + Q_opt) directly, so report C and alpha only.

  list(
    n          = sum(keep),
    C_nls      = C_hat,
    alpha_nls  = alpha_hat,
    R2_nls     = R2_nls,
    rmse_nls   = rmse_nls,
    C_ols      = C_ols,
    alpha_ols  = alpha_ols,
    R2_ols_log = R2_ols_log,
    R2_ols_lin = R2_ols_linear,
    z_fit      = z_fit,
    R_fit      = R_fit,
    R_fitted_nls = R_fitted_nls,
    R_fitted_ols = R_fitted_ols
  )
}

# ---------------------------------------------------------------------------
# 2. Process each by_z file
# ---------------------------------------------------------------------------

files <- sort(list.files(INDIR, pattern = "by_z\\.csv$", full.names = TRUE))
cat("Found", length(files), "by_z files.\n\n")

results <- vector("list", length(files))

for (fi in seq_along(files)) {
  f   <- files[fi]
  dat <- read.csv(f)
  # m may be missing in scaling files — infer from max(z)
  m   <- if ("m" %in% names(dat)) dat$m[1] else max(dat$z)
  T_  <- dat$T_len[1]
  stype <- dat$sigma_type[1]
  bdes  <- if ("beta_design" %in% names(dat)) dat$beta_design[1] else "B4"
  lab   <- paste0(stype, "_T", T_, "_", bdes)

  cat("Processing:", lab, "\n")

  # --- (A) In-sample Q_rs_bag: monotonically decreasing, floor = Q at k=m ---
  Q_floor   <- dat$Q_rs_bag[dat$z == m]
  if (length(Q_floor) == 0) Q_floor <- min(dat$Q_rs_bag, na.rm = TRUE)
  R_Q       <- dat$Q_rs_bag - Q_floor   # positive, decreasing toward 0 at k=m

  # --- (B) Out-of-sample RMSFE on the DECREASING branch (k <= k*) ----------
  k_star_idx  <- which.min(dat$rmsfe_rs)
  k_star_val  <- dat$z[k_star_idx]
  rmsfe_star  <- dat$rmsfe_rs[k_star_idx]
  # Keep only k <= k* for the power-law fit on the decreasing branch
  left_idx    <- dat$z <= k_star_val
  R_rmsfe_sq  <- dat$rmsfe_rs^2 - rmsfe_star^2  # excess MSPE vs optimum

  fit_Q     <- fit_powerlaw(dat$z,           R_Q,          m, label = paste(lab, "Q_inSample"))
  fit_rmsfe <- fit_powerlaw(dat$z[left_idx], R_rmsfe_sq[left_idx], m,
                            label = paste(lab, "RMSFE_left_branch"))

  results[[fi]] <- list(
    label   = lab,
    sigma   = stype,
    T_len   = T_,
    beta    = bdes,
    m       = m,
    k_star  = k_star_val,
    sqrt_T  = floor(sqrt(T_)),
    Q_fit   = fit_Q,
    rmsfe   = fit_rmsfe
  )

  # --- Plot: two panels (Q in-sample + RMSFE left branch) -------------------
  png(file.path(OUTDIR, paste0(lab, "_powerlaw_fit.png")),
      width = 1300, height = 520)
  par(mfrow = c(1, 2), mar = c(4, 4.5, 3, 1))

  ## Panel 1: in-sample Q
  if (!is.null(fit_Q)) {
    z_seq     <- seq(min(fit_Q$z_fit) * 0.5, max(fit_Q$z_fit) * 1.02, length.out = 300)
    R_nls_seq <- fit_Q$C_nls * z_seq^fit_Q$alpha_nls
    R_ols_seq <- fit_Q$C_ols * z_seq^fit_Q$alpha_ols
    ylim_up   <- max(fit_Q$R_fit) * 1.08

    plot(fit_Q$z_fit, fit_Q$R_fit,
         pch = 19, col = "grey30", cex = 1.1,
         xlab = expression(z[k] == 1/k - 1/m),
         ylab = expression(Q[rs]^bag * "(k)" - Q[rs]^bag * "(m)"),
         main = paste0(lab, "\nIn-sample excess Q"),
         ylim = c(0, ylim_up))
    lines(z_seq, R_nls_seq, col = "#0072B2", lwd = 2.5)
    lines(z_seq, R_ols_seq, col = "#D55E00", lwd = 2, lty = 2)
    legend("topleft", bty = "n",
           legend = c(
             sprintf("NLS:    C=%.3f  α=%.3f  R²=%.4f", fit_Q$C_nls, fit_Q$alpha_nls, fit_Q$R2_nls),
             sprintf("LogOLS: C=%.3f  α=%.3f  R²(log)=%.4f", fit_Q$C_ols, fit_Q$alpha_ols, fit_Q$R2_ols_log)
           ),
           col = c("#0072B2", "#D55E00"), lwd = 2, lty = c(1, 2), cex = 0.85)
  } else {
    plot.new(); title(paste0(lab, "\nQ fit failed"))
  }

  ## Panel 2: RMSFE left branch (k <= k*)
  if (!is.null(fit_rmsfe) && length(fit_rmsfe$R_fit) >= 3) {
    z_seq2    <- seq(min(fit_rmsfe$z_fit) * 0.8, max(fit_rmsfe$z_fit) * 1.02, length.out = 300)
    R_nls2    <- fit_rmsfe$C_nls * z_seq2^fit_rmsfe$alpha_nls
    R_ols2    <- fit_rmsfe$C_ols * z_seq2^fit_rmsfe$alpha_ols
    ylim_up2  <- max(fit_rmsfe$R_fit) * 1.08

    plot(fit_rmsfe$z_fit, fit_rmsfe$R_fit,
         pch = 19, col = "grey30", cex = 1.1,
         xlab = expression(z[k] == 1/k - 1/m),
         ylab = expression(RMSFE(k)^2 - RMSFE(k^"*")^2),
         main = paste0("RMSFE excess (k ≤ k*=", k_star_val, ")"),
         ylim = c(0, ylim_up2))
    lines(z_seq2, R_nls2, col = "#0072B2", lwd = 2.5)
    lines(z_seq2, R_ols2, col = "#D55E00", lwd = 2, lty = 2)
    legend("topleft", bty = "n",
           legend = c(
             sprintf("NLS:    C=%.3f  α=%.3f  R²=%.4f", fit_rmsfe$C_nls, fit_rmsfe$alpha_nls, fit_rmsfe$R2_nls),
             sprintf("LogOLS: C=%.3f  α=%.3f  R²(log)=%.4f", fit_rmsfe$C_ols, fit_rmsfe$alpha_ols, fit_rmsfe$R2_ols_log)
           ),
           col = c("#0072B2", "#D55E00"), lwd = 2, lty = c(1, 2), cex = 0.85)
  } else {
    plot.new()
    title(paste0("RMSFE left branch: only ", if(!is.null(fit_rmsfe)) fit_rmsfe$n else 0, " points"))
  }

  dev.off()
}

# ---------------------------------------------------------------------------
# 3. Summary table
# ---------------------------------------------------------------------------

safe <- function(x, digits = 4) if (!is.null(x) && is.finite(x)) round(x, digits) else NA_real_

rows <- lapply(results, function(r) {
  if (is.null(r)) return(NULL)
  Q_ <- r$Q_fit;  rf <- r$rmsfe
  data.frame(
    label            = r$label,
    sigma            = r$sigma,
    T_len            = r$T_len,
    beta             = r$beta,
    k_star           = r$k_star,
    sqrt_T           = r$sqrt_T,
    C_Q_nls          = safe(if (!is.null(Q_)) Q_$C_nls),
    alpha_Q_nls      = safe(if (!is.null(Q_)) Q_$alpha_nls, 3),
    R2_Q_nls         = safe(if (!is.null(Q_)) Q_$R2_nls),
    alpha_Q_ols      = safe(if (!is.null(Q_)) Q_$alpha_ols, 3),
    R2_Q_ols_log     = safe(if (!is.null(Q_)) Q_$R2_ols_log),
    R2_Q_ols_lin     = safe(if (!is.null(Q_)) Q_$R2_ols_lin),
    n_rmsfe_pts      = if (!is.null(rf)) rf$n     else NA_integer_,
    alpha_rmsfe_nls  = safe(if (!is.null(rf)) rf$alpha_nls, 3),
    R2_rmsfe_nls     = safe(if (!is.null(rf)) rf$R2_nls),
    alpha_rmsfe_ols  = safe(if (!is.null(rf)) rf$alpha_ols, 3),
    R2_rmsfe_ols_log = safe(if (!is.null(rf)) rf$R2_ols_log),
    stringsAsFactors = FALSE
  )
})

summary_df <- do.call(rbind, rows)
write.csv(summary_df, file.path(OUTDIR, "powerlaw_fit_summary.csv"), row.names = FALSE)

cat("\n=== Power-law fit: in-sample Q_rs_bag (full k range) ===\n")
cat(sprintf("%-32s %5s %4s %5s | %6s %6s %6s | %6s %6s\n",
  "Label", "T", "k*", "sqT", "C_NLS", "a_NLS", "R2_NLS", "a_OLS", "R2_log"))
for (i in seq_len(nrow(summary_df))) {
  r <- summary_df[i, ]
  cat(sprintf("%-32s %5d %4d %5d | %6.3f %6.3f %6.4f | %6.3f %6.4f\n",
    r$label, r$T_len, r$k_star, r$sqrt_T,
    r$C_Q_nls, r$alpha_Q_nls, r$R2_Q_nls,
    r$alpha_Q_ols, r$R2_Q_ols_log))
}

cat("\n=== Power-law fit: RMSFE excess on left branch (k <= k*) ===\n")
cat(sprintf("%-32s %5s %4s %4s %4s | %6s %6s | %6s %6s\n",
  "Label", "T", "k*", "sqT", "npts", "a_NLS", "R2_NLS", "a_OLS", "R2_log"))
for (i in seq_len(nrow(summary_df))) {
  r <- summary_df[i, ]
  if (is.na(r$alpha_rmsfe_nls)) {
    cat(sprintf("%-32s %5d %4d  -- too few points\n", r$label, r$T_len, r$k_star))
  } else {
    cat(sprintf("%-32s %5d %4d %4d %4d | %6.3f %6.4f | %6.3f %6.4f\n",
      r$label, r$T_len, r$k_star, r$sqrt_T, r$n_rmsfe_pts,
      r$alpha_rmsfe_nls, r$R2_rmsfe_nls,
      r$alpha_rmsfe_ols, r$R2_rmsfe_ols_log))
  }
}

# ---------------------------------------------------------------------------
# 4. Multi-panel summary plot: NLS alpha vs T for each sigma type
# ---------------------------------------------------------------------------
png(file.path(OUTDIR, "alpha_nls_vs_T.png"), width = 900, height = 500)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
for (stype in c("toeplitz", "factor")) {
  sub <- summary_df[summary_df$sigma == stype & grepl("B4", summary_df$beta), ]
  sub <- sub[order(sub$T_len), ]
  if (nrow(sub) == 0) next

  ylim <- range(c(sub$alpha_Q_nls, sub$alpha_Q_ols,
                  sub$alpha_rmsfe_nls, sub$alpha_rmsfe_ols, 1, 2), na.rm = TRUE)
  plot(sub$T_len, sub$alpha_Q_nls,
       type = "b", pch = 19, lwd = 2, col = "#0072B2",
       xlab = "T", ylab = expression(hat(alpha)),
       main = paste0(stype, ": fitted alpha vs T (B4)"),
       ylim = ylim)
  lines(sub$T_len, sub$alpha_Q_ols,    type = "b", pch = 17, lwd = 2, lty = 2, col = "#0072B2")
  lines(sub$T_len, sub$alpha_rmsfe_nls,type = "b", pch = 19, lwd = 2, lty = 1, col = "#D55E00")
  lines(sub$T_len, sub$alpha_rmsfe_ols,type = "b", pch = 17, lwd = 2, lty = 2, col = "#D55E00")
  abline(h = 1, lty = 3, col = "grey50"); abline(h = 2, lty = 3, col = "grey50")
  legend("topright", bty = "n",
         legend = c("Q: NLS", "Q: LogOLS", "RMSFE left: NLS", "RMSFE left: LogOLS"),
         col = c("#0072B2","#0072B2","#D55E00","#D55E00"),
         lwd = 2, lty = c(1,2,1,2), pch = c(19,17,19,17))
}
dev.off()

cat("\nDone. Outputs in:", OUTDIR, "\n")
