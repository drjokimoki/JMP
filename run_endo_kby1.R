#!/usr/bin/env Rscript
###############################################################################
# run_endo_kby1.R
#
# Runs the endogenous Monte Carlo identically to the "run" profile in
# endogenous_rsm_toeplitz_factor.R, with ONE change: K_BY = 1 instead of 3.
# Everything else — N_REP, K_RS, seeds, scenarios, T-scaling — is unchanged.
#
# Results go to RSM_endogenous_kby1/ so the original RSM_endogenous_results/
# is not touched and can be used for comparison.
#
# After the simulation, fits C*(1/k - 1/m)^alpha via NLS and log-log OLS
# to both the in-sample Q_rs_bag curve and the RMSFE left-branch, and prints
# a comparison table.
###############################################################################

suppressMessages(library(minpack.lm))

KBY1_OUTDIR <- file.path(getwd(), "RSM_endogenous_kby1")
dir.create(KBY1_OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# 1. Load function definitions from the endogenous script (stop before
#    "# 9. Main run block" so the original run is not triggered).
# ---------------------------------------------------------------------------
endo_lines <- readLines("endogenous_rsm_toeplitz_factor.R")
stop_at    <- grep("^# 9\\. Main run block", endo_lines)[1] - 1L
cat(sprintf("Loading functions from lines 1-%d of endogenous_rsm_toeplitz_factor.R\n",
            stop_at))

# The endogenous script calls rm(list = ls()) near the top, which would wipe
# KBY1_OUTDIR from the environment. Save it, eval, then restore.
.kby1_saved <- KBY1_OUTDIR
eval(parse(text = endo_lines[1:stop_at]), envir = .GlobalEnv)
KBY1_OUTDIR <- .kby1_saved
rm(.kby1_saved)

# ---------------------------------------------------------------------------
# 2. Override ONLY K_BY and OUTDIR. Everything else comes from the "run"
#    profile already loaded above (N_REP=200, K_RS=200, same seeds, etc.).
# ---------------------------------------------------------------------------
K_BY   <- 1L
OUTDIR <- KBY1_OUTDIR
# Also set K_BY_scale = 1 in the T-scaling config loaded from the endogenous script
T_SCALING_CONFIG$K_BY_scale <- 1L

cat(sprintf("K_BY overridden to %d\n", K_BY))
cat(sprintf("Output directory: %s\n", OUTDIR))
cat(sprintf("N_REP=%d  K_RS=%d  OOS_LEN=%d\n", N_REP, K_RS, OOS_LEN))

# ---------------------------------------------------------------------------
# 3. Part A: main scenarios (identical to the original run)
# ---------------------------------------------------------------------------
cat("\n--- Part A: main scenarios ---\n")
scenarios <- build_scenarios()
cat("Number of scenarios:", length(scenarios), "\n")

all_summaries <- vector("list", length(scenarios))
for (ii in seq_along(scenarios)) {
  cat(sprintf("\nScenario %d / %d\n", ii, length(scenarios)))
  res              <- run_scenario(scenarios[[ii]])
  all_summaries[[ii]] <- res$summary
}

summary_df <- do.call(rbind, all_summaries)
write.csv(summary_df,
          file.path(KBY1_OUTDIR, "endogenous_rsm_summary.csv"),
          row.names = FALSE)

# ---------------------------------------------------------------------------
# 4. Part B: T-scaling (identical to original, K_BY already = 1)
# ---------------------------------------------------------------------------
cat("\n--- Part B: T-scaling ---\n")
scaling_results <- run_scaling_T_endo(T_SCALING_CONFIG, KBY1_OUTDIR)

cat("\nAll simulations completed.\n")

# ---------------------------------------------------------------------------
# 5. Power-law fitting on the dense k grid
# ---------------------------------------------------------------------------
cat("\n\n=== Power-law fitting (K_BY=1) ===\n\n")

fit_powerlaw_fn <- function(k_grid, R_k, m, label = "") {
  z    <- 1 / k_grid - 1 / m
  keep <- z > 0 & R_k > 0 & is.finite(R_k) & is.finite(z) & k_grid > 1
  if (sum(keep) < 4) {
    message(label, ": fewer than 4 usable points.")
    return(NULL)
  }
  z_fit <- z[keep];  R_fit <- R_k[keep]
  lfit  <- lm(log(R_fit) ~ log(z_fit))
  C0    <- unname(exp(coef(lfit)[1]))
  a0    <- unname(coef(lfit)[2])
  df    <- data.frame(zf = z_fit, Rf = R_fit)

  nls_fit <- tryCatch(
    nlsLM(Rf ~ Cpow * zf^alph, data = df,
          start = list(Cpow = C0, alph = a0),
          lower = c(Cpow = 1e-10, alph = 0.1),
          upper = c(Cpow = 1e3,   alph = 6),
          control = nls.lm.control(maxiter = 300, ftol = 1e-9, ptol = 1e-9)),
    error = function(e) { message("NLS [", label, "]: ", e$message); NULL }
  )
  if (is.null(nls_fit)) return(NULL)

  C_hat <- unname(coef(nls_fit)["Cpow"])
  a_hat <- unname(coef(nls_fit)["alph"])
  pred  <- C_hat * z_fit^a_hat
  ss_r  <- sum((R_fit - pred)^2)
  ss_t  <- sum((R_fit - mean(R_fit))^2)
  pred_ols <- C0 * z_fit^a0
  ss_ols   <- sum((R_fit - pred_ols)^2)

  list(n = sum(keep),
       C_nls = C_hat,  alpha_nls = a_hat,
       R2_nls = 1 - ss_r / ss_t,
       C_ols = C0,     alpha_ols = a0,
       R2_ols_log = summary(lfit)$r.squared,
       R2_ols_lin = 1 - ss_ols / ss_t,
       z_fit = z_fit, R_fit = R_fit,
       R_pred_nls = pred, R_pred_ols = pred_ols)
}

by_z_files <- sort(list.files(KBY1_OUTDIR, pattern = "by_z\\.csv$",
                               full.names = TRUE))
cat("Found", length(by_z_files), "by_z files.\n")

PLOT_DIR <- file.path(KBY1_OUTDIR, "powerlaw_plots")
dir.create(PLOT_DIR, showWarnings = FALSE)

fit_rows <- list()

for (f in by_z_files) {
  dat   <- read.csv(f)
  m_val <- if ("m" %in% names(dat)) dat$m[1] else max(dat$z)
  T_val <- dat$T_len[1]
  stype <- dat$sigma_type[1]
  bdes  <- if ("beta_design" %in% names(dat)) dat$beta_design[1] else "B4"
  lab   <- paste0(stype, "_T", T_val, "_", bdes)
  cat("Fitting:", lab, "(n =", nrow(dat), "k values)\n")

  # (A) In-sample Q excess: Q_rs_bag(k) - Q_rs_bag(m)
  Q_floor <- dat$Q_rs_bag[dat$z == m_val]
  if (length(Q_floor) == 0) Q_floor <- min(dat$Q_rs_bag, na.rm = TRUE)
  R_Q <- dat$Q_rs_bag - Q_floor

  # (B) RMSFE^2 excess on the decreasing branch (k <= k*)
  k_star    <- dat$z[which.min(dat$rmsfe_rs)]
  rmsfe_min <- min(dat$rmsfe_rs, na.rm = TRUE)
  R_mspe    <- dat$rmsfe_rs^2 - rmsfe_min^2
  left      <- dat$z <= k_star

  fitQ <- fit_powerlaw_fn(dat$z,        R_Q,         m_val, paste(lab, "Q"))
  fitR <- fit_powerlaw_fn(dat$z[left],  R_mspe[left], m_val, paste(lab, "RMSFE"))

  fit_rows[[length(fit_rows) + 1]] <- data.frame(
    label = lab, sigma = stype, T_len = T_val, beta = bdes,
    m = m_val, k_star = k_star, sqrt_T = floor(sqrt(T_val)),
    n_Q          = if (!is.null(fitQ)) fitQ$n             else NA_integer_,
    C_Q_nls      = if (!is.null(fitQ)) round(fitQ$C_nls,   4) else NA_real_,
    alpha_Q_nls  = if (!is.null(fitQ)) round(fitQ$alpha_nls,3) else NA_real_,
    R2_Q_nls     = if (!is.null(fitQ)) round(fitQ$R2_nls,  4) else NA_real_,
    alpha_Q_ols  = if (!is.null(fitQ)) round(fitQ$alpha_ols,3) else NA_real_,
    R2_Q_log     = if (!is.null(fitQ)) round(fitQ$R2_ols_log,4) else NA_real_,
    R2_Q_lin     = if (!is.null(fitQ)) round(fitQ$R2_ols_lin,4) else NA_real_,
    n_R          = if (!is.null(fitR)) fitR$n             else NA_integer_,
    alpha_R_nls  = if (!is.null(fitR)) round(fitR$alpha_nls,3) else NA_real_,
    R2_R_nls     = if (!is.null(fitR)) round(fitR$R2_nls,  4) else NA_real_,
    alpha_R_ols  = if (!is.null(fitR)) round(fitR$alpha_ols,3) else NA_real_,
    R2_R_log     = if (!is.null(fitR)) round(fitR$R2_ols_log,4) else NA_real_,
    stringsAsFactors = FALSE
  )

  # Plot
  png(file.path(PLOT_DIR, paste0(lab, "_fit.png")), width = 1300, height = 520)
  par(mfrow = c(1, 2), mar = c(4, 4.5, 3, 1))

  if (!is.null(fitQ)) {
    zs <- seq(min(fitQ$z_fit) * 0.3, max(fitQ$z_fit), length.out = 500)
    plot(fitQ$z_fit, fitQ$R_fit, pch = 19, col = "grey30", cex = 0.6,
         xlab = expression(z[k] == 1/k - 1/m),
         ylab = expression(Q[rs]^bag * "(k) - " * Q[rs]^bag * "(m)"),
         main = paste0(lab, "  (n=", fitQ$n, " k values)\nIn-sample excess Q"),
         ylim = c(0, max(fitQ$R_fit) * 1.05))
    lines(zs, fitQ$C_nls * zs^fitQ$alpha_nls, col = "#0072B2", lwd = 2.5)
    lines(zs, fitQ$C_ols * zs^fitQ$alpha_ols, col = "#D55E00", lwd = 2, lty = 2)
    legend("topleft", bty = "n", cex = 0.82,
           legend = c(
             sprintf("NLS:    C=%.3f  α=%.3f  R²=%.4f",
                     fitQ$C_nls, fitQ$alpha_nls, fitQ$R2_nls),
             sprintf("LogOLS: C=%.3f  α=%.3f  R²(log)=%.4f",
                     fitQ$C_ols, fitQ$alpha_ols, fitQ$R2_ols_log)),
           col = c("#0072B2", "#D55E00"), lwd = 2, lty = c(1, 2))
  } else { plot.new(); title("Q fit failed") }

  if (!is.null(fitR) && fitR$n >= 4) {
    zs2 <- seq(min(fitR$z_fit) * 0.3, max(fitR$z_fit), length.out = 500)
    plot(fitR$z_fit, fitR$R_fit, pch = 19, col = "grey30", cex = 0.8,
         xlab = expression(z[k]),
         ylab = expression(RMSFE(k)^2 - RMSFE(k^"*")^2),
         main = paste0("RMSFE excess  k ≤ k*=", k_star,
                       "  (n=", fitR$n, " pts)"),
         ylim = c(0, max(fitR$R_fit) * 1.05))
    lines(zs2, fitR$C_nls * zs2^fitR$alpha_nls, col = "#0072B2", lwd = 2.5)
    lines(zs2, fitR$C_ols * zs2^fitR$alpha_ols, col = "#D55E00", lwd = 2, lty = 2)
    legend("topleft", bty = "n", cex = 0.82,
           legend = c(
             sprintf("NLS:    C=%.3f  α=%.3f  R²=%.4f",
                     fitR$C_nls, fitR$alpha_nls, fitR$R2_nls),
             sprintf("LogOLS: C=%.3f  α=%.3f  R²(log)=%.4f",
                     fitR$C_ols, fitR$alpha_ols, fitR$R2_ols_log)),
           col = c("#0072B2", "#D55E00"), lwd = 2, lty = c(1, 2))
  } else {
    plot.new()
    title(paste0("RMSFE left: ", if (!is.null(fitR)) fitR$n else 0,
                 " pts (k ≤ k*=", k_star, ")"))
  }
  dev.off()
}

fit_df <- do.call(rbind, fit_rows)
write.csv(fit_df,
          file.path(KBY1_OUTDIR, "powerlaw_kby1_fits.csv"),
          row.names = FALSE)

# ---------------------------------------------------------------------------
# 6. Print summary tables
# ---------------------------------------------------------------------------
cat("\n=== In-sample Q_rs_bag fit (full k range, K_BY=1) ===\n")
cat(sprintf("%-32s %4s %4s %4s %4s | %6s %6s %6s | %6s %6s %6s\n",
  "Label","T","k*","sqT","n_Q","C_NLS","a_NLS","R2_NLS","a_OLS","R2_log","R2_lin"))
for (i in seq_len(nrow(fit_df))) {
  r <- fit_df[i, ]
  if (is.na(r$n_Q)) next
  cat(sprintf("%-32s %4d %4d %4d %4d | %6.3f %6.3f %6.4f | %6.3f %6.4f %6.4f\n",
    r$label, r$T_len, r$k_star, r$sqrt_T, r$n_Q,
    r$C_Q_nls, r$alpha_Q_nls, r$R2_Q_nls,
    r$alpha_Q_ols, r$R2_Q_log, r$R2_Q_lin))
}

cat("\n=== RMSFE excess left branch (k <= k*, K_BY=1) ===\n")
cat(sprintf("%-32s %4s %4s %4s | %6s %6s | %6s %6s\n",
  "Label","T","k*","npts","a_NLS","R2_NLS","a_OLS","R2_log"))
for (i in seq_len(nrow(fit_df))) {
  r <- fit_df[i, ]
  if (is.na(r$n_R) || r$n_R < 4) {
    cat(sprintf("%-32s %4d %4d  -- too few points (n=%s)\n",
      r$label, r$T_len, r$k_star,
      if (is.na(r$n_R)) "NA" else r$n_R))
  } else {
    cat(sprintf("%-32s %4d %4d %4d | %6.3f %6.4f | %6.3f %6.4f\n",
      r$label, r$T_len, r$k_star, r$n_R,
      r$alpha_R_nls, r$R2_R_nls, r$alpha_R_ols, r$R2_R_log))
  }
}

cat("\nDone. All outputs in:", KBY1_OUTDIR, "\n")
