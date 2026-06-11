#!/usr/bin/env Rscript
# Debug why B2 scenario crashes in run_endo_kby1.R
suppressMessages(library(minpack.lm))

endo_lines <- readLines("endogenous_rsm_toeplitz_factor.R")
stop_at <- grep("^# 9\\. Main run block", endo_lines)[1] - 1L
.saved <- "RSM_endogenous_kby1"
eval(parse(text = endo_lines[1:stop_at]), envir = .GlobalEnv)
K_BY <- 1L; OUTDIR <- .saved

sc <- build_scenarios()[[1]]
cat("Scenario:", sc$sigma_type, "beta=", sc$beta_design,
    "m=", sc$m, "T=", sc$T_len, "\n")

Sigma_e <- make_sigma_e(
  type = sc$sigma_type, m = sc$m,
  rho  = ifelse(is.na(sc$rho), 0, sc$rho),
  seed = GLOBAL_SEED + 1000L * sc$scenario_id
)
beta_r <- beta_raw_design(sc$m, sc$beta_design)
z_grid   <- seq(2L, sc$m, by = 1L)
cat("z_grid length:", length(z_grid), "\n")

set.seed(123)
for (rr in 1:60) {
  cat("rep", rr, "... ")
  result <- tryCatch(
    withCallingHandlers(
      run_one_rep_all_z(sc$T_len, sc$m, beta_r, Sigma_e,
                        sc$phi_f, sc$phi_eta, sc$phi_u,
                        z_grid = z_grid, K_RS = 200, OOS_LEN = 30),
      warning = function(w) {
        cat("\n  WARN:", conditionMessage(w), "\n")
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      cat("\n  ERROR:", conditionMessage(e), "\n")
      NULL
    }
  )
  if (!is.null(result)) cat("OK\n") else cat("FAILED\n")
}
cat("Done.\n")
