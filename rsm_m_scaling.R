#!/usr/bin/env Rscript
###############################################################################
# rsm_m_scaling.R
#
# m-scaling experiment for endogenous RSM.
# Holds T fixed, sweeps m in {20, 30, 50, 80, 120}.
# For each (m, Sigma_type), reports k*(RMSFE) and k*(L_bag).
#
# Theory prediction (Kan-Zhou + cube-root):
#   k*_endo(m,T) ~ k*_exo * [(T-2)/(T-m-2)]^{1/(alpha+1)}
# so k* should increase with m holding T fixed.
###############################################################################

rm(list = ls())

###############################################################################
# 0. Parameters
###############################################################################

GLOBAL_SEED   <- 9191
set.seed(GLOBAL_SEED)

OUTDIR <- file.path(getwd(), "RSM_m_scaling_results")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

T_FIXED    <- 200          # fixed sample size
M_GRID     <- c(20, 30, 50, 80, 120)
SIGMA_TYPES <- c("toeplitz", "factor")

N_REP      <- 500
K_RS       <- 500
OOS_LEN    <- 50
BETA_DESIGN <- "B4"
R2_TARGET  <- 0.3
PHI_F      <- 0.3
PHI_ETA    <- 0.7
PHI_U      <- 0.0
SIGMA_U    <- 1.0
RIDGE      <- 1e-8
BETA_SCALE_USING <- "innovation"

# Toeplitz rho
RHO <- 0.8

# Factor params (same as main endogenous run)
FACTOR_R       <- 2L
FACTOR_LOADING <- "random"
FACTOR_STRENGTH <- 0.5
FACTOR_IDIO_VAR <- 0.25

###############################################################################
# 1. Log helper
###############################################################################

LOG_FILE <- file.path(OUTDIR, "run_log.txt")
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ...)
  cat(msg, "\n")
  cat(msg, "\n", file = LOG_FILE, append = TRUE)
}

###############################################################################
# 2. Matrix helpers
###############################################################################

symmetrize <- function(A) 0.5 * (A + t(A))

make_pd <- function(A, ridge = RIDGE) {
  A <- symmetrize(A)
  ev <- min(eigen(A, symmetric = TRUE, only.values = TRUE)$values)
  if (!is.finite(ev) || ev <= ridge)
    A <- A + (ridge - ev + ridge) * diag(nrow(A))
  A
}

safe_solve <- function(A, b = NULL, ridge = RIDGE) {
  A <- as.matrix(A)
  if (is.null(b)) b <- diag(nrow(A))
  tryCatch(solve(A, b),
           error = function(e) solve(A + ridge * diag(nrow(A)), b))
}

optimal_w <- function(Sigma, ridge = RIDGE) {
  Sigma <- make_pd(Sigma, ridge = ridge)
  m <- nrow(Sigma)
  one <- rep(1, m)
  s1 <- as.numeric(safe_solve(Sigma, one, ridge = ridge))
  d  <- sum(s1)
  if (!is.finite(d) || abs(d) < 1e-12)
    return(list(w = rep(1/m, m), Q = as.numeric(t(rep(1/m,m)) %*% Sigma %*% rep(1/m,m))))
  w <- s1 / d
  list(w = as.numeric(w), Q = as.numeric(1 / d))
}

###############################################################################
# 3. Covariance designs
###############################################################################

toeplitz_cov <- function(m, rho) {
  idx <- seq_len(m)
  rho^abs(outer(idx, idx, "-"))
}

factor_cov <- function(m, r, strength, idio_var, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Lambda <- matrix(rnorm(m * r), m, r)
  Sigma  <- Lambda %*% (strength * diag(r)) %*% t(Lambda) + idio_var * diag(m)
  Sigma  <- make_pd(Sigma)
  d      <- sqrt(diag(Sigma))
  Sigma / tcrossprod(d, d)
}

make_sigma_e <- function(type, m, seed = NULL) {
  if (type == "toeplitz") return(toeplitz_cov(m, RHO))
  factor_cov(m, FACTOR_R, FACTOR_STRENGTH, FACTOR_IDIO_VAR, seed = seed)
}

###############################################################################
# 4. Beta design
###############################################################################

make_beta <- function(m, Sigma_e) {
  b_raw <- c(rep(1, min(10, m)), rep(0, max(0, m - 10)))  # B4 sparse
  Sigma_x <- matrix(1, m, m) + Sigma_e                    # innovation covariance
  V_raw   <- as.numeric(t(b_raw) %*% Sigma_x %*% b_raw)
  c_scale <- sqrt((R2_TARGET / (1 - R2_TARGET)) * SIGMA_U^2 / V_raw)
  c_scale * b_raw
}

###############################################################################
# 5. DGP
###############################################################################

gen_ar1 <- function(n, phi, sd = 1) {
  x <- numeric(n)
  for (t in 2:n) x[t] <- phi * x[t-1] + rnorm(1, sd = sd)
  x
}

gen_var1 <- function(n, phi, Sigma) {
  m <- nrow(Sigma)
  C <- chol(make_pd(Sigma))
  Z <- matrix(rnorm(n * m), n, m) %*% C
  E <- matrix(0, n, m)
  for (t in 2:n) E[t,] <- phi * E[t-1,] + Z[t,]
  E
}

simulate_path <- function(T_len, beta, Sigma_e) {
  m   <- length(beta)
  F_t <- gen_ar1(T_len, PHI_F)
  Eta <- gen_var1(T_len, PHI_ETA, Sigma_e)
  X   <- Eta + matrix(F_t, T_len, m)
  u   <- rnorm(T_len, sd = SIGMA_U); u[1] <- 0
  Y   <- c(NA_real_, as.numeric(X[-T_len,] %*% beta) + u[-1])
  list(X = X, Y = Y)
}

###############################################################################
# 6. RSM fast inner loop (mirrors rsm_fast_one_origin in endogenous script)
###############################################################################

rsm_one_origin <- function(Sigma_hat, f_now, z_grid, K_RS, T_tr) {
  m   <- ncol(Sigma_hat)
  one <- rep(1, m)
  n_z <- length(z_grid)

  pred_rs  <- numeric(n_z)
  Q_rs_bag <- numeric(n_z)
  L_rs_bag <- numeric(n_z)

  Sfull <- Sigma_hat + RIDGE * diag(m)
  s1    <- tryCatch(solve(Sfull, one), error = function(e) rep(1/m, m))
  d0    <- sum(s1)
  wfull <- if (is.finite(d0) && abs(d0) > 1e-12) s1/d0 else rep(1/m, m)
  Qfull <- as.numeric(t(wfull) %*% Sigma_hat %*% wfull)
  Lfull <- if (T_tr > m) Qfull * (T_tr - 1) / (T_tr - m) else NA_real_

  for (ii in seq_along(z_grid)) {
    z <- z_grid[ii]
    if (z == m) {
      pred_rs[ii]  <- sum(wfull * f_now)
      Q_rs_bag[ii] <- Qfull
      L_rs_bag[ii] <- Lfull
      next
    }
    infl  <- (T_tr - 1) / (T_tr - z)
    w_acc <- numeric(m)
    Q_acc <- 0
    for (bb in seq_len(K_RS)) {
      idx   <- sample.int(m, z, replace = FALSE)
      SS    <- Sigma_hat[idx, idx, drop = FALSE] + RIDGE * diag(z)
      s1z   <- tryCatch(solve(SS, rep(1,z)), error = function(e) rep(1/z, z))
      dz    <- sum(s1z)
      wS    <- if (is.finite(dz) && abs(dz) > 1e-12) s1z/dz else rep(1/z, z)
      w_acc[idx] <- w_acc[idx] + wS
      Q_acc <- Q_acc + as.numeric(t(wS) %*% Sigma_hat[idx, idx] %*% wS)
    }
    w_bar <- w_acc / K_RS
    sw    <- sum(w_bar)
    if (is.finite(sw) && sw > 1e-12) w_bar <- w_bar / sw
    Q_bag        <- as.numeric(t(w_bar) %*% Sigma_hat %*% w_bar)
    pred_rs[ii]  <- sum(w_bar * f_now)
    Q_rs_bag[ii] <- Q_bag
    L_rs_bag[ii] <- Q_bag * infl
  }

  list(pred_rs = pred_rs, Q_rs_bag = Q_rs_bag, L_rs_bag = L_rs_bag)
}

###############################################################################
# 7. Individual forecast helpers
###############################################################################

est_gammas <- function(X_tr, y_tr) {
  den <- colSums(X_tr^2)
  num <- as.numeric(crossprod(X_tr, y_tr))
  ifelse(abs(den) < 1e-12, 0, num / den)
}

est_sigma_hat <- function(X_tr, y_tr, gamma) {
  Uhat <- matrix(y_tr, nrow(X_tr), ncol(X_tr)) - sweep(X_tr, 2, gamma, `*`)
  make_pd(crossprod(Uhat) / nrow(Uhat))
}

###############################################################################
# 8. One replication
###############################################################################

one_rep <- function(T_len, m, beta, Sigma_e, z_grid) {
  dat   <- simulate_path(T_len, beta, Sigma_e)
  X     <- dat$X; Y <- dat$Y
  init  <- T_len - OOS_LEN
  steps <- OOS_LEN - 1
  t_vec <- init + seq_len(steps)

  bp  <- numeric(steps)
  rsp <- matrix(NA_real_, steps, length(z_grid))
  Qbm <- matrix(NA_real_, steps, length(z_grid))
  Lbm <- matrix(NA_real_, steps, length(z_grid))

  for (hh in seq_along(t_vec)) {
    T_tr <- t_vec[hh] - 1
    if (T_tr < max(z_grid) + 5L) next
    X_tr  <- X[1:T_tr,, drop = FALSE]
    y_tr  <- Y[2:(T_tr+1)]
    x_now <- X[t_vec[hh],]
    gam   <- est_gammas(X_tr, y_tr)
    f_now <- as.numeric(gam * x_now)
    bp[hh] <- mean(f_now)
    Sh <- est_sigma_hat(X_tr, y_tr, gam)
    res <- rsm_one_origin(Sh, f_now, z_grid, K_RS, T_tr)
    rsp[hh,]  <- res$pred_rs
    Qbm[hh,] <- res$Q_rs_bag
    Lbm[hh,] <- res$L_rs_bag
  }

  y_oos  <- Y[t_vec + 1]
  e_base <- y_oos - bp
  rows <- vector("list", length(z_grid))
  for (ii in seq_along(z_grid)) {
    e_rs <- y_oos - rsp[,ii]
    rows[[ii]] <- data.frame(
      z           = z_grid[ii],
      se2_base    = sum(e_base^2, na.rm = TRUE),
      cnt_base    = sum(is.finite(e_base)),
      se2_rs      = sum(e_rs^2, na.rm = TRUE),
      cnt_rs      = sum(is.finite(e_rs)),
      Q_rs_bag    = mean(Qbm[,ii], na.rm = TRUE),
      L_rs_bag    = mean(Lbm[,ii], na.rm = TRUE)
    )
  }
  do.call(rbind, rows)
}

###############################################################################
# 9. Main m-scaling loop
###############################################################################

log_msg("m-scaling experiment started.")
log_msg("T_FIXED=", T_FIXED, " N_REP=", N_REP, " K_RS=", K_RS,
        " OOS_LEN=", OOS_LEN)
log_msg("M_GRID=", paste(M_GRID, collapse=","))

summary_all <- data.frame()

for (stype in SIGMA_TYPES) {
  log_msg("=== Sigma type: ", stype, " ===")
  seed_sigma <- if (stype == "toeplitz") 42001L else 42002L

  for (m in M_GRID) {
    log_msg("  m=", m, " T=", T_FIXED, " starting...")

    set.seed(seed_sigma + m)
    Sigma_e <- make_sigma_e(stype, m, seed = seed_sigma + m)
    beta    <- make_beta(m, Sigma_e)

    # Adaptive grid step: keep ~25-30 z values regardless of m
    k_by  <- max(2L, as.integer(floor(m / 25)))
    z_grid <- seq(2L, m, by = k_by)
    if (tail(z_grid, 1) != m) z_grid <- c(z_grid, m)

    # Kan-Zhou inflation factor (T-2)/(T-m-2) for reference
    kz_factor <- if (T_FIXED - m - 2 > 0) (T_FIXED - 2) / (T_FIXED - m - 2) else Inf

    rep_list <- vector("list", N_REP)
    for (rr in seq_len(N_REP)) {
      if (rr %% 100 == 0)
        log_msg("    m=", m, " stype=", stype, " rep ", rr, "/", N_REP)
      set.seed(GLOBAL_SEED + 1e6 * (match(stype, SIGMA_TYPES)-1) +
                 1e4 * match(m, M_GRID) + rr)
      rep_list[[rr]] <- one_rep(T_FIXED, m, beta, Sigma_e, z_grid)
    }

    all_r <- do.call(rbind, rep_list)
    agg <- do.call(rbind, lapply(z_grid, function(z0) {
      dd <- all_r[all_r$z == z0,]
      rmsfe_b  <- sqrt(sum(dd$se2_base) / sum(dd$cnt_base))
      rmsfe_rs <- sqrt(sum(dd$se2_rs)   / sum(dd$cnt_rs))
      data.frame(
        sigma_type = stype, m = m, T_len = T_FIXED, z = z0,
        rmsfe_base = rmsfe_b, rmsfe_rs = rmsfe_rs,
        rel_rmsfe  = rmsfe_rs / rmsfe_b,
        Q_rs_bag   = mean(dd$Q_rs_bag, na.rm = TRUE),
        L_rs_bag   = mean(dd$L_rs_bag, na.rm = TRUE)
      )
    }))

    # Save per-(m, stype) curve
    fname <- paste0("m_scaling_", stype, "_m", m, "_T", T_FIXED, "_by_z.csv")
    write.csv(agg, file.path(OUTDIR, fname), row.names = FALSE)

    best_z_rmsfe <- agg$z[which.min(agg$rmsfe_rs)]
    best_z_L     <- agg$z[which.min(agg$L_rs_bag)]
    rel_best     <- min(agg$rel_rmsfe, na.rm = TRUE)

    log_msg("  m=", m, " stype=", stype,
            " k*(RMSFE)=", best_z_rmsfe,
            " k*(L)=", best_z_L,
            " gain=", round((1 - rel_best)*100, 1), "%",
            " KZ_factor=", round(kz_factor, 3))

    summary_all <- rbind(summary_all, data.frame(
      sigma_type    = stype,
      m             = m,
      T_len         = T_FIXED,
      kz_factor     = round(kz_factor, 4),
      best_z_rmsfe  = best_z_rmsfe,
      best_z_L_bag  = best_z_L,
      rel_rmsfe_best = rel_best,
      min_rmsfe_rs  = min(agg$rmsfe_rs, na.rm = TRUE),
      rmsfe_base    = agg$rmsfe_base[1]
    ))

    write.csv(summary_all,
              file.path(OUTDIR, "m_scaling_summary.csv"),
              row.names = FALSE)
  }
}

###############################################################################
# 10. Plots and log-log regression
###############################################################################

log_msg("Generating plots and log-log regressions...")

for (stype in SIGMA_TYPES) {
  sub <- summary_all[summary_all$sigma_type == stype,]
  sub <- sub[order(sub$m),]

  # --- Plot A: k* vs m ---
  png(file.path(OUTDIR, paste0("kstar_vs_m_", stype, ".png")),
      width = 1100, height = 600)
  par(mfrow = c(1,2))

  plot(sub$m, sub$best_z_rmsfe, type = "b", pch = 19, lwd = 2,
       col = "#0072B2", xlab = "m (pool size)", ylab = "k*",
       main = paste("k* vs m [", stype, "], T=", T_FIXED))
  lines(sub$m, sub$best_z_L_bag, type = "b", pch = 17, lty = 2,
        lwd = 2, col = "#D55E00")
  legend("topleft", c("k*(RMSFE)", "k*(L_bag)"),
         col = c("#0072B2","#D55E00"), lty = c(1,2),
         pch = c(19,17), lwd = 2, bty = "n")

  # Panel B: log-log with Kan-Zhou adjusted fit
  valid <- sub$m < T_FIXED - 2 & sub$best_z_rmsfe > 0
  if (sum(valid) >= 3) {
    lm1 <- lm(log(best_z_rmsfe) ~ log(m), data = sub[valid,])
    beta_m <- round(coef(lm1)[2], 3)
    plot(log(sub$m[valid]), log(sub$best_z_rmsfe[valid]),
         type = "b", pch = 19, lwd = 2, col = "#0072B2",
         xlab = "log m", ylab = "log k*",
         main = paste0("log-log [", stype, "]  slope=", beta_m))
    abline(lm1, lty = 2, col = "#0072B2")
  } else {
    plot.new()
    title("Insufficient data for log-log")
  }
  par(mfrow = c(1,1))
  dev.off()

  # --- Log-log regression to file ---
  sink(file.path(OUTDIR, paste0("m_scaling_loglog_", stype, ".txt")))
  cat("Sigma type:", stype, "  T=", T_FIXED, "\n\n")
  if (sum(valid) >= 3) {
    lm_rmsfe <- lm(log(best_z_rmsfe) ~ log(m), data = sub[valid,])
    lm_L     <- lm(log(best_z_L_bag) ~ log(m),
                   data = sub[valid & sub$best_z_L_bag > 0,])
    cat("k*(RMSFE) ~ m^beta:\n"); print(summary(lm_rmsfe))
    cat("k*(L_bag) ~ m^beta:\n"); print(summary(lm_L))
    cat("\nKan-Zhou prediction: k*_endo/k*_exo = ((T-2)/(T-m-2))^{1/(alpha+1)}\n")
    cat("Observed kz_factor by m:\n")
    print(sub[, c("m","kz_factor","best_z_rmsfe","best_z_L_bag")])
  }
  sink()
}

log_msg("m-scaling experiment complete. Results in: ", OUTDIR)
cat("Done.\n")
