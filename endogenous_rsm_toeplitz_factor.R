#!/usr/bin/env Rscript

###############################################################################
# Endogenous-forecast Monte Carlo for random-subset forecast combination
#
# This is an R recreation of the uploaded C program FInal_sum_up_one.c.
# It keeps the endogenous-predictor timing:
#
#   x_t = common_factor_t * 1_m + eta_t,
#   y_{t+1} = beta' x_t + u_{t+1},
#
# and each individual forecast is estimated from the univariate predictive
# regression:
#
#   y_{t+1} = gamma_i x_{i,t} + error_{i,t+1}.
#
# The final forecast pool consists of the m fitted individual forecasts
# f_{i,t+1|t} = gamma_i x_{i,t}. Random-subset forecast combination then
# computes sum-to-one minimum-variance weights inside random subsets using the
# rolling training-window covariance matrix of individual forecast errors.
#
# New relative to the C file:
#   1. The innovation covariance Sigma_e of eta_t can be either
#        - Toeplitz: Sigma_e[i,j] = sigma_eta^2 * rho^|i-j|;
#        - Factor:   Sigma_e = Lambda_e Omega_e Lambda_e' + D_e.
#   2. You can iterate over m, T, R2, beta designs, phi values, rho values,
#      factor dimensions, and loading scenarios.
#   3. The script reports both bagged RSM and average-subset RSM, plus the
#      covariance-based Q and inflation-adjusted L diagnostics analogous to
#      the C output.
#
# The script uses only base R.
###############################################################################

rm(list = ls())

###############################################################################
# 0. User controls
###############################################################################

RUN_PROFILE <- "run"
# Options:
#   "fast"  : quick pilot run (logic / plot checks only).
#   "run"   : production run saved to RSM_endogenous_results.
#   "paper" : full paper Monte Carlo (very slow).

GLOBAL_SEED <- 778
set.seed(GLOBAL_SEED)

OUTDIR <- Sys.getenv("RSM_OUTDIR",
                     unset = file.path(getwd(), "RSM_endogenous_results"))
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

if (RUN_PROFILE == "fast") {
  N_REP        <- 10
  K_RS         <- 100
  K_BY         <- 2
  OOS_LEN      <- 20
  SAVE_PLOTS   <- TRUE

  M_LIST       <- c(60)
  T_LEN_LIST   <- c(120)
  R2_LIST      <- c(0.30, 0.1)
  BETA_DESIGNS <- c("sparse10")

  SIGMA_TYPES  <- c("toeplitz", "factor")
  RHO_LIST     <- c(0.0, 0.8)

  FACTOR_R_LIST          <- c(2)
  FACTOR_LOADING_SCEN    <- c("random", "nearly_common")
  FACTOR_STRENGTH_LIST   <- c(1.0)
  FACTOR_IDIO_VAR_LIST   <- c(0.5)

  PHI_GRID <- data.frame(
    phi_f   = c(0.3, 0.7),
    phi_eta = c(0.7, 0.3),
    phi_u   = c(0.0, 0.0)
  )

} else if (RUN_PROFILE == "run") {
  # Production run: enough reps for stable empirical results, feasible runtime.
  N_REP        <- 200
  K_RS         <- 200
  K_BY         <- 3
  OOS_LEN      <- 30
  SAVE_PLOTS   <- TRUE

  M_LIST       <- c(50)
  T_LEN_LIST   <- c(200)
  R2_LIST      <- c(0.3)
  BETA_DESIGNS <- c("B2", "B4")

  SIGMA_TYPES  <- c("toeplitz", "factor")
  RHO_LIST     <- c(0.8)

  FACTOR_R_LIST          <- c(2)
  FACTOR_LOADING_SCEN    <- c("random")
  FACTOR_STRENGTH_LIST   <- c(0.5)
  FACTOR_IDIO_VAR_LIST   <- c(0.25)

  PHI_GRID <- data.frame(
    phi_f   = c(0.3),
    phi_eta = c(0.7),
    phi_u   = c(0.0)
  )

} else if (RUN_PROFILE == "paper") {
  N_REP        <- 500
  K_RS         <- 500
  K_BY         <- 3
  OOS_LEN      <- 50
  SAVE_PLOTS   <- TRUE

  M_LIST       <- c(50)
  T_LEN_LIST   <- c(200)
  R2_LIST      <- c(0.3)
  BETA_DESIGNS <- c("B2", "B4")

  SIGMA_TYPES  <- c("toeplitz", "factor")
  RHO_LIST     <- c(0.8)

  FACTOR_R_LIST          <- c(2)
  FACTOR_LOADING_SCEN    <- c("random")
  FACTOR_STRENGTH_LIST   <- c(0.5)
  FACTOR_IDIO_VAR_LIST   <- c(0.25)

  PHI_GRID <- data.frame(
    phi_f   = c(0.3),
    phi_eta = c(0.7),
    phi_u   = c(0.0)
  )
} else {
  stop("RUN_PROFILE must be 'fast', 'run', or 'paper'.")
}

###############################################################################
# T-scaling experiment: fix Sigma_e, vary T.
# Purpose: observe how empirical optimal k changes with T while holding
# everything else constant.  This is the endogenous analogue of the
# exogenous scaling_T_m experiment.
###############################################################################

T_SCALING_CONFIG <- list(
  m             = 50,
  R2_target     = 0.3,
  beta_design   = "B4",          # sparse beta: cleaner signal
  phi_f         = 0.3,
  phi_eta       = 0.7,
  phi_u         = 0.0,
  T_len_grid    = c(100, 150, 200, 300, 500, 800),
  OOS_LEN_scale = 20,            # fixed OOS length so comparisons are clean
  N_REP_scale   = if (RUN_PROFILE == "fast") 10L else 120L,
  K_RS_scale    = if (RUN_PROFILE == "fast") 80L  else 150L,
  K_BY_scale    = 4L,
  # Two Sigma_e types, each with a fixed seed so Sigma_e is identical
  # across all T values within a type.
  sigma_specs   = list(
    list(type = "toeplitz", rho = 0.8,
         factor_r = NA, loading = NA, strength = NA, idio = NA,
         seed_sigma = 42001L),
    list(type = "factor",   rho = NA,
         factor_r = 2L, loading = "random", strength = 0.5, idio = 0.25,
         seed_sigma = 42002L)
  )
)

# Target-noise variance in y_{t+1} = beta' x_t + u_{t+1}.
SIGMA_U <- 1.0

# The uploaded C file scales beta using the innovation covariance of x_t, not the
# stationary covariance of the AR(1) process. Keep this default to reproduce it.
# Use "stationary" if you want target R2 defined from unconditional variances.
BETA_SCALE_USING <- "innovation"   # choices: "innovation", "stationary"

RIDGE <- 1e-8

PLOT_COL <- c(
  bagged        = "#0072B2",
  subset        = "#D55E00",
  base          = "#000000",
  full          = "#E69F00",
  q_bag         = "#56B4E9",
  q_sub         = "#CC79A7",
  q_ols         = "#009E73"
)

###############################################################################
# 1. Numerical helpers
###############################################################################

# Population MC for the general Sigma_e case.
# Estimates Q_sub(k), Q_bag_inf(k), Q_bag_G(k) for an arbitrary positive-
# definite forecast-error covariance matrix.  No common-loading or small-
# heterogeneity assumptions are needed.
q_population_mc_general <- function(Sigma_e, k_grid, n_subsets = 2000,
                                    G_bag = NULL) {
  Sigma_e <- as.matrix(Sigma_e)
  m   <- ncol(Sigma_e)
  one <- rep(1, m)

  # Full optimal weights and risk.
  Sinv_one <- as.numeric(tryCatch(solve(Sigma_e, one),
                                  error = function(e) rep(1/m, m)))
  denom    <- sum(Sinv_one)
  w_opt    <- if (is.finite(denom) && abs(denom) > 1e-12)
                Sinv_one / denom else rep(1/m, m)
  w_equal  <- rep(1/m, m)
  Q_opt    <- as.numeric(t(w_opt)   %*% Sigma_e %*% w_opt)
  Q_ave    <- as.numeric(t(w_equal) %*% Sigma_e %*% w_equal)

  out <- data.frame(
    k                    = k_grid,
    Q_subset_mc          = NA_real_,
    Q_bag_inf_mc         = NA_real_,
    Q_bag_G_mc           = NA_real_,
    Q_opt                = Q_opt,
    Q_ave                = Q_ave,
    R_weight_exact_mc    = NA_real_,
    R_bag_inf_mc         = NA_real_,
    R_weight_identity_gap = NA_real_,
    w_l2_mc              = NA_real_,
    w_max_abs_mc         = NA_real_
  )

  for (ii in seq_along(k_grid)) {
    k     <- as.integer(k_grid[ii])
    w_sum <- rep(0, m)
    q_sub_vals <- numeric(n_subsets)

    for (bb in seq_len(n_subsets)) {
      idx       <- sample.int(m, k, replace = FALSE)
      Sigma_S   <- Sigma_e[idx, idx, drop = FALSE]
      one_S     <- rep(1, k)
      Sinv1_S   <- as.numeric(tryCatch(solve(Sigma_S, one_S),
                                       error = function(e) {
                                         solve(Sigma_S + 1e-8*diag(k), one_S)
                                       }))
      ds        <- sum(Sinv1_S)
      w_S       <- if (is.finite(ds) && abs(ds) > 1e-12)
                     Sinv1_S / ds else rep(1/k, k)
      w_emb     <- rep(0, m); w_emb[idx] <- w_S
      w_sum     <- w_sum + w_emb
      q_sub_vals[bb] <- as.numeric(t(w_emb) %*% Sigma_e %*% w_emb)
    }

    w_bar    <- w_sum / n_subsets
    diff     <- w_bar - w_opt
    Q_sub    <- mean(q_sub_vals)
    Q_bag_inf <- as.numeric(t(w_bar) %*% Sigma_e %*% w_bar)
    R_weight <- as.numeric(t(diff)  %*% Sigma_e %*% diff)
    Q_bag_G  <- if (!is.null(G_bag) && is.finite(G_bag) && G_bag > 0)
                  Q_bag_inf + (Q_sub - Q_bag_inf) / G_bag else NA_real_

    out$Q_subset_mc[ii]           <- Q_sub
    out$Q_bag_inf_mc[ii]          <- Q_bag_inf
    out$Q_bag_G_mc[ii]            <- Q_bag_G
    out$R_weight_exact_mc[ii]     <- R_weight
    out$R_bag_inf_mc[ii]          <- Q_bag_inf - Q_opt
    out$R_weight_identity_gap[ii] <- (Q_bag_inf - Q_opt) - R_weight
    out$w_l2_mc[ii]               <- sqrt(sum(diff^2))
    out$w_max_abs_mc[ii]          <- max(abs(diff))
  }
  out
}

safe_solve <- function(A, b = NULL, ridge = RIDGE) {
  A <- as.matrix(A)
  if (is.null(b)) b <- diag(nrow(A))
  out <- tryCatch(
    solve(A, b),
    error = function(e) solve(A + ridge * diag(nrow(A)), b)
  )
  out
}

symmetrize <- function(A) {
  0.5 * (A + t(A))
}

make_pd <- function(A, ridge = RIDGE) {
  A <- symmetrize(A)
  ev <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
  min_ev <- min(ev)
  if (!is.finite(min_ev) || min_ev <= ridge) {
    A <- A + (ridge - min_ev + ridge) * diag(nrow(A))
  }
  A
}

optimal_sum_to_one_weights <- function(Sigma, ridge = RIDGE) {
  Sigma <- make_pd(Sigma, ridge = ridge)
  m <- nrow(Sigma)
  one <- rep(1, m)
  Sinv_one <- as.numeric(safe_solve(Sigma, one, ridge = ridge))
  denom <- sum(Sinv_one)
  if (!is.finite(denom) || abs(denom) < 1e-12) {
    w <- rep(1 / m, m)
    Q <- as.numeric(t(w) %*% Sigma %*% w)
  } else {
    w <- Sinv_one / denom
    Q <- as.numeric(1 / denom)
  }
  list(w = as.numeric(w), Q = as.numeric(Q))
}

mspe <- function(err) mean(err^2, na.rm = TRUE)

safe_label <- function(...) {
  x <- paste(..., sep = "_")
  x <- gsub("[^A-Za-z0-9_.-]", "_", x)
  x <- gsub("_+", "_", x)
  substr(x, 1, 180)
}

###############################################################################
# 2. Covariance designs: Toeplitz and factor Sigma_e
###############################################################################

toeplitz_cov <- function(m, rho, sigma_eta2 = 1.0) {
  idx <- seq_len(m)
  sigma_eta2 * rho^abs(outer(idx, idx, "-"))
}

make_factor_loadings <- function(m, r, scenario = "random", seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  if (scenario == "random") {
    Lambda <- matrix(rnorm(m * r, sd = 1.0), nrow = m, ncol = r)
  } else if (scenario == "nearly_common") {
    common_vec <- rep(1 / sqrt(r), r)
    Lambda <- matrix(rep(common_vec, each = m), nrow = m, ncol = r)
    Lambda <- Lambda + matrix(rnorm(m * r, sd = 0.15), nrow = m, ncol = r)
  } else if (scenario == "block") {
    Lambda <- matrix(0, nrow = m, ncol = r)
    groups <- split(seq_len(m), cut(seq_len(m), breaks = r, labels = FALSE))
    for (j in seq_len(r)) Lambda[groups[[j]], j] <- 1
    Lambda <- Lambda + matrix(rnorm(m * r, sd = 0.10), nrow = m, ncol = r)
  } else if (scenario == "decay") {
    Lambda <- outer(seq_len(m), seq_len(r), function(i, j) exp(-abs(i - j * m / r) / (m / r)))
    Lambda <- scale(Lambda, center = TRUE, scale = FALSE)
  } else {
    stop("Unknown factor loading scenario: ", scenario)
  }

  as.matrix(Lambda)
}

factor_cov <- function(m, r = 2, scenario = "random", factor_strength = 1.0,
                       idio_var = 0.5, scale_to_unit_diag = TRUE, seed = NULL) {
  Lambda <- make_factor_loadings(m, r, scenario = scenario, seed = seed)
  Omega <- diag(factor_strength, r)
  D <- diag(idio_var, m)
  Sigma <- Lambda %*% Omega %*% t(Lambda) + D
  Sigma <- make_pd(Sigma)

  if (scale_to_unit_diag) {
    d <- sqrt(diag(Sigma))
    Sigma <- Sigma / tcrossprod(d, d)
  }

  attr(Sigma, "Lambda") <- Lambda
  attr(Sigma, "Omega") <- Omega
  attr(Sigma, "idio_var") <- idio_var
  Sigma
}

make_sigma_e <- function(type, m, rho = 0.0, factor_r = 2,
                         factor_loading_scenario = "random",
                         factor_strength = 1.0, factor_idio_var = 0.5,
                         seed = NULL) {
  if (type == "toeplitz") {
    toeplitz_cov(m, rho = rho, sigma_eta2 = 1.0)
  } else if (type == "factor") {
    factor_cov(
      m = m,
      r = factor_r,
      scenario = factor_loading_scenario,
      factor_strength = factor_strength,
      idio_var = factor_idio_var,
      scale_to_unit_diag = TRUE,
      seed = seed
    )
  } else {
    stop("Unknown Sigma_e type: ", type)
  }
}

###############################################################################
# 3. Beta designs and R2 scaling
###############################################################################

beta_raw_design <- function(m, design = "sparse10") {
  i <- seq_len(m)

  if (design == "B1") {
    b <- rep(1, m)
  } else if (design == "B2") {
    b <- rev(i)
  } else if (design == "B3") {
    b <- 1 / i
  } else if (design == "B4" || design == "sparse10") {
    b <- c(rep(1, min(10, m)), rep(0, max(0, m - 10)))
  } else if (design == "B5") {
    b <- exp(-i)
  } else if (design == "linear_sparse") {
    b <- pmax(0, 1 - (i - 1) / max(1, floor(m / 2)))
  } else {
    stop("Unknown beta design: ", design)
  }

  as.numeric(b)
}

scale_beta_to_R2 <- function(beta_raw, Sigma_e, R2_target, sigma_u = SIGMA_U,
                             phi_f = 0.0, phi_eta = 0.0,
                             scale_using = BETA_SCALE_USING) {
  m <- length(beta_raw)
  one_cov <- matrix(1, m, m)

  if (scale_using == "innovation") {
    Sigma_x <- one_cov + Sigma_e
  } else if (scale_using == "stationary") {
    Sigma_x <- one_cov / (1 - phi_f^2) + Sigma_e / (1 - phi_eta^2)
  } else {
    stop("scale_using must be 'innovation' or 'stationary'.")
  }

  V_raw <- as.numeric(t(beta_raw) %*% Sigma_x %*% beta_raw)
  if (!is.finite(V_raw) || V_raw <= 0) stop("Nonpositive raw signal variance.")

  c_scale <- sqrt((R2_target / (1 - R2_target)) * sigma_u^2 / V_raw)
  as.numeric(c_scale * beta_raw)
}

###############################################################################
# 4. Time-series DGP
###############################################################################

gen_ar1_scalar <- function(T_len, phi, innov_var = 1.0) {
  x <- numeric(T_len)
  if (T_len >= 2) {
    for (tt in 2:T_len) {
      x[tt] <- phi * x[tt - 1] + rnorm(1, sd = sqrt(innov_var))
    }
  }
  x
}

gen_vector_ar1 <- function(T_len, phi, Sigma_innov) {
  m <- nrow(Sigma_innov)
  C <- tryCatch(
    chol(make_pd(Sigma_innov)),
    error = function(e) chol(make_pd(Sigma_innov, ridge = 1e-6))
  )
  Z <- matrix(rnorm(T_len * m), nrow = T_len, ncol = m) %*% C
  E <- matrix(0, nrow = T_len, ncol = m)
  if (T_len >= 2) {
    for (tt in 2:T_len) {
      E[tt, ] <- phi * E[tt - 1, ] + Z[tt, ]
    }
  }
  E
}

simulate_endogenous_path <- function(T_len, beta_true, Sigma_e,
                                     phi_f = 0.3, phi_eta = 0.7,
                                     phi_u = 0.0, sigma_u = SIGMA_U) {
  m <- length(beta_true)

  F_common <- gen_ar1_scalar(T_len, phi = phi_f, innov_var = 1.0)
  Eta <- gen_vector_ar1(T_len, phi = phi_eta, Sigma_innov = Sigma_e)
  X <- Eta + matrix(F_common, nrow = T_len, ncol = m)

  if (abs(phi_u) > 1e-12) {
    u <- gen_ar1_scalar(T_len, phi = phi_u, innov_var = sigma_u^2)
  } else {
    u <- rnorm(T_len, sd = sigma_u)
    u[1] <- 0
  }

  Y <- rep(NA_real_, T_len)
  if (T_len >= 2) {
    for (tt in 2:T_len) {
      Y[tt] <- sum(beta_true * X[tt - 1, ]) + u[tt]
    }
  }

  list(X = X, Y = Y, F_common = F_common, Eta = Eta, u = u)
}

###############################################################################
# 5. Forecast estimation and RSM evaluation
###############################################################################

estimate_univariate_gammas <- function(X_train, y_train) {
  den <- colSums(X_train^2)
  num <- as.numeric(crossprod(X_train, y_train))
  gamma <- ifelse(abs(den) < 1e-12, 0, num / den)
  as.numeric(gamma)
}

forecast_error_covariance <- function(X_train, y_train, gamma, ridge = RIDGE) {
  Uhat <- matrix(y_train, nrow = nrow(X_train), ncol = ncol(X_train)) -
    sweep(X_train, 2, gamma, `*`)
  Sigma_u_hat <- crossprod(Uhat) / nrow(Uhat)
  make_pd(Sigma_u_hat, ridge = ridge)
}

evaluate_one_origin_all_z <- function(X, Y, t_now0, z_grid, K_RS, ridge = RIDGE) {
  # t_now0 is the 0-based forecasted target time, matching the C code.
  # In R indexes:
  #   x_now = X[t_now0, ] predicts Y[t_now0 + 1].
  #   training uses X[1:T_tr, ] with y_train = Y[2:(T_tr+1)].
  T_tr <- t_now0 - 1
  if (T_tr < 5) stop("Training window too short.")

  X_train <- X[1:T_tr, , drop = FALSE]
  y_train <- Y[2:(T_tr + 1)]
  x_now <- X[t_now0, ]

  gamma <- estimate_univariate_gammas(X_train, y_train)
  f_now <- as.numeric(gamma * x_now)
  base_pred <- mean(f_now)

  Sigma_u_hat <- forecast_error_covariance(X_train, y_train, gamma, ridge = ridge)
  m <- ncol(X)

  full <- optimal_sum_to_one_weights(Sigma_u_hat, ridge = ridge)
  full_pred <- sum(full$w * f_now)
  Q_ols <- full$Q
  # Inflation factor consistent with the paper: (T-1)/(T-m).
  # Requires T_tr > m; guard against singularity.
  L_ols <- if (T_tr > m) Q_ols * (T_tr - 1) / (T_tr - m) else NA_real_

  out <- data.frame(
    z = z_grid,
    pred_rs = NA_real_,
    Q_ols = Q_ols,
    Q_rs_bag = NA_real_,
    Q_rs_subavg = NA_real_,
    L_ols = L_ols,
    L_rs_bag = NA_real_,
    L_rs_subavg = NA_real_
  )

  for (ii in seq_along(z_grid)) {
    z <- z_grid[ii]

    if (z == m) {
      out$pred_rs[ii] <- full_pred
      out$Q_rs_bag[ii] <- Q_ols
      out$Q_rs_subavg[ii] <- Q_ols
      out$L_rs_bag[ii] <- L_ols
      out$L_rs_subavg[ii] <- L_ols
      next
    }

    acc_pred <- 0
    acc_Q_sub <- 0
    w_bar <- numeric(m)
    succ <- 0L

    for (bb in seq_len(K_RS)) {
      idx <- sample.int(m, z, replace = FALSE)
      Sigma_S <- Sigma_u_hat[idx, idx, drop = FALSE]
      ow <- optimal_sum_to_one_weights(Sigma_S, ridge = ridge)
      wS <- ow$w

      acc_pred <- acc_pred + sum(wS * f_now[idx])
      acc_Q_sub <- acc_Q_sub + ow$Q
      w_bar[idx] <- w_bar[idx] + wS
      succ <- succ + 1L
    }

    if (succ > 0) {
      w_bar <- w_bar / succ
      sw <- sum(w_bar)
      if (is.finite(sw) && abs(sw) > 1e-12) w_bar <- w_bar / sw

      Q_bag <- as.numeric(t(w_bar) %*% Sigma_u_hat %*% w_bar)
      Q_sub <- acc_Q_sub / succ
      # Inflation factor consistent with the paper and exogenous code: (T-1)/(T-k).
      infl <- (T_tr - 1) / (T_tr - z)

      out$pred_rs[ii] <- acc_pred / succ
      out$Q_rs_bag[ii] <- Q_bag
      out$Q_rs_subavg[ii] <- Q_sub
      out$L_rs_bag[ii] <- Q_bag * infl
      out$L_rs_subavg[ii] <- Q_sub * infl
    }
  }

  list(base_pred = base_pred, out = out)
}

run_one_rep_all_z <- function(T_len, m, beta_true, Sigma_e, phi_f, phi_eta, phi_u,
                              z_grid, K_RS, OOS_LEN, sigma_u = SIGMA_U,
                              ridge = RIDGE) {
  dat <- simulate_endogenous_path(
    T_len = T_len,
    beta_true = beta_true,
    Sigma_e = Sigma_e,
    phi_f = phi_f,
    phi_eta = phi_eta,
    phi_u = phi_u,
    sigma_u = sigma_u
  )

  X <- dat$X
  Y <- dat$Y

  init_window <- T_len - OOS_LEN
  steps <- OOS_LEN - 1
  t_now0_vec <- init_window + seq_len(steps)

  base_pred <- numeric(steps)
  rs_pred <- matrix(NA_real_, nrow = steps, ncol = length(z_grid))
  colnames(rs_pred) <- paste0("z", z_grid)

  Q_ols_vec <- L_ols_vec <- numeric(steps)
  Q_bag_mat <- Q_sub_mat <- L_bag_mat <- L_sub_mat <- matrix(NA_real_, steps, length(z_grid))

  for (hh in seq_along(t_now0_vec)) {
    oo <- evaluate_one_origin_all_z(
      X = X,
      Y = Y,
      t_now0 = t_now0_vec[hh],
      z_grid = z_grid,
      K_RS = K_RS,
      ridge = ridge
    )

    base_pred[hh] <- oo$base_pred
    rs_pred[hh, ] <- oo$out$pred_rs
    Q_ols_vec[hh] <- oo$out$Q_ols[1]
    L_ols_vec[hh] <- oo$out$L_ols[1]
    Q_bag_mat[hh, ] <- oo$out$Q_rs_bag
    Q_sub_mat[hh, ] <- oo$out$Q_rs_subavg
    L_bag_mat[hh, ] <- oo$out$L_rs_bag
    L_sub_mat[hh, ] <- oo$out$L_rs_subavg
  }

  y_oos <- Y[t_now0_vec + 1]
  err_base <- y_oos - base_pred

  rows <- vector("list", length(z_grid))
  for (ii in seq_along(z_grid)) {
    err_rs <- y_oos - rs_pred[, ii]

    rows[[ii]] <- data.frame(
      z = z_grid[ii],
      rmsfe_base_rep = sqrt(mean(err_base^2, na.rm = TRUE)),
      rmsfe_rs_rep = sqrt(mean(err_rs^2, na.rm = TRUE)),
      sum_e_base = sum(err_base, na.rm = TRUE),
      sum_e2_base = sum(err_base^2, na.rm = TRUE),
      cnt_base = sum(is.finite(err_base)),
      sum_e_rs = sum(err_rs, na.rm = TRUE),
      sum_e2_rs = sum(err_rs^2, na.rm = TRUE),
      cnt_rs = sum(is.finite(err_rs)),
      Q_ols = mean(Q_ols_vec, na.rm = TRUE),
      Q_rs_bag = mean(Q_bag_mat[, ii], na.rm = TRUE),
      Q_rs_subavg = mean(Q_sub_mat[, ii], na.rm = TRUE),
      L_ols = mean(L_ols_vec, na.rm = TRUE),
      L_rs_bag = mean(L_bag_mat[, ii], na.rm = TRUE),
      L_rs_subavg = mean(L_sub_mat[, ii], na.rm = TRUE)
    )
  }

  do.call(rbind, rows)
}

aggregate_replications <- function(rep_dfs) {
  all <- do.call(rbind, rep_dfs)
  z_vals <- sort(unique(all$z))

  out <- lapply(z_vals, function(z0) {
    dd <- all[all$z == z0, ]

    mse_base <- sum(dd$sum_e2_base) / sum(dd$cnt_base)
    bias_base <- (sum(dd$sum_e_base) / sum(dd$cnt_base))^2
    var_base <- max(0, mse_base - bias_base)

    mse_rs <- sum(dd$sum_e2_rs) / sum(dd$cnt_rs)
    bias_rs <- (sum(dd$sum_e_rs) / sum(dd$cnt_rs))^2
    var_rs <- max(0, mse_rs - bias_rs)

    data.frame(
      z = z0,
      rmsfe_base = mean(dd$rmsfe_base_rep, na.rm = TRUE),
      rmsfe_rs = mean(dd$rmsfe_rs_rep, na.rm = TRUE),
      rel_rmsfe = mean(dd$rmsfe_rs_rep, na.rm = TRUE) / mean(dd$rmsfe_base_rep, na.rm = TRUE),
      rmse_pooled_base = sqrt(mse_base),
      rmse_pooled_rs = sqrt(mse_rs),
      bias_base = bias_base,
      var_base = var_base,
      bias_rs = bias_rs,
      var_rs = var_rs,
      Q_ols = mean(dd$Q_ols, na.rm = TRUE),
      Q_rs_bag = mean(dd$Q_rs_bag, na.rm = TRUE),
      Q_rs_subavg = mean(dd$Q_rs_subavg, na.rm = TRUE),
      L_ols = mean(dd$L_ols, na.rm = TRUE),
      L_rs_bag = mean(dd$L_rs_bag, na.rm = TRUE),
      L_rs_subavg = mean(dd$L_rs_subavg, na.rm = TRUE)
    )
  })

  do.call(rbind, out)
}

###############################################################################
# 6. Scenario grids and plotting
###############################################################################

build_scenarios <- function() {
  scenarios <- list()
  ss <- 1L

  for (m in M_LIST) {
    for (T_len in T_LEN_LIST) {
      for (R2_target in R2_LIST) {
        for (beta_design in BETA_DESIGNS) {
          for (pp in seq_len(nrow(PHI_GRID))) {
            phi_f <- PHI_GRID$phi_f[pp]
            phi_eta <- PHI_GRID$phi_eta[pp]
            phi_u <- PHI_GRID$phi_u[pp]

            if ("toeplitz" %in% SIGMA_TYPES) {
              for (rho in RHO_LIST) {
                scenarios[[ss]] <- list(
                  scenario_id = ss,
                  sigma_type = "toeplitz",
                  m = m,
                  T_len = T_len,
                  R2_target = R2_target,
                  beta_design = beta_design,
                  phi_f = phi_f,
                  phi_eta = phi_eta,
                  phi_u = phi_u,
                  rho = rho,
                  factor_r = NA_integer_,
                  factor_loading_scenario = NA_character_,
                  factor_strength = NA_real_,
                  factor_idio_var = NA_real_
                )
                ss <- ss + 1L
              }
            }

            if ("factor" %in% SIGMA_TYPES) {
              for (factor_r in FACTOR_R_LIST) {
                for (loading_scen in FACTOR_LOADING_SCEN) {
                  for (factor_strength in FACTOR_STRENGTH_LIST) {
                    for (factor_idio_var in FACTOR_IDIO_VAR_LIST) {
                      scenarios[[ss]] <- list(
                        scenario_id = ss,
                        sigma_type = "factor",
                        m = m,
                        T_len = T_len,
                        R2_target = R2_target,
                        beta_design = beta_design,
                        phi_f = phi_f,
                        phi_eta = phi_eta,
                        phi_u = phi_u,
                        rho = NA_real_,
                        factor_r = factor_r,
                        factor_loading_scenario = loading_scen,
                        factor_strength = factor_strength,
                        factor_idio_var = factor_idio_var
                      )
                      ss <- ss + 1L
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  scenarios
}

scenario_label <- function(sc) {
  if (sc$sigma_type == "toeplitz") {
    safe_label("endo", sc$sigma_type, "m", sc$m, "T", sc$T_len,
               "R2", sc$R2_target, "rho", sc$rho,
               "beta", sc$beta_design, "phif", sc$phi_f, "phie", sc$phi_eta)
  } else {
    safe_label("endo", sc$sigma_type, "m", sc$m, "T", sc$T_len,
               "R2", sc$R2_target, "r", sc$factor_r,
               "load", sc$factor_loading_scenario,
               "fs", sc$factor_strength, "id", sc$factor_idio_var,
               "beta", sc$beta_design, "phif", sc$phi_f, "phie", sc$phi_eta)
  }
}

plot_scenario <- function(sum_df, sc, outdir) {
  lab <- scenario_label(sc)

  png(file.path(outdir, paste0(lab, "_rmsfe_by_z.png")), width = 1000, height = 680)
  ylim <- range(c(sum_df$rmsfe_base, sum_df$rmsfe_rs), na.rm = TRUE)
  plot(sum_df$z, sum_df$rmsfe_rs, type = "b", pch = 19, lwd = 2,
       col = PLOT_COL["bagged"], ylim = ylim,
       xlab = "subset size z", ylab = "RMSFE",
       main = paste("Endogenous forecasts:", sc$sigma_type, "Sigma_e"))
  abline(h = mean(sum_df$rmsfe_base, na.rm = TRUE), lty = 2, lwd = 2,
         col = PLOT_COL["base"])
  legend("topright",
         legend = c("RS bagged", "equal-weight baseline"),
         col = c(PLOT_COL["bagged"], PLOT_COL["base"]),
         lty = c(1, 2), pch = c(19, NA), lwd = 2, bty = "n")
  dev.off()

  png(file.path(outdir, paste0(lab, "_Q_by_z.png")), width = 1000, height = 680)
  ylim <- range(c(sum_df$Q_ols, sum_df$Q_rs_bag, sum_df$Q_rs_subavg), na.rm = TRUE)
  plot(sum_df$z, sum_df$Q_rs_bag, type = "b", pch = 19, lwd = 2,
       col = PLOT_COL["q_bag"], ylim = ylim,
       xlab = "subset size z", ylab = "estimated population component Q",
       main = paste("Q diagnostics:", sc$sigma_type, "Sigma_e"))
  lines(sum_df$z, sum_df$Q_rs_subavg, type = "b", pch = 17, lty = 2, lwd = 2,
        col = PLOT_COL["q_sub"])
  lines(sum_df$z, sum_df$Q_ols, lty = 3, lwd = 2,
        col = PLOT_COL["q_ols"])
  legend("topright",
         legend = c("Q RS bagged", "Q RS subset avg", "Q full OLS weights"),
         col = c(PLOT_COL["q_bag"], PLOT_COL["q_sub"], PLOT_COL["q_ols"]),
         lty = c(1, 2, 3), pch = c(19, 17, NA), lwd = 2, bty = "n")
  dev.off()

  png(file.path(outdir, paste0(lab, "_L_by_z.png")), width = 1000, height = 680)
  ylim <- range(c(sum_df$L_ols, sum_df$L_rs_bag, sum_df$L_rs_subavg), na.rm = TRUE)
  plot(sum_df$z, sum_df$L_rs_bag, type = "b", pch = 19, lwd = 2,
       col = PLOT_COL["q_bag"], ylim = ylim,
       xlab = "subset size z", ylab = "inflation-adjusted L",
       main = paste("L diagnostics:", sc$sigma_type, "Sigma_e"))
  lines(sum_df$z, sum_df$L_rs_subavg, type = "b", pch = 17, lty = 2, lwd = 2,
        col = PLOT_COL["q_sub"])
  lines(sum_df$z, sum_df$L_ols, lty = 3, lwd = 2,
        col = PLOT_COL["q_ols"])
  legend("topright",
         legend = c("L RS bagged", "L RS subset avg", "L full OLS weights"),
         col = c(PLOT_COL["q_bag"], PLOT_COL["q_sub"], PLOT_COL["q_ols"]),
         lty = c(1, 2, 3), pch = c(19, 17, NA), lwd = 2, bty = "n")
  dev.off()
}

###############################################################################
# 7. Main scenario runner
###############################################################################

run_scenario <- function(sc) {
  cat("\nScenario", sc$scenario_id, ":", sc$sigma_type,
      "m=", sc$m, "T=", sc$T_len,
      "R2=", sc$R2_target, "beta=", sc$beta_design, "\n")

  Sigma_e <- make_sigma_e(
    type = sc$sigma_type,
    m = sc$m,
    rho = ifelse(is.na(sc$rho), 0, sc$rho),
    factor_r = ifelse(is.na(sc$factor_r), 2, sc$factor_r),
    factor_loading_scenario = ifelse(is.na(sc$factor_loading_scenario), "random", sc$factor_loading_scenario),
    factor_strength = ifelse(is.na(sc$factor_strength), 1.0, sc$factor_strength),
    factor_idio_var = ifelse(is.na(sc$factor_idio_var), 0.5, sc$factor_idio_var),
    seed = GLOBAL_SEED + 1000L * sc$scenario_id
  )

  beta_raw <- beta_raw_design(sc$m, sc$beta_design)
  beta_true <- scale_beta_to_R2(
    beta_raw = beta_raw,
    Sigma_e = Sigma_e,
    R2_target = sc$R2_target,
    sigma_u = SIGMA_U,
    phi_f = sc$phi_f,
    phi_eta = sc$phi_eta,
    scale_using = BETA_SCALE_USING
  )

  z_grid <- seq(2, sc$m, by = K_BY)
  if (tail(z_grid, 1) != sc$m) z_grid <- c(z_grid, sc$m)

  # -------------------------------------------------------------------
  # Population theory benchmark (Fix 3).
  # Compute the true population Q_bag_inf(k) and Q_sub(k) curves using
  # q_population_mc_general on the true Sigma_e.  These are the
  # infeasible population objects that the MC RMSFE curve should track.
  # They also provide the exact dominance threshold kappa^exact from the
  # appendix: the smallest k where R_k^{bag,inf}*(T-m) < A*(m-k).
  # -------------------------------------------------------------------
  N_POP_MC_ENDO <- if (RUN_PROFILE == "fast") 1000L else 2000L

  cat("  Computing population Q curves via MC on true Sigma_e...\n")
  pop_mc <- q_population_mc_general(
    Sigma_e  = Sigma_e,
    k_grid   = z_grid,
    n_subsets = N_POP_MC_ENDO
  )

  # Inflation-adjusted population loss curves (using paper formula).
  # T_ref is the median training window length over the OOS period.
  T_ref <- as.integer(sc$T_len - OOS_LEN + floor((OOS_LEN - 1) / 2))
  A_pop  <- SIGMA_U^2 + pop_mc$Q_opt[1]   # s + Q_opt

  pop_mc$L_pop_bag_inf <- ifelse(
    z_grid < T_ref,
    (SIGMA_U^2 + pop_mc$Q_bag_inf_mc) * (T_ref - 1) / (T_ref - z_grid),
    NA_real_
  )
  pop_mc$L_pop_opt <- if (T_ref > sc$m) {
    (SIGMA_U^2 + pop_mc$Q_opt[1]) * (T_ref - 1) / (T_ref - sc$m)
  } else {
    NA_real_
  }

  # Exact dominance threshold from the appendix eq (172):
  #   R_k^{bag,inf} * (T - m) < A * (m - k).
  R_bag_inf_mc <- pop_mc$R_bag_inf_mc
  exact_dom_cond <- !is.na(R_bag_inf_mc) &
    R_bag_inf_mc * max(0, T_ref - sc$m) < A_pop * pmax(0, sc$m - z_grid)
  kappa_exact <- if (any(exact_dom_cond)) {
    z_grid[which(exact_dom_cond)[1]]
  } else {
    NA_integer_
  }

  lab_pop <- scenario_label(sc)
  write.csv(
    pop_mc,
    file.path(OUTDIR, paste0(lab_pop, "_pop_mc_general.csv")),
    row.names = FALSE
  )
  cat("  Population Q_opt =", round(pop_mc$Q_opt[1], 4),
      " | kappa_exact (dominance threshold) =", kappa_exact, "\n")

  rep_dfs <- vector("list", N_REP)
  for (rr in seq_len(N_REP)) {
    if (rr %% max(1, floor(N_REP / 5)) == 0) {
      cat("  replication", rr, "of", N_REP, "\n")
    }
    set.seed(GLOBAL_SEED + 100000L * sc$scenario_id + rr)
    rep_dfs[[rr]] <- run_one_rep_all_z(
      T_len = sc$T_len,
      m = sc$m,
      beta_true = beta_true,
      Sigma_e = Sigma_e,
      phi_f = sc$phi_f,
      phi_eta = sc$phi_eta,
      phi_u = sc$phi_u,
      z_grid = z_grid,
      K_RS = K_RS,
      OOS_LEN = OOS_LEN,
      sigma_u = SIGMA_U,
      ridge = RIDGE
    )
  }

  sum_df <- aggregate_replications(rep_dfs)

  best_idx <- which.min(sum_df$rmsfe_rs)
  best_z_rmsfe <- sum_df$z[best_idx]
  best_idx_L_bag <- which.min(sum_df$L_rs_bag)
  best_idx_L_sub <- which.min(sum_df$L_rs_subavg)

  scenario_cols <- data.frame(
    scenario_id = sc$scenario_id,
    sigma_type = sc$sigma_type,
    m = sc$m,
    T_len = sc$T_len,
    R2_target = sc$R2_target,
    beta_design = sc$beta_design,
    phi_f = sc$phi_f,
    phi_eta = sc$phi_eta,
    phi_u = sc$phi_u,
    rho = sc$rho,
    factor_r = sc$factor_r,
    factor_loading_scenario = sc$factor_loading_scenario,
    factor_strength = sc$factor_strength,
    factor_idio_var = sc$factor_idio_var,
    beta_scale_using = BETA_SCALE_USING
  )

  sum_df <- cbind(scenario_cols[rep(1, nrow(sum_df)), ], sum_df)

  lab <- scenario_label(sc)
  write.csv(sum_df, file.path(OUTDIR, paste0(lab, "_by_z.csv")), row.names = FALSE)
  if (SAVE_PLOTS) {
    plot_scenario(sum_df, sc, OUTDIR)
    # Additional plot: population Q_bag_inf vs empirical Q_rs_bag (Fix 3).
    png(file.path(OUTDIR, paste0(lab, "_pop_vs_empirical_Q.png")),
        width = 1000, height = 680)
    ylim_q <- range(c(pop_mc$Q_bag_inf_mc, pop_mc$Q_opt,
                      sum_df$Q_rs_bag, sum_df$Q_ols), na.rm = TRUE)
    plot(pop_mc$k, pop_mc$Q_bag_inf_mc, type = "l", lwd = 2.5,
         col = PLOT_COL["q_bag"], ylim = ylim_q,
         xlab = "subset size z", ylab = "population component Q",
         main = paste("Population vs empirical Q:", sc$sigma_type, "Sigma_e"))
    lines(sum_df$z, sum_df$Q_rs_bag, lty = 2, lwd = 2,
          col = PLOT_COL["bagged"])
    abline(h = pop_mc$Q_opt[1], lty = 3, lwd = 2,
           col = PLOT_COL["q_ols"])
    if (!is.na(kappa_exact))
      abline(v = kappa_exact, lty = 4, lwd = 1.5, col = "grey40")
    legend("topright",
           legend = c("Pop Q_bag,inf (true Sigma_e)", "Empirical Q_rs_bag",
                      "Pop Q_opt", "Dominance threshold kappa"),
           col = c(PLOT_COL["q_bag"], PLOT_COL["bagged"],
                   PLOT_COL["q_ols"], "grey40"),
           lty = c(1, 2, 3, 4), lwd = 2, bty = "n")
    dev.off()
  }

  summary_row <- cbind(
    scenario_cols,
    data.frame(
      avg_abs_beta = mean(abs(beta_true)),
      max_abs_beta = max(abs(beta_true)),
      avg_corr_Sigma_e = mean(cov2cor(Sigma_e)[upper.tri(Sigma_e)]),
      min_eigen_Sigma_e = min(eigen(Sigma_e, symmetric = TRUE, only.values = TRUE)$values),
      # Population theory benchmarks (Fix 3).
      Q_pop_opt = pop_mc$Q_opt[1],
      kappa_exact_dom = kappa_exact,
      best_z_pop_bag_inf = if (all(is.na(pop_mc$L_pop_bag_inf))) NA_integer_ else
        z_grid[which.min(pop_mc$L_pop_bag_inf)],
      max_R_bag_inf_mc = max(R_bag_inf_mc, na.rm = TRUE),
      # Empirical results.
      best_z_rmsfe = best_z_rmsfe,
      best_z_L_bag = sum_df$z[best_idx_L_bag],
      best_z_L_subavg = sum_df$z[best_idx_L_sub],
      min_rmsfe_rs = min(sum_df$rmsfe_rs, na.rm = TRUE),
      rmsfe_base = mean(sum_df$rmsfe_base, na.rm = TRUE),
      rel_rmsfe_at_best = min(sum_df$rel_rmsfe, na.rm = TRUE),
      min_L_rs_bag = min(sum_df$L_rs_bag, na.rm = TRUE),
      min_L_rs_subavg = min(sum_df$L_rs_subavg, na.rm = TRUE),
      mean_L_ols = mean(sum_df$L_ols, na.rm = TRUE)
    )
  )

  print(summary_row[, c("scenario_id", "sigma_type", "m", "T_len", "R2_target",
                       "beta_design", "best_z_rmsfe", "rel_rmsfe_at_best")])

  list(summary = summary_row, by_z = sum_df)
}

###############################################################################
# 8. T-scaling experiment functions
###############################################################################

# Fast inner bag loop used only by the T-scaling experiment.
# Avoids calling eigen() on every submatrix by adding the ridge directly.
rsm_fast_one_origin <- function(Sigma_est, f_now, z_grid, K_RS, T_tr,
                                ridge = RIDGE) {
  m   <- ncol(Sigma_est)
  one <- rep(1, m)
  n_z <- length(z_grid)

  pred_rs  <- numeric(n_z)
  Q_rs_bag <- numeric(n_z)
  L_rs_bag <- numeric(n_z)

  # Full estimated OLS for reference (z == m case)
  Sigma_full <- Sigma_est + ridge * diag(m)
  Sinv1_full <- tryCatch(solve(Sigma_full, one),
                         error = function(e) rep(1 / m, m))
  denom_full <- sum(Sinv1_full)
  w_full <- if (is.finite(denom_full) && abs(denom_full) > 1e-12)
              Sinv1_full / denom_full else rep(1 / m, m)
  Q_full <- as.numeric(t(w_full) %*% Sigma_est %*% w_full)
  L_full <- if (T_tr > m) Q_full * (T_tr - 1) / (T_tr - m) else NA_real_

  for (ii in seq_along(z_grid)) {
    z <- z_grid[ii]
    if (z == m) {
      pred_rs[ii]  <- sum(w_full * f_now)
      Q_rs_bag[ii] <- Q_full
      L_rs_bag[ii] <- L_full
      next
    }
    infl   <- (T_tr - 1) / (T_tr - z)
    w_acc  <- numeric(m)
    Q_acc  <- 0
    for (bb in seq_len(K_RS)) {
      idx    <- sample.int(m, z, replace = FALSE)
      Sig_S  <- Sigma_est[idx, idx, drop = FALSE] + ridge * diag(z)
      s1     <- tryCatch(solve(Sig_S, rep(1, z)),
                         error = function(e) rep(1 / z, z))
      ds     <- sum(s1)
      wS     <- if (is.finite(ds) && abs(ds) > 1e-12) s1 / ds else rep(1/z, z)
      w_acc[idx] <- w_acc[idx] + wS
      Q_acc  <- Q_acc + as.numeric(t(wS) %*% Sigma_est[idx, idx] %*% wS)
    }
    w_bar        <- w_acc / K_RS
    sw           <- sum(w_bar)
    if (is.finite(sw) && sw > 1e-12) w_bar <- w_bar / sw
    Q_bag        <- as.numeric(t(w_bar) %*% Sigma_est %*% w_bar)
    pred_rs[ii]  <- sum(w_bar * f_now)
    Q_rs_bag[ii] <- Q_bag
    L_rs_bag[ii] <- Q_bag * infl
  }

  list(pred_rs  = pred_rs,
       Q_rs_bag = Q_rs_bag,
       L_rs_bag = L_rs_bag,
       Q_full   = Q_full,
       L_full   = L_full,
       w_full   = w_full)
}

run_one_rep_scaling <- function(T_len, m, beta_true, Sigma_e,
                                phi_f, phi_eta, phi_u,
                                z_grid, K_RS, OOS_LEN,
                                sigma_u = SIGMA_U, ridge = RIDGE) {
  dat <- simulate_endogenous_path(T_len = T_len, beta_true = beta_true,
                                  Sigma_e = Sigma_e, phi_f = phi_f,
                                  phi_eta = phi_eta, phi_u = phi_u,
                                  sigma_u = sigma_u)
  X <- dat$X; Y <- dat$Y

  init_window <- T_len - OOS_LEN
  steps       <- OOS_LEN - 1
  t_now0_vec  <- init_window + seq_len(steps)

  base_pred <- numeric(steps)
  rs_pred   <- matrix(NA_real_, steps, length(z_grid))
  Q_bag_mat <- matrix(NA_real_, steps, length(z_grid))
  L_bag_mat <- matrix(NA_real_, steps, length(z_grid))

  for (hh in seq_along(t_now0_vec)) {
    T_tr     <- t_now0_vec[hh] - 1
    if (T_tr < max(z_grid) + 5L) next
    X_train  <- X[1:T_tr, , drop = FALSE]
    y_train  <- Y[2:(T_tr + 1)]
    x_now    <- X[t_now0_vec[hh], ]
    gamma    <- estimate_univariate_gammas(X_train, y_train)
    f_now    <- as.numeric(gamma * x_now)
    base_pred[hh] <- mean(f_now)
    Sigma_hat <- forecast_error_covariance(X_train, y_train, gamma, ridge)
    res <- rsm_fast_one_origin(Sigma_hat, f_now, z_grid, K_RS, T_tr, ridge)
    rs_pred[hh, ]   <- res$pred_rs
    Q_bag_mat[hh, ] <- res$Q_rs_bag
    L_bag_mat[hh, ] <- res$L_rs_bag
  }

  y_oos     <- Y[t_now0_vec + 1]
  err_base  <- y_oos - base_pred
  rows <- vector("list", length(z_grid))
  for (ii in seq_along(z_grid)) {
    err_rs <- y_oos - rs_pred[, ii]
    rows[[ii]] <- data.frame(
      z              = z_grid[ii],
      rmsfe_base_rep = sqrt(mean(err_base^2, na.rm = TRUE)),
      rmsfe_rs_rep   = sqrt(mean(err_rs^2, na.rm = TRUE)),
      sum_e2_base    = sum(err_base^2, na.rm = TRUE),
      cnt_base       = sum(is.finite(err_base)),
      sum_e2_rs      = sum(err_rs^2, na.rm = TRUE),
      cnt_rs         = sum(is.finite(err_rs)),
      Q_rs_bag       = mean(Q_bag_mat[, ii], na.rm = TRUE),
      L_rs_bag       = mean(L_bag_mat[, ii], na.rm = TRUE)
    )
  }
  do.call(rbind, rows)
}

run_scaling_T_endo <- function(cfg, outdir) {
  cat("\n=== T-scaling experiment (fixed Sigma_e, varying T) ===\n")
  all_rows <- data.frame()

  for (ss in seq_along(cfg$sigma_specs)) {
    spec <- cfg$sigma_specs[[ss]]
    cat("\n  Sigma type:", spec$type, "\n")

    # Build Sigma_e once with fixed seed so it is identical across all T.
    Sigma_e <- make_sigma_e(
      type                    = spec$type,
      m                       = cfg$m,
      rho                     = ifelse(is.na(spec$rho), 0, spec$rho),
      factor_r                = ifelse(is.na(spec$factor_r), 2L, spec$factor_r),
      factor_loading_scenario = ifelse(is.na(spec$loading), "random", spec$loading),
      factor_strength         = ifelse(is.na(spec$strength), 1.0, spec$strength),
      factor_idio_var         = ifelse(is.na(spec$idio), 0.5, spec$idio),
      seed                    = spec$seed_sigma
    )

    # Build beta once (same Sigma_e -> same beta for all T).
    beta_raw  <- beta_raw_design(cfg$m, cfg$beta_design)
    beta_true <- scale_beta_to_R2(
      beta_raw   = beta_raw,
      Sigma_e    = Sigma_e,
      R2_target  = cfg$R2_target,
      sigma_u    = SIGMA_U,
      phi_f      = cfg$phi_f,
      phi_eta    = cfg$phi_eta,
      scale_using = BETA_SCALE_USING
    )

    for (T_len in cfg$T_len_grid) {
      cat("    T =", T_len, ": ")
      if (T_len - cfg$OOS_LEN_scale <= cfg$m + 5L) {
        cat("skipped (training window too short)\n")
        next
      }
      z_grid <- seq(2L, cfg$m, by = cfg$K_BY_scale)
      if (tail(z_grid, 1L) != cfg$m) z_grid <- c(z_grid, cfg$m)

      rep_dfs <- vector("list", cfg$N_REP_scale)
      for (rr in seq_len(cfg$N_REP_scale)) {
        set.seed(GLOBAL_SEED + 900000L + spec$seed_sigma + T_len * 100L + rr)
        rep_dfs[[rr]] <- run_one_rep_scaling(
          T_len   = T_len, m = cfg$m, beta_true = beta_true,
          Sigma_e = Sigma_e, phi_f = cfg$phi_f,
          phi_eta = cfg$phi_eta, phi_u = cfg$phi_u,
          z_grid  = z_grid, K_RS = cfg$K_RS_scale,
          OOS_LEN = cfg$OOS_LEN_scale
        )
      }
      cat(cfg$N_REP_scale, "reps done\n")

      # Aggregate over reps.
      all_r <- do.call(rbind, rep_dfs)
      agg <- do.call(rbind, lapply(z_grid, function(z0) {
        dd <- all_r[all_r$z == z0, ]
        mse_b  <- sum(dd$sum_e2_base) / sum(dd$cnt_base)
        mse_rs <- sum(dd$sum_e2_rs)   / sum(dd$cnt_rs)
        data.frame(
          sigma_type  = spec$type,
          T_len       = T_len,
          z           = z0,
          rmsfe_base  = sqrt(mse_b),
          rmsfe_rs    = sqrt(mse_rs),
          rel_rmsfe   = sqrt(mse_rs / mse_b),
          Q_rs_bag    = mean(dd$Q_rs_bag, na.rm = TRUE),
          L_rs_bag    = mean(dd$L_rs_bag, na.rm = TRUE)
        )
      }))

      # Save per-(T, sigma_type) RMSFE curve.
      lab_T <- safe_label("scaling", spec$type, "T", T_len, "m", cfg$m,
                          "R2", cfg$R2_target, "beta", cfg$beta_design)
      write.csv(agg, file.path(outdir, paste0(lab_T, "_by_z.csv")),
                row.names = FALSE)

      # Summary row: best z by RMSFE and by L_rs_bag.
      best_z_rmsfe <- agg$z[which.min(agg$rmsfe_rs)]
      best_z_L     <- agg$z[which.min(agg$L_rs_bag)]
      all_rows <- rbind(all_rows, data.frame(
        sigma_type      = spec$type,
        T_len           = T_len,
        m               = cfg$m,
        R2_target       = cfg$R2_target,
        beta_design     = cfg$beta_design,
        N_REP_scale     = cfg$N_REP_scale,
        K_RS_scale      = cfg$K_RS_scale,
        best_z_rmsfe    = best_z_rmsfe,
        best_z_L_bag    = best_z_L,
        min_rmsfe_rs    = min(agg$rmsfe_rs, na.rm = TRUE),
        rmsfe_base      = agg$rmsfe_base[1],
        rel_rmsfe_best  = min(agg$rel_rmsfe, na.rm = TRUE),
        min_L_rs_bag    = min(agg$L_rs_bag, na.rm = TRUE)
      ))
    }
  }

  write.csv(all_rows, file.path(outdir, "scaling_T_summary.csv"),
            row.names = FALSE)

  # ---- Plots ----

  sigma_types <- unique(all_rows$sigma_type)

  # 1. RMSFE curves overlaid by T for each sigma type.
  for (stype in sigma_types) {
    sub <- all_rows[all_rows$sigma_type == stype, ]
    T_vals <- sort(unique(sub$T_len))
    col_pal <- hcl.colors(length(T_vals), "Zissou 1")

    png(file.path(outdir, paste0("scaling_T_rmsfe_curves_", stype, ".png")),
        width = 1100, height = 720)
    first <- TRUE
    ylim_all <- NULL
    for (ii in seq_along(T_vals)) {
      lab_T <- safe_label("scaling", stype, "T", T_vals[ii], "m", cfg$m,
                          "R2", cfg$R2_target, "beta", cfg$beta_design)
      csv_f <- file.path(outdir, paste0(lab_T, "_by_z.csv"))
      if (!file.exists(csv_f)) next
      dd <- read.csv(csv_f)
      ylim_all <- range(c(ylim_all, dd$rmsfe_rs, dd$rmsfe_base), na.rm = TRUE)
    }
    for (ii in seq_along(T_vals)) {
      lab_T <- safe_label("scaling", stype, "T", T_vals[ii], "m", cfg$m,
                          "R2", cfg$R2_target, "beta", cfg$beta_design)
      csv_f <- file.path(outdir, paste0(lab_T, "_by_z.csv"))
      if (!file.exists(csv_f)) next
      dd <- read.csv(csv_f)
      if (first) {
        plot(dd$z, dd$rmsfe_rs, type = "l", lwd = 2, col = col_pal[ii],
             ylim = ylim_all, xlab = "subset size k", ylab = "RMSFE",
             main = paste0("RMSFE by k for varying T  [", stype, ", m=",
                           cfg$m, ", beta=", cfg$beta_design, "]"))
        first <- FALSE
      } else {
        lines(dd$z, dd$rmsfe_rs, lwd = 2, col = col_pal[ii])
      }
      abline(v = dd$z[which.min(dd$rmsfe_rs)], lty = 3,
             col = col_pal[ii], lwd = 1.2)
    }
    legend("topright", legend = paste("T =", T_vals),
           col = col_pal, lty = 1, lwd = 2, bty = "n")
    dev.off()
  }

  # 2. Optimal k* vs T  (linear scale + log-log).
  for (stype in sigma_types) {
    sub <- all_rows[all_rows$sigma_type == stype, ]
    sub <- sub[order(sub$T_len), ]
    if (nrow(sub) < 2) next

    png(file.path(outdir, paste0("scaling_T_kstar_", stype, ".png")),
        width = 1100, height = 600)
    par(mfrow = c(1, 2))

    # Panel A: k* vs T^(1/3)
    plot(sub$T_len^(1/3), sub$best_z_rmsfe, type = "b", pch = 19, lwd = 2,
         col = PLOT_COL["bagged"],
         xlab = expression(T^{1/3}),
         ylab = "optimal k*  (by RMSFE)",
         main = paste("k* vs T^(1/3) [", stype, "]"))
    lines(sub$T_len^(1/3), sub$best_z_L_bag, type = "b", pch = 17, lty = 2,
          lwd = 2, col = PLOT_COL["q_bag"])
    # Fit OLS line through origin (cube-root scaling)
    if (sum(!is.na(sub$best_z_rmsfe)) >= 2) {
      lm_r <- lm(best_z_rmsfe ~ 0 + I(T_len^(1/3)), data = sub)
      lm_L <- lm(best_z_L_bag  ~ 0 + I(T_len^(1/3)), data = sub)
      x_seq <- seq(min(sub$T_len)^(1/3), max(sub$T_len)^(1/3), length.out = 50)
      lines(x_seq, coef(lm_r) * x_seq, lty = 4, col = PLOT_COL["bagged"])
      lines(x_seq, coef(lm_L) * x_seq, lty = 4, col = PLOT_COL["q_bag"])
    }
    legend("topleft",
           legend = c("RMSFE-optimal k*", "L_bag-optimal k*",
                      "cube-root fit (origin)"),
           col = c(PLOT_COL["bagged"], PLOT_COL["q_bag"], "black"),
           lty = c(1, 2, 4), pch = c(19, 17, NA), lwd = 2, bty = "n")

    # Panel B: log-log with OLS slope
    valid <- sub$T_len > 0 & sub$best_z_rmsfe > 0
    if (sum(valid) >= 2) {
      lm_loglog <- lm(log(best_z_rmsfe) ~ log(T_len), data = sub[valid, ])
      beta_T    <- round(coef(lm_loglog)[2], 3)
      plot(log(sub$T_len[valid]), log(sub$best_z_rmsfe[valid]),
           type = "b", pch = 19, lwd = 2, col = PLOT_COL["bagged"],
           xlab = "log T", ylab = "log k*",
           main = paste0("log-log  [", stype, "]  slope=", beta_T,
                         "  (theory 1/3)"))
      abline(lm_loglog, lty = 2, col = PLOT_COL["bagged"])
    } else {
      plot.new()
      title(paste("log-log [", stype, "] — insufficient data"))
    }
    par(mfrow = c(1, 1))
    dev.off()
  }

  # 3. Table: k* and relative RMSFE by T.
  cat("\n--- T-scaling summary table ---\n")
  print(all_rows[, c("sigma_type", "T_len", "best_z_rmsfe",
                     "best_z_L_bag", "rel_rmsfe_best", "min_rmsfe_rs",
                     "rmsfe_base")])
  cat("\n")

  # 4. Log-log regression summary saved to file.
  sink(file.path(outdir, "scaling_T_loglog_regression.txt"))
  for (stype in sigma_types) {
    sub <- all_rows[all_rows$sigma_type == stype, ]
    sub <- sub[order(sub$T_len) & sub$best_z_rmsfe > 0, ]
    cat("Sigma type:", stype, "\n")
    if (nrow(sub) >= 3) {
      lm_r <- lm(log(best_z_rmsfe) ~ log(T_len), data = sub)
      lm_L <- lm(log(best_z_L_bag) ~ log(T_len),
                  data = sub[sub$best_z_L_bag > 0, ])
      cat("  RMSFE-optimal k*: log(k*) ~ log(T)\n")
      print(summary(lm_r))
      cat("  L_bag-optimal k*: log(k*) ~ log(T)\n")
      print(summary(lm_L))
      cat("  Theoretical slope: 1/3 ≈ 0.333\n\n")
    } else {
      cat("  Not enough T values for regression.\n\n")
    }
  }
  sink()

  invisible(all_rows)
}

###############################################################################
# 9. Main run block
###############################################################################

cat("Endogenous RSM Monte Carlo started.\n")
cat("Output directory:", OUTDIR, "\n")
cat("Profile:", RUN_PROFILE, "\n")
cat("N_REP:", N_REP, " K_RS:", K_RS, " OOS_LEN:", OOS_LEN, "\n")

# ---- Part A: main scenarios ----
scenarios <- build_scenarios()
cat("Number of main scenarios:", length(scenarios), "\n")

all_summaries <- vector("list", length(scenarios))
for (ii in seq_along(scenarios)) {
  res <- run_scenario(scenarios[[ii]])
  all_summaries[[ii]] <- res$summary
}

summary_df <- do.call(rbind, all_summaries)
write.csv(summary_df,
          file.path(OUTDIR, "endogenous_rsm_summary.csv"),
          row.names = FALSE)

# ---- Part B: T-scaling experiment ----
scaling_results <- run_scaling_T_endo(T_SCALING_CONFIG, OUTDIR)

cat("\nAll simulations completed.\n")
cat("Results saved in:", OUTDIR, "\n")
