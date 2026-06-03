#!/usr/bin/env Rscript

###############################################################################
# Clean common-loading RSM Monte Carlo for the small-heterogeneity proposal
#
# Purpose
# -------
# This script keeps only the common-loading heteroskedastic design. It separates:
#   1. Population objects: Q_sub(k), Q_bag_inf(k), Q_bag_G(k).
#   2. Finite-sample Monte Carlo MSPE curves.
#   3. Benchmark comparisons: equal weights, oracle optimal weights, estimated
#      full optimal weights.
#   4. Scaling checks for optimal k as T and m vary.
#   5. Joint B_scale x H grid checks and weight-distance diagnostics.
#
# Notation
# --------
# B_scale = 1 / mean(tau_i), the proposal's idiosyncratic scale B.
# G_bag   = number of subset draws used for bagging. This avoids confusing the
#           proposal's B with the bagging count.
# H       = mean(sigma_i^2) * mean(tau_i). H close to 1 means small heterogeneity.
# V_delta = mean((tau_i / mean(tau_i) - 1)^2). Under small heterogeneity,
#           H - 1 is approximately V_delta.
###############################################################################

rm(list = ls())

###############################################################################
# 0. User controls
###############################################################################

RUN_PROFILE <- Sys.getenv("RSM_PROFILE", unset = "paper")
# Options: "fast", "paper".

OUTDIR <- Sys.getenv("RSM_OUTDIR", unset = file.path(getwd(), "rsm_common_loading_clean_output"))
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

GLOBAL_SEED <- as.integer(Sys.getenv("RSM_SEED", unset = "12350"))
set.seed(GLOBAL_SEED)

parse_num_vec <- function(env_name, default) {
  raw <- Sys.getenv(env_name, unset = "")
  if (!nzchar(raw)) return(default)
  vals <- suppressWarnings(as.numeric(trimws(unlist(strsplit(raw, ",")))))
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0) default else vals
}

parse_int_vec <- function(env_name, default) {
  vals <- parse_num_vec(env_name, default)
  vals <- as.integer(round(vals))
  vals[!is.na(vals)]
}

make_scenario_id <- function(prefix, m, T_train, B_scale, H_target) {
  paste0(
    prefix,
    "_m", m,
    "_T", T_train,
    "_B", gsub("\\.", "p", format(B_scale, trim = TRUE, scientific = FALSE)),
    "_H", gsub("\\.", "p", format(H_target, trim = TRUE, scientific = FALSE))
  )
}

RUN_EXPERIMENTS <- c(
  "baseline",
  "scaling_T_m",
  "joint_B_H_grid"
)

SAVE_PLOTS <- TRUE
SAVE_CSV <- TRUE

if (RUN_PROFILE == "fast") {
  # Quick pilot. Use this to check logic and plots.
  N_REP_BASE      <- 80
  N_REP_SCALING   <- 40
  N_REP_GRID      <- 40
  G_BAG_MC        <- 250
  N_SUBSETS_POP   <- 2500
  T_TEST          <- 600
  K_BY            <- 2
} else if (RUN_PROFILE == "paper") {
  # Slower, smoother figures.
  N_REP_BASE      <- 500
  N_REP_SCALING   <- 500
  N_REP_GRID      <- 500
  G_BAG_MC        <- 800
  N_SUBSETS_POP   <- 3000
  T_TEST          <- 100
  K_BY            <- 1
} else {
  stop("RUN_PROFILE must be 'fast' or 'paper'.")
}

parse_int_scalar <- function(env_name, default) {
  raw <- Sys.getenv(env_name, unset = "")
  if (!nzchar(raw)) return(default)
  val <- suppressWarnings(as.integer(raw))
  if (!is.finite(val) || val <= 0) default else val
}

# Population Monte Carlo size for the joint B x H grid. By default it equals
# N_SUBSETS_POP. Set RSM_N_SUBSETS_POP_JOINT explicitly if you want a cheaper
# but noisier joint-grid diagnostic run.
N_SUBSETS_POP_JOINT <- parse_int_scalar("RSM_N_SUBSETS_POP_JOINT", N_SUBSETS_POP)

# Correlations between the exact weight error and the small-heterogeneity
# approximation are meaningful only when both vectors have non-negligible
# dispersion. Near k = m the approximation is exactly zero, so we suppress those
# correlations rather than letting Monte Carlo noise create artificial values.
WEIGHT_CORR_TOL <- 1e-10

# Main baseline design. B_scale and H_target may be scalar or grids.
# You can override them without editing the file, for example:
#   RSM_BASELINE_B_GRID="1,3,6" RSM_BASELINE_H_GRID="1.05,1.5,2" Rscript this_file.R
BASELINE <- list(
  m = 50,
  T_train = 240,
  s = 1.0,
  q = 1.0,
  B_scale = parse_num_vec("RSM_BASELINE_B_GRID", c(3.0, 8)),
  H_target = parse_num_vec("RSM_BASELINE_H_GRID", c(1.05, 6)),
  k_min = 1,
  k_max = NA_integer_,
  k_by = K_BY
)

# Selected k values for the per-forecast weight comparison
# between w_opt and the averaged bagged weights \bar v_k.
WEIGHT_DIAG_K <- parse_int_vec("RSM_WEIGHT_DIAG_K", c(1, 3, 5, 9, 21, 41, 79))

# Scaling design. Keep B_scale and H fixed if you want the clean proposal check:
# k*_bag,inf should be close to T^(1/3) m^(-1/3), up to constants.
GRID_SCALING <- list(
  m = c(50),
  T_train = c(80, 240, 640, 1000),
  s = 1.0,
  q = 1.0,
  B_scale = parse_num_vec("RSM_SCALING_B", c(3.0))[1],
  H_target = parse_num_vec("RSM_SCALING_H", c(1.05))[1],
  k_min = 1,
  k_by = K_BY
)

# Joint B x H grid. This replaces the old one-at-a-time B and H sensitivities.
# It lets you study the interaction between scale B and heterogeneity H.
GRID_JOINT_B_H <- list(
  m = 50,
  T_train = 240,
  s = 1.0,
  q = 1.0,
  B_scale = parse_num_vec("RSM_BH_B_GRID", c(3.0, 8)),
  H_target = parse_num_vec("RSM_BH_H_GRID", c(1.05, 6)),
  k_min = 1,
  k_by = K_BY
)

PLOT_COL <- c(
  bagged       = "#0072B2",
  subset       = "#D55E00",
  equal        = "#000000",
  full         = "#E69F00",
  oracle       = "#882255",
  theory       = "#009E73",
  finite       = "#332288",
  grey         = "#666666"
)

###############################################################################
# 1. Small helpers
###############################################################################

safe_solve <- function(A, b = NULL, ridge = 1e-10) {
  A <- as.matrix(A)
  if (is.null(b)) b <- diag(nrow(A))
  tryCatch(
    solve(A, b),
    error = function(e) solve(A + ridge * diag(nrow(A)), b)
  )
}

safe_write_csv <- function(x, filename) {
  if (isTRUE(SAVE_CSV)) write.csv(x, filename, row.names = FALSE)
}

mspe <- function(pred, y) {
  mean((as.numeric(pred) - as.numeric(y))^2)
}

nearest_feasible_integer <- function(x, lower, upper) {
  if (!is.finite(x)) return(NA_integer_)
  cand <- unique(c(floor(x), ceiling(x), round(x)))
  cand <- cand[cand >= lower & cand <= upper]
  if (length(cand) == 0) return(NA_integer_)
  cand[which.min(abs(cand - x))]
}

make_k_grid <- function(m, T_train, k_min = 1, k_max = NA_integer_, k_by = 1) {
  max_allowed <- min(m, T_train - 1)
  if (is.na(k_max)) k_max <- max_allowed
  k_max <- min(k_max, max_allowed)
  if (k_min > k_max) return(integer(0))
  unique(as.integer(seq(k_min, k_max, by = k_by)))
}

expand_grid_list <- function(grid) {
  grid <- lapply(grid, function(x) if (length(x) == 1) x else as.vector(x))
  expand.grid(grid, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
}

find_kstar_grid <- function(df, column, k_col = "k") {
  vals <- df[[column]]
  if (all(!is.finite(vals))) return(NA_integer_)
  df[[k_col]][which.min(vals)]
}

inflation <- function(T_train, k) {
  ifelse(T_train > k, (T_train - 1) / (T_train - k), NA_real_)
}

###############################################################################
# 2. Precision design and common-loading data
###############################################################################

make_tau_lognormal <- function(m, B_scale, H_target, seed = NULL) {
  # tau_i = 1 / sigma_i^2.
  # For lognormal tau in population, CV(tau)^2 = H - 1.
  # The sample is rescaled so mean(tau) = 1 / B_scale exactly.
  if (!is.null(seed)) set.seed(seed)
  if (H_target < 1) stop("H_target must be >= 1.")

  if (abs(H_target - 1) < 1e-14) {
    tau_raw <- rep(1, m)
  } else {
    omega <- sqrt(log(H_target))
    mu <- -0.5 * omega^2
    tau_raw <- exp(rnorm(m, mean = mu, sd = omega))
  }

  tau <- tau_raw / mean(tau_raw) * (1 / B_scale)
  sigma2 <- 1 / tau

  list(
    tau = tau,
    sigma2 = sigma2,
    B_scale = 1 / mean(tau),
    H = mean(sigma2) * mean(tau),
    V_delta = mean((tau / mean(tau) - 1)^2)
  )
}

precision_summary <- function(tau, q) {
  m <- length(tau)
  tau_bar <- mean(tau)
  B_scale <- 1 / tau_bar
  sigma2 <- 1 / tau
  delta <- tau / tau_bar - 1
  V_delta <- mean(delta^2)
  H <- mean(sigma2) * mean(tau)

  Q_opt <- q + B_scale / m
  Q_ave <- q + mean(sigma2) / m
  C_delta <- B_scale * V_delta / m * (m / (m - 1))^2

  list(
    m = m,
    B_scale = B_scale,
    H = H,
    V_delta = V_delta,
    Q_opt = Q_opt,
    Q_ave = Q_ave,
    C_delta = C_delta
  )
}

simulate_common_loading <- function(T_train, T_test, m, s, q, sigma2,
                                    phi_mu = 0.6, sd_mu = 1.0) {
  # Forecasts: f_it = mu_t + u_it.
  # Target:    y_t  = mu_t + eps_t.
  # Common-loading forecast errors: u_it = common_t + xi_it.
  # Because combination weights sum to one, the predictable component mu_t
  # cancels in pairwise forecast differences used by constrained OLS.
  n <- T_train + T_test

  mu <- numeric(n)
  mu[1] <- rnorm(1, sd = sd_mu / sqrt(1 - phi_mu^2))
  for (tt in 2:n) {
    mu[tt] <- phi_mu * mu[tt - 1] + rnorm(1, sd = sd_mu)
  }

  common <- rnorm(n, mean = 0, sd = sqrt(q))
  Xi <- matrix(rnorm(n * m), nrow = n, ncol = m)
  Xi <- sweep(Xi, 2, sqrt(sigma2), "*")
  eps <- rnorm(n, mean = 0, sd = sqrt(s))

  Fmat <- mu + common + Xi
  y <- mu + eps

  list(
    F_train = Fmat[seq_len(T_train), , drop = FALSE],
    y_train = y[seq_len(T_train)],
    F_test = Fmat[(T_train + 1):n, , drop = FALSE],
    y_test = y[(T_train + 1):n]
  )
}

###############################################################################
# 3. Weight estimation and finite-sample evaluation
###############################################################################

constrained_ols_weights <- function(y, Fmat, ridge = 1e-8) {
  # Sum-to-one constrained least squares without intercept.
  Fmat <- as.matrix(Fmat)
  k <- ncol(Fmat)
  if (k == 1) return(1)

  base <- Fmat[, k]
  Z <- Fmat[, 1:(k - 1), drop = FALSE] - base
  r <- y - base

  beta <- as.numeric(safe_solve(crossprod(Z) + ridge * diag(k - 1), crossprod(Z, r)))
  w <- c(beta, 1 - sum(beta))

  if (any(!is.finite(w))) w <- rep(1 / k, k)
  w
}

oracle_weights_common <- function(tau) {
  tau / sum(tau)
}

evaluate_rsm_for_k <- function(F_train, y_train, F_test, y_test, k,
                               G_bag = 500, ridge = 1e-8) {
  m <- ncol(F_train)
  T_test <- nrow(F_test)

  bagged_pred <- rep(0, T_test)
  subset_mspes <- numeric(G_bag)

  for (gg in seq_len(G_bag)) {
    idx <- sample.int(m, k, replace = FALSE)
    w <- constrained_ols_weights(y_train, F_train[, idx, drop = FALSE], ridge = ridge)
    pred <- as.numeric(F_test[, idx, drop = FALSE] %*% w)
    bagged_pred <- bagged_pred + pred
    subset_mspes[gg] <- mspe(pred, y_test)
  }

  bagged_pred <- bagged_pred / G_bag
  c(
    rsm_bagged = mspe(bagged_pred, y_test),
    rsm_subset_avg = mean(subset_mspes)
  )
}

evaluate_one_rep_common <- function(T_train, T_test, m, s, q, sigma2, tau,
                                    k_grid, G_bag) {
  dat <- simulate_common_loading(
    T_train = T_train,
    T_test = T_test,
    m = m,
    s = s,
    q = q,
    sigma2 = sigma2
  )

  F_train <- dat$F_train
  y_train <- dat$y_train
  F_test <- dat$F_test
  y_test <- dat$y_test

  w_equal <- rep(1 / m, m)
  mspe_equal <- mspe(F_test %*% w_equal, y_test)

  w_oracle <- oracle_weights_common(tau)
  mspe_oracle <- mspe(F_test %*% w_oracle, y_test)

  mspe_full_est <- NA_real_
  if (T_train > m) {
    w_full <- constrained_ols_weights(y_train, F_train)
    mspe_full_est <- mspe(F_test %*% w_full, y_test)
  }

  rsm <- matrix(NA_real_, nrow = length(k_grid), ncol = 2)
  colnames(rsm) <- c("rsm_bagged", "rsm_subset_avg")
  for (ii in seq_along(k_grid)) {
    rsm[ii, ] <- evaluate_rsm_for_k(
      F_train = F_train,
      y_train = y_train,
      F_test = F_test,
      y_test = y_test,
      k = k_grid[ii],
      G_bag = G_bag
    )
  }

  data.frame(
    k = k_grid,
    rsm_bagged = rsm[, "rsm_bagged"],
    rsm_subset_avg = rsm[, "rsm_subset_avg"],
    equal_weight = mspe_equal,
    oracle_optimal = mspe_oracle,
    full_estimated = mspe_full_est
  )
}

summarise_reps_by_k <- function(rep_list) {
  all <- do.call(rbind, Map(function(x, rr) {
    x$rep <- rr
    x
  }, rep_list, seq_along(rep_list)))

  numeric_cols <- setdiff(names(all), c("k", "rep"))
  aggregate(all[, numeric_cols, drop = FALSE], by = list(k = all$k), FUN = mean, na.rm = TRUE)
}

run_reps_common <- function(T_train, T_test, m, s, q, sigma2, tau, k_grid,
                            n_rep, G_bag, progress_label = "") {
  reps <- vector("list", n_rep)
  for (rr in seq_len(n_rep)) {
    if (rr %% max(1, floor(n_rep / 5)) == 0) {
      cat("    ", progress_label, " replication ", rr, " / ", n_rep, "\n", sep = "")
    }
    reps[[rr]] <- evaluate_one_rep_common(
      T_train = T_train,
      T_test = T_test,
      m = m,
      s = s,
      q = q,
      sigma2 = sigma2,
      tau = tau,
      k_grid = k_grid,
      G_bag = G_bag
    )
  }
  summarise_reps_by_k(reps)
}

###############################################################################
# 4. Population theory and k* calculations
###############################################################################

q_subset_delta_approx <- function(tau, q, k_grid) {
  m <- length(tau)
  tau_bar <- mean(tau)
  delta <- tau / tau_bar - 1
  V_delta <- mean(delta^2)
  B_scale <- 1 / tau_bar

  q + B_scale / k_grid + B_scale * (m - k_grid) / (k_grid^2 * (m - 1)) * V_delta
}

q_subset_leading <- function(tau, q, k_grid) {
  q + (1 / mean(tau)) / k_grid
}

q_bag_inf_smallhet <- function(tau, q, k_grid) {
  ps <- precision_summary(tau, q)
  z <- 1 / k_grid - 1 / length(tau)
  ps$Q_opt + ps$C_delta * z^2
}

q_bag_G_smallhet <- function(tau, q, k_grid, G_bag) {
  ps <- precision_summary(tau, q)
  z <- 1 / k_grid - 1 / length(tau)
  ps$Q_opt + (ps$B_scale / G_bag) * z + (1 - 1 / G_bag) * ps$C_delta * z^2
}

q_population_mc_common <- function(tau, q, k_grid, n_subsets = 5000, G_bag = NULL) {
  # Estimates:
  #   Q_sub(k)      = E[v_S' Sigma v_S]
  #   Q_bag_inf(k) = E[v_S]' Sigma E[v_S]
  #   Q_bag_G(k)   = Q_bag_inf(k) + (Q_sub(k)-Q_bag_inf(k))/G_bag.
  # Also reports diagnostics for the distance between the averaged bagged
  # weights \bar v_k and the full oracle weights w_opt.
  m <- length(tau)
  ps <- precision_summary(tau, q)
  w_opt <- oracle_weights_common(tau)
  delta <- tau / mean(tau) - 1

  out <- data.frame(
    k = k_grid,
    Q_subset_mc = NA_real_,
    Q_bag_inf_mc = NA_real_,
    Q_bag_G_mc = NA_real_,
    R_bag_inf_mc = NA_real_,
    R_weight_exact_mc = NA_real_,
    R_weight_identity_gap = NA_real_,
    R_bag_inf_smallhet = NA_real_,        # legacy name: excess-risk approximation
    R_bag_excess_smallhet = NA_real_,     # preferred name for the same excess risk
    Q_bag_inf_smallhet_level = NA_real_,  # level: Q_opt + R_bag_excess_smallhet
    R_smallhet_abs_error = NA_real_,
    R_smallhet_rel_error = NA_real_,
    R_ratio_mc_to_smallhet = NA_real_,
    w_l2_mc = NA_real_,
    w_l2_smallhet = NA_real_,
    w_max_abs_mc = NA_real_,
    w_max_abs_smallhet = NA_real_,
    w_rmse_approx_error = NA_real_,
    w_max_approx_error = NA_real_,
    w_corr_with_approx = NA_real_
  )

  for (ii in seq_along(k_grid)) {
    k <- k_grid[ii]
    w_sum <- rep(0, m)
    q_sub_vals <- numeric(n_subsets)

    for (bb in seq_len(n_subsets)) {
      idx <- sample.int(m, k, replace = FALSE)
      denom <- sum(tau[idx])
      w <- rep(0, m)
      w[idx] <- tau[idx] / denom
      w_sum <- w_sum + w
      q_sub_vals[bb] <- q + 1 / denom
    }

    w_bar <- w_sum / n_subsets
    diff <- w_bar - w_opt

    # First-order small-heterogeneity approximation from the proposal:
    # \bar v_{i,k} - w_i^opt approximately
    #   - delta_i / m * (m-k)/(k(m-1)).
    diff_approx <- -delta / m * (m - k) / (k * (m - 1))

    Q_sub <- mean(q_sub_vals)
    Q_bag_inf <- q + sum(w_bar^2 / tau)
    Q_bag_G <- if (!is.null(G_bag) && is.finite(G_bag) && G_bag > 0) {
      Q_bag_inf + (Q_sub - Q_bag_inf) / G_bag
    } else {
      NA_real_
    }

    # Exact weight-distance risk. This identity is valid for every level of
    # heterogeneity in the common-loading model, not only for small H:
    #   Q_bag_inf(k) - Q_opt = sum_i (vbar_i(k) - w_i^opt)^2 / tau_i.
    # The small-heterogeneity expression below is only an approximation to this
    # exact quantity.
    R_mc <- Q_bag_inf - ps$Q_opt
    R_weight_exact <- sum(diff^2 / tau)
    z <- 1 / k - 1 / m
    R_smallhet <- ps$C_delta * z^2

    out$Q_subset_mc[ii] <- Q_sub
    out$Q_bag_inf_mc[ii] <- Q_bag_inf
    out$Q_bag_G_mc[ii] <- Q_bag_G
    out$R_bag_inf_mc[ii] <- R_mc
    out$R_weight_exact_mc[ii] <- R_weight_exact
    out$R_weight_identity_gap[ii] <- R_mc - R_weight_exact
    out$R_bag_inf_smallhet[ii] <- R_smallhet
    out$R_bag_excess_smallhet[ii] <- R_smallhet
    out$Q_bag_inf_smallhet_level[ii] <- ps$Q_opt + R_smallhet
    out$R_smallhet_abs_error[ii] <- R_mc - R_smallhet
    out$R_smallhet_rel_error[ii] <- ifelse(abs(R_mc) > 1e-14, (R_mc - R_smallhet) / R_mc, NA_real_)
    out$R_ratio_mc_to_smallhet[ii] <- ifelse(abs(R_smallhet) > 1e-14, R_mc / R_smallhet, NA_real_)
    out$w_l2_mc[ii] <- sqrt(sum(diff^2))
    out$w_l2_smallhet[ii] <- sqrt(sum(diff_approx^2))
    out$w_max_abs_mc[ii] <- max(abs(diff))
    out$w_max_abs_smallhet[ii] <- max(abs(diff_approx))
    out$w_rmse_approx_error[ii] <- sqrt(mean((diff - diff_approx)^2))
    out$w_max_approx_error[ii] <- max(abs(diff - diff_approx))
    diff_sd <- sd(diff)
    diff_approx_sd <- sd(diff_approx)
    diff_norm <- sqrt(sum(diff^2))
    diff_approx_norm <- sqrt(sum(diff_approx^2))
    out$w_corr_with_approx[ii] <- if (
      k < m &&
      is.finite(diff_sd) && is.finite(diff_approx_sd) &&
      is.finite(diff_norm) && is.finite(diff_approx_norm) &&
      diff_sd > WEIGHT_CORR_TOL &&
      diff_approx_sd > WEIGHT_CORR_TOL &&
      diff_norm > WEIGHT_CORR_TOL &&
      diff_approx_norm > WEIGHT_CORR_TOL
    ) {
      suppressWarnings(cor(diff, diff_approx))
    } else {
      NA_real_
    }
  }

  out
}


# Generic version for arbitrary covariance matrices. This is not needed for the
# common-loading baseline, but it lets you estimate the same object for Toeplitz,
# heterogeneous-loading factor, random covariance, or any other positive-definite
# Sigma_e design:
#   R_weight(k) = (vbar_k - w_opt)' Sigma_e (vbar_k - w_opt)
#               = Q_bag_inf(k) - Q_opt.
# No small-heterogeneity assumption is used here.
q_population_mc_general <- function(Sigma_e, k_grid, n_subsets = 5000, G_bag = NULL) {
  Sigma_e <- as.matrix(Sigma_e)
  m <- ncol(Sigma_e)
  one <- rep(1, m)

  Sinv_one <- as.numeric(safe_solve(Sigma_e, one))
  w_opt <- Sinv_one / sum(Sinv_one)
  w_equal <- rep(1 / m, m)
  Q_opt <- as.numeric(t(w_opt) %*% Sigma_e %*% w_opt)
  Q_ave <- as.numeric(t(w_equal) %*% Sigma_e %*% w_equal)

  out <- data.frame(
    k = k_grid,
    Q_subset_mc = NA_real_,
    Q_bag_inf_mc = NA_real_,
    Q_bag_G_mc = NA_real_,
    Q_opt = Q_opt,
    Q_ave = Q_ave,
    R_weight_exact_mc = NA_real_,
    R_bag_inf_mc = NA_real_,
    R_weight_identity_gap = NA_real_,
    w_l2_mc = NA_real_,
    w_max_abs_mc = NA_real_
  )

  for (ii in seq_along(k_grid)) {
    k <- as.integer(k_grid[ii])
    w_sum <- rep(0, m)
    q_sub_vals <- numeric(n_subsets)

    for (bb in seq_len(n_subsets)) {
      idx <- sample.int(m, k, replace = FALSE)
      Sigma_S <- Sigma_e[idx, idx, drop = FALSE]
      one_S <- rep(1, k)
      Sinv1_S <- as.numeric(safe_solve(Sigma_S, one_S))
      w_S <- Sinv1_S / sum(Sinv1_S)
      w_emb <- rep(0, m)
      w_emb[idx] <- w_S
      w_sum <- w_sum + w_emb
      q_sub_vals[bb] <- as.numeric(t(w_emb) %*% Sigma_e %*% w_emb)
    }

    w_bar <- w_sum / n_subsets
    diff <- w_bar - w_opt
    Q_sub <- mean(q_sub_vals)
    Q_bag_inf <- as.numeric(t(w_bar) %*% Sigma_e %*% w_bar)
    R_weight <- as.numeric(t(diff) %*% Sigma_e %*% diff)
    Q_bag_G <- if (!is.null(G_bag) && is.finite(G_bag) && G_bag > 0) {
      Q_bag_inf + (Q_sub - Q_bag_inf) / G_bag
    } else {
      NA_real_
    }

    out$Q_subset_mc[ii] <- Q_sub
    out$Q_bag_inf_mc[ii] <- Q_bag_inf
    out$Q_bag_G_mc[ii] <- Q_bag_G
    out$R_weight_exact_mc[ii] <- R_weight
    out$R_bag_inf_mc[ii] <- Q_bag_inf - Q_opt
    # Identity check. Since w_bar and w_opt both sum to one and w_opt solves
    # min_w w' Sigma_e w subject to 1'w = 1, the linear term vanishes:
    #   Q(w_bar) - Q(w_opt) = (w_bar - w_opt)' Sigma_e (w_bar - w_opt).
    # Hence R_weight_identity_gap should be close to zero; nonzero values mostly
    # indicate Monte Carlo or numerical error.
    out$R_weight_identity_gap[ii] <- (Q_bag_inf - Q_opt) - R_weight
    out$w_l2_mc[ii] <- sqrt(sum(diff^2))
    out$w_max_abs_mc[ii] <- max(abs(diff))
  }

  out
}

bagged_weights_mc_common <- function(tau, k, n_subsets = 10000) {
  m <- length(tau)
  w_sum <- rep(0, m)
  for (bb in seq_len(n_subsets)) {
    idx <- sample.int(m, k, replace = FALSE)
    denom <- sum(tau[idx])
    w <- rep(0, m)
    w[idx] <- tau[idx] / denom
    w_sum <- w_sum + w
  }
  w_sum / n_subsets
}

weight_details_common <- function(tau, k_grid, n_subsets = 10000) {
  # Per-forecast diagnostic table for selected k values. This is the most direct
  # way to inspect whether the small-heterogeneity weight approximation is good.
  m <- length(tau)
  tau_bar <- mean(tau)
  delta <- tau / tau_bar - 1
  sigma2 <- 1 / tau
  w_opt <- oracle_weights_common(tau)

  rows <- vector("list", length(k_grid))
  for (jj in seq_along(k_grid)) {
    k <- as.integer(k_grid[jj])
    w_bag <- bagged_weights_mc_common(tau, k, n_subsets = n_subsets)
    diff_approx <- -delta / m * (m - k) / (k * (m - 1))
    w_bag_approx <- w_opt + diff_approx

    rows[[jj]] <- data.frame(
      k = k,
      i = seq_len(m),
      tau = tau,
      sigma2 = sigma2,
      delta = delta,
      w_opt = w_opt,
      w_bag_mc = w_bag,
      w_bag_smallhet_approx = w_bag_approx,
      diff_mc = w_bag - w_opt,
      diff_smallhet_approx = diff_approx,
      diff_approx_error = (w_bag - w_opt) - diff_approx
    )
  }

  do.call(rbind, rows)
}

loss_from_Q <- function(s, Q, T_train, k) {
  (s + Q) * inflation(T_train, k)
}

kstar_subset_leading <- function(s, q, tau, T_train, lower, upper) {
  B_scale <- 1 / mean(tau)
  A0 <- s + q
  k_cont <- (-B_scale + sqrt(B_scale^2 + A0 * B_scale * T_train)) / A0
  list(
    k_cont = k_cont,
    k_int = nearest_feasible_integer(k_cont, lower, upper)
  )
}

kstar_bag_cube_root_local <- function(s, q, tau, T_train, lower, upper) {
  ps <- precision_summary(tau, q)
  A <- s + ps$Q_opt
  if (!is.finite(ps$C_delta) || ps$C_delta <= 0) {
    return(list(k_cont = lower, k_int = lower))
  }
  k_cont <- (2 * ps$C_delta * T_train / A)^(1 / 3)
  list(
    k_cont = k_cont,
    k_int = nearest_feasible_integer(k_cont, lower, upper)
  )
}

make_theory_curves <- function(tau, q, s, T_train, k_grid, G_bag) {
  ps <- precision_summary(tau, q)
  Q_bag_inf_sh <- q_bag_inf_smallhet(tau, q, k_grid)
  Q_bag_G_sh <- q_bag_G_smallhet(tau, q, k_grid, G_bag)
  Q_sub_lead <- q_subset_leading(tau, q, k_grid)
  Q_sub_delta <- q_subset_delta_approx(tau, q, k_grid)

  out <- data.frame(
    k = k_grid,
    Q_bag_inf_smallhet = Q_bag_inf_sh,
    Q_bag_G_smallhet = Q_bag_G_sh,
    R_bag_excess_smallhet = Q_bag_inf_sh - ps$Q_opt,
    R_bag_G_excess_smallhet = Q_bag_G_sh - ps$Q_opt,
    Q_subset_leading = Q_sub_lead,
    Q_subset_delta = Q_sub_delta,
    L_bag_inf_smallhet = loss_from_Q(s, Q_bag_inf_sh, T_train, k_grid),
    L_bag_G_smallhet = loss_from_Q(s, Q_bag_G_sh, T_train, k_grid),
    L_subset_leading = loss_from_Q(s, Q_sub_lead, T_train, k_grid),
    L_subset_delta = loss_from_Q(s, Q_sub_delta, T_train, k_grid),
    L_equal = s + ps$Q_ave,
    L_oracle_optimal = s + ps$Q_opt,
    L_full_estimated_theory = if (T_train > length(tau)) {
      (s + ps$Q_opt) * (T_train - 1) / (T_train - length(tau))
    } else {
      NA_real_
    }
  )

  out
}

theory_summary_one <- function(tau, q, s, T_train, k_grid, G_bag) {
  m <- length(tau)
  ps <- precision_summary(tau, q)
  curves <- make_theory_curves(tau, q, s, T_train, k_grid, G_bag)
  ks_sub <- kstar_subset_leading(s, q, tau, T_train, min(k_grid), max(k_grid))
  ks_cube <- kstar_bag_cube_root_local(s, q, tau, T_train, min(k_grid), max(k_grid))

  data.frame(
    m = m,
    T_train = T_train,
    s = s,
    q = q,
    B_scale = ps$B_scale,
    H = ps$H,
    V_delta = ps$V_delta,
    C_delta = ps$C_delta,
    Q_ave = ps$Q_ave,
    Q_opt = ps$Q_opt,
    L_equal = s + ps$Q_ave,
    L_oracle_optimal = s + ps$Q_opt,
    L_full_estimated_theory = curves$L_full_estimated_theory[1],
    kstar_subset_sqrt_cont = ks_sub$k_cont,
    kstar_subset_sqrt_int = ks_sub$k_int,
    kstar_bag_cube_local_cont = ks_cube$k_cont,
    kstar_bag_cube_local_int = ks_cube$k_int,
    kstar_bag_inf_grid = find_kstar_grid(curves, "L_bag_inf_smallhet"),
    kstar_bag_G_grid = find_kstar_grid(curves, "L_bag_G_smallhet"),
    kstar_subset_grid = find_kstar_grid(curves, "L_subset_leading")
  )
}

###############################################################################
# 5. Plotting functions with fewer curves
###############################################################################

plot_clean_mspe_bagged_vs_benchmarks <- function(sum_df, title, filename) {
  if (!isTRUE(SAVE_PLOTS)) return(invisible(NULL))
  png(filename, width = 1000, height = 680)
  ylim <- range(c(
    sum_df$rsm_bagged,
    sum_df$equal_weight,
    sum_df$full_estimated,
    sum_df$oracle_optimal
  ), na.rm = TRUE)

  plot(sum_df$k, sum_df$rsm_bagged, type = "l", lwd = 2.5,
       col = PLOT_COL["bagged"], ylim = ylim,
       xlab = "subset size k", ylab = "out-of-sample MSPE", main = title)
  abline(h = mean(sum_df$equal_weight, na.rm = TRUE), lty = 2, lwd = 2,
         col = PLOT_COL["equal"])
  abline(h = mean(sum_df$oracle_optimal, na.rm = TRUE), lty = 3, lwd = 2,
         col = PLOT_COL["oracle"])
  abline(h = mean(sum_df$full_estimated, na.rm = TRUE), lty = 4, lwd = 2,
         col = PLOT_COL["full"])
  abline(v = find_kstar_grid(sum_df, "rsm_bagged"), lty = 5, lwd = 1.5,
         col = PLOT_COL["bagged"])

  legend("topright",
         legend = c("bagged RSM", "equal weights", "oracle optimal", "estimated full optimal", "bagged empirical k*"),
         lty = c(1, 2, 3, 4, 5),
         col = c(PLOT_COL["bagged"], PLOT_COL["equal"], PLOT_COL["oracle"], PLOT_COL["full"], PLOT_COL["bagged"]),
         lwd = c(2.5, 2, 2, 2, 1.5), bty = "n")
  dev.off()
}

plot_subset_vs_bagged <- function(sum_df, title, filename) {
  if (!isTRUE(SAVE_PLOTS)) return(invisible(NULL))
  png(filename, width = 1000, height = 680)
  ylim <- range(c(sum_df$rsm_bagged, sum_df$rsm_subset_avg), na.rm = TRUE)
  plot(sum_df$k, sum_df$rsm_bagged, type = "l", lwd = 2.5,
       col = PLOT_COL["bagged"], ylim = ylim,
       xlab = "subset size k", ylab = "out-of-sample MSPE", main = title)
  lines(sum_df$k, sum_df$rsm_subset_avg, lty = 2, lwd = 2.5,
        col = PLOT_COL["subset"])
  legend("topright", legend = c("bagged RSM", "average subset RSM"),
         lty = c(1, 2), col = c(PLOT_COL["bagged"], PLOT_COL["subset"]),
         lwd = 2.5, bty = "n")
  dev.off()
}

plot_population_bagged <- function(pop_df, theory_df, ps, title, filename) {
  if (!isTRUE(SAVE_PLOTS)) return(invisible(NULL))
  png(filename, width = 1000, height = 680)
  ylim <- range(c(
    pop_df$Q_bag_inf_mc,
    pop_df$Q_bag_G_mc,
    theory_df$Q_bag_inf_smallhet,
    ps$Q_ave,
    ps$Q_opt
  ), na.rm = TRUE)

  plot(pop_df$k, pop_df$Q_bag_inf_mc, type = "l", lwd = 2.5,
       col = PLOT_COL["bagged"], ylim = ylim,
       xlab = "subset size k", ylab = "population component Q", main = title)
  lines(pop_df$k, pop_df$Q_bag_G_mc, lty = 2, lwd = 2,
        col = PLOT_COL["finite"])
  lines(theory_df$k, theory_df$Q_bag_inf_smallhet, lty = 3, lwd = 2,
        col = PLOT_COL["theory"])
  abline(h = ps$Q_ave, lty = 4, lwd = 1.5, col = PLOT_COL["equal"])
  abline(h = ps$Q_opt, lty = 5, lwd = 1.5, col = PLOT_COL["oracle"])

  legend("topright",
         legend = c("Q_bag,inf population MC", "Q_bag,G finite-bag bridge", "small-heterogeneity approximation", "Q_ave", "Q_opt"),
         lty = c(1, 2, 3, 4, 5),
         col = c(PLOT_COL["bagged"], PLOT_COL["finite"], PLOT_COL["theory"], PLOT_COL["equal"], PLOT_COL["oracle"]),
         lwd = c(2.5, 2, 2, 1.5, 1.5), bty = "n")
  dev.off()
}

plot_weight_diagnostics <- function(pop_df, title, filename) {
  if (!isTRUE(SAVE_PLOTS)) return(invisible(NULL))
  if (!all(c("k", "R_weight_exact_mc", "R_bag_excess_smallhet") %in% names(pop_df))) return(invisible(NULL))

  png(filename, width = 1000, height = 680)
  ylim <- range(c(pop_df$R_weight_exact_mc, pop_df$R_bag_excess_smallhet), na.rm = TRUE)
  plot(pop_df$k, pop_df$R_weight_exact_mc, type = "l", lwd = 2.5,
       col = PLOT_COL["bagged"], ylim = ylim,
       xlab = "subset size k",
       ylab = "weight-distance excess risk",
       main = title)
  lines(pop_df$k, pop_df$R_bag_excess_smallhet, lty = 2, lwd = 2.5,
        col = PLOT_COL["theory"])
  legend("topright",
         legend = c("MC weight-distance risk", "small-heterogeneity approximation"),
         lty = c(1, 2), col = c(PLOT_COL["bagged"], PLOT_COL["theory"]),
         lwd = 2.5, bty = "n")
  dev.off()
}

plot_weight_approx_error <- function(pop_df, title, filename) {
  if (!isTRUE(SAVE_PLOTS)) return(invisible(NULL))
  if (!all(c("k", "w_rmse_approx_error", "w_corr_with_approx") %in% names(pop_df))) return(invisible(NULL))

  png(filename, width = 1000, height = 680)
  plot(pop_df$k, pop_df$w_rmse_approx_error, type = "l", lwd = 2.5,
       col = PLOT_COL["bagged"],
       xlab = "subset size k",
       ylab = "RMSE of weight-difference approximation",
       main = title)
  dev.off()
}

plot_joint_heatmap <- function(df, value_col, filename, title, zlab = value_col) {
  if (!isTRUE(SAVE_PLOTS)) return(invisible(NULL))
  x_col <- if ("B_target" %in% names(df)) "B_target" else "B_scale"
  y_col <- if ("H_target" %in% names(df)) "H_target" else "H"
  if (!all(c(x_col, y_col, value_col) %in% names(df))) return(invisible(NULL))

  xs <- sort(unique(df[[x_col]]))
  ys <- sort(unique(df[[y_col]]))
  z <- matrix(NA_real_, nrow = length(xs), ncol = length(ys))
  for (ii in seq_along(xs)) {
    for (jj in seq_along(ys)) {
      hit <- df[df[[x_col]] == xs[ii] & df[[y_col]] == ys[jj], value_col]
      if (length(hit) > 0) z[ii, jj] <- hit[1]
    }
  }

  png(filename, width = 1000, height = 760)
  image(xs, ys, z, xlab = "target B_scale", ylab = "target H", main = title)
  contour(xs, ys, z, add = TRUE, drawlabels = TRUE)
  grid_xy <- expand.grid(xs, ys)
  text(grid_xy[, 1], grid_xy[, 2], labels = round(as.vector(z), 3), cex = 0.75)
  mtext(zlab, side = 4, line = 2.5)
  dev.off()
}

plot_scaling_T_by_m <- function(df, filename) {
  if (!isTRUE(SAVE_PLOTS)) return(invisible(NULL))
  png(filename, width = 1000, height = 680)
  ms <- sort(unique(df$m))
  ylim <- range(c(df$kstar_bag_inf_grid, df$kstar_empirical_bagged), na.rm = TRUE)
  xlim <- range(df$T_train^(1 / 3), na.rm = TRUE)

  plot(NA, NA, xlim = xlim, ylim = ylim,
       xlab = expression(T^{1/3}), ylab = "optimal subset size k*",
       main = "Bagged k*: scaling in T, separated by m")

  ltys <- seq_along(ms)
  for (ii in seq_along(ms)) {
    dd <- df[df$m == ms[ii], ]
    dd <- dd[order(dd$T_train), ]
    lines(dd$T_train^(1 / 3), dd$kstar_bag_inf_grid, lty = ltys[ii], lwd = 2,
          col = PLOT_COL["theory"])
    points(dd$T_train^(1 / 3), dd$kstar_empirical_bagged, pch = 16 + ii,
           col = PLOT_COL["bagged"])
  }

  legend("topleft",
         legend = paste("m =", ms), lty = ltys, lwd = 2,
         col = PLOT_COL["theory"], bty = "n")
  legend("bottomright",
         legend = c("line: theory grid k*", "points: empirical bagged k*"),
         lty = c(1, NA), pch = c(NA, 16),
         col = c(PLOT_COL["theory"], PLOT_COL["bagged"]),
         lwd = c(2, NA), bty = "n")
  dev.off()
}

plot_sensitivity <- function(df, xvar, yvar, filename, title, xlab) {
  if (!isTRUE(SAVE_PLOTS)) return(invisible(NULL))
  png(filename, width = 1000, height = 680)
  plot(df[[xvar]], df[[yvar]], type = "b", pch = 19, lwd = 2,
       col = PLOT_COL["bagged"], xlab = xlab, ylab = yvar, main = title)
  dev.off()
}


value_at_k <- function(df, k_value, column) {
  if (!column %in% names(df) || !is.finite(k_value)) return(NA_real_)
  hit <- df[df$k == k_value, column]
  if (length(hit) == 0) return(NA_real_)
  hit[1]
}

###############################################################################
# 6. Experiments
###############################################################################

run_baseline <- function() {
  cat("\nRunning baseline common-loading experiment...\n")

  base_grid <- expand.grid(
    B_scale = BASELINE$B_scale,
    H_target = BASELINE$H_target,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  all_summary <- data.frame()

  for (cc in seq_len(nrow(base_grid))) {
    B_scale <- as.numeric(base_grid$B_scale[cc])
    H_target <- as.numeric(base_grid$H_target[cc])
    m <- BASELINE$m
    T_train <- BASELINE$T_train
    s <- BASELINE$s
    q <- BASELINE$q
    scenario_id <- make_scenario_id("baseline", m, T_train, B_scale, H_target)

    cat("  baseline combo ", cc, " / ", nrow(base_grid),
        ": B=", B_scale, ", H=", H_target, "\n", sep = "")

    design <- make_tau_lognormal(
      m, B_scale = B_scale, H_target = H_target,
      seed = GLOBAL_SEED + 101 + round(100 * B_scale) + round(1000 * H_target)
    )
    tau <- design$tau
    sigma2 <- design$sigma2
    k_grid <- make_k_grid(m, T_train, BASELINE$k_min, BASELINE$k_max, BASELINE$k_by)

    ps <- precision_summary(tau, q)
    cat("    Actual B_scale = ", round(ps$B_scale, 4), "\n", sep = "")
    cat("    Actual H       = ", round(ps$H, 4), "\n", sep = "")
    cat("    Actual V_delta = ", round(ps$V_delta, 4), "\n", sep = "")

    theory_df <- make_theory_curves(tau, q, s, T_train, k_grid, G_BAG_MC)
    theory_sum <- theory_summary_one(tau, q, s, T_train, k_grid, G_BAG_MC)

    pop_df <- q_population_mc_common(
      tau = tau,
      q = q,
      k_grid = k_grid,
      n_subsets = N_SUBSETS_POP,
      G_bag = G_BAG_MC
    )
    pop_df <- merge(pop_df, theory_df, by = "k", all.x = TRUE)

    # Per-forecast weight differences for selected k values.
    weight_k <- sort(unique(WEIGHT_DIAG_K[WEIGHT_DIAG_K %in% k_grid]))
    if (length(weight_k) > 0) {
      weight_details <- weight_details_common(
        tau = tau,
        k_grid = weight_k,
        n_subsets = N_SUBSETS_POP
      )
    } else {
      weight_details <- data.frame()
    }

    sum_df <- run_reps_common(
      T_train = T_train,
      T_test = T_TEST,
      m = m,
      s = s,
      q = q,
      sigma2 = sigma2,
      tau = tau,
      k_grid = k_grid,
      n_rep = N_REP_BASE,
      G_bag = G_BAG_MC,
      progress_label = scenario_id
    )
    sum_df <- merge(sum_df, theory_df, by = "k", all.x = TRUE)
    emp_k_bag <- find_kstar_grid(sum_df, "rsm_bagged")
    emp_k_sub <- find_kstar_grid(sum_df, "rsm_subset_avg")

    summary <- cbind(
      data.frame(
        scenario_id = scenario_id,
        B_target = B_scale,
        H_target = H_target
      ),
      theory_sum,
      data.frame(
        n_rep = N_REP_BASE,
        G_bag = G_BAG_MC,
        N_subsets_pop = N_SUBSETS_POP,
        kstar_empirical_bagged = emp_k_bag,
        kstar_empirical_subset_avg = emp_k_sub,
        min_bagged_mspe = min(sum_df$rsm_bagged, na.rm = TRUE),
        min_subset_avg_mspe = min(sum_df$rsm_subset_avg, na.rm = TRUE),
        equal_weight_mspe = mean(sum_df$equal_weight, na.rm = TRUE),
        oracle_optimal_mspe = mean(sum_df$oracle_optimal, na.rm = TRUE),
        full_estimated_mspe = mean(sum_df$full_estimated, na.rm = TRUE),
        bagged_gain_vs_equal = mean(sum_df$equal_weight, na.rm = TRUE) - min(sum_df$rsm_bagged, na.rm = TRUE),
        bagged_gap_vs_oracle = min(sum_df$rsm_bagged, na.rm = TRUE) - mean(sum_df$oracle_optimal, na.rm = TRUE),
        min_weight_R_mc = min(pop_df$R_weight_exact_mc, na.rm = TRUE),
        max_abs_weight_identity_gap = max(abs(pop_df$R_weight_identity_gap), na.rm = TRUE),
        mean_weight_R_ratio = mean(pop_df$R_ratio_mc_to_smallhet, na.rm = TRUE),
        mean_smallhet_rel_error = mean(pop_df$R_smallhet_rel_error, na.rm = TRUE),
        median_weight_corr_with_approx = median(pop_df$w_corr_with_approx, na.rm = TRUE)
      )
    )

    all_summary <- rbind(all_summary, summary)

    safe_write_csv(sum_df, file.path(OUTDIR, paste0(scenario_id, "_finite_sample_by_k.csv")))
    safe_write_csv(pop_df, file.path(OUTDIR, paste0(scenario_id, "_population_by_k.csv")))
    if (nrow(weight_details) > 0) {
      safe_write_csv(weight_details, file.path(OUTDIR, paste0(scenario_id, "_weight_details_by_k.csv")))
    }

    # Backward-compatible names for the first baseline combo.
    if (cc == 1) {
      safe_write_csv(sum_df, file.path(OUTDIR, "baseline_finite_sample_by_k.csv"))
      safe_write_csv(pop_df, file.path(OUTDIR, "baseline_population_by_k.csv"))
      if (nrow(weight_details) > 0) {
        safe_write_csv(weight_details, file.path(OUTDIR, "baseline_weight_details_by_k.csv"))
      }
    }

    plot_clean_mspe_bagged_vs_benchmarks(
      sum_df,
      title = paste0("Common loading: bagged RSM versus benchmarks\nB=", round(ps$B_scale, 3), ", H=", round(ps$H, 3)),
      filename = file.path(OUTDIR, paste0(scenario_id, "_bagged_vs_benchmarks.png"))
    )
    plot_subset_vs_bagged(
      sum_df,
      title = paste0("Common loading: bagged RSM versus average subset RSM\nB=", round(ps$B_scale, 3), ", H=", round(ps$H, 3)),
      filename = file.path(OUTDIR, paste0(scenario_id, "_bagged_vs_subset.png"))
    )
    plot_population_bagged(
      pop_df,
      theory_df,
      ps,
      title = paste0("Common loading: population bagged Q objects\nB=", round(ps$B_scale, 3), ", H=", round(ps$H, 3)),
      filename = file.path(OUTDIR, paste0(scenario_id, "_population_bagged_Q.png"))
    )
    plot_weight_diagnostics(
      pop_df,
      title = paste0("Weight-distance diagnostic: bagged average versus optimal weights\nB=", round(ps$B_scale, 3), ", H=", round(ps$H, 3)),
      filename = file.path(OUTDIR, paste0(scenario_id, "_weight_distance_Q.png"))
    )
    plot_weight_approx_error(
      pop_df,
      title = paste0("Weight approximation error by k\nB=", round(ps$B_scale, 3), ", H=", round(ps$H, 3)),
      filename = file.path(OUTDIR, paste0(scenario_id, "_weight_approx_error.png"))
    )

    if (cc == 1) {
      plot_clean_mspe_bagged_vs_benchmarks(
        sum_df,
        title = "Common loading: bagged RSM versus benchmarks",
        filename = file.path(OUTDIR, "baseline_bagged_vs_benchmarks.png")
      )
      plot_subset_vs_bagged(
        sum_df,
        title = "Common loading: bagged RSM versus average subset RSM",
        filename = file.path(OUTDIR, "baseline_bagged_vs_subset.png")
      )
      plot_population_bagged(
        pop_df,
        theory_df,
        ps,
        title = "Common loading: population bagged Q objects",
        filename = file.path(OUTDIR, "baseline_population_bagged_Q.png")
      )
      plot_weight_diagnostics(
        pop_df,
        title = "Weight-distance diagnostic: bagged average versus optimal weights",
        filename = file.path(OUTDIR, "baseline_weight_distance_Q.png")
      )
    }
  }

  safe_write_csv(all_summary, file.path(OUTDIR, "baseline_grid_summary.csv"))
  # Backward-compatible name if only one baseline combo is run.
  if (nrow(all_summary) == 1) {
    safe_write_csv(all_summary, file.path(OUTDIR, "baseline_summary.csv"))
  }

  all_summary
}

run_scaling_T_m <- function() {
  cat("\nRunning scaling experiment over T and m...\n")

  combos <- expand_grid_list(GRID_SCALING)
  all_summary <- data.frame()

  for (cc in seq_len(nrow(combos))) {
    combo <- combos[cc, , drop = FALSE]
    m <- as.integer(combo$m)
    T_train <- as.integer(combo$T_train)
    s <- as.numeric(combo$s)
    q <- as.numeric(combo$q)
    B_scale <- as.numeric(combo$B_scale)
    H_target <- as.numeric(combo$H_target)
    k_min <- as.integer(combo$k_min)
    k_by <- as.integer(combo$k_by)

    cat("  combo ", cc, " / ", nrow(combos), ": m=", m, ", T=", T_train, "\n", sep = "")

    # The same cross-sectional design is used for all T for a given m, B, H.
    seed_design <- GLOBAL_SEED + 2000 + m + round(100 * B_scale) + round(1000 * H_target)
    design <- make_tau_lognormal(m, B_scale = B_scale, H_target = H_target,
                                 seed = seed_design)
    tau <- design$tau
    sigma2 <- design$sigma2
    k_grid <- make_k_grid(m, T_train, k_min = k_min, k_by = k_by)
    if (length(k_grid) == 0) next

    theory_sum <- theory_summary_one(tau, q, s, T_train, k_grid, G_BAG_MC)

    sum_df <- run_reps_common(
      T_train = T_train,
      T_test = T_TEST,
      m = m,
      s = s,
      q = q,
      sigma2 = sigma2,
      tau = tau,
      k_grid = k_grid,
      n_rep = N_REP_SCALING,
      G_bag = G_BAG_MC,
      progress_label = paste0("scaling m=", m, " T=", T_train)
    )

    one <- cbind(
      theory_sum,
      data.frame(
        n_rep = N_REP_SCALING,
        G_bag = G_BAG_MC,
        kstar_empirical_bagged = find_kstar_grid(sum_df, "rsm_bagged"),
        kstar_empirical_subset_avg = find_kstar_grid(sum_df, "rsm_subset_avg"),
        min_bagged_mspe = min(sum_df$rsm_bagged, na.rm = TRUE),
        equal_weight_mspe = mean(sum_df$equal_weight, na.rm = TRUE),
        oracle_optimal_mspe = mean(sum_df$oracle_optimal, na.rm = TRUE),
        full_estimated_mspe = mean(sum_df$full_estimated, na.rm = TRUE)
      )
    )
    all_summary <- rbind(all_summary, one)

    safe_write_csv(
      sum_df,
      file.path(OUTDIR, paste0("scaling_finite_sample_m", m, "_T", T_train, ".csv"))
    )
  }

  safe_write_csv(all_summary, file.path(OUTDIR, "scaling_T_m_summary.csv"))

  # Log-log checks. The proposal predicts beta_T about 1/3 and beta_m about -1/3
  # for infinite-bagged small-heterogeneity theory when B, V_delta, s, and q are stable.
  reg_theory <- lm(log(kstar_bag_inf_grid) ~ log(T_train) + log(m), data = all_summary)
  reg_emp <- NULL
  if (all(is.finite(all_summary$kstar_empirical_bagged)) && nrow(all_summary) >= 4) {
    reg_emp <- lm(log(kstar_empirical_bagged) ~ log(T_train) + log(m), data = all_summary)
  }

  sink(file.path(OUTDIR, "scaling_loglog_regressions.txt"))
  cat("Theory-grid regression: log k*_bag,inf ~ log T + log m\n")
  print(summary(reg_theory))
  cat("\nTarget from proposal under stable B, V_delta, s, q: beta_T = 1/3, beta_m = -1/3\n")
  if (!is.null(reg_emp)) {
    cat("\nEmpirical bagged regression: log k*_empirical ~ log T + log m\n")
    print(summary(reg_emp))
  } else {
    cat("\nEmpirical regression skipped because empirical k* has missing or non-finite values.\n")
  }
  sink()

  plot_scaling_T_by_m(all_summary, file.path(OUTDIR, "scaling_kstar_T_by_m.png"))
  all_summary
}

run_joint_B_H_grid <- function() {
  cat("\nRunning joint B_scale x H grid experiment...\n")

  combos <- expand_grid_list(GRID_JOINT_B_H)
  all_summary <- data.frame()
  all_pop <- data.frame()

  for (cc in seq_len(nrow(combos))) {
    combo <- combos[cc, , drop = FALSE]
    m <- as.integer(combo$m)
    T_train <- as.integer(combo$T_train)
    s <- as.numeric(combo$s)
    q <- as.numeric(combo$q)
    B_scale <- as.numeric(combo$B_scale)
    H_target <- as.numeric(combo$H_target)
    k_min <- as.integer(combo$k_min)
    k_by <- as.integer(combo$k_by)

    cat("  joint grid combo ", cc, " / ", nrow(combos),
        ": B=", B_scale, ", H=", H_target, "\n", sep = "")

    design <- make_tau_lognormal(
      m, B_scale = B_scale, H_target = H_target,
      seed = GLOBAL_SEED + 5000 + round(100 * B_scale) + round(1000 * H_target)
    )
    tau <- design$tau
    sigma2 <- design$sigma2
    k_grid <- make_k_grid(m, T_train, k_min = k_min, k_by = k_by)
    if (length(k_grid) == 0) next

    theory_sum <- theory_summary_one(tau, q, s, T_train, k_grid, G_BAG_MC)

    # Population + weight-distance diagnostics are cheap relative to the
    # finite-sample simulation and help explain when the approximation is good.
    pop_df <- q_population_mc_common(
      tau = tau,
      q = q,
      k_grid = k_grid,
      n_subsets = N_SUBSETS_POP_JOINT,
      G_bag = G_BAG_MC
    )
    scenario_id <- make_scenario_id("joint_BH", m, T_train, B_scale, H_target)
    pop_df_labeled <- cbind(
      data.frame(
        scenario_id = scenario_id,
        B_target = B_scale,
        H_target = H_target,
        B_actual = theory_sum$B_scale[1],
        H_actual = theory_sum$H[1],
        V_delta = theory_sum$V_delta[1]
      ),
      pop_df
    )
    all_pop <- rbind(all_pop, pop_df_labeled)

    sum_df <- run_reps_common(
      T_train = T_train,
      T_test = T_TEST,
      m = m,
      s = s,
      q = q,
      sigma2 = sigma2,
      tau = tau,
      k_grid = k_grid,
      n_rep = N_REP_GRID,
      G_bag = G_BAG_MC,
      progress_label = paste0("joint B=", B_scale, " H=", H_target)
    )

    emp_k_bag <- find_kstar_grid(sum_df, "rsm_bagged")
    emp_k_sub <- find_kstar_grid(sum_df, "rsm_subset_avg")

    one <- cbind(
      data.frame(
        B_target = B_scale,
        H_target = H_target
      ),
      theory_sum,
      data.frame(
        n_rep = N_REP_GRID,
        G_bag = G_BAG_MC,
        N_subsets_pop_joint = N_SUBSETS_POP_JOINT,
        kstar_empirical_bagged = emp_k_bag,
        kstar_empirical_subset_avg = emp_k_sub,
        min_bagged_mspe = min(sum_df$rsm_bagged, na.rm = TRUE),
        min_subset_avg_mspe = min(sum_df$rsm_subset_avg, na.rm = TRUE),
        equal_weight_mspe = mean(sum_df$equal_weight, na.rm = TRUE),
        oracle_optimal_mspe = mean(sum_df$oracle_optimal, na.rm = TRUE),
        full_estimated_mspe = mean(sum_df$full_estimated, na.rm = TRUE),
        bagged_gain_vs_equal = mean(sum_df$equal_weight, na.rm = TRUE) - min(sum_df$rsm_bagged, na.rm = TRUE),
        bagged_gap_vs_oracle = min(sum_df$rsm_bagged, na.rm = TRUE) - mean(sum_df$oracle_optimal, na.rm = TRUE),
        R_weight_at_empirical_bagged_k = value_at_k(pop_df, emp_k_bag, "R_weight_exact_mc"),
        R_smallhet_at_empirical_bagged_k = value_at_k(pop_df, emp_k_bag, "R_bag_excess_smallhet"),
        R_ratio_at_empirical_bagged_k = value_at_k(pop_df, emp_k_bag, "R_ratio_mc_to_smallhet"),
        R_smallhet_rel_error_at_empirical_bagged_k = value_at_k(pop_df, emp_k_bag, "R_smallhet_rel_error"),
        max_abs_weight_identity_gap = max(abs(pop_df$R_weight_identity_gap), na.rm = TRUE),
        mean_weight_R_ratio = mean(pop_df$R_ratio_mc_to_smallhet, na.rm = TRUE),
        mean_smallhet_rel_error = mean(pop_df$R_smallhet_rel_error, na.rm = TRUE),
        median_weight_corr_with_approx = median(pop_df$w_corr_with_approx, na.rm = TRUE),
        max_weight_rmse_approx_error = max(pop_df$w_rmse_approx_error, na.rm = TRUE)
      )
    )
    all_summary <- rbind(all_summary, one)

    safe_write_csv(
      sum_df,
      file.path(OUTDIR, paste0(scenario_id, "_finite_sample_by_k.csv"))
    )
    safe_write_csv(
      pop_df,
      file.path(OUTDIR, paste0(scenario_id, "_population_weight_distance_by_k.csv"))
    )
  }

  safe_write_csv(all_summary, file.path(OUTDIR, "joint_B_H_summary.csv"))
  safe_write_csv(all_pop, file.path(OUTDIR, "joint_B_H_weight_distance_by_k.csv"))

  plot_joint_heatmap(
    all_summary,
    "kstar_empirical_bagged",
    file.path(OUTDIR, "joint_B_H_kstar_empirical_bagged.png"),
    "Empirical bagged k* over joint B x H grid",
    "empirical bagged k*"
  )
  plot_joint_heatmap(
    all_summary,
    "bagged_gain_vs_equal",
    file.path(OUTDIR, "joint_B_H_gain_vs_equal.png"),
    "Bagged gain versus equal weights over joint B x H grid",
    "equal MSPE - min bagged MSPE"
  )
  plot_joint_heatmap(
    all_summary,
    "median_weight_corr_with_approx",
    file.path(OUTDIR, "joint_B_H_weight_approx_corr.png"),
    "Weight-approximation correlation over joint B x H grid",
    "median correlation"
  )
  plot_joint_heatmap(
    all_summary,
    "R_weight_at_empirical_bagged_k",
    file.path(OUTDIR, "joint_B_H_weight_distance_at_empirical_k.png"),
    "Exact weight-distance risk at empirical bagged k*",
    "exact R_weight at empirical k*"
  )
  plot_joint_heatmap(
    all_summary,
    "R_smallhet_rel_error_at_empirical_bagged_k",
    file.path(OUTDIR, "joint_B_H_smallhet_rel_error_at_empirical_k.png"),
    "Relative error of small-heterogeneity approximation at empirical k*",
    "relative approximation error"
  )

  all_summary
}

###############################################################################
# 7. Main run block
###############################################################################

run_all <- function() {
  cat("Clean common-loading RSM script started.\n")
  cat("Output directory: ", OUTDIR, "\n", sep = "")
  cat("Run profile: ", RUN_PROFILE, "\n", sep = "")
  cat("Experiments: ", paste(RUN_EXPERIMENTS, collapse = ", "), "\n", sep = "")

  results <- list()

  if ("baseline" %in% RUN_EXPERIMENTS) {
    results$baseline <- run_baseline()
    print(results$baseline)
  }

  if ("scaling_T_m" %in% RUN_EXPERIMENTS) {
    results$scaling_T_m <- run_scaling_T_m()
    print(results$scaling_T_m)
  }

  if ("joint_B_H_grid" %in% RUN_EXPERIMENTS) {
    results$joint_B_H_grid <- run_joint_B_H_grid()
    print(results$joint_B_H_grid)
  }

  cat("\nDone. CSV files and PNG figures saved in: ", OUTDIR, "\n", sep = "")
  invisible(results)
}

if (!exists("AUTO_RUN")) AUTO_RUN <- TRUE
if (isTRUE(AUTO_RUN)) {
  results <- run_all()
}
