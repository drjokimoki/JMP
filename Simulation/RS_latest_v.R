# =========================
# RS vs Mean, Time-Series DGP (AR(1))
# Three specific phi triplets; no duplicates; OLS when z == N
# =========================

# --- Clear & libs ---
rm(list = ls())
library(MASS)
library(foreach)
library(doParallel)
library(doRNG)
library(ggplot2)

# --- Global seed (reproducibility for %dorng%) ---
set.seed(778)
main_seed <- 778

# =========================
# Parameters
# =========================
n_rep      <- 500              # Monte Carlo replications
T_len_list <- c(100, 200, 1000)      # sample sizes
N          <- 50               # number of base predictors
sigma_e    <- 1                # target disturbance sd (uncond.)
K_RS       <- 1000             # RS Monte Carlo per forecast origin
rho_list   <- c(0.8)        # cross-sec Toeplitz rho for contemporaneous idio cov
R2_list    <- c(0.3, 0.5, 0.1) # target R^2 for signal strength
beta_list_raw <- list(
  first_ten = c(rep(1, 10), rep(0, N - 10))
)

# fraction reserved for out-of-sample (expanding window)
oos_prop <- 0.2

# =========================
# EXACT phi triplets (no expand.grid — avoids duplicates)
# =========================
phi_grid <- data.frame(
  phi_f   = c(0.7, 0.3, 0.3),
  phi_eta = c(0.3, 0.7, 0.3),
  phi_u   = c(0.3, 0.3, 0.7)
)

# Safety: ensure no duplicates
stopifnot(nrow(phi_grid) == nrow(unique(phi_grid)))

# =========================
# Helpers
# =========================

# Toeplitz covariance for contemporaneous idiosyncratic correlation across predictors
toeplitz_cov <- function(rho, p) toeplitz(rho^(0:(p - 1)))

# Simulate T x N matrix following vector AR(1): E_t = phi * E_{t-1} + z_t,  z_t ~ N(0, Sigma)
gen_E_AR1 <- function(T_len, N, phi, Sigma) {
  Z <- MASS::mvrnorm(n = T_len, mu = rep(0, N), Sigma = Sigma)
  E <- matrix(0, T_len, N)
  for (t in 2:T_len) E[t, ] <- phi * E[t - 1, ] + Z[t, ]
  E
}

# Scalar AR(1) with chosen unconditional variance s2 (so innovations keep Var(x_t) = s2)
gen_AR1_scalar <- function(T_len, phi, s2 = 1) {
  sd_eps <- sqrt((1 - phi^2) * s2)
  e <- rnorm(T_len, 0, sd_eps)
  x <- numeric(T_len)
  for (t in 2:T_len) x[t] <- phi * x[t - 1] + e[t]
  x
}

# =========================
# Build unique scenario grid for BASELINE
# =========================
base_grid <- expand.grid(
  T_len     = T_len_list,
  rho       = rho_list,
  R2_target = R2_list,
  beta_name = names(beta_list_raw),
  stringsAsFactors = FALSE
)

scenario_grid <- merge(base_grid, phi_grid, all = TRUE)
scenario_grid <- unique(scenario_grid)  # just in case
scenario_grid$rmsfe_base <- NA_real_

# Assert uniqueness by scenario ID
id_cols <- c("T_len","rho","R2_target","beta_name","phi_f","phi_eta","phi_u")
stopifnot(nrow(scenario_grid) == nrow(unique(scenario_grid[id_cols])))

# =========================
# Parallel backend
# =========================
num_cores <- max(parallel::detectCores() - 1, 1)
cl <- makeCluster(num_cores)
on.exit({ try(stopCluster(cl), silent = TRUE) }, add = TRUE)
registerDoParallel(cl)
cat(sprintf("Computing baseline RMSFE on %d cores...\n", num_cores))

# =========================
# Baseline: simple mean combiner
# =========================
for (j in seq_len(nrow(scenario_grid))) {
  
  # --- Extract scenario ---
  T_len     <- scenario_grid$T_len[j]
  rho       <- scenario_grid$rho[j]
  R2_target <- scenario_grid$R2_target[j]
  beta_raw  <- beta_list_raw[[scenario_grid$beta_name[j]]]
  
  phi_f   <- scenario_grid$phi_f[j]
  phi_eta <- scenario_grid$phi_eta[j]
  phi_u   <- scenario_grid$phi_u[j]
  
  # OOS configuration (expanding window)
  oos_length  <- max(1, floor(T_len * oos_prop))
  init_window <- T_len - oos_length
  
  # Contemporaneous structures
  Qvec   <- rep(1, N)
  Sigma  <- toeplitz_cov(rho, N)
  V_mat  <- Qvec %*% t(Qvec) + Sigma
  
  # Scale beta to hit target in-sample R^2
  V_signal  <- as.numeric(t(beta_raw) %*% V_mat %*% beta_raw)
  C         <- sqrt((R2_target / (1 - R2_target)) / V_signal)
  beta_true <- C * beta_raw
  
  # --- Monte Carlo for baseline ---
  base_results <- foreach(repl = 1:n_rep, .combine = c, .packages = 'MASS', .options.RNG = main_seed) %dorng% {
    
    # DGP
    F_t  <- gen_AR1_scalar(T_len, phi = phi_f,  s2 = 1)
    E_t  <- gen_E_AR1(T_len, N, phi = phi_eta, Sigma = Sigma)
    X    <- outer(F_t, Qvec) + E_t
    u_t  <- gen_AR1_scalar(T_len, phi = phi_u, s2 = sigma_e^2)
    
    # y_t = x_{t-1} %*% beta + u_t
    Y <- c(NA, as.vector(X[1:(T_len - 1), ] %*% beta_true + u_t[2:T_len]))
    
    # Mean combiner via bivariate slopes
    f_base <- numeric(oos_length - 1)
    for (k in seq_len(oos_length - 1)) {
      t_now <- init_window + k
      X_tr  <- X[1:(t_now - 1), , drop = FALSE]
      y_tr  <- Y[2:t_now]
      x_now <- X[t_now, , drop = FALSE]
      
      csq       <- colSums(X_tr^2)
      gamma_hat <- ifelse(csq == 0, 0, colSums(X_tr * y_tr) / csq)
      
      f_base[k] <- mean(gamma_hat * as.vector(x_now), na.rm = TRUE)
    }
    
    sqrt(mean((Y[(init_window + 2):T_len] - f_base)^2, na.rm = TRUE))
  }
  
  scenario_grid$rmsfe_base[j] <- mean(base_results)
  cat(sprintf("Baseline %d/%d (T=%d, rho=%.2f, R2=%.2f, beta=%s, phis=(%.1f,%.1f,%.1f)): RMSFE_base=%.4f\n",
              j, nrow(scenario_grid), T_len, rho, R2_target, scenario_grid$beta_name[j],
              phi_f, phi_eta, phi_u, scenario_grid$rmsfe_base[j]))
}

# Keep only what we need and guarantee uniqueness
scenario_grid_unique <- unique(scenario_grid[, c(id_cols, "rmsfe_base")])
stopifnot(nrow(scenario_grid_unique) == nrow(scenario_grid))  # should match

# =========================
# RS param grid (clean, no duplicates)
# =========================
param_grid <- merge(base_grid, phi_grid, all = TRUE)
param_grid <- merge(param_grid, scenario_grid_unique,
                    by = id_cols, all.x = TRUE)
param_grid <- unique(param_grid)
stopifnot(nrow(param_grid) == nrow(unique(param_grid[id_cols])))

cat("\nComputing RS forecasts for z in {2, z_suggest, N}...\n")

# =========================
# RS over z in {2, round(sqrt(T)), N}; OLS when z == N
# =========================
results_rows <- list()
row_out <- 0

for (ii in seq_len(nrow(param_grid))) {
  
  T_len     <- param_grid$T_len[ii]
  rho       <- param_grid$rho[ii]
  R2_target <- param_grid$R2_target[ii]
  beta_raw  <- beta_list_raw[[param_grid$beta_name[ii]]]
  base_val  <- param_grid$rmsfe_base[ii]
  
  phi_f   <- param_grid$phi_f[ii]
  phi_eta <- param_grid$phi_eta[ii]
  phi_u   <- param_grid$phi_u[ii]
  
  oos_length  <- max(1, floor(T_len * oos_prop))
  init_window <- T_len - oos_length
  
  Qvec   <- rep(1, N)
  Sigma  <- toeplitz_cov(rho, N)
  V_mat  <- Qvec %*% t(Qvec) + Sigma
  
  V_signal  <- as.numeric(t(beta_raw) %*% V_mat %*% beta_raw)
  C         <- sqrt((R2_target / (1 - R2_target)) / V_signal)
  beta_true <- C * beta_raw
  
  # z_try = {2, sqrt(T), N}
  z_suggest <- max(1, round(sqrt(T_len)))
  offsets <- seq(-1, 6, by = 1)
  z_try <- sort(unique(pmin(N, pmax(1, c(2, z_suggest+offsets, N)))))
  
  for (current_z in z_try) {
    row_out <- row_out + 1
    
    rs_results <- foreach(repl = 1:n_rep, .combine = c, .packages = 'MASS', .options.RNG = main_seed) %dorng% {
      
      # DGP
      F_t  <- gen_AR1_scalar(T_len, phi = phi_f,  s2 = 1)
      E_t  <- gen_E_AR1(T_len, N, phi = phi_eta, Sigma = Sigma)
      X    <- outer(F_t, Qvec) + E_t
      u_t  <- gen_AR1_scalar(T_len, phi = phi_u, s2 = sigma_e^2)
      Y    <- c(NA, as.vector(X[1:(T_len - 1), ] %*% beta_true + u_t[2:T_len]))
      
      f_rs <- numeric(oos_length - 1)
      
      for (k in seq_len(oos_length - 1)) {
        t_now <- init_window + k
        X_tr  <- X[1:(t_now - 1), , drop = FALSE]
        y_tr  <- Y[2:t_now]
        x_now <- X[t_now, , drop = FALSE]
        
        csq       <- colSums(X_tr^2)
        gamma_hat <- ifelse(csq == 0, 0, colSums(X_tr * y_tr) / csq)
        
        if (current_z == N) {
          # ---- FULL OLS once (no RS replication) ----
          M_all    <- cbind(1, sweep(X_tr, 2, gamma_hat, `*`))
          coef_all <- lm.fit(M_all, y_tr)$coefficients
          f_rs[k]  <- sum(coef_all * c(1, gamma_hat * as.vector(x_now)), na.rm = TRUE)
        } else {
          # ---- Random Subspace with k = current_z ----
          preds <- replicate(K_RS, {
            sub    <- sample.int(N, current_z)
            M_sub  <- cbind(1, sweep(X_tr[, sub, drop = FALSE], 2, gamma_hat[sub], `*`))
            coef_s <- lm.fit(M_sub, y_tr)$coefficients
            sum(coef_s * c(1, gamma_hat[sub] * as.vector(x_now)[sub]), na.rm = TRUE)
          })
          f_rs[k] <- mean(preds)
        }
      }
      
      sqrt(mean((Y[(init_window + 2):T_len] - f_rs)^2, na.rm = TRUE))
    }
    
    avg_rs <- mean(rs_results)
    
    results_rows[[row_out]] <- data.frame(
      T_len      = T_len,
      rho        = rho,
      R2_target  = R2_target,
      beta_name  = param_grid$beta_name[ii],
      phi_f      = phi_f,
      phi_eta    = phi_eta,
      phi_u      = phi_u,
      rmsfe_base = base_val,
      z          = current_z,
      rmsfe_rs   = avg_rs,
      rel_rmsfe  = avg_rs / base_val,
      stringsAsFactors = FALSE
    )
    
    cat(sprintf("T=%d rho=%.2f R2=%.2f beta=%s z=%2d | phis=(%.1f,%.1f,%.1f) | RS=%.4f base=%.4f rel=%.4f%s\n",
                T_len, rho, R2_target, param_grid$beta_name[ii], current_z,
                phi_f, phi_eta, phi_u, avg_rs, base_val, avg_rs / base_val,
                if (current_z == N) " [full OLS]" else ""))
  }
}

# Collect RS results
rs_results_df <- do.call(rbind, results_rows)

# Quick peek
cat("\n--- Head of RS results ---\n")
print(head(rs_results_df))

# Optional: save
# Optional: save
write.csv(rs_results_df, "/Users/boris/Desktop/RSShapley/Combination_Monte/Final/rs_results.csv", row.names = FALSE)