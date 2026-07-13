suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})

panel <- read.csv(
  "/Users/boris/Desktop/RSShapley/Marco/marco_empirical_results_unemployment_all_eu/eurostat_unemployment_panel_all_eu.csv",
  check.names = FALSE
)
out_dir <- file.path(getwd(), "marco_empirical_results_unemployment_k_profile")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

k_grid <- 4:10
q_max <- 5
n_boot <- 40
n_calib_subsets <- 1200
n_random_subsets <- 5000
p_req <- function(q) ifelse(q == 1, 0.75, 0.55)
set.seed(20260708)

X <- as.matrix(panel[, setdiff(names(panel), c("date", "period")), drop = FALSE])
X <- scale(X, center = TRUE, scale = TRUE)
storage.mode(X) <- "double"

random_k_subsets <- function(N, k, M) replicate(M, sort(sample.int(N, k)), simplify = FALSE)

subset_pca <- function(X, subset, qbar) {
  Xs <- scale(X[, subset, drop = FALSE], center = TRUE, scale = FALSE)
  s <- svd(Xs, nu = qbar, nv = 0)
  list(evals = (s$d[seq_len(qbar)]^2) / nrow(Xs), Q = s$u[, seq_len(qbar), drop = FALSE])
}

nested_pass <- function(evals, eta) cumprod(as.integer(evals > eta))

compute_dispersion_fast <- function(Qlist, qbar) {
  M <- length(Qlist)
  Tn <- nrow(Qlist[[1]])
  out <- numeric(qbar)
  for (q in seq_len(qbar)) {
    Psum <- matrix(0, Tn, Tn)
    for (Q in Qlist) {
      Qq <- Q[, seq_len(q), drop = FALSE]
      Psum <- Psum + tcrossprod(Qq)
    }
    total_ordered_offdiag <- sum(Psum * Psum) - M * q
    avg_trace <- total_ordered_offdiag / (M * (M - 1))
    out[q] <- (q - avg_trace) / q
  }
  out
}

run_subsets <- function(X, subsets, qbar, need_dispersion = TRUE) {
  evals_mat <- matrix(NA_real_, nrow = length(subsets), ncol = qbar)
  Qlist <- if (need_dispersion) vector("list", length(subsets)) else NULL
  for (m in seq_along(subsets)) {
    pc <- subset_pca(X, subsets[[m]], qbar)
    evals_mat[m, ] <- pc$evals
    if (need_dispersion) Qlist[[m]] <- pc$Q
  }
  list(evals_mat = evals_mat, Qlist = Qlist)
}

circular_shift_columns <- function(X) {
  Tn <- nrow(X)
  Xn <- X
  for (j in seq_len(ncol(X))) {
    shift <- sample.int(Tn, 1L) - 1L
    Xn[, j] <- X[((seq_len(Tn) + shift - 1L) %% Tn) + 1L, j]
  }
  scale(Xn, center = TRUE, scale = TRUE)
}

rank_q_reference <- function(X, q) {
  s <- svd(scale(X, center = TRUE, scale = FALSE), nu = q, nv = q)
  common <- s$u[, seq_len(q), drop = FALSE] %*%
    diag(s$d[seq_len(q)], q, q) %*%
    t(s$v[, seq_len(q), drop = FALSE])
  common + circular_shift_columns(X - common)
}

calibrate_for_k <- function(X, k, qbar) {
  null_evals <- matrix(NA_real_, nrow = n_boot * n_calib_subsets, ncol = qbar)
  delta_ref <- matrix(NA_real_, nrow = n_boot, ncol = qbar)
  row <- 1L
  for (b in seq_len(n_boot)) {
    subs <- random_k_subsets(ncol(X), k, n_calib_subsets)
    res_null <- run_subsets(circular_shift_columns(X), subs, qbar, need_dispersion = FALSE)
    null_evals[row:(row + n_calib_subsets - 1L), ] <- res_null$evals_mat
    row <- row + n_calib_subsets
    res_ref <- run_subsets(rank_q_reference(X, qbar), subs, qbar, need_dispersion = TRUE)
    delta_ref[b, ] <- compute_dispersion_fast(res_ref$Qlist, qbar)
  }
  list(
    eta = apply(null_evals, 2, quantile, probs = 0.95, na.rm = TRUE),
    delta = apply(delta_ref, 2, quantile, probs = 0.95, na.rm = TRUE)
  )
}

profile_for_k <- function(k) {
  qbar <- min(q_max, k)
  cat(sprintf("\n=== unemployment k=%d, qbar=%d ===\n", k, qbar))
  thresholds <- calibrate_for_k(X, k, qbar)
  subsets <- random_k_subsets(ncol(X), k, n_random_subsets)
  res <- run_subsets(X, subsets, qbar, need_dispersion = TRUE)
  pass_surv <- t(apply(res$evals_mat, 1, nested_pass, eta = thresholds$eta))
  pi_q <- colMeans(pass_surv)
  Delta_q <- compute_dispersion_fast(res$Qlist, qbar)
  p_q <- p_req(seq_len(qbar))
  passes <- pi_q >= p_q & Delta_q <= thresholds$delta
  selected <- if (any(passes)) max(which(passes)) else 0L
  tibble(
    k = k,
    q = seq_len(qbar),
    n_subsets = n_random_subsets,
    pi_q = pi_q,
    p_required = p_q,
    survival_pass = pi_q >= p_q,
    Delta_q = Delta_q,
    delta_threshold = thresholds$delta,
    dispersion_pass = Delta_q <= thresholds$delta,
    eta_threshold = thresholds$eta,
    selected_rank = selected
  )
}

profile_list <- list()
for (kk in k_grid) {
  profile_list[[as.character(kk)]] <- profile_for_k(kk)
  write_csv(bind_rows(profile_list), file.path(out_dir, "unemployment_k_profile_partial.csv"))
}

profiles <- bind_rows(profile_list)
summary_tbl <- profiles |>
  group_by(k) |>
  summarize(
    selected_rank = max(selected_rank),
    n_subsets = max(n_subsets),
    survival_q1 = pi_q[q == 1],
    survival_q2 = pi_q[q == 2],
    dispersion_q1 = Delta_q[q == 1],
    delta_q1 = delta_threshold[q == 1],
    dispersion_q2 = Delta_q[q == 2],
    delta_q2 = delta_threshold[q == 2],
    .groups = "drop"
  )

write_csv(profiles, file.path(out_dir, "unemployment_k_profile.csv"))
write_csv(summary_tbl, file.path(out_dir, "unemployment_k_profile_summary.csv"))

p_surv <- ggplot(profiles, aes(k, pi_q, color = factor(q), group = factor(q))) +
  geom_line() + geom_point() +
  geom_hline(yintercept = 0.75, linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = 0.55, linetype = "dotted", color = "gray40") +
  scale_x_continuous(breaks = k_grid) +
  labs(x = "Subset size k", y = "Survival frequency", color = "q")
p_rank <- ggplot(summary_tbl, aes(k, selected_rank)) +
  geom_line() + geom_point() +
  scale_x_continuous(breaks = k_grid) +
  scale_y_continuous(breaks = 0:q_max, limits = c(0, q_max)) +
  labs(x = "Subset size k", y = expression(hat(r)[global](k)))
p_disp <- ggplot(profiles, aes(k, Delta_q, color = factor(q), group = factor(q))) +
  geom_line() + geom_point() +
  geom_line(aes(y = delta_threshold), linetype = "dashed") +
  scale_x_continuous(breaks = k_grid) +
  labs(x = "Subset size k", y = "Dispersion", color = "q")

ggsave(file.path(out_dir, "unemployment_survival_by_k.png"), p_surv, width = 7, height = 4.5, dpi = 160)
ggsave(file.path(out_dir, "unemployment_rank_by_k.png"), p_rank, width = 7, height = 4.5, dpi = 160)
ggsave(file.path(out_dir, "unemployment_dispersion_by_k.png"), p_disp, width = 7, height = 4.5, dpi = 160)

print(summary_tbl)
