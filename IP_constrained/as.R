# ================================================================
# Backward pseudo real-time forecasting across multiple vintages
#   + RS/RSM forecast combination with SUM-TO-ONE weights (can be negative):
#       (i)  no intercept
#       (ii) weights can be negative
#       (iii) sum(w)=1
#   + OUTPUTS per vintage:
#       - w_bar (unconditional mean sum-to-one weight across successful RS draws)
#       - "contributions": w_bar * (Prediction - Hist_mean)
#       - sum checks printed per vintage and stored in summary
# ================================================================

rm(list = ls())

# ---------------- Packages ----------------
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(lubridate)
  library(purrr); library(tibble); library(stringr)
  library(forecast)
})

# ---------------- Settings ----------------
monthly_dir   <- "/Users/boris/Downloads/Historical FRED-MD Vintages Final"
quarterly_dir <- "/Users/boris/Desktop/RSShapley/Combination_Monte/FRED_QD"

uni_window_type <- "expanding"
uni_window_size <- 40
uni_min_start   <- 40

rs_draws <- 1000
rs_seed  <- 12893

suppressWarnings(RNGkind(kind = "Mersenne-Twister",
                         normal.kind = "Inversion",
                         sample.kind = "Rounding"))
set.seed(rs_seed)

# ---------------- Helpers -----------------
lag1  <- function(x) c(NA, head(x, -1))
d1    <- function(x) x - lag1(x)
d2    <- function(x) d1(d1(x))
lg    <- function(x) { y <- x; y[x <= 0] <- NA_real_; log(y) }
d1lg  <- function(x) d1(lg(x))
d2lg  <- function(x) d2(lg(x))
d_ret <- function(x) {
  prev <- lag1(x)
  r <- ifelse(prev <= 0 | is.na(prev), NA_real_, x/prev - 1)
  d1(r)
}
parse_date <- function(x) if (inherits(x,"Date")) x else coalesce(ymd(x, quiet=TRUE), mdy(x, quiet=TRUE))

transform_by_first_row_df <- function(df) {
  codes <- suppressWarnings(as.numeric(df[1, ]))
  out <- df[-1, , drop = FALSE]
  for (i in seq_along(out)) {
    if (!is.numeric(out[[i]])) next
    out[[i]] <- switch(as.character(codes[i]),
                       "1" = out[[i]],
                       "2" = d1(out[[i]]),
                       "3" = d2(out[[i]]),
                       "4" = lg(out[[i]]),
                       "5" = d1lg(out[[i]]),
                       "6" = d2lg(out[[i]]),
                       "7" = d_ret(out[[i]]),
                       out[[i]])
  }
  rownames(out) <- NULL
  out
}

safe_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)

# ================================================================
# Sum-to-one weights (CAN BE NEGATIVE): closed-form ridge OLS
#   min_w ||y - X w||^2 + ridge * ||w||^2   s.t. 1'w = 1
# ================================================================
ols_sum1 <- function(X, y, ridge = 1e-8) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  k <- ncol(X)
  if (k == 1L) return(1)
  
  A  <- crossprod(X) + diag(ridge, k)
  b  <- crossprod(X, y)
  
  w0 <- tryCatch(solve(A, b), error = function(e) rep(NA_real_, k))
  if (anyNA(w0) || any(!is.finite(w0))) return(rep(NA_real_, k))
  
  cvec <- tryCatch(solve(A, rep(1, k)), error = function(e) rep(NA_real_, k))
  if (anyNA(cvec) || any(!is.finite(cvec))) return(rep(NA_real_, k))
  
  denom <- sum(cvec)
  if (!is.finite(denom) || abs(denom) <= 1e-12) return(rep(NA_real_, k))
  
  adj <- (sum(w0) - 1) / denom
  as.numeric(w0 - cvec * adj)
}

# ---------------- Univariate forecasts (split across time) -----------------
univariate_forecasts_split <- function(y, X, y_levels, n_start, expanding = TRUE) {
  stopifnot(length(y) == length(y_levels))
  T  <- length(y)
  TX <- nrow(X)
  p  <- ncol(X)
  if (n_start < 2) stop("n_start must be >= 2")
  
  last_t <- if (TX == T) T - 1 else T
  
  z_hat      <- matrix(NA_real_, nrow = TX, ncol = p, dimnames = list(NULL, colnames(X)))
  lvl_med    <- matrix(NA_real_, nrow = TX, ncol = p, dimnames = list(NULL, colnames(X)))
  growth_med <- matrix(NA_real_, nrow = TX, ncol = p, dimnames = list(NULL, colnames(X)))
  
  for (t in seq(n_start, last_t)) {
    idx <- if (expanding) 1:t else (t - n_start + 1):t
    base_level <- y_levels[t]
    for (j in seq_len(p)) {
      x_in <- X[idx, j]
      if (all(is.na(x_in)) || length(unique(stats::na.omit(x_in))) < 2) next
      fit <- lm(y[idx] ~ x, data = data.frame(y = y[idx], x = x_in))
      x_next <- X[t + 1, j]
      if (is.na(x_next)) next
      z <- as.numeric(predict(fit, newdata = data.frame(x = x_next)))
      z_hat[t + 1, j] <- z
      y_next_med <- base_level * exp(z)
      lvl_med[t + 1, j] <- y_next_med
      growth_med[t + 1, j] <- ((y_next_med / base_level) - 1) * 100
    }
  }
  
  list(
    z_hat      = as.data.frame(z_hat),
    lvl_med    = as.data.frame(lvl_med),
    growth_med = as.data.frame(growth_med)
  )
}

# ---------------- Random subsampling RS with sum-to-one weights -----------------
rand_subsample_from_df <- function(df, y_col, x_next, k, B = 1000,
                                   scale_predictors = FALSE,
                                   weight_scheme = c("uniform", "abs_t"),
                                   alpha = 1,
                                   t_no_intercept = FALSE,
                                   ridge = 1e-8) {
  weight_scheme <- match.arg(weight_scheme)
  
  if (!y_col %in% names(df)) stop("y_col not found in df.")
  y <- as.numeric(df[[y_col]])
  Xnames_all <- setdiff(names(df), y_col)
  
  # accept named vector or 1-row data frame for x_next
  if (is.vector(x_next)) {
    if (is.null(names(x_next))) stop("If x_next is a vector, it must be named.")
    x_next <- as.data.frame(as.list(x_next))
  } else {
    x_next <- as.data.frame(x_next)
  }
  if (nrow(x_next) != 1) stop("x_next must be a single row.")
  
  # candidates with a non-NA next value
  cand <- intersect(Xnames_all, names(x_next))
  cand <- cand[!is.na(as.numeric(x_next[1, cand]))]
  if (length(cand) < k) stop("Not enough overlapping predictors between df and x_next to sample k of them.")
  
  # ---- (A) selection probabilities ----
  prob_vec <- rep(1/length(cand), length(cand))
  if (weight_scheme == "abs_t") {
    t_stat <- sapply(cand, function(v) {
      subdat <- df[, c(y_col, v), drop = FALSE]
      cc <- stats::complete.cases(subdat)
      if (sum(cc) < 3) return(NA_real_)
      yy <- as.numeric(subdat[[y_col]][cc])
      xx <- as.numeric(subdat[[v]][cc])
      if (all(xx == xx[1L])) return(NA_real_)
      if (t_no_intercept) {
        n   <- length(yy)
        xy  <- sum(xx * yy); xx2 <- sum(xx^2)
        beta <- xy / xx2
        SSE  <- sum(yy^2) - 2 * beta * xy + beta^2 * xx2
        SSE  <- max(SSE, 1e-12)
        beta * sqrt((n - 1) * xx2 / SSE)
      } else {
        X <- cbind(1, xx)
        fit <- lm.fit(X, yy)
        if (anyNA(fit$coefficients)) return(NA_real_)
        n <- length(yy); p <- 2
        res <- yy - drop(X %*% fit$coefficients)
        s2  <- sum(res^2) / max(n - p, 1)
        XtX_inv <- tryCatch(solve(crossprod(X)), error = function(e) NA)
        if (anyNA(XtX_inv)) return(NA_real_)
        se_beta <- sqrt(s2 * XtX_inv[2, 2])
        if (!is.finite(se_beta) || se_beta <= 0) return(NA_real_)
        fit$coefficients[2] / se_beta
      }
    })
    w <- exp(alpha * abs(t_stat))
    w[!is.finite(w)] <- 0
    sw <- sum(w)
    prob_vec <- if (sw > 0) as.numeric(w / sw) else rep(1/length(cand), length(cand))
  }
  names(prob_vec) <- cand
  
  # ---- (B) RS draws ----
  forecasts     <- rep(NA_real_, B)
  selected_vars <- matrix(NA_character_, nrow = B, ncol = k)
  coef_entries  <- vector("list", B)
  n_failed <- 0L
  
  for (b in seq_len(B)) {
    vars <- sample(cand, k, replace = FALSE, prob = prob_vec)
    selected_vars[b, ] <- vars
    
    subdat <- df[, c(y_col, vars), drop = FALSE]
    cc <- stats::complete.cases(subdat)
    if (!any(cc)) { n_failed <- n_failed + 1L; next }
    
    y_train <- as.numeric(subdat[[y_col]][cc])
    X_train <- as.matrix(subdat[cc, vars, drop = FALSE])
    
    x_row <- as.numeric(as.matrix(x_next[1, vars, drop = FALSE]))
    if (any(is.na(x_row))) { n_failed <- n_failed + 1L; next }
    
    if (scale_predictors) {
      mu  <- colMeans(X_train)
      sdv <- apply(X_train, 2, sd)
      sdv[sdv == 0 | is.na(sdv)] <- 1
      X_train <- scale(X_train, center = mu, scale = sdv)
      x_row   <- (x_row - mu) / sdv
    }
    
    w_hat <- ols_sum1(X_train, y_train, ridge = ridge)
    if (anyNA(w_hat)) { n_failed <- n_failed + 1L; next }
    
    forecasts[b] <- drop(x_row %*% w_hat)
    coef_entries[[b]] <- tibble::tibble(draw = b, variable = vars, coef = as.numeric(w_hat))
  }
  
  coef_long <- dplyr::bind_rows(coef_entries)
  if (is.null(coef_long) || nrow(coef_long) == 0) {
    coef_long <- tibble::tibble(draw = integer(), variable = character(), coef = numeric())
  }
  
  # Successful draws
  B_succ <- dplyr::n_distinct(coef_long$draw)
  
  # Per-draw sum(w) check among successful draws
  draw_sum_check <- coef_long %>%
    dplyr::group_by(draw) %>%
    dplyr::summarise(sumw = sum(coef, na.rm = TRUE), .groups = "drop")
  
  max_abs_sumw_dev <- if (nrow(draw_sum_check) > 0) max(abs(draw_sum_check$sumw - 1)) else NA_real_
  
  # selection frequency (sampled into subset, regardless of success)
  sel_tab <- table(factor(as.vector(selected_vars), levels = cand))
  selection_counts <- as.integer(sel_tab)
  names(selection_counts) <- cand
  
  # mean weight conditional on nonzero + mean over included (selected) draws
  coef_stats <- coef_long |>
    dplyr::group_by(variable) |>
    dplyr::summarise(
      n_with_coef = dplyr::n(),
      nz = sum(coef != 0, na.rm = TRUE),
      mean_coef_nonzero = ifelse(nz > 0, mean(coef[coef != 0], na.rm = TRUE), NA_real_),
      mean_coef_all = mean(coef, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::right_join(tibble::tibble(variable = cand), by = "variable") |>
    dplyr::mutate(
      n_with_coef = dplyr::coalesce(n_with_coef, 0L),
      nz          = dplyr::coalesce(nz, 0L)
    ) |>
    dplyr::arrange(variable)
  
  mean_coef_nonzero <- stats::setNames(coef_stats$mean_coef_nonzero, coef_stats$variable)
  selection_counts_success <- stats::setNames(coef_stats$n_with_coef, coef_stats$variable)
  
  # unconditional mean weight across successful draws (zeros for unselected are implicit)
  if (B_succ > 0) {
    w_bar <- coef_long %>%
      dplyr::group_by(variable) %>%
      dplyr::summarise(w_bar = sum(coef, na.rm = TRUE) / B_succ, .groups = "drop")
  } else {
    w_bar <- tibble::tibble(variable = cand, w_bar = 0)
  }
  w_bar_vec <- stats::setNames(w_bar$w_bar, w_bar$variable)
  
  out <- list(
    forecasts = forecasts,
    selected_vars = selected_vars,
    selected_vars_matrix = selected_vars,
    coef_long = coef_long,
    failed = n_failed,
    probs = prob_vec,
    selection_counts = selection_counts,
    selection_counts_success = selection_counts_success,
    coef_stats = coef_stats,
    mean_coef_nonzero = mean_coef_nonzero,
    B_succ = B_succ,
    max_abs_sumw_dev = max_abs_sumw_dev,
    w_bar = w_bar_vec
  )
  class(out) <- c("rand_subsample_result", class(out))
  out
}

# ---------------- One-vintage pipeline -----------------
run_forecast <- function(monthly_path, quarterly_path, n_start = 40, ridge = 1e-8, verbose_checks = TRUE) {
  
  Monthly_FRED   <- read.csv(monthly_path, stringsAsFactors = FALSE, check.names = FALSE)
  Quarterly_FRED <- read.csv(quarterly_path, stringsAsFactors = FALSE, check.names = FALSE)
  
  Monthly_FRED_clean   <- Monthly_FRED
  Quarterly_FRED_clean <- Quarterly_FRED[-1,]
  Quarterly_raw_for_levels <- Quarterly_FRED_clean
  
  # Monthly (apply codes)
  Monthly_without_time <- Monthly_FRED_clean[, -1, drop = FALSE] %>%
    mutate(across(everything(), ~ suppressWarnings(as.numeric(as.character(.)))))
  Final_monthly <- transform_by_first_row_df(Monthly_without_time)
  Final_monthly$Date <- parse_date(Monthly_FRED_clean[-1, 1, drop = TRUE])
  
  # Quarterly (apply codes)
  Quarterly_without_time <- Quarterly_FRED_clean[, -1, drop = FALSE] %>%
    mutate(across(everything(), ~ suppressWarnings(as.numeric(as.character(.)))))
  Final_quarterly <- transform_by_first_row_df(Quarterly_without_time)
  Final_quarterly$Date <- parse_date(Quarterly_FRED_clean[-1, 1, drop = TRUE]) %>% as_date()
  
  # GDP levels (for back-transform)
  y_levels <- tibble(
    Date = parse_date(Quarterly_raw_for_levels[-1, 1, drop = TRUE]) %>% as_date(),
    GDPC1_level = suppressWarnings(as.numeric(as.character(
      Quarterly_raw_for_levels[-1, which(names(Quarterly_raw_for_levels)=="GDPC1"), drop = TRUE]
    )))
  ) %>% arrange(Date) %>%
    drop_na()
  
  # Aggregate monthly → quarterly
  m_q <- Final_monthly %>%
    mutate(
      Date = parse_date(Date),
      q_date = make_date(year(Date), c(3,6,9,12)[quarter(Date)], 1L)
    ) %>%
    group_by(q_date) %>%
    summarise(across(where(is.numeric), safe_mean), .groups = "drop") %>%
    rename(Date = q_date) %>%
    arrange(Date)
  
  Final_quarterly <- Final_quarterly[c("GDPC1","Date")] %>% arrange(Date)
  Final_quarterly <- Final_quarterly[-1,]
  m_q <- m_q[-1,]
  
  estimation0 <- Final_quarterly %>% inner_join(m_q, by = "Date") %>% arrange(Date)
  last_x <- m_q %>% filter(Date == max(Date)) %>% select(-Date)
  
  y_last_level <- y_levels %>%
    semi_join(estimation0, by = "Date") %>%
    arrange(Date) %>%
    pull(GDPC1_level) %>%
    tail(1)
  
  if (length(y_last_level) == 0 || is.na(y_last_level)) {
    stop("Could not determine the last GDP level aligned with the estimation sample.")
  }
  
  # ========= MEAN comb step (your original) =========
  y <- estimation0$GDPC1
  X <- estimation0 %>% select(-Date, -GDPC1)
  
  keep_vars <- names(X)[!is.na(as.numeric(last_x[1, names(X)]))]
  X <- X[, keep_vars, drop = FALSE]
  last_x <- last_x[, keep_vars, drop = FALSE]
  
  results <- map_dfr(names(X), function(varname) {
    x <- estimation0[[varname]]
    if (all(is.na(x)) || length(unique(na.omit(x))) < 2) {
      return(tibble(variable = varname, alpha = NA_real_, beta = NA_real_,
                    r2 = NA_real_, z_hat_next = NA_real_,
                    y_forecast_next = NA_real_, y_forecast_next_mean_bc = NA_real_))
    }
    fit <- lm(y ~ x)
    x_future <- as.numeric(last_x[[varname]])
    if (is.na(x_future)) {
      return(tibble(variable = varname,
                    alpha = coef(fit)[1], beta = coef(fit)[2],
                    r2 = summary(fit)$r.squared, z_hat_next = NA_real_,
                    y_forecast_next = NA_real_, y_forecast_next_mean_bc = NA_real_))
    }
    
    z_hat_next <- as.numeric(predict(fit, newdata = data.frame(x = x_future)))
    sigma2 <- summary(fit)$sigma^2
    
    y_next_med  <- as.numeric(y_last_level * exp(z_hat_next))
    y_next_mean <- as.numeric(y_last_level * exp(z_hat_next + 0.5 * sigma2))
    
    tibble(
      variable = varname,
      alpha = unname(coef(fit)[1]),
      beta  = unname(coef(fit)[2]),
      r2 = summary(fit)$r.squared,
      z_hat_next = z_hat_next,
      y_forecast_next = y_next_med,
      y_forecast_next_mean_bc = y_next_mean
    )
  }) %>% arrange(desc(r2))
  
  results <- results %>%
    mutate(
      y_forecast_growth = ((y_forecast_next / y_last_level) - 1)*100,
      y_forecast_growth_mean_bc = ((y_forecast_next_mean_bc / y_last_level) - 1)*100
    )
  
  MEAN_mean <- mean(results$y_forecast_growth, na.rm = TRUE)
  
  # ========= Build RS regression dataset =========
  X_plus <- rbind(X, last_x)
  fcasts <- univariate_forecasts_split(
    y, X_plus, y_levels = y_levels$GDPC1_level[-1], n_start = n_start
  )
  
  growth_table <- fcasts$growth_med
  growth_table$Date <- (Final_monthly %>%
                          mutate(q_date = make_date(year(Date), c(3,6,9,12)[quarter(Date)], 1L)) %>%
                          distinct(q_date) %>% arrange(q_date) %>% pull(q_date))[-1]
  
  growth_y <- (y_levels$GDPC1_level / dplyr::lag(y_levels$GDPC1_level) - 1)*100
  growth_y <- growth_y[-1]
  growth_y <- as.data.frame(growth_y)
  growth_y$Date <- y_levels$Date[-1]
  
  estimation <- suppressMessages(
    growth_y %>%
      inner_join(growth_table, by = "Date") %>%
      arrange(Date) %>%
      drop_na()
  )
  
  last_x_rs <- growth_table %>% filter(Date == max(Date)) %>% select(-Date)
  
  k <- max(1L, round(sqrt(nrow(estimation)))+5)
  
  # ========= RS/RSM sum-to-one combination (weights can be negative) =========
  RSSS <- rand_subsample_from_df(
    estimation, "growth_y", last_x_rs, k,
    B = rs_draws,
    scale_predictors = FALSE,        # set TRUE if you want extra stability
    weight_scheme = "uniform",
    alpha = 1,
    t_no_intercept = FALSE,
    ridge = ridge
  )
  
  RS_mean <- mean(RSSS$forecasts, na.rm = TRUE)
  RS_fore <- RSSS$forecasts
  
  # ========= Build output table with w_bar and contributions =========
  pred_vars <- names(last_x_rs)
  base_tbl <- tibble(
    Var = pred_vars,
    Selected = as.numeric(RSSS$selection_counts[pred_vars]) / rs_draws,
    mean_weight_nonzero = as.numeric(RSSS$mean_coef_nonzero[pred_vars]),
    w_bar = as.numeric(RSSS$w_bar[pred_vars])
  )
  
  pred_tbl <- tibble(
    Var = pred_vars,
    Prediction = as.numeric(last_x_rs[1, pred_vars])
  )
  
  historic_mean <- estimation %>% summarise_if(is.numeric, mean, na.rm = TRUE)
  hist_tbl <- as.data.frame(t(historic_mean)) %>%
    rownames_to_column("Var") %>%
    setNames(c("Var", "Hist_mean"))
  
  Fin_table2 <- base_tbl %>%
    left_join(pred_tbl, by = "Var") %>%
    left_join(hist_tbl, by = "Var") %>%
    mutate(
      w_bar = ifelse(is.na(w_bar), 0, w_bar),
      Contribution = w_bar * (Prediction - Hist_mean)
    )
  
  # GDP historical mean (of realized growth_y)
  mean_hist_gdp <- hist_tbl %>% filter(Var == "growth_y") %>% pull(Hist_mean)
  if (length(mean_hist_gdp) == 0) mean_hist_gdp <- NA_real_
  
  # ========= Sum checks per vintage =========
  sum_w_bar <- sum(Fin_table2$w_bar, na.rm = TRUE)
  
  # check: RS forecast implied by w_bar on last_x_rs (should match RS_mean closely)
  RS_from_w_bar <- sum(Fin_table2$w_bar * Fin_table2$Prediction, na.rm = TRUE)
  
  # baseline: weighted historic mean of each predictor forecast
  base_weighted <- sum(Fin_table2$w_bar * Fin_table2$Hist_mean, na.rm = TRUE)
  
  sum_contrib <- sum(Fin_table2$Contribution, na.rm = TRUE)
  
  # should be near 0:
  contrib_residual <- (RS_from_w_bar - base_weighted) - sum_contrib
  
  # compare to RS - mean(growth_y)
  phi0 <- base_weighted - mean_hist_gdp
  total_vs_hist_gdp <- sum_contrib + phi0
  rs_vs_hist_gdp <- RS_from_w_bar - mean_hist_gdp
  
  # Vintage tags
  vintage_monthly   <- basename(monthly_path)
  vintage_quarterly <- basename(quarterly_path)
  vintage_key <- tryCatch({
    vm <- gsub("\\.csv$", "", vintage_monthly)
    vq <- gsub("\\.csv$", "", vintage_quarterly)
    if (grepl("^\\d{4}-\\d{2}$", vm)) {
      vm
    } else if (grepl("_(\\d{4})m(\\d{2})$", vq)) {
      sub(".*_(\\d{4})m(\\d{2})$", "\\1-\\2", vq)
    } else {
      vm
    }
  }, error = function(e) NA_character_)
  
  if (isTRUE(verbose_checks)) {
    cat("\n",
        "================ Vintage ", vintage_key, " ================\n", sep = "")
    cat("B_succ (successful RS draws): ", RSSS$B_succ, " / ", rs_draws, "\n", sep = "")
    cat("max |sum_w(draw)-1| (successful draws): ", format(RSSS$max_abs_sumw_dev, digits = 6), "\n", sep = "")
    cat("sum(w_bar): ", format(sum_w_bar, digits = 8), "\n", sep = "")
    cat("RS_mean (avg over draws): ", format(RS_mean, digits = 8), "\n", sep = "")
    cat("RS_from_w_bar (w_bar' * Prediction): ", format(RS_from_w_bar, digits = 8),
        "  |diff|=", format(abs(RS_from_w_bar - RS_mean), digits = 6), "\n", sep = "")
    cat("sum(Contribution): ", format(sum_contrib, digits = 8), "\n", sep = "")
    cat("Check: (RS_from_w_bar - base_weighted) - sum(Contribution) = ",
        format(contrib_residual, digits = 8), "\n", sep = "")
    cat("Compare to RS - mean(growth_y):\n")
    cat("  RS_from_w_bar - mean_hist_gdp = ", format(rs_vs_hist_gdp, digits = 8), "\n", sep = "")
    cat("  sum(Contribution) + (base_weighted - mean_hist_gdp) = ", format(total_vs_hist_gdp, digits = 8),
        "  |diff|=", format(abs(total_vs_hist_gdp - rs_vs_hist_gdp), digits = 6), "\n", sep = "")
  }
  
  list(
    summary = tibble(
      vintage = vintage_key,
      vintage_monthly = vintage_monthly,
      vintage_quarterly = vintage_quarterly,
      MEAN_mean = MEAN_mean,
      RS_mean = RS_mean,
      B_succ = RSSS$B_succ,
      max_abs_sumw_dev = RSSS$max_abs_sumw_dev,
      sum_w_bar = sum_w_bar,
      RS_from_w_bar = RS_from_w_bar,
      abs_diff_RS = abs(RS_from_w_bar - RS_mean),
      base_weighted = base_weighted,
      sum_contrib = sum_contrib,
      contrib_residual = contrib_residual,
      mean_hist_gdp = mean_hist_gdp,
      rs_vs_hist_gdp = rs_vs_hist_gdp,
      total_vs_hist_gdp = total_vs_hist_gdp,
      abs_diff_vs_hist_gdp = abs(total_vs_hist_gdp - rs_vs_hist_gdp)
    ),
    Fin_table2,
    RS_fore,
    RSSS
  )
}

# ---------------- Build vintage pairs -----------------
monthly_files   <- list.files(monthly_dir,   pattern = "\\.csv$", full.names = TRUE)
quarterly_files <- list.files(quarterly_dir, pattern = "\\.csv$", full.names = TRUE)

monthly_tbl <- tibble(
  monthly_path = monthly_files,
  monthly_file = basename(monthly_files),
  vintage = gsub("\\.csv$", "", basename(monthly_files))
)

quarterly_tbl <- tibble(
  quarterly_path = quarterly_files,
  quarterly_file = basename(quarterly_files),
  vintage = str_replace(basename(quarterly_files), ".*_(\\d{4})m(\\d{2})\\.csv$", "\\1-\\2")
)

pairs <- inner_join(monthly_tbl, quarterly_tbl, by = "vintage") %>%
  arrange(vintage)
pairs <- pairs[-18, ]

if (nrow(pairs) == 0) stop("No matched monthly/quarterly vintages found. Check file names and directories.")

# ---------------- Run across all matched vintages -----------------
res_list <- map2(
  pairs$monthly_path, pairs$quarterly_path,
  ~ run_forecast(.x, .y, n_start = uni_min_start, ridge = 1e-8, verbose_checks = TRUE)
)

summary_results <- map_dfr(res_list, "summary")
#saveRDS(res_list, file="/Users/boris/Desktop/RSShapley/Combination_Monte/Final/codess/JMP/GDP_constrained/final_list.RData")
ak <- read.csv("/Users/boris/Downloads/Historical vintages of FRED-QD 2018-05 to 2024-12-2/FRED-QD_2025m06.csv",
               stringsAsFactors = FALSE, check.names = FALSE)
ak <- ak[-c(1:2),]

growth_y <- (ak$GDPC1 / dplyr::lag(ak$GDPC1) - 1)*100
growth_y <- growth_y[-1]
growth_y <- as.data.frame(growth_y)
growth_y$Date <- ak$sasdate[-1]

growth_y_sa <- tail(growth_y, n = nrow(pairs)+1)[-18,]

aa <- rep(NA_real_, length(res_list))
bb <- rep(NA_real_, length(res_list))

exclude <- c(1845)  # keep your exclusion
for (i in seq_along(res_list)) {
  if (i %in% exclude) next
  aa[i] <- unique(as.numeric(res_list[[i]][["summary"]][["RS_mean"]]))
  bb[i] <- unique(as.numeric(res_list[[i]][["summary"]][["MEAN_mean"]]))
}

summary_results$actual <- growth_y_sa$growth_y
summary_results <- data.frame(RS_mean = aa, MEAN_mean = bb, actual = growth_y_sa$growth_y)

# drop excluded rows (if they exist)
if (any(exclude <= nrow(summary_results))) summary_results <- summary_results[-exclude, , drop = FALSE]

mean(abs(summary_results$MEAN_mean - summary_results$actual), na.rm = TRUE)
mean(abs(summary_results$RS_mean   - summary_results$actual), na.rm = TRUE)
sqrt(mean((summary_results$MEAN_mean - summary_results$actual)^2, na.rm = TRUE))
sqrt(mean((summary_results$RS_mean   - summary_results$actual)^2, na.rm = TRUE))

e2 <- summary_results$MEAN_mean - summary_results$actual
e1 <- summary_results$RS_mean   - summary_results$actual

dm.test(e1, e2, h = 1, alternative = "less")

summary_results$MEAN_mean[c(83,84)]

mean(abs(summary_results$MEAN_mean[-c(83,84)] - summary_results$actual[-c(83,84)]), na.rm = TRUE)
mean(abs(summary_results$RS_mean[-c(83,84)]   - summary_results$actual[-c(83,84)]), na.rm = TRUE)
sqrt(mean((summary_results$MEAN_mean[-c(83,84)] - summary_results$actual[-c(83,84)])^2, na.rm = TRUE))
sqrt(mean((summary_results$RS_mean[-c(83,84)]   - summary_results$actual[-c(83,84)])^2, na.rm = TRUE))

e2 <- summary_results$MEAN_mean[-c(83,84)] - summary_results$actual[-c(83,84)]
e1 <- summary_results$RS_mean[-c(83,84)]  - summary_results$actual[-c(83,84)]

dm.test(e1, e2, h = 1, alternative = "less")
