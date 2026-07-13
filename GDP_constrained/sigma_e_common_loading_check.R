#!/usr/bin/env Rscript
# =============================================================================
# GDP_ar_bivariate_allk.R
#
# Like GDP_all_k_results.R but individual forecasts use AR(P)+X bivariate
# regressions (following Elliott et al. 2015 / Maung 2023): each predictor's
# individual forecast is produced from  y ~ ylag_1 + ... + ylag_P + x_j
# where P is selected by AIC at each point in time.
#
# The RSM combination step is identical to the baseline script.
# Produces the same tables/plots plus a comparison figure overlaying the
# RMSFE-by-k curves from both setups.
# =============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(lubridate)
  library(purrr); library(tibble); library(stringr)
  library(forecast); library(glmnet); library(ranger)
})

# ---- Paths ------------------------------------------------------------------
WORKDIR       <- "/Users/boris/Desktop/RSShapley/Combination_Monte/Final/codess/JMP"
monthly_dir   <- file.path(WORKDIR, "Historical FRED-MD Vintages Final")
quarterly_dir <- file.path(WORKDIR, "FRED_QD")
OUTDIR        <- file.path(WORKDIR, "GDP_constrained", "ar_bivariate_allk_results")
ORIG_OUTDIR   <- file.path(WORKDIR, "GDP_constrained", "all_k_results")  # for comparison
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ---- Parameters -------------------------------------------------------------
rs_draws_main <- 2L
rs_draws_allk <- 2L
rs_seed       <- 12893L
k_grid_app    <- 2:4
uni_min_start <- 40L
p_max_ar      <- 4L          # max AR lag order for AIC selection
covid_outlier_dates <- as.Date(c("2020-06-01", "2020-09-01"))

suppressWarnings(RNGkind(kind = "Mersenne-Twister",
                         normal.kind = "Inversion",
                         sample.kind = "Rounding"))
set.seed(rs_seed)

cat("GDP AR-bivariate all-k analysis started.\n")
cat("  Individual forecasts: AR(P)+X bivariate, P chosen by AIC (p_max =", p_max_ar, ")\n")
cat("  Output directory:", OUTDIR, "\n\n")

# ---- Helper functions (identical to GDP_all_k_results.R) --------------------
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
parse_date <- function(x) {
  if (inherits(x, "Date")) x else coalesce(ymd(x, quiet = TRUE), mdy(x, quiet = TRUE))
}

transform_by_first_row_df <- function(df) {
  codes <- suppressWarnings(as.numeric(df[1, ]))
  out   <- df[-1, , drop = FALSE]
  for (i in seq_along(out)) {
    if (!is.numeric(out[[i]])) next
    out[[i]] <- switch(as.character(codes[i]),
      "1" = out[[i]], "2" = d1(out[[i]]),  "3" = d2(out[[i]]),
      "4" = lg(out[[i]]), "5" = d1lg(out[[i]]), "6" = d2lg(out[[i]]),
      "7" = d_ret(out[[i]]), out[[i]])
  }
  rownames(out) <- NULL
  out
}

safe_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)

ols_sum1 <- function(X, y, ridge = 1e-8) {
  X <- as.matrix(X); y <- as.numeric(y); k <- ncol(X)
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

rand_subsample_from_df <- function(df, y_col, x_next, k, B = 1000,
                                   scale_predictors = FALSE,
                                   weight_scheme = c("uniform", "abs_t"),
                                   alpha = 1, t_no_intercept = FALSE,
                                   ridge = 1e-8) {
  weight_scheme <- match.arg(weight_scheme)
  if (!y_col %in% names(df)) stop("y_col not found in df.")
  y <- as.numeric(df[[y_col]])
  Xnames_all <- setdiff(names(df), y_col)
  if (is.vector(x_next)) {
    if (is.null(names(x_next))) stop("x_next must be named if vector.")
    x_next <- as.data.frame(as.list(x_next))
  } else { x_next <- as.data.frame(x_next) }
  if (nrow(x_next) != 1) stop("x_next must be a single row.")
  cand <- intersect(Xnames_all, names(x_next))
  cand <- cand[!is.na(as.numeric(x_next[1, cand]))]
  if (length(cand) < k) stop("Not enough predictors.")
  prob_vec <- rep(1/length(cand), length(cand))
  names(prob_vec) <- cand
  forecasts     <- rep(NA_real_, B)
  coef_entries  <- vector("list", B)
  n_failed <- 0L
  for (b in seq_len(B)) {
    vars   <- sample(cand, k, replace = FALSE, prob = prob_vec)
    subdat <- df[, c(y_col, vars), drop = FALSE]
    cc     <- stats::complete.cases(subdat)
    if (!any(cc)) { n_failed <- n_failed + 1L; next }
    y_train <- as.numeric(subdat[[y_col]][cc])
    X_train <- as.matrix(subdat[cc, vars, drop = FALSE])
    x_row   <- as.numeric(as.matrix(x_next[1, vars, drop = FALSE]))
    if (any(is.na(x_row))) { n_failed <- n_failed + 1L; next }
    if (scale_predictors) {
      mu <- colMeans(X_train); sdv <- apply(X_train, 2, sd)
      sdv[sdv == 0 | is.na(sdv)] <- 1
      X_train <- scale(X_train, center = mu, scale = sdv)
      x_row   <- (x_row - mu) / sdv
    }
    w_hat <- ols_sum1(X_train, y_train, ridge = ridge)
    if (anyNA(w_hat)) { n_failed <- n_failed + 1L; next }
    forecasts[b]    <- drop(x_row %*% w_hat)
    coef_entries[[b]] <- tibble(draw = b, variable = vars, coef = as.numeric(w_hat))
  }
  coef_long <- dplyr::bind_rows(coef_entries)
  if (is.null(coef_long) || nrow(coef_long) == 0)
    coef_long <- tibble(draw = integer(), variable = character(), coef = numeric())
  B_succ <- dplyr::n_distinct(coef_long$draw)
  w_bar  <- if (B_succ > 0) {
    coef_long %>% group_by(variable) %>%
      summarise(w_bar = sum(coef, na.rm = TRUE) / B_succ, .groups = "drop")
  } else { tibble(variable = cand, w_bar = 0) }
  w_bar_vec <- stats::setNames(w_bar$w_bar, w_bar$variable)
  list(forecasts = forecasts, B_succ = B_succ, w_bar = w_bar_vec,
       failed = n_failed, coef_long = coef_long)
}

prep_combo_data <- function(df, y_col, x_next) {
  df <- as.data.frame(df)
  if (is.vector(x_next)) {
    if (is.null(names(x_next))) stop("x_next must be named.")
    x_next <- as.data.frame(as.list(x_next))
  } else { x_next <- as.data.frame(x_next) }
  if ("Date" %in% names(df))     df$Date     <- NULL
  if ("Date" %in% names(x_next)) x_next$Date <- NULL
  y    <- as.numeric(df[[y_col]])
  X_all <- df[, setdiff(names(df), y_col), drop = FALSE]
  cand <- intersect(names(X_all), names(x_next))
  cand <- cand[!is.na(as.numeric(x_next[1, cand]))]
  if (length(cand) < 1) stop("No usable predictors.")
  X    <- X_all[, cand, drop = FALSE]
  xrow <- x_next[, cand, drop = FALSE]
  nzv  <- sapply(X, function(z) { z <- z[is.finite(z)]; length(z) >= 2 && sd(z, na.rm=TRUE) > 0 })
  cand <- cand[nzv]; X <- X[, cand, drop=FALSE]; xrow <- xrow[, cand, drop=FALSE]
  cc <- stats::complete.cases(y, X)
  list(y = y[cc], X = X[cc, , drop=FALSE], xrow = xrow, vars = cand)
}

glmnet_combo_forecast <- function(df, y_col, x_next, alpha = 1,
                                   lambda_choice = c("lambda.1se","lambda.min"),
                                   nfolds = 10, standardize = TRUE) {
  lambda_choice <- match.arg(lambda_choice)
  dat  <- prep_combo_data(df, y_col, x_next)
  y    <- dat$y; X <- as.matrix(dat$X); xrow <- as.matrix(dat$xrow); n <- length(y)
  if (n < 8) {
    fit  <- glmnet(X, y, alpha = alpha, standardize = standardize)
    lam  <- fit$lambda[length(fit$lambda)]
    pred <- as.numeric(predict(fit, newx = xrow, s = lam))
    co   <- as.matrix(coef(fit, s = lam))[, 1]
    names(co) <- rownames(as.matrix(coef(fit, s = lam)))
    return(list(forecast = pred, lambda = lam, coefs = co, vars = dat$vars))
  }
  nfolds_eff <- min(nfolds, n)
  cvfit <- cv.glmnet(X, y, alpha = alpha, nfolds = nfolds_eff, standardize = standardize)
  lam  <- if (lambda_choice == "lambda.1se") cvfit$lambda.1se else cvfit$lambda.min
  pred <- as.numeric(predict(cvfit, newx = xrow, s = lam))
  co   <- as.matrix(coef(cvfit, s = lam))[, 1]
  names(co) <- rownames(as.matrix(coef(cvfit, s = lam)))
  list(forecast = pred, lambda = lam, coefs = co, vars = dat$vars)
}

rf_combo_forecast <- function(df, y_col, x_next, num_trees = 500,
                               mtry = NULL, min_node_size = 5, seed = 123) {
  dat   <- prep_combo_data(df, y_col, x_next)
  train <- data.frame(y = dat$y, dat$X)
  test  <- data.frame(dat$xrow)
  p     <- ncol(dat$X)
  if (is.null(mtry)) mtry <- max(1, floor(sqrt(p)))
  set.seed(seed)
  rf   <- ranger::ranger(y ~ ., data = train, num.trees = num_trees, mtry = mtry,
                         min.node.size = min_node_size, importance = "none", seed = seed)
  pred <- as.numeric(predict(rf, data = test)$predictions)
  list(forecast = pred, vars = dat$vars)
}

# =============================================================================
# NEW: AR(P)+X bivariate individual forecasts (Maung / Elliott et al. style)
# =============================================================================
# For each predictor j and each time t, fits:
#   y[s] = a0 + a1*y[s-1] + ... + aP*y[s-P] + b*x_j[s] + e[s]
# where P is chosen by AIC from {0,...,p_max} using the pure AR model on y.
# Prediction for y[t+1]: uses y[t],...,y[t-P+1] and x_j[t+1].
# =============================================================================
univariate_forecasts_ar_split <- function(y, X, y_levels, n_start, p_max = 4L,
                                          expanding = TRUE, dates = NULL,
                                          exclude_dates = as.Date(character(0))) {
  T_y  <- length(y)
  TX   <- nrow(X)
  pcol <- ncol(X)
  if (n_start < 2) stop("n_start must be >= 2")
  excl_idx <- if (!is.null(dates)) which(as.Date(dates) %in% exclude_dates) else integer(0)
  last_t   <- if (TX == T_y) T_y - 1L else T_y

  z_hat      <- matrix(NA_real_, nrow = TX, ncol = pcol, dimnames = list(NULL, colnames(X)))
  lvl_med    <- matrix(NA_real_, nrow = TX, ncol = pcol, dimnames = list(NULL, colnames(X)))
  growth_med <- matrix(NA_real_, nrow = TX, ncol = pcol, dimnames = list(NULL, colnames(X)))

  for (t in seq(n_start, last_t)) {
    idx <- if (expanding) seq_len(t) else seq(t - n_start + 1L, t)
    if (length(excl_idx) > 0L) idx <- setdiff(idx, excl_idx)
    n_idx <- length(idx)
    if (n_idx < max(3L, p_max + 2L)) next

    y_in <- y[idx]

    # ---- AIC-based AR lag order selection (once per t) ----------------------
    # Build AR lag matrix for 1:p_max within idx
    # Y_lag_list[[l]] = lag-l of y evaluated at idx
    Y_lag_list <- lapply(seq_len(p_max), function(l) {
      lag_positions <- idx - l
      out <- rep(NA_real_, length(idx))
      valid <- lag_positions >= 1L
      if (any(valid)) out[valid] <- y[lag_positions[valid]]
      out
    })

    # Baseline AIC: intercept-only model (compute directly from RSS)
    cc0 <- !is.na(y_in)
    if (sum(cc0) < 2L) next
    n0 <- sum(cc0)
    rss0 <- sum((y_in[cc0] - mean(y_in[cc0]))^2)
    # AIC = n*log(RSS/n) + 2*k, k=2 (intercept + sigma)
    aic_fn <- function(rss, n, k) n * log(rss / n) + 2 * k
    best_aic <- aic_fn(rss0, n0, 2L)
    best_p   <- 0L

    for (pp in seq_len(p_max)) {
      Y_lag_pp <- do.call(cbind, Y_lag_list[seq_len(pp)])
      cc_pp    <- complete.cases(y_in, Y_lag_pp)
      if (sum(cc_pp) < pp + 2L) next
      Xm_pp  <- cbind(1, Y_lag_pp[cc_pp, , drop = FALSE])
      beta_pp <- tryCatch(
        solve(crossprod(Xm_pp), crossprod(Xm_pp, y_in[cc_pp])),
        error = function(e) NULL)
      if (is.null(beta_pp)) next
      rss_pp  <- sum((y_in[cc_pp] - Xm_pp %*% beta_pp)^2)
      aic_pp  <- aic_fn(rss_pp, sum(cc_pp), pp + 2L)
      if (is.finite(aic_pp) && aic_pp < best_aic) {
        best_aic <- aic_pp
        best_p   <- pp
      }
    }
    p_star <- best_p

    # ---- AR prediction values for time t+1 ----------------------------------
    # lag l for prediction = y[t - l + 1]
    ar_pred <- if (p_star > 0L) {
      vals <- vapply(seq_len(p_star), function(l) {
        pos <- t - l + 1L
        if (pos >= 1L) y[pos] else NA_real_
      }, numeric(1))
      vals
    } else {
      numeric(0)
    }

    # If any AR prediction value is NA, fall back to p=0 for this t
    if (p_star > 0L && any(is.na(ar_pred))) p_star <- 0L

    # Pre-build lag matrix used for all j (avoids repeated construction)
    if (p_star > 0L) {
      Y_lag_star <- do.call(cbind, Y_lag_list[seq_len(p_star)])
    }

    base_level <- y_levels[t]

    # ---- Per-predictor AR(p_star)+X regression ------------------------------
    for (j in seq_len(pcol)) {
      x_next <- X[t + 1L, j]
      if (is.na(x_next)) next

      x_in <- X[idx, j]

      if (p_star == 0L) {
        # Pure bivariate via matrix algebra (avoids lm/predict scoping issues)
        cc <- complete.cases(y_in, x_in)
        if (sum(cc) < 3L || length(unique(x_in[cc])) < 2L) next
        Xm <- cbind(1, x_in[cc])
        beta <- tryCatch(
          solve(crossprod(Xm), crossprod(Xm, y_in[cc])),
          error = function(e) rep(NA_real_, 2L))
        if (anyNA(beta) || any(!is.finite(beta))) next
        z <- as.numeric(c(1, x_next) %*% beta)
      } else {
        # AR(p_star) + X: use matrix algebra to avoid formula/newdata issues
        cc <- complete.cases(y_in, Y_lag_star, x_in)
        if (sum(cc) < p_star + 3L) next
        if (length(unique(x_in[cc])) < 2L) next

        y_tr  <- y_in[cc]
        # Design matrix: [intercept, ylag1, ..., ylabP, x_j]
        Xm_tr <- cbind(1, Y_lag_star[cc, , drop = FALSE], x_in[cc])

        beta <- tryCatch(
          solve(crossprod(Xm_tr), crossprod(Xm_tr, y_tr)),
          error = function(e) rep(NA_real_, ncol(Xm_tr)))
        if (anyNA(beta) || any(!is.finite(beta))) next

        # Prediction: [1, y[t], y[t-1], ..., y[t-p+1], x_next]
        x_pred <- c(1, ar_pred, x_next)
        z <- as.numeric(x_pred %*% beta)
      }

      if (length(z) != 1L || is.na(z) || !is.finite(z)) next

      z_hat[t + 1L, j]      <- z
      y_next_med             <- base_level * exp(z)
      lvl_med[t + 1L, j]    <- y_next_med
      growth_med[t + 1L, j] <- ((y_next_med / base_level) - 1) * 100
    }
  }

  list(z_hat      = as.data.frame(z_hat),
       lvl_med    = as.data.frame(lvl_med),
       growth_med = as.data.frame(growth_med))
}

# ---- Main per-vintage pipeline (AR version) ---------------------------------
run_forecast_allk_ar <- function(monthly_path, quarterly_path,
                                  n_start = 40L, ridge = 1e-8,
                                  k_grid = 2:55,
                                  B_main = 1000L, B_allk = 200L,
                                  p_max = 4L) {

  Monthly_FRED   <- read.csv(monthly_path,   stringsAsFactors = FALSE, check.names = FALSE)
  Quarterly_FRED <- read.csv(quarterly_path, stringsAsFactors = FALSE, check.names = FALSE)

  Monthly_FRED_clean        <- Monthly_FRED
  Quarterly_FRED_clean      <- Quarterly_FRED[-1, ]
  Quarterly_raw_for_levels  <- Quarterly_FRED_clean

  Monthly_without_time <- Monthly_FRED_clean[, -1, drop = FALSE] %>%
    mutate(across(everything(), ~ suppressWarnings(as.numeric(as.character(.)))))
  Final_monthly <- transform_by_first_row_df(Monthly_without_time)
  Final_monthly$Date <- parse_date(Monthly_FRED_clean[-1, 1, drop = TRUE])

  Quarterly_without_time <- Quarterly_FRED_clean[, -1, drop = FALSE] %>%
    mutate(across(everything(), ~ suppressWarnings(as.numeric(as.character(.)))))
  Final_quarterly <- transform_by_first_row_df(Quarterly_without_time)
  Final_quarterly$Date <- parse_date(Quarterly_FRED_clean[-1, 1, drop = TRUE]) %>% as_date()

  y_levels <- tibble(
    Date = parse_date(Quarterly_raw_for_levels[-1, 1, drop = TRUE]) %>% as_date(),
    GDPC1_level = suppressWarnings(as.numeric(as.character(
      Quarterly_raw_for_levels[-1,
        which(names(Quarterly_raw_for_levels) == "GDPC1"), drop = TRUE])))
  ) %>% arrange(Date) %>% drop_na()

  m_q <- Final_monthly %>%
    mutate(Date = parse_date(Date),
           q_date = make_date(year(Date), c(3,6,9,12)[quarter(Date)], 1L)) %>%
    group_by(q_date) %>%
    summarise(across(where(is.numeric), safe_mean), .groups = "drop") %>%
    rename(Date = q_date) %>% arrange(Date)

  Final_quarterly <- Final_quarterly[c("GDPC1","Date")] %>% arrange(Date)
  Final_quarterly <- Final_quarterly[-1, ]
  m_q <- m_q[-1, ]

  estimation0 <- Final_quarterly %>% inner_join(m_q, by = "Date") %>% arrange(Date)
  last_x      <- m_q %>% filter(Date == max(Date)) %>% select(-Date)

  y_last_level <- y_levels %>% semi_join(estimation0, by = "Date") %>%
    arrange(Date) %>% pull(GDPC1_level) %>% tail(1)
  if (length(y_last_level) == 0 || is.na(y_last_level))
    stop("Could not determine last GDP level.")

  apply_covid_excl <- max(estimation0$Date, na.rm = TRUE) > max(covid_outlier_dates)
  keep_rows0 <- if (apply_covid_excl) !(estimation0$Date %in% covid_outlier_dates) else
    rep(TRUE, nrow(estimation0))

  y      <- estimation0$GDPC1
  X      <- estimation0 %>% select(-Date, -GDPC1)
  keep_vars <- names(X)[!is.na(as.numeric(last_x[1, names(X)]))]
  X      <- X[, keep_vars, drop = FALSE]
  last_x <- last_x[, keep_vars, drop = FALSE]

  # Simple bivariate-OLS mean benchmark (unchanged from baseline)
  results <- map_dfr(names(X), function(varname) {
    x <- estimation0[[varname]]
    if (all(is.na(x)) || length(unique(na.omit(x))) < 2)
      return(tibble(variable=varname, z_hat_next=NA_real_, y_forecast_next=NA_real_))
    fit     <- lm(y[keep_rows0] ~ x[keep_rows0])
    x_future <- as.numeric(last_x[[varname]])
    if (is.na(x_future))
      return(tibble(variable=varname, z_hat_next=NA_real_, y_forecast_next=NA_real_))
    z_hat_next <- as.numeric(predict(fit, newdata = data.frame(x = x_future)))
    tibble(variable=varname, z_hat_next=z_hat_next,
           y_forecast_next = y_last_level * exp(z_hat_next))
  })
  MEAN_bivariate <- mean(((results$y_forecast_next / y_last_level) - 1) * 100, na.rm = TRUE)

  # ---- AR(P)+X individual forecasts ----------------------------------------
  X_plus <- rbind(X, last_x)
  fcasts  <- univariate_forecasts_ar_split(
    y, X_plus,
    y_levels      = y_levels$GDPC1_level[-1],
    n_start       = n_start,
    p_max         = p_max,
    dates         = estimation0$Date,
    exclude_dates = if (apply_covid_excl) covid_outlier_dates else as.Date(character(0))
  )

  growth_table <- fcasts$growth_med
  growth_table$Date <- (Final_monthly %>%
    mutate(q_date = make_date(year(Date), c(3,6,9,12)[quarter(Date)], 1L)) %>%
    distinct(q_date) %>% arrange(q_date) %>% pull(q_date))[-1]

  growth_y_rs <- (y_levels$GDPC1_level / dplyr::lag(y_levels$GDPC1_level) - 1) * 100
  growth_y_rs <- growth_y_rs[-1]
  growth_y_df <- data.frame(growth_y = growth_y_rs, Date = y_levels$Date[-1])

  estimation <- suppressMessages(
    growth_y_df %>% inner_join(growth_table, by = "Date") %>% arrange(Date) %>% drop_na()
  )
  if (apply_covid_excl)
    estimation <- estimation %>% filter(!(Date %in% covid_outlier_dates))

  last_x_rs   <- growth_table %>% filter(Date == max(Date)) %>% select(-Date)
  MEAN_fore_rs <- mean(as.numeric(last_x_rs[1, ]), na.rm = TRUE)

  # ---- Main RSM at k = round(sqrt(T)) --------------------------------------
  k_main <- max(1L, round(sqrt(nrow(estimation))))
  RSSS_main <- rand_subsample_from_df(estimation, "growth_y", last_x_rs, k_main,
                                      B = B_main, ridge = ridge)
  RS_fore_main <- mean(RSSS_main$forecasts, na.rm = TRUE)

  # LASSO, Ridge, RF
  LASSO_fore <- tryCatch(
    glmnet_combo_forecast(estimation, "growth_y", last_x_rs, alpha=1)$forecast,
    error = function(e) NA_real_)
  RIDGE_fore <- tryCatch(
    glmnet_combo_forecast(estimation, "growth_y", last_x_rs, alpha=0)$forecast,
    error = function(e) NA_real_)
  RF_fore <- tryCatch(
    rf_combo_forecast(estimation, "growth_y", last_x_rs, num_trees=500, seed=rs_seed)$forecast,
    error = function(e) NA_real_)

  # ---- All-k sweep ---------------------------------------------------------
  T_est <- nrow(estimation)
  allk_fore <- setNames(
    vapply(k_grid, function(kk) {
      if (kk >= T_est - 2L || kk < 2L) return(NA_real_)
      tryCatch({
        r <- rand_subsample_from_df(estimation, "growth_y", last_x_rs, kk,
                                    B = B_allk, ridge = ridge)
        mean(r$forecasts, na.rm = TRUE)
      }, error = function(e) NA_real_)
    }, numeric(1)),
    paste0("k", k_grid)
  )

  vintage_key <- tryCatch({
    vm <- gsub("\\.csv$", "", basename(monthly_path))
    if (grepl("^\\d{4}-\\d{2}$", vm)) vm else
      sub(".*_(\\d{4})m(\\d{2})$", "\\1-\\2",
          gsub("\\.csv$", "", basename(quarterly_path)))
  }, error = function(e) NA_character_)

  cat("  Vintage", vintage_key, "| p*_mode ~", "AIC-sel", "| k_main =", k_main,
      "| T_est =", T_est, "| RS_AR =", round(RS_fore_main, 3), "\n")

  list(vintage        = vintage_key,
       target_date    = max(growth_table$Date),
       target_year    = year(max(growth_table$Date)),
       target_q       = quarter(max(growth_table$Date)),
       k_main         = k_main,
       MEAN_bivariate = MEAN_bivariate,
       MEAN_fore_rs   = MEAN_fore_rs,
       RS_fore        = RS_fore_main,
       LASSO_fore     = LASSO_fore,
       RIDGE_fore     = RIDGE_fore,
       RF_fore        = RF_fore,
       covid_excl     = apply_covid_excl,
       allk_fore      = allk_fore,
       indiv_estimation = estimation)
}

# ---- Build vintage pairs ----------------------------------------------------
monthly_files   <- sort(list.files(monthly_dir,   pattern = "\\.csv$", full.names = TRUE))
quarterly_files <- sort(list.files(quarterly_dir, pattern = "\\.csv$", full.names = TRUE))

monthly_tbl <- tibble(
  monthly_path = monthly_files,
  vintage = gsub("\\.csv$", "", basename(monthly_files))
)
quarterly_tbl <- tibble(
  quarterly_path = quarterly_files,
  vintage = str_replace(basename(quarterly_files),
                        ".*_(\\d{4})m(\\d{2})\\.csv$", "\\1-\\2")
)
pairs <- inner_join(monthly_tbl, quarterly_tbl, by = "vintage") %>% arrange(vintage)

exclude_pair_idx <- 18L   # Q4 2003 vintage (missing values)
pairs <- pairs[-exclude_pair_idx, , drop = FALSE]

cat("Running", nrow(pairs), "vintages (AR+X individual forecasts)...\n")
cat("k_grid:", min(k_grid_app), "to", max(k_grid_app),
    "| B_allk =", rs_draws_allk, "\n\n")

# Log file
log_file <- file.path(OUTDIR, "run_log.txt")
sink(log_file, append = FALSE, split = TRUE)
cat("GDP AR-bivariate all-k | started:", as.character(Sys.time()), "\n")
cat("p_max =", p_max_ar, "| B_main =", rs_draws_main, "| B_allk =", rs_draws_allk, "\n\n")

# ---- Run all vintages -------------------------------------------------------
pairs <- tail(pairs, 1)
res_list <- map2(
  pairs$monthly_path, pairs$quarterly_path,
  ~ run_forecast_allk_ar(.x, .y, n_start = uni_min_start, ridge = 1e-8,
                         k_grid = k_grid_app, p_max = p_max_ar,
                         B_main = rs_draws_main, B_allk = rs_draws_allk)
)

saveRDS(res_list, file.path(OUTDIR, "res_list_allk_ar.rds"))

# ===== Sigma_e common-loading diagnostics =====
est <- res_list[[1]]$indiv_estimation
dt  <- est$Date
y   <- est$growth_y
F   <- as.matrix(est[, setdiff(names(est), c("growth_y","Date"))])
E   <- y - F                                   # error panel: T x m  (e_{t,j}=actual-forecast)
ok  <- colSums(is.na(E)) == 0
E   <- E[, ok, drop=FALSE]
cat(sprintf("\nPanel: T=%d periods, m=%d predictors (complete-case)\n", nrow(E), ncol(E)))
S   <- cov(E)                                  # m x m error covariance
m   <- ncol(S); offd <- S[upper.tri(S)]
dg  <- diag(S)
cat(sprintf("mean diag var = %.4f ; mean off-diag cov = %.4f\n", mean(dg), mean(offd)))
# common-loading fit: q = mean off-diag, D = diag - q
q  <- mean(offd); Dhat <- dg - q
cat(sprintf("common-loading q (common-shock var) = %.4f\n", q))
cat(sprintf("off-diag cov: sd = %.4f, CV = %.2f, min = %.4f, max = %.4f\n",
            sd(offd), sd(offd)/abs(mean(offd)), min(offd), max(offd)))
# correlations
R  <- cov2cor(S); offr <- R[upper.tri(R)]
cat(sprintf("off-diag CORR: mean = %.3f, sd = %.3f, min = %.3f, max = %.3f\n",
            mean(offr), sd(offr), min(offr), max(offr)))
# fraction of pairs with corr>0.9
cat(sprintf("share of pairs |corr|>0.9: %.1f%%\n", 100*mean(abs(offr)>0.9)))
# eigen / leading factor
ev <- eigen(S, symmetric=TRUE)
lam <- ev$values; v1 <- ev$vectors[,1]
cat(sprintf("top-5 eigvals: %s\n", paste(round(lam[1:min(5,m)],3), collapse=", ")))
cat(sprintf("var explained by factor 1: %.1f%% ; factor 2: %.1f%%\n",
            100*lam[1]/sum(lam), 100*lam[2]/sum(lam)))
# is leading eigenvector ~ iota (equal/common-loading signature)?
iota <- rep(1/sqrt(m), m)
cos_iota <- abs(sum(v1*iota))                  # |cos angle| with equal-weights direction
cat(sprintf("|cos(angle(eigvec1, iota))| = %.3f  (1=common loading)\n", cos_iota))
cat(sprintf("eigvec1 sign-consistent (all same sign): %s ; loadings CV = %.2f\n",
            all(sign(v1)==sign(v1[1])), sd(abs(v1))/mean(abs(v1))))
# best common-loading reconstruction error
CL <- matrix(q, m, m); diag(CL) <- dg          # q*iota*iota' + D
relF_all  <- norm(S-CL,"F")/norm(S,"F")
relF_off  <- norm(S-CL,"F")/norm(S-diag(dg),"F")  # vs off-diagonal mass only
cat(sprintf("rel Frobenius ||S-CL||/||S|| = %.3f ; vs off-diag mass = %.3f\n", relF_all, relF_off))
# residual after removing rank-1 common shock: how diagonal is D really?
S_minus_q <- S - matrix(q,m,m)                 # remove q*iota*iota'
res_off   <- S_minus_q[upper.tri(S_minus_q)]
cat(sprintf("after removing q*iota*iota': residual off-diag cov mean=%.4f sd=%.4f (D diagonal => ~0)\n",
            mean(res_off), sd(res_off)))
cat(sprintf("residual off-diag / diag ratio (mean|offd|/mean diag) = %.3f\n",
            mean(abs(res_off))/mean(dg)))
saveRDS(list(S=S, dates=dt, y=y), "_sigma_e_panel.rds")
cat("\nDONE\n")
