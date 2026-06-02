# ================================================================
# Backward pseudo real-time forecasting across multiple vintages
# ================================================================

rm(list = ls())

# ---------------- Packages ----------------
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(lubridate)
  library(purrr); library(tibble); library(stringr); library(forecast)
})

# ---------------- Settings ----------------
monthly_dir <- "/Users/boris/Downloads/Historical FRED-MD Vintages Final"
quarterly_dir <- "/Users/boris/Desktop/RSShapley/Combination_Monte/FRED_QD"

# Univariate estimation window (for the x->y regressions)
uni_window_type <- "expanding" # "expanding" or "rolling" (kept for reference)
uni_window_size <- 40 # used if rolling
uni_min_start <- 40 # first forecasts after at least this many quarters

# RS combination settings
rs_draws <- 1000
rs_seed <- 12893
comb_window_type <- "expanding" # "expanding" or "rolling"
comb_window_size <- 40 # used if rolling

suppressWarnings(RNGkind(kind = "Mersenne-Twister",
                         normal.kind = "Inversion",
                         sample.kind = "Rounding"))


set.seed(rs_seed)
# ---------------- Helpers -----------------
lag1 <- function(x) c(NA, head(x, -1))
d1 <- function(x) x - lag1(x)
d2 <- function(x) d1(d1(x))
lg <- function(x) { y <- x; y[x <= 0] <- NA_real_; log(y) }
d1lg <- function(x) d1(lg(x))
d2lg <- function(x) d2(lg(x))
d_ret <- function(x) { prev <- lag1(x); r <- ifelse(prev <= 0 | is.na(prev), NA_real_, x/prev - 1); d1(r) }
parse_date <- function(x) if (inherits(x,"Date")) x else coalesce(ymd(x, quiet=TRUE), mdy(x, quiet=TRUE))

transform_by_first_row_df <- function(df) {
  codes <- suppressWarnings(as.numeric(df[1, ]))
  out <- df[-1, , drop = FALSE]
  for (i in seq_along(out)) {
    if (!is.numeric(out[[i]])) next
    out[[i]] <- switch(as.character(codes[i]),
                       "1" = out[[i]], "2" = d1(out[[i]]), "3" = d2(out[[i]]),
                       "4" = lg(out[[i]]), "5" = d1lg(out[[i]]), "6" = d2lg(out[[i]]),
                       "7" = d_ret(out[[i]]), out[[i]]
    )
  }
  rownames(out) <- NULL; out
}

# Safe mean over months in a quarter
safe_mean <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

# ---------------- Core building blocks from your script -----------------

univariate_forecasts_split <- function(y, X, y_levels, n_start, expanding = TRUE) {
  stopifnot(length(y) == length(y_levels))
  T <- length(y)
  TX <- nrow(X)
  p <- ncol(X)
  if (n_start < 2) stop("n_start must be >= 2")
  
  last_t <- if (TX == T) T - 1 else T
  
  z_hat <- matrix(NA_real_, nrow = TX, ncol = p, dimnames = list(NULL, colnames(X)))
  lvl_med <- matrix(NA_real_, nrow = TX, ncol = p, dimnames = list(NULL, colnames(X)))
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
      growth_med[t + 1, j] <- ((y_next_med / base_level) - 1)*100
    }
  }
  
  list(
    z_hat = as.data.frame(z_hat),
    lvl_med = as.data.frame(lvl_med),
    growth_med = as.data.frame(growth_med)
  )
}

rand_subsample_from_df <- function(df, y_col, x_next, k, B = 1000,
                                   scale_predictors = FALSE,
                                   weight_scheme = c("uniform", "abs_t"),
                                   alpha = 1,
                                   t_no_intercept = FALSE) {
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
  
  # (only candidates with a non-NA next value)
  cand <- intersect(Xnames_all, names(x_next))
  cand <- cand[!is.na(as.numeric(x_next[1, cand]))]
  if (length(cand) < k) stop("Not enough overlapping predictors between df and x_next to sample k of them.")
  
  ## ---------- (A) selection probabilities ----------
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
  
  ## ---------- (B) RS draws ----------
  forecasts <- numeric(B)
  selected_vars <- matrix(NA_character_, nrow = B, ncol = k)
  coef_entries <- vector("list", B)  # holds tibble(draw, variable, coef)
  n_failed <- 0L
  
  for (b in seq_len(B)) {
    vars <- sample(cand, k, replace = FALSE, prob = prob_vec)
    selected_vars[b, ] <- vars
    
    subdat <- df[, c(y_col, vars), drop = FALSE]
    cc <- stats::complete.cases(subdat)
    if (!any(cc)) { forecasts[b] <- NA_real_; n_failed <- n_failed + 1L; next }
    
    y_train <- as.numeric(subdat[[y_col]][cc])
    X_train <- as.matrix(subdat[cc, vars, drop = FALSE])
    
    x_row <- as.numeric(as.matrix(x_next[1, vars, drop = FALSE]))
    if (any(is.na(x_row))) { forecasts[b] <- NA_real_; n_failed <- n_failed + 1L; next }
    
    if (scale_predictors) {
      mu <- colMeans(X_train)
      sdv <- apply(X_train, 2, sd)
      sdv[sdv == 0 | is.na(sdv)] <- 1
      X_train <- scale(X_train, center = mu, scale = sdv)
      x_row <- (x_row - mu) / sdv
    }
    
    X_train_i <- cbind(`(Intercept)` = 1, X_train)
    x_i <- c(1, x_row)
    
    fit <- try(lm.fit(X_train_i, y_train), silent = TRUE)
    if (inherits(fit, "try-error") || anyNA(fit$coefficients)) {
      forecasts[b] <- NA_real_; n_failed <- n_failed + 1L; next
    }
    
    forecasts[b] <- drop(x_i %*% fit$coefficients)
    
    # store coefficients (exclude intercept)
    cf <- as.numeric(fit$coefficients[-1])
    coef_entries[[b]] <- tibble::tibble(draw = b, variable = vars, coef = cf)
  }
  
  coef_long <- dplyr::bind_rows(coef_entries)
  # be robust if all draws failed
  if (is.null(coef_long) || nrow(coef_long) == 0) {
    coef_long <- tibble::tibble(draw = integer(), variable = character(), coef = numeric())
  }
  
  ## ---------- NEW: selection counts & mean(non-zero coefs) ----------
  # how many times each var was sampled (regardless of fit success)
  sel_tab <- table(factor(as.vector(selected_vars), levels = cand))
  selection_counts <- as.integer(sel_tab)
  names(selection_counts) <- cand
  
  # per-variable stats across successful fits
  coef_stats <- coef_long |>
    dplyr::group_by(variable) |>
    dplyr::summarise(
      n_with_coef = dplyr::n(),
      nz = sum(coef != 0, na.rm = TRUE),
      mean_coef_nonzero = ifelse(nz > 0, mean(coef[coef != 0], na.rm = TRUE), NA_real_),
      mean_coef_all = mean(coef, na.rm = TRUE),
      .groups = "drop"
    )
  # ensure all candidates appear in the table
  coef_stats <- dplyr::right_join(
    coef_stats,
    tibble::tibble(variable = cand),
    by = "variable"
  ) |>
    dplyr::mutate(
      n_with_coef = dplyr::coalesce(n_with_coef, 0L),
      nz = dplyr::coalesce(nz, 0L),
      mean_coef_nonzero = mean_coef_nonzero,   # keep NA if none
      mean_coef_all = mean_coef_all            # NA if never had a coef
    ) |>
    dplyr::arrange(variable)
  
  # convenience named vectors
  mean_coef_nonzero <- stats::setNames(coef_stats$mean_coef_nonzero, coef_stats$variable)
  selection_counts_success <- stats::setNames(coef_stats$n_with_coef, coef_stats$variable)
  
  out <- list(
    forecasts = forecasts,
    selected_vars = selected_vars,          # raw selections (B x k)
    selected_vars_matrix = selected_vars,   # alias for clarity
    coef_long = coef_long,                  # per-draw coefficients
    failed = n_failed,
    probs = prob_vec,
    # NEW:
    selection_counts = selection_counts,                    # sampled counts (all draws)
    selection_counts_success = selection_counts_success,    # had a coef (successful fits only)
    coef_stats = coef_stats,                                # tibble with means & counts
    mean_coef_nonzero = mean_coef_nonzero                   # named vector for convenience
  )
  class(out) <- c("rand_subsample_result", class(out))
  out
}


# Per-quarter RS panel with:
# 1) times selected; 2) mean OLS coef across RS submodels;
# 3) one-step growth prediction from each univariate model;
# 4) historical mean of that univariate prediction up to t.

# ---------------- Settings ----------------
monthly_path <- "/Users/boris/Downloads/Historical FRED-MD Vintages Final/1999-12.csv"
quarterly_path <- "/Users/boris/Desktop/RSShapley/Combination_Monte/FRED_QD/FRED-QD_1999m12.csv"

# ---------------- One-vintage pipeline -----------------
run_forecast <- function(monthly_path, quarterly_path, n_start = 40) {
  # 1) READ & TRANSFORM (row 1 already = transform codes)
  Monthly_FRED <- read.csv(monthly_path, stringsAsFactors = FALSE, check.names = FALSE)
  Quarterly_FRED <- read.csv(quarterly_path, stringsAsFactors = FALSE, check.names = FALSE)
  
  Monthly_FRED_clean <- Monthly_FRED
  Quarterly_FRED_clean <- Quarterly_FRED[-1,]
  Quarterly_raw_for_levels <- Quarterly_FRED_clean # for GDP levels
  
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
  
  # GDP Levels (for realized growth & back-transform)
  y_levels <- tibble(
    Date = parse_date(Quarterly_raw_for_levels[-1, 1, drop = TRUE]) %>% as_date(),
    GDPC1_level = suppressWarnings(as.numeric(as.character(
      Quarterly_raw_for_levels[-1, which(names(Quarterly_raw_for_levels)=="GDPC1"), drop = TRUE]
    )))
  ) %>% arrange(Date) %>%
    drop_na()
  
  # (Optional) Sanity check: code for GDPC1 expected to be 5 (Δlog)
  gdp_code <- suppressWarnings(as.integer(Quarterly_FRED_clean[1, "GDPC1"]))
  if (!is.na(gdp_code) && gdp_code != 5) {
    warning(sprintf("GDPC1 transform code is %s, not 5. Back-transform may be inappropriate.", gdp_code))
  }
  
  # 4) Aggregate monthly → quarterly (mean across months)
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
  
  # 5) Align estimation sample and prepare last X
  estimation0 <- Final_quarterly %>% inner_join(m_q, by = "Date") %>% arrange(Date)
  
  last_x <- m_q %>% filter(Date == max(Date)) %>% select(-Date)
  
  # 6) Run bivariate OLS and back-transform Δlog(GDP) → levels
  y_last_level <- y_levels %>%
    semi_join(estimation0, by = "Date") %>%
    arrange(Date) %>%
    pull(GDPC1_level) %>%
    tail(1)
  
  if (length(y_last_level) == 0 || is.na(y_last_level)) {
    stop("Could not determine the last GDP level aligned with the estimation sample.")
  }
  
  y <- estimation0$GDPC1
  X <- estimation0 %>% select(-Date, -GDPC1)
  
  # Drop variables that have NA in the future x (can't predict)
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
    
    z_hat_next <- as.numeric(predict(fit, newdata = data.frame(x = x_future))) # Δlog forecast
    sigma2 <- summary(fit)$sigma^2
    
    y_next_med <- as.numeric(y_last_level * exp(z_hat_next))
    y_next_mean <- as.numeric(y_last_level * exp(z_hat_next + 0.5 * sigma2)) # bias-corrected mean
    
    tibble(
      variable = varname,
      alpha = unname(coef(fit)[1]),
      beta = unname(coef(fit)[2]),
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
  init_forecast<-results$y_forecast_growth
  
  # Split univariate forecasts across time (for RS combo step)
  X_plus <- rbind(X, last_x)
  fcasts <- univariate_forecasts_split(
    y, X_plus, y_levels = y_levels$GDPC1_level[-1], n_start = n_start
  )
  growth_table <- fcasts$growth_med
  growth_table$Date <- (Final_monthly %>%
                          mutate(q_date = make_date(year(Date), c(3,6,9,12)[quarter(Date)], 1L)) %>%
                          distinct(q_date) %>% arrange(q_date) %>% pull(q_date))[-1]
  # Align growth_y
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
  
  # last available indicators for RS
  last_x_rs <- growth_table %>% filter(Date == max(Date)) %>% select(-Date)
  
  # k = sqrt(# predictors); ensure >=1
  #k <- max(1L, round(sqrt(nrow(estimation) - 1L)))
  k <- max(1L, round(sqrt(nrow(estimation))))
  RSSS <- rand_subsample_from_df(
    estimation, "growth_y", last_x_rs, k,
    B = rs_draws,
    weight_scheme = "uniform",   # <-- turn on |t|-based weighting
    alpha = 1,                 # <-- softmax temperature
    t_no_intercept = FALSE     # <-- TRUE to use no-intercept t-stats
  )
  
  RS_mean <- mean(RSSS[["forecasts"]], na.rm = TRUE)
  RS_forecast<-
  MM<-as.data.frame(RSSS[["selection_counts"]])
  MM<-rownames_to_column(MM, "Var")
  colnames(MM)<-c("Var", "Selected")
  MM$Selected<-MM$Selected/rs_draws
  betas<-as.data.frame(RSSS[["mean_coef_nonzero"]])
  betas<-rownames_to_column(betas, "Var")
  Fin_table<-merge(MM, betas, by=c("Var"))
  Fin_table<-Fin_table %>%
    drop_na()
  colnames(Fin_table)<-c("Var", "Selected", "Mean_OLS")
  
  
  last_x_rsq<-as.data.frame(t(last_x_rs))
  last_x_rsq<-rownames_to_column(last_x_rsq, "Var")
  colnames(last_x_rsq)<-c("Var", "Prediction")
  
  Fin_table1<-merge(Fin_table, last_x_rsq, by=c("Var"))
  
  historic_mean<-estimation %>%
    summarise_if(is.numeric, mean, na.rm = TRUE)
  
  historic_mean<-as.data.frame(t(historic_mean))
  historic_mean<-rownames_to_column(historic_mean, "Var")
  colnames(historic_mean)<-c("Var", "Hist_mean")
  
  Fin_table2<-merge(Fin_table1, historic_mean, by=c("Var"))
  
  aqa<-historic_mean[(historic_mean$Var== "growth_y"),]
  mean_hist_gdp<-aqa$Hist_mean
  
  #RS_mean-mean_hist_gdp
  
  Fin_table2$Shapley<- Fin_table2$Selected*Fin_table2$Mean_OLS*(Fin_table2$Prediction-Fin_table2$Hist_mean)
  #sum(Fin_table2$Shapley)
  RS_fore<-RSSS[["forecasts"]]
  
  # Vintage tags
  vintage_monthly <- basename(monthly_path)
  vintage_quarterly <- basename(quarterly_path)
  vintage_key <- tryCatch({
    # unify to YYYY-MM
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
  
  
  # Return both summary and detailed tables
  list(
    summary = tibble(
      vintage = vintage_key,
      vintage_monthly = vintage_monthly,
      vintage_quarterly = vintage_quarterly,
      MEAN_mean = MEAN_mean,
      RS_mean = RS_mean
    ), 
    Fin_table2, RS_fore, init_forecast
  )
}  # <-- close run_forecast() here


# ---------------- Build vintage pairs -----------------

monthly_files <- list.files(monthly_dir, pattern = "\\.csv$", full.names = TRUE)
quarterly_files <- list.files(quarterly_dir, pattern = "\\.csv$", full.names = TRUE)

monthly_tbl <- tibble(
  monthly_path = monthly_files,
  monthly_file = basename(monthly_files),
  vintage = gsub("\\.csv$", "", basename(monthly_files)) # e.g., "2018-06"
)

quarterly_tbl <- tibble(
  quarterly_path = quarterly_files,
  quarterly_file = basename(quarterly_files),
  # Extract "YYYY-MM" from "..._YYYYmMM.csv"
  vintage = str_replace(basename(quarterly_files), ".*_(\\d{4})m(\\d{2})\\.csv$", "\\1-\\2")
)

pairs <- inner_join(monthly_tbl, quarterly_tbl, by = "vintage") %>%
  arrange(vintage)

if (nrow(pairs) == 0) stop("No matched monthly/quarterly vintages found. Check file names and directories.")

# ---------------- Run across all matched vintages -----------------

res_list <- map2(pairs$monthly_path, pairs$quarterly_path, ~ run_forecast(.x, .y))

summary_results <- map_dfr(res_list, "summary")


ak<-read.csv("/Users/boris/Downloads/Historical vintages of FRED-QD 2018-05 to 2024-12-2/FRED-QD_2025m06.csv", stringsAsFactors = FALSE, check.names = FALSE)

ak<-ak[-c(1:2),]

growth_y <- (ak$GDPC1/ dplyr::lag(ak$GDPC1) - 1)*100
growth_y <- growth_y[-1]
growth_y <- as.data.frame(growth_y)
growth_y$Date <- ak$sasdate[-1]

growth_y_sa<-tail(growth_y, n = nrow(pairs))
#summary_results<-summary_results[-nrow(summary_results),]

aa<-vector("numeric", 103L)
bb<-vector("numeric", 103L)
exclude <- c(18)
for(i in 1:length(res_list)){
  if(i == 18) {
    next
  }
  aa[i]<-unique(as.numeric(res_list[[i]][["summary"]][["RS_mean"]]))
  bb[i]<-unique(res_list[[i]][["summary"]][["MEAN_mean"]])
}


summary_results$actual<-growth_y_sa$growth_y
summary_results<-data.frame(RS_mean=aa, MEAN_mean=bb, actual=growth_y_sa$growth_y)
summary_results<-summary_results[-18,]

mean(abs(summary_results$MEAN_mean-summary_results$actual))
mean(abs(summary_results$RS_mean-summary_results$actual))
sqrt(mean((summary_results$MEAN_mean-summary_results$actual)^2))
sqrt(mean((summary_results$RS_mean-summary_results$actual)^2))

e2<-summary_results$MEAN_mean-summary_results$actual
e1<-summary_results$RS_mean-summary_results$actual

dm.test(e1,e2,h=1, alternative = "less")

histogram <- map_dfr(res_list, "[[3]]")


out <- matrix(data = NA_integer_, nrow = length(res_list), ncol = rs_draws)
for(i in 1:length(res_list)){
  if(i == 18) {
    next
  }
  out[i, ]<-res_list[[i]][[3]]
}

out <- matrix(data = NA_integer_, nrow = length(res_list), ncol = rs_draws)
for(i in 1:length(res_list)){
  if(i == 18) {
    next
  }
  ae<-res_list[[i]][["summary"]][["vintage"]]
  #unpack a 3D array
  mat <-res_list[[i]][[2]]
  form = sprintf('Shapley_%s.csv', ae)
  write.csv(mat, file = form)
  
}

out <- matrix(data = NA_integer_, nrow = length(res_list), ncol = length(bb))
for(i in 1:length(res_list)){
  if(i == 18) {
    next
  }
  out[i, ]<-res_list[[i]][[4]]
}
getwd()
out<-as.data.frame(out[-18,])
write.csv(summary_results, "/Users/boris/Desktop/RSShapley/Combination_Monte/Final/codess/JMP/uniform_forecast_error.csv")
write.csv(out, "/Users/boris/Desktop/RSShapley/Combination_Monte/Final/codess/JMP/RS_predi_for_histograms.csv")

saveRDS(res_list, file="/Users/boris/Desktop/RSShapley/Combination_Monte/Final/codess/JMP/Uniform_shapley_per_quarter/final_list.RData")

