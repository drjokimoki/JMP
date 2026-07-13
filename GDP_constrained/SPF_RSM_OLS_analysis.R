#!/usr/bin/env Rscript
# =============================================================================
# SPF_RSM_OLS_analysis.R
#
# Full RSM for SPF with within-subset sum-to-one OLS weights.
# This differs from SPF_RSM_analysis.R, which uses simple subset averages.
#
# For each survey round:
#   1. Draw k currently available respondent IDs.
#   2. Estimate sum-to-one OLS weights using past rounds where all selected
#      respondents and the actual outcome are observed.
#   3. Forecast the current target and bag over B successful subset draws.
#
# Draws without enough complete historical observations are skipped. No equal-
# weight fallback is used, so the reported forecasts are genuinely OLS-weighted.
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(lubridate)
})

setwd("/Users/boris/Desktop/RSShapley/Combination_Monte/Final/codess/JMP")

DATADIR <- "GDP_constrained"
OUTDIR  <- "GDP_constrained/SPF_ols_results"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

RSM_SEED  <- 20260704L
RSM_DRAWS <- as.integer(Sys.getenv("SPF_OLS_DRAWS", unset = "500"))
K_GRID    <- 2:20
RUN_ONLY  <- strsplit(Sys.getenv("SPF_OLS_SERIES", unset = ""), ",")[[1]]
RUN_ONLY  <- trimws(RUN_ONLY[nzchar(trimws(RUN_ONLY))])
MIN_N     <- 5L
MIN_TRAIN_BASE <- as.integer(Sys.getenv("SPF_OLS_MIN_TRAIN", unset = "10"))
OOS_START <- as.Date("1985-01-01")
COVID_QS  <- as.Date(c("2020-04-01", "2020-07-01"))

sink(file.path(OUTDIR, "run_log.txt"), split = TRUE)
cat("SPF OLS-weighted RSM analysis --", format(Sys.time()), "\n")
cat("RSM_DRAWS =", RSM_DRAWS, "| K_GRID =", paste(range(K_GRID), collapse = ":"), "\n\n")

num_safe <- function(x) suppressWarnings(as.numeric(x))

qdate <- function(year, quarter) {
  as.Date(paste(year, (quarter - 1L) * 3L + 1L, "01", sep = "-"))
}

safe_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)

ols_sum1 <- function(X, y, ridge = 1e-8) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  k <- ncol(X)
  if (k == 1L) return(1)

  A <- crossprod(X) + diag(ridge, k)
  b <- crossprod(X, y)
  w0 <- tryCatch(solve(A, b), error = function(e) rep(NA_real_, k))
  cvec <- tryCatch(solve(A, rep(1, k)), error = function(e) rep(NA_real_, k))
  if (anyNA(w0) || anyNA(cvec) || any(!is.finite(w0)) || any(!is.finite(cvec))) {
    return(rep(NA_real_, k))
  }
  denom <- sum(cvec)
  if (!is.finite(denom) || abs(denom) < 1e-12) return(rep(NA_real_, k))
  as.numeric(w0 - cvec * ((sum(w0) - 1) / denom))
}

# ---- Actuals ----------------------------------------------------------------
gdpc1 <- read.csv(file.path(DATADIR, "fred_GDPC1.csv"),
                  col.names = c("date", "gdpc1")) |>
  mutate(date = as.Date(date)) |>
  filter(!is.na(gdpc1)) |>
  arrange(date) |>
  mutate(
    yoy_gdp = 100 * (gdpc1 / lag(gdpc1, 4) - 1),
    qoq_gdp = 100 * (gdpc1 / lag(gdpc1, 1) - 1) * 4
  )

cpi_m <- read.csv(file.path(DATADIR, "fred_CPIAUCSL.csv"),
                  col.names = c("date", "cpi")) |>
  mutate(date = as.Date(date)) |>
  filter(!is.na(cpi))

cpi_q <- cpi_m |>
  mutate(qdate = as.Date(paste(year(date),
                               c(1,1,1,4,4,4,7,7,7,10,10,10)[month(date)],
                               "01", sep = "-"))) |>
  group_by(qdate) |>
  summarise(cpi = mean(cpi), .groups = "drop") |>
  arrange(qdate) |>
  mutate(yoy_cpi = 100 * (cpi / lag(cpi, 4) - 1))

# ---- SPF forecast panels -----------------------------------------------------
spf_rgdp <- read_excel(file.path(DATADIR, "spf_rgdp.xlsx"), sheet = "RGDP") |>
  mutate(across(starts_with("RGDP"), num_safe),
         survey_date = qdate(YEAR, QUARTER))

spf_cpi <- read_excel(file.path(DATADIR, "spf_cpi.xlsx"), sheet = "CPI") |>
  mutate(across(starts_with("CPI"), num_safe),
         survey_date = qdate(YEAR, QUARTER))

rgdp_h1 <- spf_rgdp |>
  filter(!is.na(RGDP1), !is.na(RGDP2), RGDP1 > 0, RGDP2 > 0) |>
  mutate(forecast = (RGDP2 / RGDP1 - 1) * 400,
         target_q = QUARTER %% 4L + 1L,
         target_y = YEAR + as.integer(QUARTER == 4L),
         target_date = qdate(target_y, target_q)) |>
  select(survey_date, ID, forecast, target_date) |>
  left_join(gdpc1 |> select(date, actual = qoq_gdp),
            by = c("target_date" = "date")) |>
  filter(!is.na(actual))

rgdp_h4 <- spf_rgdp |>
  filter(!is.na(RGDP1), !is.na(RGDP5), RGDP1 > 0, RGDP5 > 0) |>
  mutate(forecast = (RGDP5 / RGDP1 - 1) * 100,
         tq = ((QUARTER - 1L + 4L) %% 4L) + 1L,
         ty = YEAR + ((QUARTER + 3L) %/% 4L),
         target_date = qdate(ty, tq)) |>
  select(survey_date, ID, forecast, target_date) |>
  left_join(gdpc1 |> select(date, actual = yoy_gdp),
            by = c("target_date" = "date")) |>
  filter(!is.na(actual))

cpi_h1 <- spf_cpi |>
  filter(!is.na(CPI2)) |>
  mutate(forecast = num_safe(CPI2),
         target_q = QUARTER %% 4L + 1L,
         target_y = YEAR + as.integer(QUARTER == 4L),
         target_date = qdate(target_y, target_q)) |>
  select(survey_date, ID, forecast, target_date) |>
  left_join(cpi_q |> select(qdate, actual = yoy_cpi),
            by = c("target_date" = "qdate")) |>
  filter(!is.na(forecast), !is.na(actual))

cpi_h4 <- spf_cpi |>
  filter(!is.na(CPI5)) |>
  mutate(forecast = num_safe(CPI5),
         tq = ((QUARTER - 1L + 4L) %% 4L) + 1L,
         ty = YEAR + ((QUARTER + 3L) %/% 4L),
         target_date = qdate(ty, tq)) |>
  select(survey_date, ID, forecast, target_date) |>
  left_join(cpi_q |> select(qdate, actual = yoy_cpi),
            by = c("target_date" = "qdate")) |>
  filter(!is.na(forecast), !is.na(actual))

make_wide <- function(data) {
  data |>
    mutate(ID = paste0("ID_", ID)) |>
    group_by(survey_date, target_date, ID) |>
    summarise(forecast = mean(forecast, na.rm = TRUE),
              actual = first(actual),
              .groups = "drop") |>
    pivot_wider(names_from = ID, values_from = forecast) |>
    arrange(survey_date, target_date)
}

rsm_ols_one_round <- function(wide, row_idx, k, B) {
  current <- wide[row_idx, , drop = FALSE]
  id_cols <- setdiff(names(wide), c("survey_date", "target_date", "actual"))
  available <- id_cols[!is.na(as.numeric(current[1, id_cols]))]
  if (length(available) < k) {
    return(c(forecast = NA_real_, success = 0, attempts = B, avg_train_n = NA_real_))
  }

  hist <- wide |>
    filter(target_date < current$target_date[1])
  if (nrow(hist) == 0) {
    return(c(forecast = NA_real_, success = 0, attempts = B, avg_train_n = NA_real_))
  }

  preds <- numeric(0)
  train_ns <- numeric(0)
  min_train <- max(MIN_TRAIN_BASE, k + 2L)

  for (bb in seq_len(B)) {
    ids <- sample(available, k, replace = FALSE)
    sub <- hist[, c("actual", ids), drop = FALSE]
    cc <- complete.cases(sub)
    if (sum(cc) < min_train) next
    y <- sub$actual[cc]
    X <- as.matrix(sub[cc, ids, drop = FALSE])
    w <- ols_sum1(X, y)
    if (anyNA(w) || any(!is.finite(w))) next
    x_now <- as.numeric(current[1, ids])
    if (anyNA(x_now) || any(!is.finite(x_now))) next
    preds <- c(preds, sum(x_now * w))
    train_ns <- c(train_ns, sum(cc))
  }

  if (length(preds) == 0) {
    c(forecast = NA_real_, success = 0, attempts = B, avg_train_n = NA_real_)
  } else {
    c(forecast = mean(preds), success = length(preds), attempts = B,
      avg_train_n = mean(train_ns))
  }
}

run_series <- function(data, label) {
  cat("\n===", label, "===\n")
  wide_all <- make_wide(data)

  # Evaluate only post-1985, but let the OLS weights train on all earlier SPF
  # history that is complete for the selected respondents.
  id_cols <- setdiff(names(wide_all), c("survey_date", "target_date", "actual"))
  n_avail_all <- rowSums(!is.na(wide_all[, id_cols, drop = FALSE]))
  eval_idx <- which(wide_all$survey_date >= OOS_START & n_avail_all >= MIN_N)
  wide_eval <- wide_all[eval_idx, , drop = FALSE]
  n_avail <- n_avail_all[eval_idx]

  k_sqrtN <- as.integer(round(sqrt(median(n_avail))))
  k_sqrtT <- as.integer(round(sqrt(nrow(wide_eval))))
  k_eval <- sort(unique(c(K_GRID, k_sqrtN, k_sqrtT)))

  cat("Rounds:", nrow(wide_eval),
      "| N range:", min(n_avail), "-", max(n_avail),
      "| median N:", median(n_avail),
      "| k_sqrtN:", k_sqrtN,
      "| k_sqrtT:", k_sqrtT, "\n")

  mean_pred <- apply(wide_eval[, id_cols, drop = FALSE], 1, safe_mean)
  mean_fe <- mean_pred - wide_eval$actual
  mean_mafe <- mean(abs(mean_fe))
  mean_rmsfe <- sqrt(mean(mean_fe^2))

  set.seed(RSM_SEED + match(label, c("RGDP_h1", "RGDP_h4", "CPI_h1", "CPI_h4")))
  rows <- vector("list", length(k_eval))

  for (jj in seq_along(k_eval)) {
    k <- k_eval[jj]
    cat("  k =", k, "\n")
    fc <- rep(NA_real_, nrow(wide_eval))
    succ <- rep(0, nrow(wide_eval))
    avg_train <- rep(NA_real_, nrow(wide_eval))

    for (rr in seq_len(nrow(wide_eval))) {
      out <- rsm_ols_one_round(wide_all, eval_idx[rr], k, RSM_DRAWS)
      fc[rr] <- out[["forecast"]]
      succ[rr] <- out[["success"]]
      avg_train[rr] <- out[["avg_train_n"]]
    }

    ok <- !is.na(fc)
    fe <- fc[ok] - wide_eval$actual[ok]
    mean_fe_same <- mean_pred[ok] - wide_eval$actual[ok]
    covid_mask <- wide_eval$target_date %in% COVID_QS
    ok_nc <- ok & !covid_mask
    fe_nc <- fc[ok_nc] - wide_eval$actual[ok_nc]
    mean_fe_same_nc <- mean_pred[ok_nc] - wide_eval$actual[ok_nc]

    rows[[jj]] <- data.frame(
      k = k,
      n_eval = sum(ok),
      coverage = mean(ok),
      mean_success_draws = mean(succ[ok]),
      mean_success_rate = mean(succ[ok] / RSM_DRAWS),
      mean_train_n = mean(avg_train[ok], na.rm = TRUE),
      mafe = mean(abs(fe)),
      rmsfe = sqrt(mean(fe^2)),
      mean_same_mafe = mean(abs(mean_fe_same)),
      mean_same_rmsfe = sqrt(mean(mean_fe_same^2)),
      rmsfe_gain_vs_mean_same_pct =
        100 * (sqrt(mean(mean_fe_same^2)) - sqrt(mean(fe^2))) /
        sqrt(mean(mean_fe_same^2)),
      mafe_nc = mean(abs(fe_nc)),
      rmsfe_nc = sqrt(mean(fe_nc^2)),
      mean_same_mafe_nc = mean(abs(mean_fe_same_nc)),
      mean_same_rmsfe_nc = sqrt(mean(mean_fe_same_nc^2))
    )
  }

  by_k <- bind_rows(rows)
  write.csv(by_k, file.path(OUTDIR, paste0("rmsfe_by_k_ols_", label, ".csv")), row.names = FALSE)

  # Oracle only among k with at least 80% forecast coverage.
  eligible <- by_k$coverage >= 0.80
  k_star <- if (any(eligible)) by_k$k[eligible][which.min(by_k$rmsfe[eligible])] else NA_integer_
  row_sqrtN <- by_k[by_k$k == k_sqrtN, , drop = FALSE]
  row_sqrtT <- by_k[by_k$k == k_sqrtT, , drop = FALSE]
  row_star <- by_k[by_k$k == k_star, , drop = FALSE]

  tab <- data.frame(
    Model = c("Mean", paste0("OLS-RSM(k*=", k_star, ")"),
              paste0("OLS-RSM(k=sqrtN=", k_sqrtN, ")"),
              paste0("OLS-RSM(k=sqrtT=", k_sqrtT, ")")),
    k = c(NA, k_star, k_sqrtN, k_sqrtT),
    n_eval = c(nrow(wide_eval), row_star$n_eval, row_sqrtN$n_eval, row_sqrtT$n_eval),
    MAFE = c(mean_mafe, row_star$mafe, row_sqrtN$mafe, row_sqrtT$mafe),
    RMSFE = c(mean_rmsfe, row_star$rmsfe, row_sqrtN$rmsfe, row_sqrtT$rmsfe),
    comparable_mean_MAFE = c(mean_mafe, row_star$mean_same_mafe,
                             row_sqrtN$mean_same_mafe, row_sqrtT$mean_same_mafe),
    comparable_mean_RMSFE = c(mean_rmsfe, row_star$mean_same_rmsfe,
                              row_sqrtN$mean_same_rmsfe, row_sqrtT$mean_same_rmsfe),
    rmsfe_gain_vs_comparable_mean_pct =
      c(NA, row_star$rmsfe_gain_vs_mean_same_pct,
        row_sqrtN$rmsfe_gain_vs_mean_same_pct,
        row_sqrtT$rmsfe_gain_vs_mean_same_pct),
    coverage = c(1, row_star$coverage, row_sqrtN$coverage, row_sqrtT$coverage),
    success_rate = c(NA, row_star$mean_success_rate,
                     row_sqrtN$mean_success_rate, row_sqrtT$mean_success_rate)
  )
  write.csv(tab, file.path(OUTDIR, paste0("table_ols_", label, ".csv")), row.names = FALSE)
  print(tab)
  invisible(list(by_k = by_k, tab = tab))
}

all_series <- list(
  RGDP_h1 = rgdp_h1,
  RGDP_h4 = rgdp_h4,
  CPI_h1  = cpi_h1,
  CPI_h4  = cpi_h4
)
if (length(RUN_ONLY) > 0) {
  all_series <- all_series[names(all_series) %in% RUN_ONLY]
}

res <- lapply(names(all_series), function(nm) run_series(all_series[[nm]], nm))
names(res) <- names(all_series)

master <- bind_rows(lapply(names(res), function(nm) {
  cbind(series = nm, res[[nm]]$tab)
}))
write.csv(master, file.path(OUTDIR, "master_table_ols.csv"), row.names = FALSE)

cat("\nMaster table:\n")
print(master)
cat("\nDone. Outputs saved in", OUTDIR, "\n")
sink()
