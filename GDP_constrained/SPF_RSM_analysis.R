#!/usr/bin/env Rscript
# =============================================================================
# SPF_RSM_analysis.R
#
# Applies RSM to Philadelphia Fed Survey of Professional Forecasters.
# Replicates Conflitti, De Mol & Giannone (2015) setup using US SPF:
#   - Real GDP growth h=1 (one quarter ahead, annualised)
#   - Real GDP growth h=4 (four quarters ahead, year-on-year)
#   - CPI inflation    h=4 (four quarters ahead, year-on-year)
#
# RSM vs equal-weighted mean; missing data handled by using only
# available respondents at each round (no imputation needed for RSM).
#
# Outputs:  SPF_results/
# =============================================================================

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(lubridate)
})

setwd("/Users/boris/Desktop/RSShapley/Combination_Monte/Final/codess/JMP")

DATADIR <- "GDP_constrained"
OUTDIR  <- "GDP_constrained/SPF_results"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

RSM_SEED  <- 20250604L
RSM_DRAWS     <- 1000L      # subsets per k per round
MIN_N         <- 5L         # minimum respondents to include a round
MIN_N_K_COV   <- 0.80       # min fraction of rounds that must have N_t >= k
OOS_START     <- as.Date("1985-01-01")  # first OOS quarter
COVID_QS  <- as.Date(c("2020-04-01","2020-07-01"))  # exclude these quarters

sink(file.path(OUTDIR, "run_log.txt"), split = TRUE)
cat("SPF RSM Analysis —", format(Sys.time()), "\n\n")

# =============================================================================
# 1.  Load FRED actuals
# =============================================================================
gdpc1 <- read.csv(file.path(DATADIR, "fred_GDPC1.csv"),
                  col.names = c("date","gdpc1")) |>
  mutate(date = as.Date(date)) |>
  filter(!is.na(gdpc1))

cpi_m <- read.csv(file.path(DATADIR, "fred_CPIAUCSL.csv"),
                  col.names = c("date","cpi")) |>
  mutate(date = as.Date(date)) |>
  filter(!is.na(cpi))

# Aggregate monthly CPI to quarterly (average of 3 months)
cpi_q <- cpi_m |>
  mutate(qdate = as.Date(paste(year(date),
                               c(1,1,1,4,4,4,7,7,7,10,10,10)[month(date)],
                               "01", sep="-"))) |>
  group_by(qdate) |>
  summarise(cpi = mean(cpi), .groups = "drop")

# Helper: quarter-start date from YEAR + QUARTER
qdate <- function(year, quarter) {
  as.Date(paste(year, (quarter-1L)*3L+1L, "01", sep="-"))
}

# GDPC1 year-on-year growth at quarter d: 100*(GDPC1_d / GDPC1_{d-4q} - 1)
gdpc1 <- gdpc1 |> arrange(date) |>
  mutate(yoy_gdp = 100 * (gdpc1 / lag(gdpc1, 4) - 1),
         qoq_gdp = 100 * (gdpc1 / lag(gdpc1, 1) - 1) * 4)  # annualised

# =============================================================================
# 2.  Load SPF individual responses
# =============================================================================
num_safe <- function(x) suppressWarnings(as.numeric(x))

spf_rgdp <- read_excel(file.path(DATADIR, "spf_rgdp.xlsx"), sheet = "RGDP") |>
  mutate(across(starts_with("RGDP"), num_safe),
         survey_date = qdate(YEAR, QUARTER))

spf_cpi <- read_excel(file.path(DATADIR, "spf_cpi.xlsx"), sheet = "CPI") |>
  mutate(across(starts_with("CPI"), num_safe),
         survey_date = qdate(YEAR, QUARTER))

cat("SPF RGDP: rows =", nrow(spf_rgdp),
    "  years:", range(spf_rgdp$YEAR), "\n")
cat("SPF CPI:  rows =", nrow(spf_cpi),
    "  years:", range(spf_cpi$YEAR, na.rm=TRUE), "\n\n")

# =============================================================================
# 3.  Construct individual forecasts per round
# =============================================================================

# ---- RGDP h=1 (QoQ annualised) ---------------------------------------------
#   Survey (Y, q): RGDP2/RGDP1 → forecast for Q(q+1)
#   Actual: QoQ annualised growth of GDPC1 at Q(q+1)
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

# ---- RGDP h=4 (year-on-year) -----------------------------------------------
#   Survey (Y, q): RGDP5/RGDP1 → 4-quarter growth ending at Q(q+4)
#   Actual: YoY GDPC1 growth at Q(q+4)
rgdp_h4 <- spf_rgdp |>
  filter(!is.na(RGDP1), !is.na(RGDP5), RGDP1 > 0, RGDP5 > 0) |>
  mutate(forecast = (RGDP5 / RGDP1 - 1) * 100,
         tq = ((QUARTER - 1L + 4L) %% 4L) + 1L,
         ty = YEAR + as.integer(QUARTER + 4L > 4L &
                                 (QUARTER + 4L) %/% 5L > 0L),
         ty = YEAR + ((QUARTER + 3L) %/% 4L),
         target_date = qdate(ty, tq)) |>
  select(survey_date, ID, forecast, target_date) |>
  left_join(gdpc1 |> select(date, actual = yoy_gdp),
            by = c("target_date" = "date")) |>
  filter(!is.na(actual))

# ---- CPI h=1 and h=4 -------------------------------------------------------
# IMPORTANT: SPF CPI columns (CPI1..CPI5) are ALREADY percentage rates
# (year-on-year % change), NOT index levels.  Do NOT divide CPI5/CPI1.
# CPI1 = current-quarter YoY inflation (observed/nowcast)
# CPI2..CPI5 = forecasted YoY inflation for horizons h=1..4
# Actual: CPIAUCSL quarterly average, YoY percent change from FRED.

cpi_q <- cpi_q |> arrange(qdate) |>
  mutate(yoy_cpi = 100 * (cpi / lag(cpi, 4) - 1))

# h=1: CPI2 = YoY CPI forecast one quarter ahead
cpi_h1 <- spf_cpi |>
  filter(!is.na(CPI1), !is.na(CPI2)) |>
  mutate(CPI2n = num_safe(CPI2)) |>
  filter(!is.na(CPI2n)) |>
  mutate(forecast    = CPI2n,
         target_q    = QUARTER %% 4L + 1L,
         target_y    = YEAR + as.integer(QUARTER == 4L),
         target_date = qdate(target_y, target_q)) |>
  select(survey_date, ID, forecast, target_date) |>
  left_join(cpi_q |> select(qdate, actual = yoy_cpi),
            by = c("target_date" = "qdate")) |>
  filter(!is.na(actual))

# h=4: CPI5 = YoY CPI forecast four quarters ahead
cpi_h4 <- spf_cpi |>
  filter(!is.na(CPI1), !is.na(CPI5)) |>
  mutate(CPI5n = num_safe(CPI5)) |>
  filter(!is.na(CPI5n)) |>
  mutate(forecast    = CPI5n,
         tq          = ((QUARTER - 1L + 4L) %% 4L) + 1L,
         ty          = YEAR + ((QUARTER + 3L) %/% 4L),
         target_date = qdate(ty, tq)) |>
  select(survey_date, ID, forecast, target_date) |>
  left_join(cpi_q |> select(qdate, actual = yoy_cpi),
            by = c("target_date" = "qdate")) |>
  filter(!is.na(actual))

cat("Usable observations:\n")
cat(sprintf("  RGDP h=1: %d forecaster-rounds, %d survey rounds\n",
            nrow(rgdp_h1), n_distinct(rgdp_h1$survey_date)))
cat(sprintf("  RGDP h=4: %d forecaster-rounds, %d survey rounds\n",
            nrow(rgdp_h4), n_distinct(rgdp_h4$survey_date)))
cat(sprintf("  CPI  h=1: %d forecaster-rounds, %d survey rounds\n",
            nrow(cpi_h1), n_distinct(cpi_h1$survey_date)))
cat(sprintf("  CPI  h=4: %d forecaster-rounds, %d survey rounds\n\n",
            nrow(cpi_h4), n_distinct(cpi_h4$survey_date)))

# =============================================================================
# 4.  RSM engine: for one survey round, one k, return B subset-mean forecasts
# =============================================================================
rsm_round <- function(forecasts, k, B) {
  n <- length(forecasts)
  if (k > n) return(NA_real_)
  draws <- replicate(B, mean(sample(forecasts, k, replace = FALSE)))
  mean(draws)
}

# =============================================================================
# 5.  Full OOS evaluation: RSM(k) and Mean across all rounds
# =============================================================================
run_spf_rsm <- function(data, label, k_max_global = 55L) {
  cat("\n===", label, "===\n")

  rounds <- data |>
    filter(survey_date >= OOS_START) |>
    group_by(survey_date, target_date) |>
    summarise(n      = n(),
              mean_f = mean(forecast, na.rm = TRUE),
              actual = first(actual),
              fvec   = list(forecast),
              .groups = "drop") |>
    filter(n >= MIN_N) |>
    arrange(survey_date)

  n_rounds <- nrow(rounds)
  cat(sprintf("  OOS rounds: %d  |  N range: %d-%d  |  median N: %.0f\n",
              n_rounds, min(rounds$n), max(rounds$n), median(rounds$n)))

  # coverage(k) = fraction of rounds with N_t >= k
  # k_max_cov = largest k with coverage >= MIN_N_K_COV (unbiased RMSFE range)
  k_all <- seq(2L, min(max(rounds$n), k_max_global))
  coverage <- sapply(k_all, function(k) mean(rounds$n >= k))
  k_max_cov <- max(k_all[coverage >= MIN_N_K_COV], na.rm = TRUE)
  cat(sprintf("  k_max at >=%.0f%% coverage: %d\n",
              MIN_N_K_COV*100, k_max_cov))

  # ---- per-round errors for mean and RSM(k) --------------------------------
  set.seed(RSM_SEED)
  # Use full k_max_global for the curve but mark unbiased region
  k_grid <- seq(2L, min(max(rounds$n), k_max_global))

  # Matrix: rows = rounds, cols = k values
  rsm_fe <- matrix(NA_real_, nrow(rounds), length(k_grid))
  colnames(rsm_fe) <- paste0("k", k_grid)

  for (i in seq_len(nrow(rounds))) {
    fv <- rounds$fvec[[i]]
    act <- rounds$actual[i]
    ni  <- rounds$n[i]
    for (j in seq_along(k_grid)) {
      kk <- k_grid[j]
      if (kk > ni) next
      rsm_fe[i, j] <- rsm_round(fv, kk, RSM_DRAWS) - act
    }
    if (i %% 20 == 0) cat(sprintf("    round %d / %d\n", i, nrow(rounds)))
  }

  mean_fe <- rounds$mean_f - rounds$actual
  covid_mask <- rounds$target_date %in% COVID_QS

  # ---- summarise RMSFE and MAFE by k ---------------------------------------
  rmsfe_k <- apply(rsm_fe, 2, function(e) sqrt(mean(e^2, na.rm=TRUE)))
  mafe_k  <- apply(rsm_fe, 2, function(e) mean(abs(e), na.rm=TRUE))
  rmsfe_k_nc <- apply(rsm_fe[!covid_mask,,drop=FALSE], 2,
                      function(e) sqrt(mean(e^2, na.rm=TRUE)))
  mafe_k_nc  <- apply(rsm_fe[!covid_mask,,drop=FALSE], 2,
                      function(e) mean(abs(e), na.rm=TRUE))

  rmsfe_mean    <- sqrt(mean(mean_fe^2))
  mafe_mean     <- mean(abs(mean_fe))
  rmsfe_mean_nc <- sqrt(mean(mean_fe[!covid_mask]^2))
  mafe_mean_nc  <- mean(abs(mean_fe[!covid_mask]))

  # k* restricted to coverage-unbiased region (all rounds contribute)
  cov_ok      <- k_grid <= k_max_cov
  kstar_rmsfe    <- k_grid[which.min(ifelse(cov_ok, rmsfe_k, NA))]
  kstar_rmsfe_nc <- k_grid[which.min(ifelse(cov_ok, rmsfe_k_nc, NA))]
  kstar_mafe     <- k_grid[which.min(ifelse(cov_ok, mafe_k, NA))]

  avg_n  <- round(median(rounds$n))
  k_sqrt <- round(sqrt(avg_n))

  cat(sprintf("  Median N = %d  =>  k=sqrt(N) ~ %d\n", avg_n, k_sqrt))
  cat(sprintf("  k*(RMSFE) full = %d  |  k*(RMSFE) ex-COVID = %d\n",
              kstar_rmsfe, kstar_rmsfe_nc))

  # DM-style t-test: RSM(k*) vs Mean (simple, no HAC)
  rms_best_fe <- rsm_fe[, paste0("k", kstar_rmsfe)]
  dm_ok <- !is.na(rms_best_fe)
  loss_diff_sq  <- mean_fe[dm_ok]^2 - rms_best_fe[dm_ok]^2
  loss_diff_abs <- abs(mean_fe[dm_ok]) - abs(rms_best_fe[dm_ok])
  t_sq  <- mean(loss_diff_sq)  / (sd(loss_diff_sq)  / sqrt(sum(dm_ok)))
  t_abs <- mean(loss_diff_abs) / (sd(loss_diff_abs) / sqrt(sum(dm_ok)))
  sig_sq  <- ifelse(abs(t_sq) > 2.58, "**", ifelse(abs(t_sq) > 1.96, "*", ""))
  sig_abs <- ifelse(abs(t_abs) > 2.58, "**", ifelse(abs(t_abs) > 1.96, "*", ""))

  # k=sqrt(N) performance
  rmsfe_sqrtN <- if (paste0("k",k_sqrt) %in% names(rmsfe_k)) rmsfe_k[paste0("k",k_sqrt)] else NA
  mafe_sqrtN  <- if (paste0("k",k_sqrt) %in% names(mafe_k))  mafe_k[paste0("k",k_sqrt)]  else NA

  # ---- Table printout -------------------------------------------------------
  cat("\n  === Full sample ===\n")
  cat(sprintf("  %-22s  MAFE=%.4f  RMSFE=%.4f\n","Mean",mafe_mean,rmsfe_mean))
  cat(sprintf("  %-22s  MAFE=%.4f%s  RMSFE=%.4f%s   (k*=%d)\n",
              "RSM(k*)",
              mafe_k[paste0("k",kstar_mafe)], sig_abs,
              rmsfe_k[paste0("k",kstar_rmsfe)], sig_sq,
              kstar_rmsfe))
  cat(sprintf("  %-22s  MAFE=%.4f  RMSFE=%.4f   (k=sqrt(N)~%d)\n",
              "RSM(sqrt(N))", mafe_sqrtN, rmsfe_sqrtN, k_sqrt))

  cat("\n  === Excluding COVID-19 ===\n")
  cat(sprintf("  %-22s  MAFE=%.4f  RMSFE=%.4f\n","Mean",mafe_mean_nc,rmsfe_mean_nc))
  cat(sprintf("  %-22s  MAFE=%.4f  RMSFE=%.4f   (k*=%d)\n",
              "RSM(k*_nc)", mafe_k_nc[paste0("k",kstar_rmsfe_nc)],
              rmsfe_k_nc[paste0("k",kstar_rmsfe_nc)], kstar_rmsfe_nc))

  # ---- Save CSV summaries ---------------------------------------------------
  res_by_k <- data.frame(
    k        = k_grid,
    rmsfe    = rmsfe_k,
    mafe     = mafe_k,
    rmsfe_nc = rmsfe_k_nc,
    mafe_nc  = mafe_k_nc
  )
  fname_csv <- file.path(OUTDIR, paste0("rmsfe_by_k_", label, ".csv"))
  write.csv(res_by_k, fname_csv, row.names = FALSE)

  tab <- data.frame(
    Model     = c("Mean", sprintf("RSM(k*=%d)",kstar_rmsfe),
                  sprintf("RSM(k=sqrt(N)~%d)",k_sqrt)),
    MAFE      = c(mafe_mean,  mafe_k[paste0("k",kstar_mafe)],  mafe_sqrtN),
    MAFE_sig  = c("", sig_abs, ""),
    RMSFE     = c(rmsfe_mean, rmsfe_k[paste0("k",kstar_rmsfe)], rmsfe_sqrtN),
    RMSFE_sig = c("", sig_sq, ""),
    MAFE_nc   = c(mafe_mean_nc,  mafe_k_nc[paste0("k",kstar_rmsfe_nc)], NA),
    RMSFE_nc  = c(rmsfe_mean_nc, rmsfe_k_nc[paste0("k",kstar_rmsfe_nc)], NA)
  )
  write.csv(tab, file.path(OUTDIR, paste0("table_", label, ".csv")),
            row.names = FALSE)

  # ---- RMSFE-by-k figure ----------------------------------------------------
  png(file.path(OUTDIR, paste0("fig_rmsfe_by_k_", label, ".png")),
      width = 820, height = 500, res = 110)
  par(mar=c(4.2,4.2,3,1.5), mgp=c(2.5,0.7,0))

  valid_k <- !is.na(rmsfe_k)
  yk <- rmsfe_k[valid_k]; xk <- k_grid[valid_k]
  yk_nc <- rmsfe_k_nc[valid_k]
  ylim <- range(c(yk[k_grid[valid_k] <= k_max_cov],
                  yk_nc[k_grid[valid_k] <= k_max_cov],
                  rmsfe_mean, rmsfe_mean_nc), na.rm=TRUE)
  ylim <- ylim + c(-0.02, 0.08)*diff(ylim)

  plot(xk, yk, type="n",
       xlab="Subset size k  (forecasters drawn per subset)",
       ylab="RMSFE",
       main=sprintf("RSM vs. SPF Equal-Weighted Mean — %s", label),
       ylim=ylim, las=1)
  # shade selection-biased region (k > k_max_cov)
  if (k_max_cov < max(xk)) {
    rect(k_max_cov, ylim[1], max(xk)+1, ylim[2],
         col=adjustcolor("grey80",0.4), border=NA)
    text(k_max_cov + (max(xk)-k_max_cov)/2, ylim[2]*0.98,
         sprintf("<%.0f%% coverage", MIN_N_K_COV*100),
         cex=0.65, col="grey50", adj=c(0.5,1))
  }
  lines(xk, yk,    lwd=2.2, col="#2166AC")
  lines(xk, yk_nc, lwd=1.8, col="#4DAC26", lty=2)
  abline(h=rmsfe_mean,    lty=1, lwd=1.4, col="#D6604D")
  abline(h=rmsfe_mean_nc, lty=2, lwd=1.2, col="#F4A582")
  abline(v=kstar_rmsfe,    lty=3, lwd=1.2, col="#2166AC")
  abline(v=kstar_rmsfe_nc, lty=3, lwd=1.2, col="#4DAC26")
  abline(v=k_sqrt,          lty=4, lwd=1.5, col="grey40")
  abline(v=k_max_cov,       lty=1, lwd=1.0, col="grey60")
  points(kstar_rmsfe,    rmsfe_k[paste0("k",kstar_rmsfe)],   pch=19, cex=1.4, col="#2166AC")
  points(kstar_rmsfe_nc, rmsfe_k_nc[paste0("k",kstar_rmsfe_nc)], pch=17, cex=1.4, col="#4DAC26")

  legend("topright", bty="n", cex=0.78,
    legend=c(
      sprintf("RSM full  (k*=%d, RMSFE=%.3f)", kstar_rmsfe,
              rmsfe_k[paste0("k",kstar_rmsfe)]),
      sprintf("RSM ex-COVID (k*=%d, RMSFE=%.3f)", kstar_rmsfe_nc,
              rmsfe_k_nc[paste0("k",kstar_rmsfe_nc)]),
      sprintf("Mean full  (RMSFE=%.3f)", rmsfe_mean),
      sprintf("Mean ex-COVID (RMSFE=%.3f)", rmsfe_mean_nc),
      sprintf("k=sqrt(N)~%d", k_sqrt),
      sprintf("80%% coverage boundary (k=%d)", k_max_cov)
    ),
    col=c("#2166AC","#4DAC26","#D6604D","#F4A582","grey40","grey60"),
    lty=c(1,2,1,2,4,1), lwd=c(2.2,1.8,1.4,1.2,1.5,1.0),
    pch=c(19,17,NA,NA,NA,NA))
  dev.off()
  cat(sprintf("\n  Saved figure: fig_rmsfe_by_k_%s.png\n", label))

  invisible(list(rounds=rounds, res_by_k=res_by_k, tab=tab,
                 kstar=kstar_rmsfe, k_sqrt=k_sqrt,
                 rmsfe_mean=rmsfe_mean, rmsfe_mean_nc=rmsfe_mean_nc))
}

# =============================================================================
# 6.  Run all three targets
# =============================================================================
res_gdp_h1 <- run_spf_rsm(rgdp_h1, "RGDP_h1", k_max_global = 55L)
res_gdp_h4 <- run_spf_rsm(rgdp_h4, "RGDP_h4", k_max_global = 55L)
res_cpi_h1 <- run_spf_rsm(cpi_h1,  "CPI_h1",  k_max_global = 55L)
res_cpi_h4 <- run_spf_rsm(cpi_h4,  "CPI_h4",  k_max_global = 55L)

# =============================================================================
# 7.  Combined two-panel figure (h=4 GDP + h=4 CPI, matching Giannone layout)
# =============================================================================
png(file.path(OUTDIR, "fig_spf_combined.png"),
    width = 1300, height = 520, res = 110)
par(mfrow=c(1,2), mar=c(4.2,4.2,3,1.5), mgp=c(2.5,0.7,0), cex.axis=0.85)

for (res in list(res_gdp_h4, res_cpi_h4)) {

  d  <- res$res_by_k
  valid <- !is.na(d$rmsfe)
  xk <- d$k[valid]; yk <- d$rmsfe[valid]; yk_nc <- d$rmsfe_nc[valid]
  ks  <- res$kstar; km <- res$k_sqrt
  rm  <- res$rmsfe_mean; rm_nc <- res$rmsfe_mean_nc
  lbl <- if (identical(res, res_gdp_h4)) "GDP h=4 (year-on-year)" else "CPI h=4 (year-on-year)"

  ylim <- range(c(yk,yk_nc,rm,rm_nc),na.rm=TRUE)
  ylim <- ylim + c(-0.03,0.07)*diff(ylim)
  plot(xk, yk, type="l", lwd=2.2, col="#2166AC",
       xlab="k", ylab="RMSFE", main=lbl, ylim=ylim, las=1)
  lines(xk, yk_nc, lwd=1.8, col="#4DAC26", lty=2)
  abline(h=rm,    lty=1, lwd=1.4, col="#D6604D")
  abline(h=rm_nc, lty=2, lwd=1.2, col="#F4A582")
  abline(v=ks, lty=3, lwd=1.2, col="#2166AC")
  abline(v=km, lty=4, lwd=1.5, col="grey40")
  points(ks, min(yk,na.rm=TRUE), pch=19, cex=1.4, col="#2166AC")
  legend("topright", bty="n", cex=0.75,
    legend=c(sprintf("RSM (k*=%d)",ks), "Mean", sprintf("k=sqrt(N)~%d",km)),
    col=c("#2166AC","#D6604D","grey40"),
    lty=c(1,1,4), lwd=c(2.2,1.4,1.5), pch=c(19,NA,NA))
}
dev.off()
cat("\nSaved: fig_spf_combined.png\n")

# =============================================================================
# 8.  Master summary table (Giannone-style)
# =============================================================================
cat("\n\n=== MASTER TABLE (SPF RSM vs Mean) ===\n")
cat(sprintf("%-28s  %8s  %8s  %8s  %8s\n",
            "Series","MAFE_f","RMSFE_f","MAFE_nc","RMSFE_nc"))
for (res in list(res_gdp_h1, res_gdp_h4, res_cpi_h1, res_cpi_h4)) {
  t <- res$tab
  for (i in seq_len(nrow(t))) {
    cat(sprintf("%-28s  %8.4f%s  %8.4f%s  %8s  %8s\n",
                t$Model[i],
                t$MAFE[i], t$MAFE_sig[i],
                t$RMSFE[i], t$RMSFE_sig[i],
                ifelse(is.na(t$MAFE_nc[i]),"",sprintf("%.4f",t$MAFE_nc[i])),
                ifelse(is.na(t$RMSFE_nc[i]),"",sprintf("%.4f",t$RMSFE_nc[i]))))
  }
  cat("\n")
}

cat("\nDone:", format(Sys.time()), "\n")
sink()
cat("All outputs saved to:", OUTDIR, "\n")
