#!/usr/bin/env Rscript
###############################################################################
# compare_kby_powerlaw.R
#
# Fits R_k = C * z_k^alpha via NLS and log-log OLS to both:
#   (A) In-sample Q_rs_bag excess:  R_k^Q  = Q_rs_bag(k) - Q_rs_bag(m)
#   (B) RMSFE excess left branch:   R_k^rf = RMSFE(k)^2 - RMSFE(k*)^2  [k<=k*]
#
# Runs on BOTH K_BY=3 (RSM_endogenous_results/) and K_BY=1 (RSM_endogenous_kby1/)
# and prints a side-by-side comparison table so we can verify the estimated
# alpha and C are stable across grid densities.
###############################################################################

suppressMessages(library(minpack.lm))

DIR3 <- file.path(getwd(), "RSM_endogenous_results")
DIR1 <- file.path(getwd(), "RSM_endogenous_kby1")

# ---------------------------------------------------------------------------
# Helper: fit C * z^alpha via NLS (primary) and log-log OLS (secondary)
# ---------------------------------------------------------------------------
fit_pw <- function(k_grid, R_k, m, label = "") {
  z    <- 1 / k_grid - 1 / m
  keep <- z > 0 & R_k > 0 & is.finite(R_k) & is.finite(z) & k_grid > 1
  if (sum(keep) < 4) {
    message(sprintf("  [%s] only %d usable pts — skip", label, sum(keep)))
    return(NULL)
  }
  zf <- z[keep];  Rf <- R_k[keep]
  lfit <- lm(log(Rf) ~ log(zf))
  C0   <- unname(exp(coef(lfit)[1]))
  a0   <- unname(coef(lfit)[2])
  df   <- data.frame(zf = zf, Rf = Rf)

  nfit <- tryCatch(
    nlsLM(Rf ~ Cpow * zf^alph, data = df,
          start = list(Cpow = C0, alph = a0),
          lower = c(Cpow = 1e-10, alph = 0.05),
          upper = c(Cpow = 1e4,   alph = 8),
          control = nls.lm.control(maxiter = 300, ftol = 1e-10, ptol = 1e-10)),
    error = function(e) { message("  NLS [", label, "]: ", e$message); NULL }
  )
  if (is.null(nfit)) return(NULL)

  Ch <- unname(coef(nfit)["Cpow"])
  ah <- unname(coef(nfit)["alph"])
  pr <- Ch * zf^ah
  ss_r <- sum((Rf - pr)^2);  ss_t <- sum((Rf - mean(Rf))^2)
  pr0  <- C0 * zf^a0
  list(n = sum(keep),
       C_nls    = Ch,   alpha_nls = ah,
       R2_nls   = 1 - ss_r / ss_t,
       C_ols    = C0,   alpha_ols = a0,
       R2_ols_log = summary(lfit)$r.squared,
       R2_ols_lin = 1 - sum((Rf - pr0)^2) / ss_t)
}

# ---------------------------------------------------------------------------
# Process one directory of by_z CSVs
# ---------------------------------------------------------------------------
process_dir <- function(dir, tag) {
  files <- sort(list.files(dir, pattern = "by_z\\.csv$", full.names = TRUE))
  if (length(files) == 0) {
    cat(sprintf("[%s] No by_z files found in %s\n", tag, dir))
    return(NULL)
  }
  cat(sprintf("\n[%s] Found %d by_z files\n", tag, length(files)))

  rows <- list()
  for (f in files) {
    dat   <- read.csv(f)
    m_val <- if ("m" %in% names(dat)) dat$m[1] else max(dat$z)
    T_val <- dat$T_len[1]
    stype <- dat$sigma_type[1]
    bdes  <- if ("beta_design" %in% names(dat)) dat$beta_design[1] else "B4"
    lab   <- paste0(stype, "_T", T_val, "_", bdes)

    n_k <- nrow(dat)
    cat(sprintf("  %-42s  n_k=%3d  k_step=%.1f\n", lab, n_k,
                if (n_k > 1) diff(sort(dat$z)[1:2]) else NA))

    # (A) Q excess
    Q_floor <- dat$Q_rs_bag[dat$z == m_val]
    if (length(Q_floor) == 0) Q_floor <- min(dat$Q_rs_bag, na.rm = TRUE)
    R_Q <- dat$Q_rs_bag - Q_floor

    # (B) RMSFE left branch
    k_star    <- dat$z[which.min(dat$rmsfe_rs)]
    rmsfe_min <- min(dat$rmsfe_rs, na.rm = TRUE)
    R_mspe    <- dat$rmsfe_rs^2 - rmsfe_min^2
    left      <- dat$z <= k_star

    fQ <- fit_pw(dat$z,       R_Q,          m_val, paste(lab, "Q"))
    fR <- fit_pw(dat$z[left], R_mspe[left], m_val, paste(lab, "RMSFE"))

    rows[[length(rows) + 1]] <- data.frame(
      tag = tag, label = lab, sigma = stype, T_len = T_val, beta = bdes,
      m = m_val, k_star = k_star, n_k = n_k,
      # Q fit
      n_Q         = if (!is.null(fQ)) fQ$n           else NA_integer_,
      C_Q_nls     = if (!is.null(fQ)) round(fQ$C_nls,   4) else NA_real_,
      a_Q_nls     = if (!is.null(fQ)) round(fQ$alpha_nls,3) else NA_real_,
      R2_Q_nls    = if (!is.null(fQ)) round(fQ$R2_nls,  4) else NA_real_,
      a_Q_ols     = if (!is.null(fQ)) round(fQ$alpha_ols,3) else NA_real_,
      R2_Q_log    = if (!is.null(fQ)) round(fQ$R2_ols_log,4) else NA_real_,
      # RMSFE fit
      n_R         = if (!is.null(fR)) fR$n           else NA_integer_,
      a_R_nls     = if (!is.null(fR)) round(fR$alpha_nls,3) else NA_real_,
      R2_R_nls    = if (!is.null(fR)) round(fR$R2_nls,  4) else NA_real_,
      a_R_ols     = if (!is.null(fR)) round(fR$alpha_ols,3) else NA_real_,
      R2_R_log    = if (!is.null(fR)) round(fR$R2_ols_log,4) else NA_real_,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

# ---------------------------------------------------------------------------
# Run both directories
# ---------------------------------------------------------------------------
df3 <- process_dir(DIR3, "K_BY=3")
df1 <- process_dir(DIR1, "K_BY=1")

all_df <- rbind(df3, df1)

# ---------------------------------------------------------------------------
# Side-by-side comparison table
# Print for each scenario that appears in BOTH grids
# ---------------------------------------------------------------------------
cat("\n\n")
cat(strrep("=", 100), "\n")
cat("COMPARISON: K_BY=3 vs K_BY=1  —  same scenario, same sigma type\n")
cat(strrep("=", 100), "\n")

## Deduplicate: if same label appears multiple times (main + scaling), keep highest n_k
dedup <- function(df) {
  do.call(rbind, lapply(split(df, df$label), function(g) g[which.max(g$n_k), ]))
}
df3d <- if (!is.null(df3)) dedup(df3) else NULL
df1d <- if (!is.null(df1)) dedup(df1) else NULL

common_labs <- intersect(df3d$label, df1d$label)
if (length(common_labs) == 0) {
  cat("No matching labels between K_BY=3 and K_BY=1 yet.\n")
} else {
  cat(sprintf("\n%-42s | %-5s | %4s %5s %6s | %4s %5s %6s | %s\n",
    "Scenario", "GRD", "n_Q", "a_Q", "R2_Q", "n_R", "a_R", "R2_R", "note"))
  cat(strrep("-", 100), "\n")

  for (lab in sort(common_labs)) {
    r3 <- df3d[df3d$label == lab, ]
    r1 <- df1d[df1d$label == lab, ]
    q_match  <- !is.na(r3$a_Q_nls) && !is.na(r1$a_Q_nls) &&
                abs(r3$a_Q_nls - r1$a_Q_nls) < 0.05
    rf_match <- !is.na(r3$a_R_nls) && !is.na(r1$a_R_nls) &&
                abs(r3$a_R_nls - r1$a_R_nls) < 0.05

    for (rx in list(r3, r1)) {
      note <- if (rx$tag == "K_BY=1") {
                paste0(if (q_match)  "Q-OK "   else "Q-DIFF ",
                       if (rf_match) "RF-OK"   else "RF-DIFF")
              } else ""
      cat(sprintf("%-42s | %-5s | %4s %5s %6s | %4s %5s %6s | %s\n",
        lab, rx$tag,
        if (is.na(rx$n_Q)) "--" else rx$n_Q,
        if (is.na(rx$a_Q_nls)) "--" else sprintf("%.3f", rx$a_Q_nls),
        if (is.na(rx$R2_Q_nls)) "--" else sprintf("%.4f", rx$R2_Q_nls),
        if (is.na(rx$n_R)) "--" else rx$n_R,
        if (is.na(rx$a_R_nls)) "--" else sprintf("%.3f", rx$a_R_nls),
        if (is.na(rx$R2_R_nls)) "--" else sprintf("%.4f", rx$R2_R_nls),
        note))
    }
    cat(strrep("-", 100), "\n")
  }
}

# ---------------------------------------------------------------------------
# Also print full tables for each grid separately
# ---------------------------------------------------------------------------
print_full <- function(df, tag) {
  cat(sprintf("\n\n=== Full results: %s ===\n", tag))
  cat(sprintf("%-42s %5s %4s %4s | %6s %5s %6s | %4s %5s %6s\n",
    "Label","T","k*","n_k","C_Q","a_Q","R2_Q","n_R","a_R","R2_R"))
  cat(strrep("-", 95), "\n")
  for (i in seq_len(nrow(df))) {
    r <- df[i, ]
    cat(sprintf("%-42s %5d %4d %4d | %6s %5s %6s | %4s %5s %6s\n",
      r$label, r$T_len, r$k_star, r$n_k,
      if (is.na(r$C_Q_nls)) "  --" else sprintf("%.3f", r$C_Q_nls),
      if (is.na(r$a_Q_nls)) "  --" else sprintf("%.3f", r$a_Q_nls),
      if (is.na(r$R2_Q_nls))"  --" else sprintf("%.4f", r$R2_Q_nls),
      if (is.na(r$n_R)) "--"         else r$n_R,
      if (is.na(r$a_R_nls)) "  --"  else sprintf("%.3f", r$a_R_nls),
      if (is.na(r$R2_R_nls))"  --"  else sprintf("%.4f", r$R2_R_nls)))
  }
}

if (!is.null(df3)) print_full(df3d, "K_BY=3 (dedup, highest n_k)")
if (!is.null(df1)) print_full(df1d, "K_BY=1 (dedup, highest n_k)")

# Save combined table
out_csv <- file.path(getwd(), "RSM_endogenous_kby1", "kby_comparison.csv")
write.csv(all_df, out_csv, row.names = FALSE)
cat(sprintf("\n\nSaved combined table to %s\n", out_csv))
cat("Done.\n")
