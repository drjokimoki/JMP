#!/usr/bin/env Rscript
# Compare L_rs_bag vs rmsfe_rs^2 in level and k* location across all K_BY=1 scenarios.

files <- sort(list.files("RSM_endogenous_kby1", pattern = "by_z\\.csv$", full.names = TRUE))
cat(sprintf("%-38s  %4s | %5s %5s %4s | %6s %6s %5s | kL  kRF  diff\n",
  "Scenario", "T", "corr", "ratio", "n_k", "L_min", "RF2min", "L/RF2"))
cat(strrep("-", 105), "\n")

rows <- list()
for (f in files) {
  d   <- read.csv(f)
  lab <- gsub("RSM_endogenous_kby1/", "", f, fixed = TRUE)
  lab <- gsub("endo_|scaling_", "", lab)
  lab <- gsub("_phif.*|_by_z\\.csv", "", lab)
  lab <- substr(lab, 1, 38)

  T_  <- d$T_len[1]
  L   <- d$L_rs_bag
  RF2 <- d$rmsfe_rs^2
  kL  <- d$z[which.min(L)]
  kRF <- d$z[which.min(RF2)]

  corr   <- round(cor(L, RF2, use = "complete.obs"), 4)
  Lmin   <- round(min(L,   na.rm = TRUE), 5)
  RF2min <- round(min(RF2, na.rm = TRUE), 5)
  ratio  <- round(Lmin / RF2min, 4)

  cat(sprintf("%-38s  %4d | %5.3f %5.3f %4d | %6.4f %6.4f %5.3f | %3d  %3d  %+4d\n",
    lab, T_, corr, ratio, nrow(d), Lmin, RF2min, ratio, kL, kRF, kL - kRF))

  rows[[length(rows)+1]] <- data.frame(
    label = lab, T_len = T_, k_L = kL, k_RF = kRF, diff_k = kL - kRF,
    corr = corr, L_min = Lmin, RF2_min = RF2min, ratio = ratio,
    stringsAsFactors = FALSE
  )
}

df <- do.call(rbind, rows)
cat("\n--- Summary ---\n")
cat(sprintf("Mean |k*(L) - k*(RMSFE)|: %.2f\n", mean(abs(df$diff_k))))
cat(sprintf("k*(L) > k*(RMSFE) in %d / %d scenarios (L predicts too large k*)\n",
  sum(df$diff_k > 0), nrow(df)))
cat(sprintf("L/RMSFE² ratio range: [%.3f, %.3f]  (1.0 = perfect level match)\n",
  min(df$ratio), max(df$ratio)))
cat(sprintf("Mean correlation L vs RMSFE²: %.4f\n", mean(df$corr)))
write.csv(df, "RSM_endogenous_kby1/L_vs_MSE_kstar.csv", row.names = FALSE)
cat("Saved to RSM_endogenous_kby1/L_vs_MSE_kstar.csv\n")
