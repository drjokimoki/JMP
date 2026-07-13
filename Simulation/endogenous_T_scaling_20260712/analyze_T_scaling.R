#!/usr/bin/env Rscript

root <- normalizePath(dirname(sub("^--file=", "", grep(
  "^--file=", commandArgs(FALSE), value = TRUE)[1])), mustWork = TRUE)
d <- read.csv(file.path(root, "rs_results.csv"))
d$msfe_rs <- d$rmse_pooled_rs^2
d$rel_pooled_rmse <- d$rmse_pooled_rs / d$rmse_pooled_base

Ts <- sort(unique(d$T_len))
summary_rows <- lapply(Ts, function(T0) {
  x <- d[d$T_len == T0, ]
  x <- x[order(x$z), ]
  best <- which.min(x$msfe_rs)
  kstar <- x$z[best]
  min_loss <- x$msfe_rs[best]
  ksqrt <- min(max(round(sqrt(T0)), min(x$z)), max(x$z))
  isqrt <- which.min(abs(x$z - ksqrt))
  kcube <- min(max(round(T0^(1/3)), min(x$z)), max(x$z))
  icube <- which.min(abs(x$z - kcube))
  near01 <- x$z[x$msfe_rs <= 1.001 * min_loss]
  near05 <- x$z[x$msfe_rs <= 1.005 * min_loss]
  data.frame(
    T = T0,
    kstar = kstar,
    k_sqrt = x$z[isqrt],
    k_cube = x$z[icube],
    signed_gap = kstar - x$z[isqrt],
    absolute_gap = abs(kstar - x$z[isqrt]),
    relative_gap = (kstar - x$z[isqrt]) / kstar,
    ratio_to_sqrt = kstar / sqrt(T0),
    cube_signed_gap = kstar - x$z[icube],
    cube_absolute_gap = abs(kstar - x$z[icube]),
    cube_relative_gap = (kstar - x$z[icube]) / kstar,
    ratio_to_cube = kstar / T0^(1/3),
    min_msfe = min_loss,
    sqrt_msfe = x$msfe_rs[isqrt],
    sqrt_excess_loss_pct = 100 * (x$msfe_rs[isqrt] / min_loss - 1),
    cube_msfe = x$msfe_rs[icube],
    cube_excess_loss_pct = 100 * (x$msfe_rs[icube] / min_loss - 1),
    near01_low = min(near01), near01_high = max(near01),
    near05_low = min(near05), near05_high = max(near05)
  )
})
s <- do.call(rbind, summary_rows)

fit <- lm(log(kstar) ~ log(T), data = s)
s$k_fit <- exp(predict(fit))
cube_c <- exp(mean(log(s$kstar) - log(s$T) / 3))
sqrt_c <- exp(mean(log(s$kstar) - log(s$T) / 2))
s$cube_reference <- cube_c * s$T^(1/3)
s$sqrt_reference <- sqrt_c * sqrt(s$T)

write.csv(s, file.path(root, "T_scaling_summary.csv"), row.names = FALSE)
capture.output(summary(fit), file = file.path(root, "T_scaling_loglog_regression.txt"))
rule_comparison <- data.frame(
  rule = c("round(sqrt(T))", "round(T^(1/3))"),
  mean_absolute_k_gap = c(mean(s$absolute_gap), mean(s$cube_absolute_gap)),
  median_absolute_k_gap = c(median(s$absolute_gap), median(s$cube_absolute_gap)),
  mean_excess_msfe_pct = c(mean(s$sqrt_excess_loss_pct),
                           mean(s$cube_excess_loss_pct)),
  max_excess_msfe_pct = c(max(s$sqrt_excess_loss_pct),
                          max(s$cube_excess_loss_pct))
)
write.csv(rule_comparison, file.path(root, "sqrt_vs_cube_rule_comparison.csv"),
          row.names = FALSE)

figdir <- file.path(root, "figures")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
cols <- c(empirical = "#0072B2", sqrt = "#009E73", cube = "#CC79A7",
          fit = "#D55E00", band = "#999999")

pdf(file.path(figdir, "kstar_vs_T_bounds.pdf"), width = 8.5, height = 5.5)
ylim <- range(c(s$near05_low, s$near05_high, s$k_sqrt,
                s$cube_reference, s$sqrt_reference))
plot(s$T, s$kstar, type = "n", ylim = ylim,
     xlab = "training sample size T", ylab = "optimal subset size k*")
polygon(c(s$T, rev(s$T)), c(s$near05_low, rev(s$near05_high)),
        col = adjustcolor(cols["band"], alpha.f = 0.22), border = NA)
lines(s$T, s$cube_reference, col = cols["cube"], lwd = 2, lty = 3)
lines(s$T, s$sqrt_reference, col = cols["sqrt"], lwd = 2, lty = 2)
lines(s$T, s$kstar, col = cols["empirical"], lwd = 2.5)
points(s$T, s$kstar, col = cols["empirical"], pch = 16, cex = 1.1)
legend("topleft",
       c("empirical k*", "scaled T^(1/3)", "scaled T^(1/2)",
         "0.5% near-optimal band"),
       col = c(cols["empirical"], cols["cube"], cols["sqrt"], cols["band"]),
       lty = c(1, 3, 2, 1), lwd = c(2.5, 2, 2, 8), bty = "n")
dev.off()

pdf(file.path(figdir, "sqrt_cube_rule_diagnostics.pdf"), width = 10, height = 8)
par(mfrow = c(2, 2), mar = c(4.2, 4.4, 2.2, 1))
plot(s$T, s$kstar, type = "b", pch = 16, lwd = 2,
     col = cols["empirical"], xlab = "T", ylab = "subset size")
lines(s$T, s$k_sqrt, type = "b", pch = 17, lwd = 2,
      col = cols["sqrt"], lty = 2)
lines(s$T, s$k_cube, type = "b", pch = 15, lwd = 2,
      col = cols["cube"], lty = 3)
legend("topleft", c("oracle k*", "round(sqrt(T))", "round(T^(1/3))"),
       col = c(cols["empirical"], cols["sqrt"], cols["cube"]),
       pch = c(16, 17, 15), lty = c(1, 2, 3), bty = "n")
plot(s$T, s$signed_gap, type = "b", pch = 16, lwd = 2,
     col = cols["sqrt"], xlab = "T", ylab = "signed k gap")
lines(s$T, s$cube_signed_gap, type = "b", pch = 15, lwd = 2,
      col = cols["cube"], lty = 3)
abline(h = 0, col = "grey55", lty = 2)
plot(s$T, s$ratio_to_sqrt, type = "b", pch = 16, lwd = 2,
     col = cols["sqrt"], xlab = "T", ylab = "oracle-to-rule ratio")
lines(s$T, s$ratio_to_cube, type = "b", pch = 15, lwd = 2,
      col = cols["cube"], lty = 3)
abline(h = 1, col = "grey55", lty = 2)
plot(s$T, s$sqrt_excess_loss_pct, type = "b", pch = 16, lwd = 2,
     col = cols["sqrt"], xlab = "T",
     ylab = "excess MSFE, percent")
lines(s$T, s$cube_excess_loss_pct, type = "b", pch = 15, lwd = 2,
      col = cols["cube"], lty = 3)
abline(h = 0, col = "grey55", lty = 2)
par(mfrow = c(1, 1))
dev.off()

pdf(file.path(figdir, "kstar_loglog_scaling.pdf"), width = 7.5, height = 5.5)
plot(s$T, s$kstar, log = "xy", pch = 16, cex = 1.1,
     col = cols["empirical"], xlab = "T (log scale)",
     ylab = "k* (log scale)")
lines(s$T, s$k_fit, col = cols["fit"], lwd = 2)
lines(s$T, s$cube_reference, col = cols["cube"], lwd = 2, lty = 3)
lines(s$T, s$sqrt_reference, col = cols["sqrt"], lwd = 2, lty = 2)
legend("topleft",
       c(sprintf("fitted slope = %.3f", coef(fit)[2]),
         "scaled T^(1/3)", "scaled T^(1/2)"),
       col = c(cols["fit"], cols["cube"], cols["sqrt"]),
       lty = c(1, 3, 2), lwd = 2, bty = "n")
dev.off()

pdf(file.path(figdir, "relative_rmsfe_curves_by_T.pdf"), width = 8.5, height = 5.5)
plot(NA, xlim = range(d$z), ylim = range(d$rel_pooled_rmse),
     xlab = "subset size k", ylab = "relative RMSFE")
palette_cols <- hcl.colors(length(Ts), "Dark 3")
for (ii in seq_along(Ts)) {
  x <- d[d$T_len == Ts[ii], ]
  x <- x[order(x$z), ]
  lines(x$z, x$rel_pooled_rmse, col = palette_cols[ii], lwd = 2)
}
abline(h = 1, col = "grey55", lty = 2)
legend("topleft", paste0("T=", Ts), col = palette_cols, lwd = 2,
       ncol = 2, bty = "n")
dev.off()

cat("Wrote T-scaling summary and figures for", length(Ts), "T values.\n")
