#!/usr/bin/env Rscript

ROOT <- "/Users/boris/Desktop/RSShapley/Combination_Monte/Final/codess/JMP"
infile <- file.path(ROOT, "RSM_endogenous_results", "scaling_T_summary.csv")
outdir <- file.path(ROOT, "RSM_endogenous_results")

df <- read.csv(infile, stringsAsFactors = FALSE)
df <- df[order(df$sigma_type, df$T_len), ]

T_grid <- sort(unique(df$T_len))
T0 <- 200
k0 <- 14
cube_ref <- k0 * (T_grid / T0)^(1 / 3)
sqrt_ref <- k0 * (T_grid / T0)^(1 / 2)

make_plot <- function(path) {
  png(path, width = 1400, height = 900, res = 180)
  par(mar = c(4.8, 5.0, 3.2, 2.0), family = "sans")

  ylim <- range(c(df$best_z_rmsfe, cube_ref, sqrt_ref), na.rm = TRUE)
  plot(
    NA, NA,
    xlim = range(T_grid),
    ylim = ylim,
    log = "x",
    xaxt = "n",
    xlab = "Training sample size T",
    ylab = "RMSFE-optimal subset size k*",
    main = "Endogenous Monte Carlo: optimal k growth"
  )
  axis(1, at = T_grid, labels = T_grid)
  grid(col = "grey88", lty = 1)

  cols <- c(toeplitz = "#2166AC", factor = "#B2182B")
  pchs <- c(toeplitz = 19, factor = 15)
  labs <- c(toeplitz = "Toeplitz endogenous MC", factor = "Factor endogenous MC")

  for (tp in c("toeplitz", "factor")) {
    sub <- df[df$sigma_type == tp, ]
    lines(sub$T_len, sub$best_z_rmsfe, col = cols[[tp]], lwd = 2.4)
    points(sub$T_len, sub$best_z_rmsfe, col = cols[[tp]], pch = pchs[[tp]], cex = 1.2)
  }

  lines(T_grid, cube_ref, col = "#4D9221", lwd = 2.2, lty = 2)
  lines(T_grid, sqrt_ref, col = "#762A83", lwd = 2.6, lty = 3)

  legend(
    "topleft",
    legend = c(labs, "Calibrated T^(1/3): k(200)=14", "Calibrated sqrt(T): k(200)=14"),
    col = c(cols, "#4D9221", "#762A83"),
    pch = c(pchs, NA, NA),
    lty = c(1, 1, 2, 3),
    lwd = c(2.4, 2.4, 2.2, 2.6),
    bty = "n",
    cex = 0.9
  )
  dev.off()
}

make_plot(file.path(outdir, "endogenous_kstar_growth.png"))

pdf(file.path(outdir, "endogenous_kstar_growth.pdf"), width = 7.2, height = 4.8)
par(mar = c(4.8, 5.0, 3.2, 2.0), family = "sans")
ylim <- range(c(df$best_z_rmsfe, cube_ref, sqrt_ref), na.rm = TRUE)
plot(
  NA, NA,
  xlim = range(T_grid),
  ylim = ylim,
  log = "x",
  xaxt = "n",
  xlab = "Training sample size T",
  ylab = "RMSFE-optimal subset size k*",
  main = "Endogenous Monte Carlo: optimal k growth"
)
axis(1, at = T_grid, labels = T_grid)
grid(col = "grey88", lty = 1)
cols <- c(toeplitz = "#2166AC", factor = "#B2182B")
pchs <- c(toeplitz = 19, factor = 15)
labs <- c(toeplitz = "Toeplitz endogenous MC", factor = "Factor endogenous MC")
for (tp in c("toeplitz", "factor")) {
  sub <- df[df$sigma_type == tp, ]
  lines(sub$T_len, sub$best_z_rmsfe, col = cols[[tp]], lwd = 2.4)
  points(sub$T_len, sub$best_z_rmsfe, col = cols[[tp]], pch = pchs[[tp]], cex = 1.2)
}
lines(T_grid, cube_ref, col = "#4D9221", lwd = 2.2, lty = 2)
lines(T_grid, sqrt_ref, col = "#762A83", lwd = 2.6, lty = 3)
legend(
  "topleft",
  legend = c(labs, "Calibrated T^(1/3): k(200)=14", "Calibrated sqrt(T): k(200)=14"),
  col = c(cols, "#4D9221", "#762A83"),
  pch = c(pchs, NA, NA),
  lty = c(1, 1, 2, 3),
  lwd = c(2.4, 2.4, 2.2, 2.6),
  bty = "n",
  cex = 0.9
)
dev.off()

cat(file.path(outdir, "endogenous_kstar_growth.png"), "\n")
