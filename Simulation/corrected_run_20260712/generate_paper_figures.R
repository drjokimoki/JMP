#!/usr/bin/env Rscript

options(warn = 1)

library(dplyr)
library(ggplot2)
library(ggrepel)

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- args_all[grep("^--file=", args_all)][1]
script_path <- sub("^--file=", "", script_arg)
root <- normalizePath(dirname(script_path), mustWork = TRUE)
infile <- file.path(root, "rs_results.csv")
rel_dir <- file.path(root, "figures", "relative_rmsfe")
q_dir <- file.path(root, "figures", "Q_rsm")
dir.create(rel_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(q_dir, recursive = TRUE, showWarnings = FALSE)

d <- read.csv(infile, stringsAsFactors = FALSE) |>
  mutate(
    rel_pooled_rmse = rmse_pooled_rs / rmse_pooled_base,
    R2_target = factor(R2_target, levels = c(0.1, 0.3, 0.5))
  )

combos <- d |>
  distinct(T_len, rho, phi_f, phi_eta) |>
  arrange(T_len, rho, phi_f, phi_eta)

tag <- function(x) format(x, trim = TRUE, scientific = FALSE)
paper_tag <- function(x) sub("\\.", "p", tag(x))

for (ii in seq_len(nrow(combos))) {
  cc <- combos[ii, ]
  dd <- d |>
    filter(T_len == cc$T_len, rho == cc$rho,
           phi_f == cc$phi_f, phi_eta == cc$phi_eta) |>
    arrange(R2_target, z)

  mins <- dd |>
    group_by(R2_target) |>
    slice_min(rel_pooled_rmse, n = 1, with_ties = FALSE) |>
    ungroup()
  sqrt_k <- round(sqrt(cc$T_len))
  sqrt_pts <- dd |>
    filter(z == sqrt_k) |>
    anti_join(select(mins, R2_target, z), by = c("R2_target", "z"))

  p_rel <- ggplot(dd, aes(z, rel_pooled_rmse,
                          color = R2_target, group = R2_target)) +
    geom_line(size = 0.9, na.rm = TRUE) +
    geom_hline(yintercept = 1, color = "grey60", show.legend = FALSE) +
    geom_vline(xintercept = sqrt_k, linetype = "dashed",
               color = "grey60", linewidth = 0.6, show.legend = FALSE) +
    geom_point(data = mins, size = 2.4, show.legend = FALSE) +
    geom_point(data = sqrt_pts, size = 2.4, show.legend = FALSE) +
    labs(title = "", x = "subset dimension", y = "rRMSFE",
         color = expression(R^2)) +
    theme_minimal(base_size = 12)

  add_labels <- function(plot_obj, dat, seed) {
    if (nrow(dat) == 0) return(plot_obj)
    plot_obj + geom_text_repel(
      data = dat,
      aes(label = round(rel_pooled_rmse, 2), x = z,
          y = rel_pooled_rmse, color = R2_target),
      show.legend = FALSE, size = 3,
      box.padding = 0.15, point.padding = 0.12,
      force = 12, max.time = 1.2, max.overlaps = Inf,
      min.segment.length = 0, segment.size = 0.25, seed = seed
    )
  }

  # This reproduces the paper script's special label treatment for panel 6.
  if (ii == 6) {
    offset_map <- data.frame(
      R2_target = levels(dd$R2_target),
      x_off = c(-0.18, 0, 0.18)
    )
    mins_lab <- left_join(mins, offset_map, by = "R2_target") |>
      mutate(label_x = z + x_off, label_y = rel_pooled_rmse)
    sqrt_lab <- left_join(sqrt_pts, offset_map, by = "R2_target") |>
      mutate(label_x = z + x_off, label_y = rel_pooled_rmse)
    p_rel <- p_rel +
      scale_x_continuous(expand = expansion(mult = c(0.02, 0.08))) +
      scale_y_continuous(expand = expansion(mult = c(0.02, 0.10))) +
      coord_cartesian(clip = "off") +
      theme(plot.margin = margin(10, 28, 10, 10))
    if (nrow(mins_lab) > 0) p_rel <- p_rel + geom_text_repel(
      data = mins_lab,
      aes(label_x, label_y, label = round(rel_pooled_rmse, 2),
          color = R2_target),
      show.legend = FALSE, size = 3, direction = "y",
      box.padding = 0.15, point.padding = 0.12,
      force = 30, max.time = 1.5, max.overlaps = Inf,
      min.segment.length = 0, segment.size = 0.25, seed = 21
    )
    if (nrow(sqrt_lab) > 0) p_rel <- p_rel + geom_text_repel(
      data = sqrt_lab,
      aes(label_x, label_y, label = round(rel_pooled_rmse, 2),
          color = R2_target),
      show.legend = FALSE, size = 3, direction = "y",
      box.padding = 0.15, point.padding = 0.12,
      force = 30, max.time = 1.5, max.overlaps = Inf,
      min.segment.length = 0, segment.size = 0.25, seed = 22
    )
  } else {
    p_rel <- add_labels(p_rel, mins, 1)
    p_rel <- add_labels(p_rel, sqrt_pts, 2)
  }

  rel_name <- sprintf("T%s_rho%s_phif%s_phieta%s.pdf",
                      tag(cc$T_len), paper_tag(cc$rho),
                      paper_tag(cc$phi_f), paper_tag(cc$phi_eta))
  ggsave(file.path(rel_dir, rel_name), p_rel,
         width = 12, height = 5, device = "pdf")

  dd_q <- dd |>
    mutate(Q_plot = ifelse(z == 50 & is.finite(Q_ols), Q_ols, Q_rs_bag))
  p_q <- ggplot(dd_q, aes(z, Q_plot,
                          color = R2_target, group = R2_target)) +
    geom_line(size = 0.9, na.rm = TRUE) +
    labs(title = "", x = "subspace dimension", y = "Q_rs",
         color = expression(R^2)) +
    theme_minimal(base_size = 12)

  q_name <- sprintf("T%s_rho%s_phif%s_phieta%s.pdf",
                    tag(cc$T_len), paper_tag(cc$rho),
                    paper_tag(cc$phi_f), paper_tag(cc$phi_eta))
  ggsave(file.path(q_dir, q_name), p_q,
         width = 12, height = 5, device = "pdf")
}

summary <- d |>
  group_by(T_len, rho, phi_f, phi_eta, R2_target) |>
  slice_min(rel_pooled_rmse, n = 1, with_ties = FALSE) |>
  transmute(T_len, rho, phi_f, phi_eta, R2_target,
            best_k = z, min_rel_pooled_rmse = rel_pooled_rmse) |>
  ungroup()
write.csv(summary, file.path(root, "figure_minima_summary.csv"), row.names = FALSE)

cat("Generated", nrow(combos), "relative-RMSFE and", nrow(combos),
    "Q_RSM figures.\n")
