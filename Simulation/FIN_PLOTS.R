rm(list = ls())

Fin_data<-read.csv("/Users/boris/Desktop/RSShapley/Combination_Monte/Final/codess/JMP/Simulation/Full_Results1.csv", sep = ",", stringsAsFactors = F)

Fin_data$rel_pooled_rmse<-Fin_data$rmse_pooled_rs/Fin_data$rmse_pooled_base

# install.packages(c("dplyr","ggplot2"))
library(dplyr)
library(ggplot2)
library(ggrepel)

# ---- OUTPUT SETTINGS ----
out_dir <- "/Users/boris/Desktop/RSShapley/Combination_Monte/Final/codess/JMP/Simulation/plots/plots_rel_pooled_rmse_by_R2"
#dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- (Optional) ensure types ----
Fin_data <- Fin_data %>%
  mutate(
    across(c(T_len, rho, phi_f, phi_eta, z, rel_pooled_rmse), as.numeric),
    # Keep R2_target as factor for clean legend (works even if numeric in file)
    R2_target = as.factor(R2_target)
  ) %>%
  filter(!is.na(z), !is.na(rel_pooled_rmse), !is.na(R2_target))


combos <- Fin_data %>%
  distinct(T_len, rho, phi_f, phi_eta) %>%
  arrange(T_len, rho, phi_f, phi_eta)

#Plots for R2
for (i in seq_len(nrow(combos))) {
  rowi <- combos[i, ]
  
  df_i <- Fin_data %>%
    dplyr::filter(
      T_len   == rowi$T_len,
      rho     == rowi$rho,
      phi_f   == rowi$phi_f,
      phi_eta == rowi$phi_eta
    ) %>%
    dplyr::arrange(R2_target, z) %>%
    dplyr::mutate(R2_target = as.factor(R2_target))
  
  if (nrow(df_i) == 0) next
  
  mins <- df_i %>%
    dplyr::group_by(R2_target) %>%
    dplyr::slice_min(rel_pooled_rmse, with_ties = FALSE) %>%
    dplyr::ungroup()
  
  z_special_val <- dplyr::case_when(
    rowi$T_len == 100 ~ 10L,
    rowi$T_len == 500 ~ 22L,
    TRUE              ~ NA_integer_
  )
  
  z_special <- df_i %>%
    dplyr::filter(!is.na(z_special_val), z == z_special_val) %>%
    dplyr::group_by(R2_target) %>% dplyr::slice_tail(n = 1) %>%
    dplyr::ungroup()
  
  # avoid double labels if min happens at special z
  z_special_unique <- z_special %>%
    dplyr::anti_join(mins %>% dplyr::select(R2_target, z),
                     by = c("R2_target","z"))
  
  p <- ggplot(df_i, aes(x = z, y = rel_pooled_rmse,
                        color = R2_target, group = R2_target)) +
    geom_line(size = 0.9, na.rm = TRUE) +
    geom_hline(yintercept = 1, color = "grey60", show.legend = FALSE) +
    labs(title = "",
         x = "subset dimension", y = "rRMSE",
         color = expression(R^2)) +
    theme_minimal(base_size = 12)
  
  if (!is.na(z_special_val)) {
    p <- p + geom_vline(xintercept = z_special_val,
                        linetype = "dashed", color = "grey60",
                        linewidth = 0.6, show.legend = FALSE)
  }
  
  # dots (same size) for mins & special z
  p <- p +
    geom_point(data = mins, aes(x = z, y = rel_pooled_rmse),
               shape = 16, size = 2.4, alpha = 0.95, show.legend = FALSE) +
    geom_point(data = z_special_unique, aes(x = z, y = rel_pooled_rmse),
               shape = 16, size = 2.4, alpha = 0.95, show.legend = FALSE)
  
  # -------------------------
  # LABELS (default behavior)
  # -------------------------
  add_labels <- function(plot_obj, data, seed = 1) {
    if (nrow(data) == 0) return(plot_obj)
    plot_obj +
      geom_text_repel(
        data = data,
        aes(label = round(rel_pooled_rmse, 2), x = z, y = rel_pooled_rmse, color = R2_target),
        show.legend = FALSE, size = 3,
        box.padding = 0.15, point.padding = 0.12,
        force = 12, max.time = 1.2, max.overlaps = Inf,
        min.segment.length = 0, segment.size = 0.25,
        seed = seed
      )
  }
  
  # -------------------------
  # SPECIAL CASE: i == 6
  #   - stagger label x by R2 level (labels only)
  #   - repel primarily in Y; add room & no clipping
  # -------------------------
  if (i == 6) {
    # offsets per R^2 level (adjust if you have different levels/order)
    # This keeps dots/lines at true x=z; labels shift a tiny bit sideways.
    offset_map <- data.frame(
      R2_target = levels(df_i$R2_target),
      x_off = c(-0.18, 0.00, 0.18)   # left / center / right
    )
    
    mins_lab <- dplyr::left_join(mins, offset_map, by = "R2_target") %>%
      dplyr::mutate(label_x = z + x_off, label_y = rel_pooled_rmse)
    
    zspec_lab <- dplyr::left_join(z_special_unique, offset_map, by = "R2_target") %>%
      dplyr::mutate(label_x = z + x_off, label_y = rel_pooled_rmse)
    
    p <- p +
      scale_x_continuous(expand = expansion(mult = c(0.02, 0.08))) +
      scale_y_continuous(expand = expansion(mult = c(0.02, 0.10))) +
      coord_cartesian(clip = "off") +
      theme(plot.margin = margin(10, 28, 10, 10))
    
    if (nrow(mins_lab) > 0) {
      p <- p + geom_text_repel(
        data = mins_lab,
        aes(x = label_x, y = label_y, label = round(rel_pooled_rmse, 2), color = R2_target),
        show.legend = FALSE, size = 3,
        direction = "y",      # spread vertically; x is mostly fixed
        box.padding = 0.15, point.padding = 0.12,
        force = 30, max.time = 1.5, max.overlaps = Inf,
        min.segment.length = 0, segment.size = 0.25,
        seed = 21
      )
    }
    if (nrow(zspec_lab) > 0) {
      p <- p + geom_text_repel(
        data = zspec_lab,
        aes(x = label_x, y = label_y, label = round(rel_pooled_rmse, 2), color = R2_target),
        show.legend = FALSE, size = 3,
        direction = "y",
        box.padding = 0.15, point.padding = 0.12,
        force = 30, max.time = 1.5, max.overlaps = Inf,
        min.segment.length = 0, segment.size = 0.25,
        seed = 22
      )
    }
  } else {
    # all other panels: your usual labels
    p <- add_labels(p, mins, seed = 1)
    p <- add_labels(p, z_special_unique, seed = 2)
  }
  
  print(p)
  
  ggsave(
    filename = sprintf("rel_rmse_byR2_T%s_rho%s_phif%s_phieta%s.pdf",
                       rowi$T_len, rowi$rho, rowi$phi_f, rowi$phi_eta),
    path = out_dir, plot = p, width = 10, height = 6, dpi = 150
  )
}



# Packages
library(dplyr)
library(ggplot2)

# ---- OUTPUT ----
out_dir <- "/Users/boris/Desktop/RSShapley/Combination_Monte/Final/codess/JMP/Simulation/plots/Q_by_R2"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# (Optional) ensure types are right
Fin_data <- Fin_data %>%
  mutate(
    R2_target  = as.factor(R2_target),
    z          = as.numeric(z),
    Q_rs_bag   = as.numeric(Q_rs_bag),
    Q_ols      = as.numeric(Q_ols),
    T_len      = as.numeric(T_len),
    rho        = as.numeric(rho),
    phi_f      = as.numeric(phi_f),
    phi_eta    = as.numeric(phi_eta)
  )

# All unique combos
combos <- Fin_data %>%
  distinct(T_len, rho, phi_f, phi_eta) %>%
  arrange(T_len, rho, phi_f, phi_eta)

for (i in seq_len(nrow(combos))) {
  rowi <- combos[i, ]
  
  # Filter one panel's data
  df_i <- Fin_data %>%
    filter(
      T_len   == rowi$T_len,
      rho     == rowi$rho,
      phi_f   == rowi$phi_f,
      phi_eta == rowi$phi_eta
    ) %>%
    arrange(R2_target, z)
  
  if (nrow(df_i) == 0) next
  
  # y value: Q_rs_bag normally, but if z==50 and Q_ols available, use Q_ols
  df_i <- df_i %>%
    mutate(
      Q_plot = ifelse(z == 50 & !is.na(Q_ols), Q_ols, Q_rs_bag)
    )
  
  title_i <- sprintf(
    "Q vs z — T=%s, rho=%s, phi_f=%s, phi_eta=%s  (y: RS-bag, at z=50 → OLS)",
    rowi$T_len, rowi$rho, rowi$phi_f, rowi$phi_eta
  )
  
  p <- ggplot(df_i, aes(x = z, y = Q_plot,
                        color = R2_target, group = R2_target)) +
    geom_line(size = 0.9, na.rm = TRUE) +
    # Mark where the y-source switches
    labs(
      title = "",
      x = "subset size",
      y = "",
      color = expression(R^2)
    ) +
    theme_minimal(base_size = 12)
  
  print(p)
  
  # Save one PDF per combo (adjust to PNG if you prefer)
  ggsave(
    filename = sprintf("Q_byR2_T%s_rho%s_phif%s_phieta%s.pdf",
                       rowi$T_len, rowi$rho, rowi$phi_f, rowi$phi_eta),
    path = out_dir, plot = p, width = 10, height = 6, dpi = 150
  )
}
