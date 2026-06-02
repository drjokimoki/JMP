# ---------------------------
# SETTINGS
# ---------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(forcats)
})

exclude <- c(8450)

# ---------- SAFE GETTERS (handle both named and positional lists) ----------
get_summary <- function(res) {
  if (!is.null(res$summary)) res$summary else res[[1]]
}
get_table <- function(res) {
  if (!is.null(res$Fin_table2)) res$Fin_table2 else res[[2]]
}
get_RSSS <- function(res) {
  if (!is.null(res$RSSS)) res$RSSS else res[[4]]
}

# (Optional) fallback vintage labels if summary$vintage is NA
# If you have 'pairs' in your environment from the main script, uncomment:
# get_vintage_fallback <- function(idx) pairs$vintage[idx]

# -------------------------------------------------------
# 1) Meta per vintage: vintage + date + k + 1/k
# -------------------------------------------------------
meta_df <- imap_dfr(res_list, function(res, idx) {
  if (idx %in% exclude) return(NULL)
  
  summ <- get_summary(res)
  vin  <- as.character(summ$vintage[1])
  if (length(vin) == 0 || is.na(vin) || vin == "") return(NULL)
  
  tab <- get_table(res)                 # Fin_table2
  p_t <- nrow(tab)                      # number of predictors
  if (!is.finite(p_t) || p_t < 1) return(NULL)
  
  vdate <- as.Date(paste0(substr(vin, 1, 7), "-01"))
  if (is.na(vdate)) return(NULL)
  
  tibble(
    unit_idx = idx,
    vintage  = vin,
    vdate    = vdate,
    p        = p_t,
    eq_w     = 1 / p_t
  )
}) %>% arrange(vdate)


stopifnot(nrow(meta_df) > 0)   # fail loudly if still empty

# -------------------------------------------------------
# 2) Long panel: Var x vintage with w_bar
# -------------------------------------------------------
weights_long <- imap_dfr(res_list, function(res, idx) {
  if (idx %in% exclude) return(NULL)
  
  summ <- get_summary(res)
  vin  <- as.character(summ$vintage[1])
  
  # fallback if needed (uncomment if you have 'pairs')
  # if (length(vin) == 0 || is.na(vin) || vin == "") vin <- get_vintage_fallback(idx)
  
  if (length(vin) == 0 || is.na(vin) || vin == "") return(NULL)
  
  tab <- get_table(res)
  if (!all(c("Var", "w_bar") %in% names(tab))) return(NULL)
  
  tibble(
    unit_idx = idx,
    vintage  = vin,
    Var      = as.character(tab$Var),
    w_bar    = as.numeric(tab$w_bar)
  )
})

stopifnot(nrow(weights_long) > 0)

panel <- weights_long %>%
  left_join(meta_df, by = c("unit_idx", "vintage")) %>%
  filter(!is.na(vdate), !is.na(eq_w)) %>%
  arrange(vdate)

stopifnot(nrow(panel) > 0)

df_list <- split(panel, panel$Var)

# -------------------------------------------------------
# 3) Plot function: w_bar vs Equal weight (1/k)
# -------------------------------------------------------
plot_weight_over_time <- function(df_var, varname) {
  df_var <- df_var %>% arrange(vdate)
  
  ggplot(df_var, aes(x = vdate)) +
    geom_line(aes(y = w_bar, linetype = "w_bar (RSM)")) +
    geom_point(aes(y = w_bar), size = 0.7) +
    geom_line(aes(y = eq_w, linetype = "Equal weight (1/n)")) +
    scale_linetype_manual(
      values = c("w_bar (RSM)" = "solid", "Equal weight (1/n)" = "dashed"),
      name = NULL
    ) +
    scale_x_date(breaks = df_var$vdate, date_labels = "%Y-%m") +
    labs(
      x = "Vintage",
      y = "Weight"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      plot.title  = element_text(hjust = 0.5)
    )
}

plots <- imap(df_list, plot_weight_over_time)

# -------------------------------------------------------
# 4) Save (PDF without Cairo)
# -------------------------------------------------------
out_dir <- "/Users/boris/Desktop/RSShapley/Combination_Monte/Final/codess/JMP/GDP_constrained/Weights_k_covid/"
out_dir <- normalizePath(out_dir, winslash = "/", mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

iwalk(plots, function(p, varname) {
  safe <- substr(gsub("[^A-Za-z0-9_]+", "_", varname), 1, 80)
  fn <- file.path(out_dir, paste0(safe, ".pdf"))
  
  grDevices::pdf(fn, width = 15, height = 4, onefile = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  print(p)
})




suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(forcats)
})

exclude <- c(8450)

# ---------- safe getters (handles both named and positional outputs) ----------
get_summary <- function(res) {
  if (!is.null(res$summary)) res$summary else res[[1]]
}
get_RSSS <- function(res) {
  if (!is.null(res$RSSS)) res$RSSS else res[[4]]
}

# -------------------------------------------------------
# 1) meta per vintage: vdate + k + 1/k
# -------------------------------------------------------
meta_df <- imap_dfr(res_list, function(res, idx) {
  if (idx %in% exclude) return(NULL)
  
  summ <- get_summary(res)
  vin  <- as.character(summ$vintage[1])
  if (length(vin) == 0 || is.na(vin) || vin == "") return(NULL)
  
  RSSS <- get_RSSS(res)
  k_t  <- tryCatch(ncol(RSSS$selected_vars), error = function(e) NA_integer_)
  if (!is.finite(k_t) || k_t < 1) return(NULL)
  
  vdate <- as.Date(paste0(substr(vin, 1, 7), "-01"))
  if (is.na(vdate)) return(NULL)
  
  tibble(
    unit_idx = idx,
    vintage  = vin,
    vdate    = vdate,
    k        = k_t,
    eq_w     = 1 / k_t
  )
}) %>% arrange(vdate)

stopifnot(nrow(meta_df) > 0)

# -------------------------------------------------------
# 2) conditional mean weights: E[w | selected] per vintage
#    computed from coef_long (successful draws only)
# -------------------------------------------------------
cond_long <- imap_dfr(res_list, function(res, idx) {
  if (idx %in% exclude) return(NULL)
  
  summ <- get_summary(res)
  vin  <- as.character(summ$vintage[1])
  if (length(vin) == 0 || is.na(vin) || vin == "") return(NULL)
  
  RSSS <- get_RSSS(res)
  cl   <- RSSS$coef_long
  
  # if no successful draws, return empty
  if (is.null(cl) || nrow(cl) == 0) return(NULL)
  
  cl %>%
    mutate(
      unit_idx = idx,
      vintage  = vin
    ) %>%
    group_by(unit_idx, vintage, variable) %>%
    summarise(
      mean_w_given_selected = mean(coef, na.rm = TRUE),  # includes zeros if QP sets them
      n_selected_success    = n_distinct(draw),
      .groups = "drop"
    ) %>%
    rename(Var = variable)
})

stopifnot(nrow(cond_long) > 0)

# join meta (adds vdate and eq_w=1/k)
panel <- cond_long %>%
  left_join(meta_df, by = c("unit_idx", "vintage")) %>%
  filter(!is.na(vdate), !is.na(eq_w)) %>%
  arrange(vdate)

stopifnot(nrow(panel) > 0)

df_list <- split(panel, panel$Var)

# -------------------------------------------------------
# 3) plot function: conditional weight vs 1/k
# -------------------------------------------------------
plot_cond_weight_over_time <- function(df_var, varname) {
  df_var <- df_var %>% arrange(vdate)
  
  ggplot(df_var, aes(x = vdate)) +
    geom_line(aes(y = mean_w_given_selected, linetype = "E[w | selected]")) +
    geom_point(aes(y = mean_w_given_selected), size = 0.7) +
    geom_line(aes(y = eq_w, linetype = "Equal weight (1/k)")) +
    scale_linetype_manual(
      values = c("E[w | selected]" = "solid", "Equal weight (1/k)" = "dashed"),
      name = NULL
    ) +
    scale_x_date(breaks = df_var$vdate, date_labels = "%Y-%m") +
    labs(
      x = "Vintage",
      y = "Weight"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      plot.title  = element_text(hjust = 0.5)
    )
}

plots <- imap(df_list, plot_cond_weight_over_time)

# -------------------------------------------------------
# 4) save PDFs (base pdf device; no Cairo dependency)
# -------------------------------------------------------
out_dir <- "/Users/boris/Desktop/RSShapley/Combination_Monte/Final/codess/JMP/GDP_constrained/Weights_k_conditional_covid/"
out_dir <- normalizePath(out_dir, winslash = "/", mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(length(plots) > 0)

iwalk(plots, function(p, varname) {
  safe <- substr(gsub("[^A-Za-z0-9_]+", "_", varname), 1, 80)
  fn <- file.path(out_dir, paste0(safe, ".pdf"))
  
  grDevices::pdf(fn, width = 15, height = 4, onefile = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  print(p)
})



suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(forcats)
  library(tidyr)
  library(lubridate)
})

# =========================
# USER SETTINGS
# =========================
exclude <- c(8450)     # indices of res_list to skip (leave integer(0) if none)
label_every <- 8       # show every 8th vintage label (use 4 for every 4th, etc.)

# Base folder where 3 subfolders will be created
out_dir_base <- "/Users/boris/Desktop/RSShapley/Combination_Monte/Final/codess/JMP/GDP_constrained/Plots_Contributions_covid/"
out_dir_base <- normalizePath(out_dir_base, winslash = "/", mustWork = FALSE)
dir.create(out_dir_base, recursive = TRUE, showWarnings = FALSE)

out_dir_1 <- file.path(out_dir_base, "01_Contribution")
out_dir_2 <- file.path(out_dir_base, "02_Contribution_plus_wbar")
out_dir_3 <- file.path(out_dir_base, "03_Contribution_vs_equalweight")

dir.create(out_dir_1, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_2, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_3, recursive = TRUE, showWarnings = FALSE)

# =========================
# SAFE GETTERS (your style)
# =========================
get_summary <- function(res) {
  if (!is.null(res$summary)) res$summary else res[[1]]
}
get_table <- function(res) {
  if (!is.null(res$Fin_table2)) res$Fin_table2 else res[[2]]
}

# =========================
# HELPERS
# =========================
safe_vdate <- function(vin) {
  vin <- as.character(vin)
  if (length(vin) == 0 || is.na(vin) || vin == "") return(NA_Date_)
  # assume "YYYY-MM"
  d <- suppressWarnings(as.Date(paste0(substr(vin, 1, 7), "-01")))
  if (!is.na(d)) return(d)
  # fallback
  suppressWarnings(ymd(paste0(vin, "-01")))
}

nth_breaks <- function(dates, n = 8) {
  d <- sort(unique(dates))
  if (length(d) == 0) return(d)
  d[seq(1, length(d), by = n)]
}

safe_filename <- function(x, maxlen = 80) {
  substr(gsub("[^A-Za-z0-9_]+", "_", x), 1, maxlen)
}

# =========================
# 1) META per vintage: p and eq_w = 1/p
# =========================
meta_df <- imap_dfr(res_list, function(res, idx) {
  if (idx %in% exclude) return(NULL)
  
  summ <- get_summary(res)
  if (is.null(summ) || !("vintage" %in% names(summ))) return(NULL)
  
  vin <- as.character(summ$vintage[1])
  vdate <- safe_vdate(vin)
  if (is.na(vdate)) return(NULL)
  
  tab <- get_table(res)
  if (is.null(tab) || !is.data.frame(tab)) return(NULL)
  
  p_t <- nrow(tab)
  if (!is.finite(p_t) || p_t < 1) return(NULL)
  
  tibble(
    unit_idx = idx,
    vintage  = vin,
    vdate    = vdate,
    p        = p_t,
    eq_w     = 1 / p_t
  )
}) %>% arrange(vdate)

stopifnot(nrow(meta_df) > 0)

# =========================
# 2) Long panel: Var x vintage with Contribution, w_bar, etc.
# =========================
panel_long <- imap_dfr(res_list, function(res, idx) {
  if (idx %in% exclude) return(NULL)
  
  summ <- get_summary(res)
  if (is.null(summ) || !("vintage" %in% names(summ))) return(NULL)
  
  vin <- as.character(summ$vintage[1])
  vdate <- safe_vdate(vin)
  if (is.na(vdate)) return(NULL)
  
  tab <- get_table(res)
  if (is.null(tab) || !is.data.frame(tab)) return(NULL)
  
  need_cols <- c("Var", "w_bar", "Contribution", "Prediction", "Hist_mean")
  if (!all(need_cols %in% names(tab))) return(NULL)
  
  tibble(
    unit_idx     = idx,
    vintage      = vin,
    vdate        = vdate,
    Var          = as.character(tab$Var),
    w_bar        = as.numeric(tab$w_bar),
    Contribution = as.numeric(tab$Contribution),
    Prediction   = as.numeric(tab$Prediction),
    Hist_mean    = as.numeric(tab$Hist_mean)
  )
})

stopifnot(nrow(panel_long) > 0)

panel <- panel_long %>%
  left_join(meta_df, by = c("unit_idx", "vintage", "vdate")) %>%
  filter(!is.na(eq_w), !is.na(vdate)) %>%
  arrange(vdate, Var) %>%
  mutate(
    Contribution_eq = eq_w * (Prediction - Hist_mean)
  )

stopifnot(nrow(panel) > 0)

df_list <- split(panel, panel$Var)

# =========================
# 3) PLOT FUNCTIONS (3 groups)
# =========================

# (1) Contribution over time
plot_contribution <- function(df_var, varname) {
  df_var <- df_var %>% arrange(vdate)
  
  ggplot(df_var, aes(x = vdate, y = Contribution)) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_line() +
    geom_point(size = 0.7) +
    scale_x_date(
      breaks = nth_breaks(df_var$vdate, n = label_every),
      date_labels = "%Y-%m"
    ) +
    labs(
      x = "Vintage",
      y = "Contribution"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      plot.title  = element_text(hjust = 0.5)
    )
}

# (2) Contribution + w_bar (stacked panels)
plot_contribution_plus_wbar <- function(df_var, varname) {
  df_var <- df_var %>% arrange(vdate) %>%
    select(vdate, Contribution, w_bar) %>%
    pivot_longer(cols = c(Contribution, w_bar),
                 names_to = "series",
                 values_to = "value") %>%
    mutate(series = factor(series, levels = c("Contribution", "w_bar")))
  
  ggplot(df_var, aes(x = vdate, y = value)) +
    geom_hline(data = subset(df_var, series == "Contribution"),
               aes(yintercept = 0), inherit.aes = FALSE, linewidth = 0.3) +
    geom_line() +
    geom_point(size = 0.6) +
    facet_wrap(~ series, ncol = 1, scales = "free_y") +
    scale_x_date(
      breaks = nth_breaks(df_var$vdate, n = label_every),
      date_labels = "%Y-%m"
    ) +
    labs(
      x = "Vintage",
      y = NULL
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      plot.title  = element_text(hjust = 0.5)
    )
}

# (3) Contribution (w_bar) vs equal-weight contribution (1/n)
plot_contribution_vs_equal <- function(df_var, varname) {
  df_var <- df_var %>% arrange(vdate) %>%
    select(vdate, Contribution, Contribution_eq) %>%
    pivot_longer(cols = c(Contribution, Contribution_eq),
                 names_to = "scheme",
                 values_to = "value") %>%
    mutate(
      scheme = recode(scheme,
                      Contribution = "w_bar contribution",
                      Contribution_eq = "equal-weight (1/n) contribution")
    )
  
  ggplot(df_var, aes(x = vdate, y = value, linetype = scheme)) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_line() +
    geom_point(size = 0.6) +
    scale_linetype_manual(
      values = c("w_bar contribution" = "solid",
                 "equal-weight (1/n) contribution" = "dashed"),
      name = NULL
    ) +
    scale_x_date(
      breaks = nth_breaks(df_var$vdate, n = label_every),
      date_labels = "%Y-%m"
    ) +
    labs(
      x = "Vintage",
      y = "Contribution"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      plot.title  = element_text(hjust = 0.5)
    )
}

# =========================
# 4) MAKE + SAVE (one pdf per variable per group)
# =========================
plots_1 <- imap(df_list, plot_contribution)
plots_2 <- imap(df_list, plot_contribution_plus_wbar)
plots_3 <- imap(df_list, plot_contribution_vs_equal)

save_plot_list <- function(plot_list, out_dir, width = 15, height = 4) {
  iwalk(plot_list, function(p, varname) {
    safe <- safe_filename(varname)
    fn <- file.path(out_dir, paste0(safe, ".pdf"))
    
    grDevices::pdf(fn, width = width, height = height, onefile = FALSE)
    on.exit(grDevices::dev.off(), add = TRUE)
    print(p)
  })
}

save_plot_list(plots_1, out_dir_1, width = 15, height = 4)
save_plot_list(plots_2, out_dir_2, width = 15, height = 6)  # taller: 2 panels
save_plot_list(plots_3, out_dir_3, width = 15, height = 4)

cat("\nSaved:\n",
    " (1) ", out_dir_1, "\n",
    " (2) ", out_dir_2, "\n",
    " (3) ", out_dir_3, "\n", sep = "")


res_list[[82]][[2]][["w_bar"]]
