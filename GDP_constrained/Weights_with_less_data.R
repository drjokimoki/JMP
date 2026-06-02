suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(forcats)
})

exclude <- c(18)

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
plot_cond_weight_over_time <- function(df_var, varname, every = 4) {
  df_var <- df_var %>% arrange(vdate)
  
  # label only every 'every'-th vintage
  brks <- df_var$vdate[seq(1, nrow(df_var), by = every)]
  
  ggplot(df_var, aes(x = vdate)) +
    geom_line(aes(y = mean_w_given_selected, linetype = "E[w | selected]")) +
    geom_point(aes(y = mean_w_given_selected), size = 0.7) +
    geom_line(aes(y = eq_w, linetype = "Equal weight (1/k)")) +
    scale_linetype_manual(
      values = c("E[w | selected]" = "solid", "Equal weight (1/k)" = "dashed"),
      name = NULL
    ) +
    scale_x_date(breaks = brks, date_labels = "%Y-%m") +
    labs(
      title = paste("Evolution of", varname, "conditional weight across vintages"),
      x = "Vintage", y = "Weight"
    ) +
    theme_minimal() +
    theme(
      legend.position = "top",
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      plot.title  = element_text(hjust = 0.5)
    )
}


plots <- imap(df_list, ~ plot_cond_weight_over_time(.x, .y, every = 4))  # or every = 8


# -------------------------------------------------------
# 4) save PDFs (base pdf device; no Cairo dependency)
# -------------------------------------------------------
out_dir <- "/Users/boris/Desktop/RSShapley/Combination_Monte/Final/codess/JMP/GDP_constrained/Weights_k_conditional/"
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
})

exclude <- c(18)   # set integer(0) if none
every   <- 4       # <- set to 8 if you want every 8th label

# ---------- safe getters (handles both named and positional outputs) ----------
get_summary <- function(res) {
  if (!is.null(res$summary)) res$summary else res[[1]]
}
get_table <- function(res) {
  if (!is.null(res$Fin_table2)) res$Fin_table2 else res[[2]]
}

# -------------------------------------------------------
# 1) Meta per vintage: vdate + p + 1/p  (equal weight 1/n)
# -------------------------------------------------------
meta_df <- imap_dfr(res_list, function(res, idx) {
  if (idx %in% exclude) return(NULL)
  
  summ <- get_summary(res)
  vin  <- as.character(summ$vintage[1])
  if (length(vin) == 0 || is.na(vin) || vin == "") return(NULL)
  
  tab <- get_table(res)
  if (!("Var" %in% names(tab))) return(NULL)
  
  # number of predictors available this vintage
  # (exclude growth_y if it happens to be in Var)
  p_t <- sum(as.character(tab$Var) != "growth_y", na.rm = TRUE)
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

stopifnot(nrow(meta_df) > 0)

# -------------------------------------------------------
# 2) Zero-bagged weights panel: Var x vintage
#    - NEW code: use w_bar
#    - OLD code: use Selected * Mean_OLS (approx zero-bagged)
# -------------------------------------------------------
zb_long <- imap_dfr(res_list, function(res, idx) {
  if (idx %in% exclude) return(NULL)
  
  summ <- get_summary(res)
  vin  <- as.character(summ$vintage[1])
  if (length(vin) == 0 || is.na(vin) || vin == "") return(NULL)
  
  tab <- get_table(res)
  
  if (!("Var" %in% names(tab))) return(NULL)
  
  # choose zero-bagged weight column:
  if ("w_bar" %in% names(tab)) {
    zb <- as.numeric(tab$w_bar)
  } else if (all(c("Selected", "Mean_OLS") %in% names(tab))) {
    zb <- as.numeric(tab$Selected) * as.numeric(tab$Mean_OLS)
  } else if (all(c("Selected", "mean_weight_nonzero") %in% names(tab))) {
    zb <- as.numeric(tab$Selected) * as.numeric(tab$mean_weight_nonzero)
  } else {
    stop("Cannot construct zero-bagged weights: need w_bar OR (Selected and Mean_OLS).")
  }
  
  tibble(
    unit_idx = idx,
    vintage  = vin,
    Var      = as.character(tab$Var),
    zero_bagged = zb
  ) %>%
    filter(Var != "growth_y")  # keep only predictors
})

stopifnot(nrow(zb_long) > 0)

panel <- zb_long %>%
  left_join(meta_df, by = c("unit_idx", "vintage")) %>%
  filter(!is.na(vdate), !is.na(eq_w)) %>%
  arrange(vdate)

stopifnot(nrow(panel) > 0)

df_list <- split(panel, panel$Var)

# -------------------------------------------------------
# 3) Plot function: zero-bagged weight vs equal 1/p
#    + label every 4th (or 8th) vintage
# -------------------------------------------------------
plot_zb_weight_over_time_1n <- function(df_var, varname, every = 4) {
  df_var <- df_var %>% arrange(vdate)
  brks <- df_var$vdate[seq(1, nrow(df_var), by = every)]
  
  ggplot(df_var, aes(x = vdate)) +
    geom_line(aes(y = zero_bagged, linetype = "Zero-bagged weight")) +
    geom_point(aes(y = zero_bagged), size = 0.7) +
    geom_line(aes(y = eq_w, linetype = "Equal weight (1/n)")) +
    scale_linetype_manual(
      values = c("Zero-bagged weight" = "solid", "Equal weight (1/n)" = "dashed"),
      name = NULL
    ) +
    scale_x_date(breaks = brks, date_labels = "%Y-%m") +
    labs(
      title = paste("Evolution of", varname, "zero-bagged weight across vintages"),
      x = "Vintage", y = "Weight"
    ) +
    theme_minimal() +
    theme(
      legend.position = "top",
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      plot.title  = element_text(hjust = 0.5)
    )
}

plots <- imap(df_list, ~ plot_zb_weight_over_time_1n(.x, .y, every = every))

# -------------------------------------------------------
# 4) Save PDFs (base pdf device; no Cairo)
# -------------------------------------------------------
out_dir <- "/Users/boris/Desktop/RSShapley/Combination_Monte/Final/codess/JMP/GDP_constrained/Weights_k/"
out_dir <- normalizePath(out_dir, winslash = "/", mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

iwalk(plots, function(p, varname) {
  safe <- substr(gsub("[^A-Za-z0-9_]+", "_", varname), 1, 80)
  fn <- file.path(out_dir, paste0(safe, ".pdf"))
  
  grDevices::pdf(fn, width = 15, height = 4, onefile = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  print(p)
})
