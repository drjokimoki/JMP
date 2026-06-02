rm(list=ls())
ffd<-readRDS("/Users/boris/Desktop/RSShapley/Combination_Monte/Final/codess/JMP/IP_Weighted/final_list.RData")
ffd1<-readRDS("/Users/boris/Desktop/RSShapley/Combination_Monte/Final/codess/JMP/IP_UNIFORM_d/final_list.RData")
ffd<-readRDS("/Users/boris/Desktop/RSShapley/Combination_Monte/Final/codess/JMP/Weighted_shapley_per_quarter/final_list.RData")

ffd<-ffd[-18]
ffd1<-ffd1[-18]


library(dplyr)
library(purrr)
library(ggplot2)
library(tidyr)

# --- Build the long data (same idea as before) ---

# --- Build the long data (same idea as before) ---
unit_labels <- c()

for(i in 1:length(ffd)){
  unit_labels[i]<- ffd[[i]][["summary"]][["vintage"]]
}

for (j in seq_along(ffd)) {
  df <- ffd[[j]][[2]]
  df$zero_bagged_ols <- df$Selected * df$Mean_OLS
  ffd[[j]][[2]] <- df
}

nn<-c()
for (j in seq_along(ffd)) {
  df <- ffd[[j]][[2]]
  nn[j]<-1/nrow(df)
}


for (j in seq_along(ffd)) {
  df <- ffd[[j]][[2]]
  df$Thr<-1/(nrow(df)*df$Selected)
  ffd[[j]][[2]] <- df
}


all_coef <- map2_dfr(ffd, seq_along(ffd), function(x, j) {
  tbl <- x[[2]]  # IMPORTANT: list indexing with [[ ]]
  if (is.null(tbl) || !all(c("Var","zero_bagged_ols") %in% names(tbl))) return(NULL)
  tibble(
    unit_idx = j,
    unit     = unit_labels[j],
    Var      = tbl$Var,
    zero_bagged_ols = tbl$zero_bagged_ols
  )
})

all_coef <- map2_dfr(ffd, seq_along(ffd), function(x, j) {
  tbl <- x[[2]]  # IMPORTANT: list indexing with [[ ]]
  if (is.null(tbl) || !all(c("Var","Mean_OLS") %in% names(tbl))) return(NULL)
  tibble(
    unit_idx = j,
    unit     = unit_labels[j],
    Var      = tbl$Var,
    Mean_OLS = tbl$Mean_OLS
  )
})

all_coef <- map2_dfr(ffd, seq_along(ffd), function(x, j) {
  tbl <- x[[2]]  # IMPORTANT: list indexing with [[ ]]
  if (is.null(tbl) || !all(c("Var","Mean_OLS") %in% names(tbl))) return(NULL)
  tibble(
    unit_idx = j,
    unit     = unit_labels[j],
    Var      = tbl$Var,
    Mean_OLS = tbl$Mean_OLS,
    Thr=tbl$Thr
  )
})

# --- Make one plot per variable ---
df_list <- split(all_coef, all_coef$Var)

plots <- imap(df_list, ~
                ggplot(.x, aes(x = unit_idx, y = Mean_OLS)) +
                geom_line() +
                geom_point(size = 0.7) +
                scale_x_continuous(
                  breaks = unique(.x$unit_idx),
                  labels = unit_labels[unique(.x$unit_idx)]
                ) +
                labs(
                  title = paste("Evolution of", .y, "OLS coefficient across time"),
                  x = "Unit",
                  y = "Mean OLS"
                ) +
                theme_minimal() +
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
                theme(plot.title = element_text(hjust = 0.5))
)

# Print one (example)
plots[[72]]





# --- Build the long data (same idea as before) ---
unit_labels <- c()

for(i in 1:length(ffd)){
  unit_labels[i]<- ffd[[i]][["summary"]][["vintage"]]
}

all_coef <- map2_dfr(ffd, seq_along(ffd), function(x, j) {
  tbl <- x[[2]]  # IMPORTANT: list indexing with [[ ]]
  if (is.null(tbl) || !all(c("Var","Shapley") %in% names(tbl))) return(NULL)
  tibble(
    unit_idx = j,
    unit     = unit_labels[j],
    Var      = tbl$Var,
    Shapley = tbl$Shapley
  )
})

# --- Make one plot per variable ---
df_list <- split(all_coef, all_coef$Var)

plots <- imap(df_list, ~
                ggplot(.x, aes(x = unit_idx, y = Shapley)) +
                geom_line() +
                geom_point(size = 0.7) +
                scale_x_continuous(
                  breaks = unique(.x$unit_idx),
                  labels = unit_labels[unique(.x$unit_idx)]
                ) +
                labs(
                  title = paste("Evolution of Shapley Value of", .y, "across time"),
                  x = "Time",
                  y = "Shapley"
                ) +
                theme_minimal() +
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
              +
                theme(plot.title = element_text(hjust = 0.5))
)

# Print one (example)
plots[[72]]

library(purrr)
library(ggplot2)

# choose a single output directory and ensure it exists
out_dir <- "/Users/boris/Desktop/RSShapley/Combination_Monte/Final/codess/JMP/Uniform_shapley_per_quarter/OLS_coef_over_time/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

iwalk(plots, ~ {
  # make a safe file name from the variable name
  safe <- make.names(.y)
  fn <- file.path(out_dir, paste0(safe, ".pdf"))
  
  ggsave(
    filename = fn,
    plot = .x,
    width = 15,
    height = 4,
    device = "pdf" # or "pdf" if Cairo isn't available
    # dpi is ignored for PDF, so omit it
  )
})

fas<-ffd[[84]][[2]]
sum(fas$Shapley)


library(dplyr)
library(purrr)
library(lubridate)

# Long table: Var | Prediction | vintage
pred_long <- imap_dfr(ffd, function(x, j) {
  tbl <- x[[2]]
  if (is.null(tbl) || !all(c("Var","Prediction") %in% names(tbl))) return(NULL)
  
  vint_raw <- tryCatch(x[["summary"]][["vintage"]], error = function(e) NA_character_)
  vint_chr <- if (length(vint_raw) == 0 || is.na(vint_raw)) NA_character_ else as.character(vint_raw)
  
  tibble(
    Var        = tbl$Var,
    Prediction = suppressWarnings(as.numeric(tbl$Prediction)),
    vintage    = vint_chr
  )
})

# (Optional) Parse vintage into a Date, then keep the 3 requested columns
pred_long <- pred_long %>%
  mutate(
    vintage = suppressWarnings(as_date(parse_date_time(
      vintage,
      orders = c("Ymd","Y-m-d","Y-m","Ym","dmy","mdy","ymd HMS","ymd HM","ymd H","ymd")
    )))
  ) %>%
  arrange(Var, vintage)

# Result: columns are Var, Prediction, vintage (date)
pred_long


for (j in seq_along(ffd)) {
  df <- ffd[[j]][[2]]
  df$zero_bagged_ols <- df$Selected * df$Mean_OLS
  dd<-ffd1[[j]][[2]]
  dd$zero_bagged_ols_uniform <- dd$Selected * dd$Mean_OLS
  fass<- df %>%
    left_join(dd, by= "Var") %>%
    select(Var, zero_bagged_ols_uniform, zero_bagged_ols)
           
  ffd[[j]][[2]] <- fass
}


fa<-ffd[[2]][[2]]
fas<-ffd1[[2]][[2]]
names(fas)[3]<-"Mean_OLS_uniform"
fass<- fa %>%
  left_join(fas, by= "Var") %>%
  select(Var, Mean_OLS, Mean_OLS_uniform)
  

all_coef <- map2_dfr(ffd, seq_along(ffd), function(x, j) {
  tbl <- x[[2]]  # IMPORTANT: list indexing with [[ ]]
  if (is.null(tbl) || !all(c("Var","zero_bagged_ols", "zero_bagged_ols_uniform") %in% names(tbl))) return(NULL)
  tibble(
    unit_idx = j,
    unit     = unit_labels[j],
    Var      = tbl$Var,
    zero_bagged_ols = tbl$zero_bagged_ols,
    zero_bagged_ols_uniform=tbl$zero_bagged_ols_uniform
  )
})

?scale_linetype_manual

# 1) Rebuild plots adding the k line (dashed) + legend
plots <- imap(df_list, ~
                ggplot(.x, aes(x = unit_idx)) +
                geom_line(aes(y = zero_bagged_ols, linetype = "zero_bagged_ols", color = "zero_bagged_ols")) +
                geom_point(aes(y = zero_bagged_ols), size = 0.7) +
                geom_line(aes(y = zero_bagged_ols_uniform, linetype = "zero_bagged_ols_uniform", color = "zero_bagged_ols_uniform")) +
                geom_point(aes(y = zero_bagged_ols_uniform), size = 0.7)+
                geom_line(data = k_df, aes(y = k, linetype = "Equal Weight")) +
                scale_linetype_manual(values = c("zero_bagged_ols" = "solid", "Equal Weight" = "dashed", "zero_bagged_ols_uniform" = "solid"), name = NULL) +
                scale_x_continuous(
                  breaks = unique(.x$unit_idx),
                  labels = unit_labels[unique(.x$unit_idx)]
                ) +
                labs(
                  title = paste("Evolution of", .y, "OLS coefficient across time"),
                  x = "Vintage",
                  y = NULL
                ) +
                theme_minimal() +
                theme(
                  legend.position = "top",
                  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                  plot.title  = element_text(hjust = 0.5)
                )
)



plots <- imap(df_list, ~
                ggplot(.x, aes(x = unit_idx)) +
                geom_line(aes(y = zero_bagged_ols,         linetype = "OLS (bagged)",  color = "OLS (bagged)")) +
                geom_point(aes(y = zero_bagged_ols,                               color = "OLS (bagged)"),  size = 0.7) +
                geom_line(aes(y = zero_bagged_ols_uniform,  linetype = "OLS (uniform)", color = "OLS (uniform)")) +
                geom_point(aes(y = zero_bagged_ols_uniform,                        color = "OLS (uniform)"), size = 0.7) +
                geom_line(data = k_df, aes(y = k,           linetype = "Equal Weight",  color = "Equal Weight")) +
                
                scale_linetype_manual(
                  name   = NULL,
                  breaks = c("OLS (bagged)", "OLS (uniform)", "Equal Weight"),
                  values = c("OLS (bagged)" = "solid", "OLS (uniform)" = "solid", "Equal Weight" = "dashed")
                ) +
                scale_color_manual(
                  name   = NULL,
                  breaks = c("OLS (bagged)", "OLS (uniform)", "Equal Weight"),
                  values = c("OLS (bagged)" = "#1f77b4", "OLS (uniform)" = "#d62728", "Equal Weight" = "grey40")
                ) +
                scale_x_continuous(
                  breaks = unique(.x$unit_idx),
                  labels = unit_labels[unique(.x$unit_idx)]
                ) +
                labs(title = paste("Evolution of", .y, "OLS coefficient across time"),
                     x = "Vintage", y = NULL) +
                theme_minimal() +
                theme(legend.position = "top",
                      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                      plot.title  = element_text(hjust = 0.5))
)
