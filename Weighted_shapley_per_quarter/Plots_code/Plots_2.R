gy<-ffd[[103]][[2]]

k <- 1/16
prop_lt_k <- mean(gy$Mean_OLS < k, na.rm = TRUE)
# or, if you prefer the denominator explicit:
prop_lt_k <- sum(df$col < k, na.rm = TRUE) / sum(!is.na(df$col))

k<-c()
for(i in 164:266){
  j<-
 k[267-i]<-1/round(sqrt(i-2),0)
  
}
k<-rev(k)
k<-k[-18]
# Align growth_y
growth_y <- (y_levels$GDPC1_level / dplyr::lag(y_levels$GDPC1_level) - 1)*100
growth_y <- growth_y[-1]
growth_y <- as.data.frame(growth_y)
growth_y$Date <- y_levels$Date[-1]

k <- max(1L, round(sqrt(nrow(estimation))))




asd<-read.csv("/Users/boris/Desktop/RSShapley/Combination_Monte/FRED_QD/FRED-QD_2025m03.csv", stringsAsFactors = FALSE, check.names = FALSE)



# 0) Build a data frame for k aligned to your units/vintages
stopifnot(length(k) == length(ffd))
k_df <- tibble(
  unit_idx = seq_along(ffd),
  unit     = unit_labels,
  k        = as.numeric(k)
)

k_df <- tibble(
  unit_idx = seq_along(ffd),
  unit     = unit_labels,
  k        = as.numeric(nn)
)

# 1) Rebuild plots adding the k line (dashed) + legend
plots <- imap(df_list, ~
                ggplot(.x, aes(x = unit_idx)) +
                geom_line(aes(y = zero_bagged_ols, linetype = "zero_bagged_ols")) +
                geom_point(aes(y = zero_bagged_ols), size = 0.7) +
                geom_line(data = k_df, aes(y = k, linetype = "Equal Weight")) +
                scale_linetype_manual(values = c("zero_bagged_ols" = "solid", "Equal Weight" = "dashed"), name = NULL) +
                scale_x_continuous(
                  breaks = unique(.x$unit_idx),
                  labels = unit_labels[unique(.x$unit_idx)]
                ) +
                labs(
                  title = paste("Evolution of", .y, "OLS coefficient across time"),
                  x = "Unit / Vintage",
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
                geom_line(aes(y = Mean_OLS, linetype = "Mean_OLS")) +
                geom_point(aes(y = Mean_OLS), size = 0.7) +
                geom_line(data = k_df, aes(y = k, linetype = "Equal Weight")) +
                scale_linetype_manual(values = c("Mean_OLS" = "solid", "Equal Weight" = "dashed"), name = NULL) +
                scale_x_continuous(
                  breaks = unique(.x$unit_idx),
                  labels = unit_labels[unique(.x$unit_idx)]
                ) +
                labs(
                  title = paste("Evolution of", .y, "OLS coefficient across time"),
                  x = "Unit / Vintage",
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
                geom_line(aes(y = Mean_OLS, linetype = "Mean_OLS")) +
                geom_point(aes(y = Mean_OLS), size = 0.7) +
                geom_line(aes(y=Thr, linetype = "Equal Weight")) +
                scale_linetype_manual(values = c("Mean_OLS" = "solid", "Equal Weight" = "dashed"), name = NULL) +
                scale_x_continuous(
                  breaks = unique(.x$unit_idx),
                  labels = unit_labels[unique(.x$unit_idx)]
                ) +
                labs(
                  title = paste("Evolution of", .y, "OLS coefficient across time"),
                  x = "Unit / Vintage",
                  y = NULL
                ) +
                theme_minimal() +
                theme(
                  legend.position = "top",
                  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                  plot.title  = element_text(hjust = 0.5)
                )
)


out_dir <- "/Users/boris/Desktop/RSShapley/Combination_Monte/Final/codess/JMP/Weighted_shapley_per_quarter/OLS_coef_over_time_with_both_lines/"
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

