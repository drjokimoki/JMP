
#clear the environment
rm(list=ls())

library("midasr")
library("readxl")
library("zoo")
library("xml2")
library("rvest")
library("tidyr")
library("tsfgrnn")
library("forecast")
library("dplyr")
library("randomForest")
library("e1071")
library("urca")
library("purrr")
library("tidyr")
library("dplyr")
library("tsfknn")
library("GMDH")
library("data.table")
library("multDM")
library(xts)
library(tidyverse)
library(magrittr)
#library("aTSA")
library(tseries)
library(forecast)
library("dfms")
library(dplyr)
library(tidyr)
library(lubridate)
library(zoo)
library(tibble)


GDP<-read_excel("/Users/boris/Downloads/OH_GRNN_Data.xlsx", sheet = 2)
IND<-read_excel("/Users/boris/Downloads/OH_GRNN_Data.xlsx", sheet = 1)



# keep only columns with at least one non-NA
IND <- IND[ , colSums(!is.na(IND)) > 0, drop = FALSE]


apply_transforms <- function(df, skip_first_col = TRUE, drop_code0 = FALSE) {
  stopifnot(nrow(df) >= 2)
  
  start_col <- if (isTRUE(skip_first_col)) 2L else 1L
  code_row  <- 1L
  
  # read codes from first row (excluding first column if requested)
  codes <- as.integer(as.vector(unlist(df[code_row, start_col:ncol(df)], use.names = FALSE)))
  names(codes) <- names(df)[start_col:ncol(df)]
  
  # data (remove the code row)
  out <- df[-code_row, , drop = FALSE]
  
  # keep/pass-through first column (e.g., date/ID)
  keep_first <- if (isTRUE(skip_first_col)) out[, 1, drop = FALSE] else NULL
  
  # transformer
  transform_series <- function(x, code) {
    x <- suppressWarnings(as.numeric(x))
    if (all(is.na(x))) return(x)
    
    if (code %in% c(0L, 4L)) {
      # level / no transform
      return(x)
    }
    
    if (code == 1L) {
      y <- diff(x)
      return(c(NA_real_, y))
    }
    
    if (code %in% c(2L, 3L)) {
      # log transforms need positive values
      x_pos <- x
      x_pos[x_pos <= 0] <- NA_real_
      lx <- log(x_pos)
      if (code == 2L) {
        y <- diff(lx)
        return(c(NA_real_, y))
      } else {
        y <- diff(diff(lx))
        return(c(NA_real_, NA_real_, y))
      }
    }
    
    stop("Unknown code: ", code)
  }
  
  # apply per column
  transformed_list <- list()
  for (j in seq_along(codes)) {
    col_name <- names(codes)[j]
    code     <- codes[j]
    
    # optionally drop code==0 columns
    if (drop_code0 && code == 0L) next
    
    res <- transform_series(out[[col_name]], code)
    transformed_list[[col_name]] <- res
  }
  
  # bind
  if (length(transformed_list) == 0L) {
    res_df <- if (is.null(keep_first)) out[FALSE, , drop = FALSE] else keep_first
  } else {
    trans_df <- as.data.frame(transformed_list, check.names = FALSE, stringsAsFactors = FALSE)
    res_df <- if (is.null(keep_first)) trans_df else cbind(keep_first, trans_df, stringsAsFactors = FALSE)
  }
  
  res_df
}

df<-apply_transforms(IND)
count_last_k_na <- function(x, k = 2) colSums(is.na(tail(x, k)))

#The number of quarter to use as observations.
numbers<-3-count_last_k_na(IND[,-1], 2)



#GDP growth rate
growth_rate <- function(x)(x/lag(x)-1)*100

GDP_GR<-growth_rate(GDP$GDP)[-1]


IND$Date<-as.Date(as.character(as.POSIXct(IND$Date)))

IND_upd<-IND %>% filter(row_number() <= n()-2)
IND_upd$Date<-as.Date(as.character(as.POSIXct(IND_upd$Date)))
IND_upd[, 2:ncol(IND_upd)] <- sapply(IND_upd[, 2:ncol(IND_upd)], as.numeric)
df<-IND_upd


# df: Date (monthly, first-of-month) + one column per variable
# numbers: either a named vector c(varA=1, varB=2, ...) or a 2-col data frame (var, k)



# numbers: named vector or 2-col data.frame (var, k)
numbers_tbl <- {
  if (is.atomic(numbers) && !is.null(names(numbers))) {
    tibble::enframe(numbers, name = "var", value = "k")
  } else if (is.data.frame(numbers)) {
    stats::setNames(numbers, c("var","k"))
  } else stop("`numbers` must be a named vector or a 2-column data.frame.")
} %>% dplyr::mutate(k = pmin(as.integer(k), 3L))

quarterly <- df %>%
  dplyr::arrange(Date) %>%
  dplyr::mutate(
    qtr   = zoo::as.yearqtr(Date),
    m_in_q = ((lubridate::month(Date) - 1L) %% 3L) + 1L
  ) %>%
  tidyr::pivot_longer(
    -c(Date, qtr, m_in_q),
    names_to = "var", values_to = "val",
    values_transform = list(val = as.numeric)
  ) %>%
  dplyr::left_join(numbers_tbl, by = "var") %>%
  dplyr::mutate(k = dplyr::coalesce(k, 3L)) %>%
  dplyr::filter(m_in_q <= k) %>%
  dplyr::group_by(qtr, var) %>%
  dplyr::summarise(
    value = if (all(is.na(val))) NA_real_ else mean(val, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(names_from = var, values_from = value) %>%
  dplyr::mutate(Date = as.Date(zoo::as.yearmon(qtr))) %>%
  dplyr::select(-qtr) %>%                     # drop helper column
  dplyr::relocate(Date, .before = 1)          # put Date first


growth_rate_safe <- function(x) {
  prev <- dplyr::lag(x)
  out <- (x / prev - 1) * 100
  out[is.na(prev) | prev == 0] <- NA_real_
  out
}

quarterly_gr <- quarterly %>%
  arrange(Date) %>%
  mutate(across(-Date, growth_rate_safe)) %>%
  slice(-1)

GDP_GR_ac<-GDP_GR[-length(GDP_GR)]

quarterly_gr_a<-quarterly_gr[ , colSums(is.na(quarterly_gr)) == 0]

(GDP_GR_ac[68] + GDP_GR_ac[11] + GDP_GR_ac[5] + GDP_GR_ac[117] + GDP_GR_ac[20] + GDP_GR_ac[127] +  GDP_GR_ac[40])/7
preds<-list()
preds_gdp<-c()
ii=2
for(ii in 2:ncol(quarterly_gr_a)){
  la<-as.numeric(unlist(quarterly_gr_a[,ii]))
  aag<-knn_forecasting(la, transform = "none", h=1, lags = 4, k=1)
  preds[[ii-1]]<-aag[["weights"]]
  preds_gdp[ii-1]<-sum(tail(GDP_GR_ac, dim(preds[[ii-1]])[2])*preds[[ii-1]])
}

?knn_forecasting


preds<-c()
ii=4
for(ii in 2:ncol(quarterly_gr_a)){
  la<-as.numeric(unlist(quarterly_gr_a[,ii]))
  aag<-grnn_forecasting(la, lags = 1, transform = "none", h=1)
  preds[ii-1]<-sum(aag[["weights"]]*GDP_GR_ac)
}

DIv <- function(x) {
  ok = x / GDP_GR_ac
  return(ok)
}

quarterly_gr_a_fordi<-quarterly_gr_a[-nrow(quarterly_gr_a),]

Fin_data <- as.data.frame(sapply(quarterly_gr_a[-nrow(quarterly_gr_a), -1], DIv))


quarterly_one<-quarterly[-1,]
quarterly_ones<-quarterly_one[ , colSums(is.na(quarterly_one)) == 0]

preds<-list()
preds_gdp<-c()
#ii=2
for(ii in 2:ncol(quarterly_ones)){
  la<-as.numeric(unlist(quarterly_ones[,ii]))
  aag<-grnn_forecasting(la, lags = "FS" , transform = "additive", h=1)
  preds[[ii-1]]<-aag[["weights"]]
  preds_gdp[ii-1]<-sum(tail(GDP_GR_ac, dim(preds[[ii-1]])[2])*preds[[ii-1]])
}
mean(preds_gdp)


?grnn_forecasting
sum(tail(GDP_GR_ac, dim(preds[[2]])[2])*preds[[2]])











GDP_GR

quarterly_ones_s<-quarterly_ones[-nrow(quarterly_ones),]
GDP_GR_a<-GDP_GR_ac[-length(GDP_GR_ac)]
preds<-list()
preds_gdp<-c()

for(ii in 2:ncol(quarterly_gr_a)){
  la<-as.numeric(unlist(quarterly_gr_a[,ii]))
  aag<-knn_forecasting(la, lags = 1:4 , transform = "none", h=1, k=5)
  preds[[ii-1]]<-as.numeric(aag[["neighbors"]])
  qas<-as.numeric(aag[["neighbors"]])
  preds_gdp[ii-1]<-mean(GDP_GR_ac[qas], na.rm = T)
}
mean(preds_gdp, na.rm = T)




qas<-as.numeric(aag[["neighbors"]])
mean(GDP_GR_ac[qas], na.rm = T)
preds<-list()
preds_gdp<-c()
ii=9
for(ii in 2:ncol(quarterly_ones)){
  la<-as.numeric(unlist(quarterly_ones[,ii]))
  aag<-knn_forecasting(la, lags = 1:4 , transform = "none", h=1, k=1)
  preds[[ii-1]]<-as.numeric(aag[["neighbors"]])
  qas<-as.numeric(aag[["neighbors"]])
  preds_gdp[ii-1]<-mean(GDP_GR_ac[qas], na.rm = T)
}

mean(preds_gdp, na.rm = T)
