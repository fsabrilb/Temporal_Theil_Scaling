# Libraries ----
library(dplyr)
library(purrr)
library(cumstats)
library(data.table)

# Estimation of Temporal Fluctuation Scaling ----
estimate_tfs_data <- function(df_time_series){
  if(
    !"index" %in% names(df_time_series) | !"value" %in% names(df_time_series)
  ) {
    cat("Data frame doesn't have the correct names 'index' and 'value'\n")
    df_tfs <- data.table()
  } else {
    df_tfs <- df_time_series %>%
      filter(!is.na(value), value != 0) %>%
      # Estimate cumulative data, mean and variance
      mutate(
        cum_value = cumsum(value),
        cum_mean = cummean(value),
        cum_variance = cumvar(value),
        cum_mean_cum = cummean(cumsum(value)),
        cum_variance_cum = cumvar(cumsum(value))
      ) %>%
      filter(
        !is.na(cum_mean_cum),
        cum_mean_cum != 0,
        !is.na(cum_variance_cum),
        cum_variance_cum != 0
      )
  }
  
  return(df_tfs)
}

# Fitting to power law of cumulative variance and mean in time series (TFS) ----
adjust_tfs_power_law <- function(df_tfs) {
  # Estimate temporal fluctuation scaling as a power law of mean
  list_tfs_parameters <- lm(
    df_tfs %>% pull(cum_variance_cum) %>% log() ~
      df_tfs %>% pull(cum_mean_cum) %>% log()
  ) %>% summary()
  
  df_tfs_parameters <- list_tfs_parameters %>%
    pluck("coefficients") %>%
    t() %>%
    data.table() %>%
    rename_all(function(x) {c("tfs_coefficient", "tfs_exponent")}) %>%
    mutate(
      statistic = c("Estimate", "Standard error", "t value", "Pr(>|t|)")
    ) %>%
    relocate(statistic) %>%
    # Coefficient of determination R2
    bind_rows(
      data.table(
        "statistic" = "R2",
        "tfs_coefficient" = list_tfs_parameters %>% pluck("r.squared"),
        "tfs_exponent" = list_tfs_parameters %>% pluck("r.squared")
      )
    )
  
  return(df_tfs_parameters)
}

