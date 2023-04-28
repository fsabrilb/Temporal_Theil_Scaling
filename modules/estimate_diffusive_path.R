# Libraries ----
library(dplyr)
library(foreach)
library(data.table)
library(doParallel)

# Generate diffusive paths through diffusive algorithm ----
estimate_diffusive_algorithm <- function(
  df_time_series,
  t_min = 2,
  t_max = 10,
  path = "/home/ASIS/Temp_Felipe",
  verbose = 1
){
  if(
    !"index" %in% names(df_time_series) | !"value" %in% names(df_time_series)
  ) {
    cat("Data frame doesn't have the correct names 'index' and 'value'\n")
    df_diffusive_path <- data.table()
  } else {
    # Generation of log-return, absolute log-return and log-volatility
    df_time_series <- df_time_series %>%
      arrange(index) %>%
      mutate(log_return = log(value) - log(lag(value))) %>%
      filter(
        !is.na(log_return),
        !is.nan(log_return),
        !is.infinite(log_return)
      ) %>%
      mutate(
        log_absolute = abs(log_return),
        log_volatility = sqrt(
          abs(log_return - mean(log_return, na.rm = TRUE)) /
            sd(log_return, na.rm = TRUE)
        ) / max(abs(log_return), na.rm = TRUE)
      )
      

    # Parallel loop over time step for construction of diffusive path at time t
    df_diffusive_path <- foreach(
      i = seq(t_min, t_max, by = 1),
      .combine = bind_rows,
      .packages = c("dplyr", "data.table")
    ) %dopar% {
      # Log output for monitoring progress
      if(verbose >= 1){
        cat(
          paste0("Estimated diffusive path: ", i, "\n"),
          file = paste0(path, "/log_diffusive_path.txt"),
          append = TRUE
        )
      }
      
      # Diffusive algorithm sub-series
      df_temp <- data.table()
      max_index <- df_time_series %>% pull(index) %>% max(na.rm = TRUE)
      for(j in c(1:i)) {
        df_temp <- bind_rows(
          df_temp,
          data.table(
            time = max_index - i + 1,
            sub_time = j - 1,
            diffusive_value = df_time_series %>%
              filter(index >= j, index <= max_index - i + j) %>%
              pull(value) %>%
              sum(na.rm = TRUE),
            diffusive_log_absolute = df_time_series %>%
              filter(index >= j, index <= max_index - i + j) %>%
              pull(log_absolute) %>%
              sum(na.rm = TRUE),
            diffusive_log_volatility = df_time_series %>%
              filter(index >= j, index <= max_index - i + j) %>%
              pull(log_volatility) %>%
              sum(na.rm = TRUE)
          )
        )
      }
      
      df_temp
    }
  }
  
  return(df_diffusive_path)
}

# Estimation of Theil index over diffusive paths ----
estimate_theil_diffusive_path <- function(df_diffusive_path) {
  df_theil_index <- data.table()
  for(i in df_diffusive_path %>% distinct(time) %>% pull()) {
    df_aux <- bind_rows(
      # Original data
      estimate_theil_index(
        df_time_series = df_diffusive_path %>%
          filter(time == i) %>%
          rename(index = sub_time, value = diffusive_value)
      ) %>%
        mutate(time_series = "original") %>%
        select(-c("diffusive_log_absolute", "diffusive_log_volatility")),
      # Absolute log-return
      estimate_theil_index(
        df_time_series = df_diffusive_path %>%
          filter(time == i) %>%
          rename(index = sub_time, value = diffusive_log_absolute)
      ) %>%
        mutate(time_series = "absolute log-return") %>%
        select(-c("diffusive_value", "diffusive_log_volatility")),
      # Volatility log-return
      estimate_theil_index(
        df_time_series = df_diffusive_path %>%
          filter(time == i) %>%
          rename(index = sub_time, value = diffusive_log_volatility)
      ) %>%
        mutate(time_series = "volatility log-return") %>%
        select(-c("diffusive_value", "diffusive_log_absolute"))
    )
    
    if(nrow(df_theil_index) == 0) {
      df_theil_index <- df_aux
    } else {
      df_theil_index <- df_theil_index %>% bind_rows(df_aux)
    }
  }
  
  df_theil_index <- df_theil_index %>% arrange(time_series, time, index)
  
  return(df_theil_index)
}

# Estimation of Temporal Theil scaling (TTS) ----
estimate_tts_data <- function(df_diffusive_path) {
  df_tts <- estimate_theil_diffusive_path(
    df_diffusive_path = df_diffusive_path
  ) %>%
    group_by(time_series, time) %>%
    summarise(
      mean_value = mean(value, na.rm = TRUE),
      theil_index = mean(theil_index, na.rm = TRUE),
      norm_sum_theil_index = mean(norm_sum_theil_index, na.rm = TRUE),
      norm_prod_theil_index = mean(norm_prod_theil_index, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    rename(index = time)
  
  return(df_tts)
}

# Fitting to power law of Temporal Theil scaling in time series (TFS) ----
adjust_tts_power_law <- function(df_tts) {
  df_tts <- df_tts %>%
    group_by(time_series) %>%
    mutate(
      mean_tts = 1 - mean_value / max(abs(mean_value), na.rm = TRUE),
      theil_tts =
        norm_sum_theil_index / max(abs(norm_sum_theil_index), na.rm = TRUE)
    ) %>%
    ungroup()
    
  df_tts_parameters <- data.table()
  for(i in df_tts %>% distinct(time_series) %>% pull()) {
    df_aux <- df_tts %>% filter(time_series == i, mean_tts > 0, theil_tts > 0)
    if(nrow(df_aux) > 0) {
      # Estimate temporal Theil scaling as a power law of mean
      list_tts_parameters <- lm(
        df_aux %>% pull(theil_tts) %>% abs() %>% log() ~
          df_aux %>% pull(mean_tts) %>% abs() %>% log()
      ) %>% summary() 
      
      df_tts_parameters_aux <- list_tts_parameters %>%
        pluck("coefficients") %>%
        t() %>%
        data.table() %>%
        rename_all(function(x) {c("tts_coefficient", "tts_exponent")}) %>%
        mutate(
          statistic = c("Estimate", "Standard error", "t value", "Pr(>|t|)")
        ) %>%
        relocate(statistic) %>%
        # Coefficient of determination R2
        bind_rows(
          data.table(
            "statistic" = "R2",
            "tts_coefficient" = list_tts_parameters %>% pluck("r.squared"),
            "tts_exponent" = list_tts_parameters %>% pluck("r.squared")
          )
        ) %>%
        mutate(time_series = i)
      
    } else {
      df_tts_parameters_aux <- data.table()
    }
    
    df_tts_parameters <- df_tts_parameters %>% bind_rows(df_tts_parameters_aux)
    
  }
  
  return(df_tts_parameters)
}
