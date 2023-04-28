# Libraries ----
library(dplyr)
library(purrr)
library(foreach)
library(data.table)
library(doParallel)

# Generate diffusive trajectory time series ----
estimate_diffusive_algorithm_raw_data <- function(
  df_raw_data,
  t_min = 2,
  t_max = 10,
  path = "/home/ASIS/Temp_Felipe",
  verbose = 1,
  saved_data = FALSE,
  input_path_processed = "./input_files/processed_data",
  input_date = "2022-12-04",
  number = 1
){
  if(saved_data == FALSE) {
    loop_index <- df_raw_data %>% distinct(label) %>% pull()
    
    # Parallel loop over raw data labels
    df_dtts <- foreach(
      i = loop_index,
      .combine = bind_rows,
      .export = c("estimate_diffusive_algorithm"),
      .packages = c("dplyr", "foreach", "data.table", "doParallel")
    ) %dopar% {
      # Log output for monitoring progress ----
      if(verbose >= 1){
        cat(
          paste0("Estimated dtts for label ", i, "\n"),
          file = paste0(path, "/log_dtts.txt"),
          append = TRUE
        )
      }
      
      # Estimate diffusive trajectory time series of label i ----
      df_temp <- estimate_diffusive_algorithm(
        df_time_series = df_raw_data %>% filter(label == i),
        t_min = t_min,
        t_max = t_max,
        path = path,
        verbose = verbose
      ) %>% mutate(time_series = i) %>% relocate(time_series)
      
      df_temp
    }
    
    # Save data in processed data for not reprocess the information ----
    write.csv(
      df_dtts,
      paste0(
        input_path_processed,
        "/df_dtts_data_",
        gsub("-", "", input_date),
        "_",
        number,
        ".csv"
      ),
      fileEncoding = "UTF-8",
      row.names = FALSE
    )
  } else {
    df_dtts <- fread(
      paste0(
        input_path_processed,
        "/df_dtts_data_",
        gsub("-", "", input_date),
        "_",
        number,
        ".csv"
      ),
      encoding = "UTF-8"
    )
  }
  
  return(df_dtts)
}

# Estimate temporal Theil scaling data ----
estimate_dtts_raw_data <- function(
  df_diffusive_path,
  path = "/home/ASIS/Temp_Felipe",
  verbose = 1
) {
  df_theil_index <- data.table()
  for(i in df_diffusive_path %>% distinct(time_series) %>% pull()) {
    df_diffusive_path_local <- df_diffusive_path %>% filter(time_series == i)
    for(j in df_diffusive_path_local %>% distinct(time) %>% pull()) {
      # Log output for monitoring progress ----
      if(verbose >= 1){
        cat(
          paste0(
            "TTS data for label ",
            i,
            " and time ",
            j,
            " columns ",
            df_theil_index %>% nrow(),
            "\n"
          ),
          file = paste0(path, "/log_tts_data.txt"),
          append = TRUE
        )
      }
      
      # Estimation of Theil index and mean per sub-series ----
      df_aux <- bind_rows(
        # Original data
        estimate_theil_index(
          df_time_series = df_diffusive_path_local %>%
            filter(time == j) %>%
            rename(index = sub_time, value = diffusive_value)
        ) %>%
          mutate(sub_time_series = "original") %>%
          select(-c("diffusive_log_absolute", "diffusive_log_volatility")),
        # Absolute log-return
        estimate_theil_index(
          df_time_series = df_diffusive_path_local %>%
            filter(time == j) %>%
            rename(index = sub_time, value = diffusive_log_absolute)
        ) %>%
          mutate(sub_time_series = "absolute log-return") %>%
          select(-c("diffusive_value", "diffusive_log_volatility")),
        # Volatility log-return
        estimate_theil_index(
          df_time_series = df_diffusive_path_local %>%
            filter(time == j) %>%
            rename(index = sub_time, value = diffusive_log_volatility)
        ) %>%
          mutate(sub_time_series = "volatility log-return") %>%
          select(-c("diffusive_value", "diffusive_log_absolute"))
      ) %>%
        group_by(time_series, sub_time_series, time) %>%
        summarise(
          mean_value = mean(value, na.rm = TRUE),
          theil_index = mean(theil_index, na.rm = TRUE),
          norm_sum_theil_index = mean(norm_sum_theil_index, na.rm = TRUE),
          norm_prod_theil_index = mean(norm_prod_theil_index, na.rm = TRUE)
        ) %>%
        ungroup() %>%
        rename(index = time)

      df_theil_index <- df_theil_index %>% bind_rows(df_aux)
    }
  }
  
  return(df_theil_index)
}

# Fitting to power law of Temporal Theil scaling in time series (TTS) ----
adjust_tts_data_power_law <- function(
  df_tts,
  path = "/home/ASIS/Temp_Felipe",
  verbose = 1
) {
  df_tts <- df_tts %>%
    group_by(time_series, sub_time_series) %>%
    mutate(
      mean_tts = 1 - mean_value / max(abs(mean_value), na.rm = TRUE),
      theil_tts =
        norm_sum_theil_index / max(abs(norm_sum_theil_index), na.rm = TRUE)
    ) %>%
    ungroup()
  
  df_tts_parameters <- data.table()
  for(i in df_tts %>% distinct(time_series) %>% pull()) {
    df_tts_local <- df_tts %>% filter(time_series == i)
    for(j in df_tts_local %>% distinct(sub_time_series) %>% pull()) {
      # Log output for monitoring progress ----
      if(verbose >= 1){
        cat(
          paste0("TTS adjust for label ", i, " and sub-series ", j, "\n"),
          file = paste0(path, "/log_tts_final.txt"),
          append = TRUE
        )
      }

      df_aux <- df_tts_local %>%
        filter(sub_time_series == j, mean_tts > 0, theil_tts > 0)
      
      if(nrow(df_aux) > 0) {
        # Estimate temporal Theil scaling as a power law of mean ----
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
          mutate(time_series = i, sub_time_series = j)
        
      } else {
        df_tts_parameters_aux <- data.table()
      }
      
      df_tts_parameters <- df_tts_parameters %>%
        bind_rows(df_tts_parameters_aux)
      
    }
  }

  return(df_tts_parameters)
}



