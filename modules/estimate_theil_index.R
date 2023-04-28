# Libraries ----
library(dplyr)
library(foreach)
library(data.table)
library(doParallel)

# Generation of exponential sample ----
generate_exponential_sample <- function(
  n_simulation,
  n_length,
  rate = 1
) {
  df_exponential_sample <- data.table()
  for(i in seq(n_simulation)) {
    df_exponential_sample <- df_exponential_sample %>%
      bind_rows(
        data.table(
          simulation = i,
          length = n_length,
          index = seq(n_length),
          value = rexp(n = n_length, rate = rate)
        )
      )
  }
  
  return(df_exponential_sample)
}

# Estimation of Generalized Entropy index for estimate Theil index ----
estimate_generalized_entropy_index <- function(df_time_series, alpha = 1){
  if(
    !"index" %in% names(df_time_series) | !"value" %in% names(df_time_series)
  ) {
    cat("Data frame doesn't have the correct names 'index' and 'value'\n")
    df_ge_index <- data.table()
  } else {
    df_ge_index <- df_time_series %>%
      # Estimate temporal variables (t_variable)
      mutate(
        t_value = value / mean(value, na.rm = TRUE),
        t_n = nrow(df_time_series),
        t_logn = log(nrow(df_time_series))
      ) %>%
      # Estimate Generalized Entropy Index (GE index)
      mutate(
        ge_index = case_when(
          alpha == 0 ~ -(1 / t_n) * sum(log(t_value), na.rm = TRUE),
          alpha == 1 ~ (1 / t_n) * sum(t_value * log(t_value), na.rm = TRUE),
          TRUE ~
            1/(t_n * alpha * (alpha - 1)) * sum(t_value^alpha - 1, na.rm = TRUE)
        )
      ) %>%
      # Normalization of Generalized entropy index with maximum value
      mutate(
        norm_sum_ge_index = case_when(
          alpha == 0 ~ -t_logn / t_n - ge_index,
          alpha == 1 ~ t_logn - ge_index,
          TRUE ~ (t_n^alpha - 1) / (t_n * alpha * (alpha - 1)) - ge_index
        ),
        norm_prod_ge_index = case_when(
          alpha == 0 ~ ge_index * (-t_n / t_logn),
          alpha == 1 ~ ge_index * (1/t_logn),
          TRUE ~ ge_index * ((t_n * alpha * (alpha - 1)) / (t_n^alpha - 1))
        )
      ) %>%
      # Delete temporal variables
      select(-c(t_value, t_n, t_logn))
  }
  
  return(df_ge_index)
}

# Estimation of Theil index ----
estimate_theil_index <- function(df_time_series) {
  df_theil_index <- estimate_generalized_entropy_index(
    df_time_series = df_time_series,
    alpha = 1
  ) %>%
    rename(
      theil_index = ge_index,
      norm_sum_theil_index = norm_sum_ge_index,
      norm_prod_theil_index = norm_prod_ge_index
    )
  
  return(df_theil_index)
}

# Estimation of Atkinson index ----
estimate_atkinson_index <- function(df_time_series, epsilon) {
  df_atkinson_index <- estimate_generalized_entropy_index(
    df_time_series = df_time_series,
    alpha = 1 - epsilon
  ) %>%
    # Estimate temporal variables (t_variable)
    mutate(
      t_a = (epsilon * (epsilon - 1)) ^ (1 / (1 - epsilon)),
      t_b = 1 / (1 - epsilon),
      t_n = nrow(df_time_series),
      t_logn = log(nrow(df_time_series))
    ) %>%
    # Estimate Generalized Atkinson Index (A index)
    mutate(
      atkinson_index = case_when(
        epsilon == 1 ~ 1 - exp(-ge_index),
        TRUE ~ t_a * ge_index^t_b
      )
    ) %>%
    # Normalization of Generalized entropy index with maximum value
    mutate(
      norm_sum_atkinson_index = case_when(
        epsilon == 1 ~ 1 - t_n^(1 / t_n) - atkinson_index,
        epsilon == 0 ~ t_a * t_logn^t_b - atkinson_index,
        TRUE ~ ((1 - t_n^(epsilon - 1)) / t_n^epsilon)^t_b - atkinson_index
      ),
      norm_prod_atkinson_index = case_when(
        epsilon == 1 ~ atkinson_index / (1 - t_n^(1 / t_n)),
        epsilon == 0 ~ atkinson_index / (t_a * t_logn^t_b),
        TRUE ~ atkinson_index / (((1 - t_n^(epsilon - 1)) / t_n^epsilon)^t_b)
      )
    ) %>%
    # Delete temporal variables
    select(
      -c(t_a, t_b, t_n, t_logn, ge_index, norm_sum_ge_index, norm_prod_ge_index)
    )
  
  return(df_atkinson_index)
}

# Estimate Theil index of exponential sample ----
estimate_theil_exponential_sample <- function(
  df_exponential_sample,
  path = "/home/ASIS/Temp_Felipe",
  verbose = 1
) {
  # Selection of simulation
  loop_index <- df_exponential_sample %>% distinct(simulation) %>% pull()
  
  # Parallel loop over loop_index for estimate Theil in every simulation
  df_exponential_sample<- foreach(
    i = loop_index,
    .combine = bind_rows,
    .export = c(
      "estimate_theil_index", "estimate_generalized_entropy_index"
    ),
    .packages = c("dplyr", "purrr", "foreach", "data.table", "doParallel")
  ) %dopar% {
    df_temp <- estimate_theil_index(
      df_time_series = df_exponential_sample %>% filter(simulation == i)
    ) %>% mutate(simulation = i)
    
    # Log output for monitoring progress
    if(verbose >= 1){
      cat(
        paste0("Estimated Theil for simulation of exponential: ", i, "\n"),
        file = paste0(path, "/log_theil_exponential.txt"),
        append = TRUE
      )
    }
    
    df_temp
  }
  
  return(df_exponential_sample)
}

# Estimate excess in the estimation of Theil in exponential sample ----
estimate_excess_theil_exponential_sample <- function(
  n_simulation_vector,
  n_length_vector,
  rate = 1,
  path = "/home/ASIS/Temp_Felipe",
  verbose = 1,
  saved_all_data = FALSE,
  input_path_processed = "./input_files/processed_data",
  input_date = "2022-12-04"
) {
  if(saved_all_data == FALSE) {
    # Selection of loop size
    loop_index <- seq(length(n_simulation_vector))
    
    # Parallel loop over loop_index for estimate Theil's excess for every length
    df_theil_exponential_sample <- foreach(
      i = loop_index,
      .combine = bind_rows,
      .export = c(
        "estimate_theil_index",
        "estimate_generalized_entropy_index",
        "estimate_theil_exponential_sample",
        "generate_exponential_sample"
      ),
      .packages = c("dplyr", "purrr", "foreach", "data.table", "doParallel")
    ) %dopar% {
      # Generate Exponential Sample
      df_exponential_sample <- generate_exponential_sample(
        n_simulation = n_simulation_vector[i],
        n_length = n_length_vector[i],
        rate = rate
      )
      
      # Log output for monitoring progress
      if(verbose >= 1){
        cat(
          paste0(
            "------------ Length of exponential sample: ",
            n_length_vector[i],
            " ------------\n"
          ),
          file = paste0(path, "/log_theil_optimum.txt"),
          append = TRUE
        )
      }
      
      # Temporal dataframe of Theil index
      df_temp <- estimate_theil_exponential_sample(
        df_exponential_sample = df_exponential_sample,
        path = path,
        verbose = verbose
      ) %>% mutate(length = n_length_vector[i])
      
      df_temp
    }
    
    # Save all Theil data in processed data for not reprocess the information
    write.csv(
      df_theil_exponential_sample,
      paste0(
        input_path_processed,
        "/df_theil_excess_data_",
        gsub("-", "", input_date),
        ".csv"
      ),
      fileEncoding = "UTF-8",
      row.names = FALSE
    )
  } else {
    # Load all data from processed data for not reprocess the information
    df_theil_exponential_sample <- fread(
      paste0(
        input_path_processed,
        "/df_theil_excess_data_",
        gsub("-", "", input_date),
        ".csv"
      )
    )
  }
  
  return(df_theil_exponential_sample)
}
