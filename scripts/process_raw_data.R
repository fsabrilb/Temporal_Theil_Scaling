# -------------------------- TEMPORAL THEIL SCALING -------------------------- #
# Author: Felipe Segundo Abril Bermúdez
rm(list = ls())
input_path <- paste0(
  "G:/PC/academico/Universidad/Tesis de Doctorado Física/physreve_202207"
)

# Local path of modules for deployment ----
setwd(input_path)
source("./modules/estimate_theil_index.R")
source("./modules/estimate_diffusive_path.R")
source("./modules/estimate_diffusive_path_raw_data.R")

# Local path of data frames and global variables ----
input_path_raw <- "./input_files/raw_data"
input_path_processed <- "./input_files/processed_data"
input_path_data_dictionary <- "./input_files/data_dictionary"
output_path <- "./output_files"
input_generation_date <- "2023-04-12"

# Load raw data ----
df_raw_data <- fread(paste0(input_path_raw, "/df_raw_data_20230411.csv"))

# Setup parallel backend to use many processors (Turn on cluster) ----
cluster <- makeCluster(detectCores() - 1)
registerDoParallel(cluster)

# Estimate diffusive trajectory time series (DTTS) ----
df_dtts <- estimate_diffusive_algorithm_raw_data(
    df_raw_data = df_raw_data,
    t_min = 2,
    t_max = 500,
    path = paste0(input_path, "/logs"),
    verbose = 1,
    saved_data = TRUE,
    input_path_processed = input_path_processed,
    input_date = input_generation_date,
    number = 1
)

# Estimate Temporal Theil scaling data ----
df_dtts_data_1 <- estimate_dtts_raw_data(
  df_diffusive_path = df_dtts %>%
    filter(time_series %in% c("^GSPC", "^GDAXI", "^FCHI", "^IXIC")),
  path = paste0(input_path, "/logs"),
  verbose = 1
)

df_dtts_data_2 <- estimate_dtts_raw_data(
  df_diffusive_path = df_dtts %>%
    filter(time_series %in% c("IMOEX.ME", "^BVSP", "^N225")),
  path = paste0(input_path, "/logs"),
  verbose = 1
)

df_dtts_data_3 <- estimate_dtts_raw_data(
  df_diffusive_path = df_dtts %>%
    filter(time_series %in% c("GBPUSD=X", "EURUSD=X", "JPY=X", "COP=X")),
  path = paste0(input_path, "/logs"),
  verbose = 1
)

df_dtts_data_4 <- estimate_dtts_raw_data(
  df_diffusive_path = df_dtts %>%
    filter(
      time_series %in% c("COVIDWLDC", "COVIDWLDD", "PRCPBOG", "TEMPAVGBOG")
    ),
  path = paste0(input_path, "/logs"),
  verbose = 1
)

# Estimate Temporal Theil scaling adjustment ----
df_dtts_adjustment_1 <- adjust_tts_data_power_law(
  df_tts = df_dtts_data_1,
  path = paste0(input_path, "/logs"),
  verbose = 1
)

df_dtts_adjustment_2 <- adjust_tts_data_power_law(
  df_tts = df_dtts_data_2,
  path = paste0(input_path, "/logs"),
  verbose = 1
)

df_dtts_adjustment_3 <- adjust_tts_data_power_law(
  df_tts = df_dtts_data_3,
  path = paste0(input_path, "/logs"),
  verbose = 1
)

df_dtts_adjustment_4 <- adjust_tts_data_power_law(
  df_tts = df_dtts_data_4,
  path = paste0(input_path, "/logs"),
  verbose = 1
)

# Save final data for plotting (DTTS data) ----
df_dtts_data <- bind_rows(
  df_dtts_data_1,
  df_dtts_data_2,
  df_dtts_data_3,
  df_dtts_data_4
)

write.csv(
  df_dtts_data,
  paste0(
    input_path_processed,
    "/df_tts_data_final_",
    gsub("-", "", input_generation_date),
    ".csv"
  ),
  fileEncoding = "UTF-8",
  row.names = FALSE
)

# Save final data for plotting (DTTS adjustments) ----
df_dtts_adjustment <- bind_rows(
  df_dtts_adjustment_1,
  df_dtts_adjustment_2,
  df_dtts_adjustment_3,
  df_dtts_adjustment_4
)

df_tts_final <- df_dtts_adjustment %>%
  filter(statistic == "Estimate") %>%
  select(-statistic) %>%
  left_join(
    df_dtts_adjustment %>%
      filter(statistic == "Standard error") %>%
      rename(
        tts_coefficient_error = tts_coefficient,
        tts_exponent_error = tts_exponent
      ) %>%
      select(-statistic),
    by = c("time_series", "sub_time_series")
  ) %>%
  left_join(
    df_dtts_adjustment %>%
      filter(statistic == "R2") %>%
      select(c(time_series, sub_time_series, tts_coefficient)) %>%
      rename(rsquared = tts_coefficient),
    by = c("time_series", "sub_time_series")
  )

write.csv(
  df_tts_final,
  paste0(
    input_path_processed,
    "/df_tts_final_",
    gsub("-", "", input_generation_date),
    ".csv"
  ),
  fileEncoding = "UTF-8",
  row.names = FALSE
)

# Uninstall parallel backend to use many processors (Turn off cluster) ----
stopCluster(cluster)

# Clear workspace and free used memory ----
rm(list = ls())
gc()
