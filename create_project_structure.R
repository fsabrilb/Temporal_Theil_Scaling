# Local path of the analysis project
local_path <- paste0(
  "G:/PC/academico/Universidad/Tesis de Doctorado Física/physreve_202207"
)
setwd(local_path)

# General structure of every project ----
dir.create("./input_files")
dir.create("./output_files")
dir.create("./logs")
dir.create("./modules")
dir.create("./scripts")

# General structure of input files ----
dir.create("./input_files/raw_data")
dir.create("./input_files/processed_data")
dir.create("./input_files/data_dictionary")
