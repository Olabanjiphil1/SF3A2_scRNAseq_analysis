suppressPackageStartupMessages({
  library(yaml)
})

config_file <- Sys.getenv("SF3A2_CONFIG", unset = "config/config.yml")
if (!file.exists(config_file)) {
  config_file <- "config/config_template.yml"
  message("config/config.yml not found. Using template config.")
}

config <- yaml::read_yaml(config_file)
source("R/helper_functions.R")

load_required_packages(c(
  "Seurat", "dplyr", "tidyr", "tibble", "ggplot2", "patchwork", "Matrix"
))

set.seed(as.integer(config$project$seed))
make_output_dirs(config)
