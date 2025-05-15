# scripts/run_model.R

#!/usr/bin/env Rscript

# ---- Optional interactive debug setup ----
if (interactive()) {
  message("Loading packages and helpers for interactive debugging...")
  library(tidyverse)
  library(yaml)
  library(fs)
  miceadds::source.all("scripts/helpers")
  
  # Set config path manually if debugging
  config_path <- "config/experiment_A.yml"
} else {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1) stop("Path to config .yml file must be provided.")
  config_path <- args[1]
}

# ---- Load config ----
experiment_config <- read_yaml(config_path)
experiment_id <- experiment_config$experiment_id

data_dir <- experiment_config$data_dir
output_base_dir <- experiment_config$output_base_dir %||% file.path("results", experiment_id)
run_type <- experiment_config$run_type %||% "individual"

# ---- Load data ----
required_data_files <- c("X_design.rds", "Y_probs.rds", "model_df.rds")
if (!all(file_exists(file.path(data_dir, required_data_files)))) {
  stop("Missing one or more required data files.")
}

data <- list(
  X_design = readRDS(file.path(data_dir, "X_design.rds")),
  Y_probs  = readRDS(file.path(data_dir, "Y_probs.rds")),
  model_df = readRDS(file.path(data_dir, "model_df.rds"))
)

# ---- Add or override runtime parameters ----
model_specs <- experiment_config$model_specs

experiment_config$threshold        <- experiment_config$threshold        %||% 1e-4
experiment_config$num_workers      <- experiment_config$num_workers      %||% 16


# ---- Run model ----
results <- run_model(model_df, model_specs)

# ---- Save results ----
run_id <- experiment_config$iter %||% format(Sys.time(), "%Y%m%d_%H%M%S")
output_dir <- file.path(output_base_dir, run_type, paste0("run_", run_id))
dir_create(output_dir)

saveRDS(results$log_integrated_likelihood, file.path(output_dir, "log_integrated_likelihood.rds"))
saveRDS(results$log_profile_likelihood, file.path(output_dir, "log_profile_likelihood.rds"))

# Optionally save updated config snapshot with runtime parameters
write_yaml(experiment_config, file.path(output_dir, "resolved_config.yml"))


