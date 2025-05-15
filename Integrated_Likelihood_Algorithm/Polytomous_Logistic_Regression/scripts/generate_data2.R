# scripts/generate_data.R

#!/usr/bin/env Rscript

if (interactive()) {
  message("Loading packages and helpers for interactive debugging...")
  library(tidyverse)
  library(yaml)
  library(fs)
  miceadds::source.all("scripts/helpers")
}

args <- commandArgs(trailingOnly = TRUE)
experiment_id <- if (length(args) > 0) args[1] else "experiment_A"

config_path <- file.path("config", paste0(experiment_id, ".yml"))
if (!file.exists(config_path)) stop("[ERROR] Config file not found: ", config_path)

experiment_config <- read_yaml(config_path)

set.seed(experiment_config$seed)

X1_levels <- experiment_config$X1_levels
ep_specs <- experiment_config$ep_specs

# Setup output directories
base_path <- file.path("results", experiment_id)
true_params_dir <- file.path(base_path, "true_params")
data_dir <- file.path(base_path, "data")
plots_dir <- file.path(base_path, "plots")

dir_create(c(true_params_dir, data_dir, plots_dir))

# Check if data already exists
data_exists <- c(file.path(data_dir, "X_design.rds"),
                 file.path(data_dir, "Y_probs.rds"),
                 file.path(data_dir, "model_df.rds")) |> 
  file_exists() |> 
  all()

if (data_exists) {
  message("[INFO] Data already exists for ", experiment_id, " — skipping generation.")
} else {
  experiment_parameters <- get_experiment_parameters(X1_levels, ep_specs)
  X1_levels <- experiment_parameters$X1_levels
  
  data <- get_data(X1_levels, ep_specs$formula, experiment_parameters$true_params$Beta_0)
  
  plots <- get_plots(X1_levels, experiment_parameters$pY_0, data$Y_probs, data$model_df)
  
  save_list_objects(experiment_parameters$true_params, true_params_dir)
  save_list_objects(data, data_dir)
  save_list_plots(plots, plots_dir)
  
  # Save resolved config separately in results directory
  resolved_config <- list(
    experiment_id = experiment_id,
    seed = experiment_config$seed,
    X1_levels = X1_levels,
    specs = experiment_parameters$specs,
    true_params_dir = true_params_dir,
    data_dir = data_dir,
    plots_dir = plots_dir
  )
  
  resolved_path <- file.path("results", experiment_id, "resolved_config.yml")
  write_yaml(resolved_config, resolved_path)
  message("[INFO] Generated data and saved resolved config to ", resolved_path)
}
