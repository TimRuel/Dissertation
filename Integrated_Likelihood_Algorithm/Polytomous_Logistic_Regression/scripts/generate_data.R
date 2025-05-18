# scripts/generate_data.R

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(yaml)
  library(fs)
})

proj_subdir <- here("Integrated_Likelihood_Algorithm", "Polytomous_Logistic_Regression")
proj_path <- function(...) here(proj_subdir, ...)
miceadds::source.all(proj_path("scripts", "helpers"), print.source = FALSE)

args <- commandArgs(trailingOnly = TRUE)
experiment_id <- if (length(args) > 0) args[1] else "experiment_A"

config_path <- proj_path("config", paste0(experiment_id, ".yml"))
if (!file.exists(config_path)) stop("[ERROR] Config file not found: ", config_path)

experiment_config <- read_yaml(config_path)

set.seed(experiment_config$seed)

X1_levels <- experiment_config$X1_levels
model_specs <- experiment_config$model_specs

# Setup output directories
base_path <- proj_path("results", experiment_id)
true_params_dir <- here(base_path, "true_params")
data_dir <- here(base_path, "data")
plots_dir <- here(base_path, "plots")

dir_create(c(true_params_dir, data_dir, plots_dir))

# Check if data already exists
required_data_files <- c("X_design.rds", "Y_probs.rds", "model_df.rds")

if (all(file_exists(here(data_dir, required_data_files)))) {
  message("[INFO] Data already exists for ", experiment_id, " — skipping generation.")
} else {
  experiment_parameters <- get_experiment_parameters(X1_levels, model_specs)
  X1_levels <- experiment_parameters$X1_levels
  
  data <- get_data(X1_levels, model_specs$formula, experiment_parameters$true_params$Beta_0)
  
  plots <- get_plots(X1_levels, experiment_parameters$pY_0, data$Y_probs, data$model_df)
  
  save_list_objects(experiment_parameters$true_params, true_params_dir)
  save_list_objects(data, data_dir)
  save_list_plots(plots, plots_dir)
  
  # Save resolved config separately in results directory
  resolved_config <- list(
    experiment_id = experiment_id,
    seed = experiment_config$seed,
    X1_levels = X1_levels,
    model_specs = experiment_parameters$model_specs,
    optimization_specs = experiment_config$optimization_specs
  )
  
  resolved_path <- here(base_path, "resolved_config.yml")
  write_yaml(resolved_config, resolved_path)
  message("[INFO] Generated data and saved resolved config to ", resolved_path)
}
