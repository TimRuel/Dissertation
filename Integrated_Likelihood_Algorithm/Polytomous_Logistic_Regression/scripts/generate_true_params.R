# scripts/generate_true_params.R

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
experiment_id <- if (length(args) > 0) args[1] else stop("Missing experiment_id")

config_path <- proj_path("config", "exps", paste0(experiment_id, ".yml"))
if (!file.exists(config_path)) stop("[ERROR] Config file not found: ", config_path)

experiment_config <- read_yaml(config_path)
X1_levels <- experiment_config$X1_levels
model_specs <- experiment_config$model_specs

# Set seed
set.seed(experiment_config$seed)

# Generate and save true parameters
experiment_parameters <- get_experiment_parameters(X1_levels, model_specs)

theoretical_entropy_plot <- get_theoretical_entropy_plot(X1_levels, experiment_parameters$pY_0)

true_params_dir <- proj_path("experiments", experiment_id, "true_params")

dir_create(true_params_dir)
save_list_objects(experiment_parameters$true_params, true_params_dir)
save_list_plots(list(theoretical_entropy_plot), true_params_dir)

message("[INFO] Saved true parameters to ", true_params_dir)
