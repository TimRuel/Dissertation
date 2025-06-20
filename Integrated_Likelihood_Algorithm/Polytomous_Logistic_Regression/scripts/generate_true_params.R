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
if (!file.exists(config_path)) stop("[ERROR] Config file not found at /", sub(".*(/?config/.*)", "\\1", config_path))

experiment_config <- read_yaml(config_path)
X1_levels <- experiment_config$X1_levels
model_specs <- experiment_config$model_specs

# Set seed
set.seed(experiment_config$optimization_specs$seed)

# Generate and save true parameters
experiment_parameters <- get_experiment_parameters(X1_levels, model_specs)
true_params_dir <- proj_path("experiments", experiment_id, "true_params")
dir_create(true_params_dir)
save_list_objects(experiment_parameters$true_params, true_params_dir)
message("[INFO] Saved true parameters to ", sub(".*(/?experiments/.*)", "\\1", true_params_dir))

# Update config file with additional model specs
experiment_config$model_specs <- experiment_parameters$model_specs
write_strict_yaml(experiment_config, config_path)
message("[INFO] Updated experiment config with additional model specs.")

