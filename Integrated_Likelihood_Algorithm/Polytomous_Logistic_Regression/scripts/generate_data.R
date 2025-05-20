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

# --- Parse arguments ---
args <- commandArgs(trailingOnly = TRUE)
experiment_id <- if (length(args) > 0) args[1] else stop("[ERROR] Missing experiment_id")
mode <- if (length(args) > 1) args[2] else stop("[ERROR] Missing experiment mode")
run_id <- if (length(args) > 2) args[3] else stop("[ERROR] Missing run id")

# --- Load config ---
config_path <- proj_path("config", "exps", paste0(experiment_id, ".yml"))
if (!file.exists(config_path)) stop("[ERROR] Config file not found: ", config_path)
experiment_config <- read_yaml(config_path)

X1_levels <- experiment_config$X1_levels
model_specs <- experiment_config$model_specs

# --- Setup directories ---
true_params_dir <- proj_path("experiments", experiment_id, "true_params")
run_dir <- proj_path("experiments", experiment_id, mode, run_id)
data_dir <- here(run_dir, "data")
plots_dir <- here(run_dir, "plots")
config_snapshot_path <- here(run_dir, "config_snapshot.yml")

dir_create(c(run_dir, data_dir, plots_dir))

# --- Step 1: Load Beta_0 ---
Beta_0_path <- here(true_params_dir, "Beta_0.rds")
if (file_exists(Beta_0_path)) {
  message("[INFO] Loading Beta_0 from: ", true_params_dir)
  Beta_0 <- readRDS(Beta_0_path)
} else {
  stop("[ERROR] Beta_0.rds not found in: ", true_params_dir)
}

# --- Step 2: Generate data and plots (always) ---
message("[INFO] Generating new data and plots for run: ", run_id)
seed <- get_seed_for_run(experiment_config$seed, run_id)
set.seed(seed)
data <- get_data(X1_levels, model_specs$formula, Beta_0)
plots <- get_observed_plots(X1_levels, data$Y_probs, data$model_df)

save_list_objects(data, data_dir)
save_list_plots(plots, plots_dir)

# --- Step 3: Write resolved config for the run ---
config_snapshot <- experiment_config
config_snapshot$optimizaton_specs <- list(seed = seed,
                                          mode = mode,
                                          run_id = run_id)

write_yaml(config_snapshot, config_snapshot_path)
message("[INFO] Saved config snapshot to ", config_snapshot)
