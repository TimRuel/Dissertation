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

# ---- Parse Arguments ----
args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) > 0) args[1] else stop("[ERROR] A config filepath was not provided.")
if (!file.exists(config_path)) stop("[ERROR] Config file not found: ", config_path)

# ---- Load Config ----
config <- read_yaml(config_path)
experiment_id <- config$experiment_id
output_dir <- config$output_dir
run_type <- config$optimization_specs$run_type %||% "individual"

# ---- Load Data ----
data_dir <- proj_path("results", experiment_id, "data")
required_files <- c("X_design.rds", "Y_probs.rds", "model_df.rds")
missing <- required_files[!file_exists(here(data_dir, required_files))]
if (length(missing) > 0) stop("Missing required data files: ", paste(missing, collapse = ", "))

X_design <- readRDS(here(data_dir, "X_design.rds"))
model_df <- readRDS(here(data_dir, "model_df.rds"))

# ---- Run Experiment ----
results <- run_experiment(config, X_design, model_df)

# ---- Save Results ----
saveRDS(results$log_integrated_likelihood, file = here(output_dir, "log_integrated_likelihood.rds"))
saveRDS(results$log_profile_likelihood,   file = here(output_dir, "log_profile_likelihood.rds"))

message("✓ Experiment results saved to: ", output_dir)
