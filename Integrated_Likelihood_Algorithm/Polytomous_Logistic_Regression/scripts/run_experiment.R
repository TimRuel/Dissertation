#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(yaml)
  library(fs)
  library(doFuture)
  library(PolytomousUtils)
})

proj_subdir <- here("Integrated_Likelihood_Algorithm", "Polytomous_Logistic_Regression")
proj_path <- function(...) here(proj_subdir, ...)
miceadds::source.all(proj_path("scripts", "helpers"), print.source = FALSE)

# ---- Parse Arguments ----
args <- commandArgs(trailingOnly = TRUE)
run_dir <- if (length(args) > 0) args[1] else stop("[ERROR] A run directory was not provided.")
if (!file.exists(run_dir)) stop("[ERROR] Run directory does not exist at /", sub(".*(/?experiments/.*)", "\\1", run_dir))

# ---- Load Config ----
config_snapshot_path <- here(run_dir, "config_snapshot.yml")
config_snapshot <- read_yaml(config_snapshot_path)

# ---- Load Data ----
data_dir <- here(run_dir, "data")
required_files <- c("X_design.rds", "model_df.rds")
missing <- required_files[!file_exists(here(data_dir, required_files))]
if (length(missing) > 0) stop("Missing required data files: ", paste(missing, collapse = ", "))

X_design <- readRDS(here(data_dir, "X_design.rds"))
model_df <- readRDS(here(data_dir, "model_df.rds"))

# ---- Run Experiment ----
results <- run_experiment(config_snapshot, X_design, model_df)

# ---- Save Results ----
results_dir <- here(run_dir, "results")
dir_create(results_dir)
saveRDS(results$integrated_LL, file = here(results_dir, "integrated_LL.rds"))
saveRDS(results$profile_LL, file = here(results_dir, "profile_LL.rds"))
message("✓ Experiment results saved to /", sub(".*(/?experiments/.*)", "\\1", results_dir))

plots_dir <- here(run_dir, "plots")
save_list_plots(results$plots, plots_dir)
message("✓ Log-likelihood plots saved to /", sub(".*(/?experiments/.*)", "\\1", plots_dir))


