# scripts/run_experiment.R

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
config <- read_yaml(config_snapshot_path)

# ---- Load Data ----
data_dir <- here(run_dir, "data")
required_files <- c("X_design.rds", "model_df.rds")
missing <- required_files[!file_exists(here(data_dir, required_files))]
if (length(missing) > 0) stop("Missing required data files: ", paste(missing, collapse = ", "))
X_design <- readRDS(here(data_dir, "X_design.rds"))
model_df <- readRDS(here(data_dir, "model_df.rds"))

# ---- Run integrated likelihood ----
num_workers <- config$optimization_specs$IL$num_workers

if (.Platform$OS.type == "unix") {
  plan(multicore, workers = I(num_workers))
} else {
  plan(multisession, workers = I(num_workers))
}

integrated_LL <- get_integrated_LL(config, X_design, model_df)

plan(sequential)

results_dir <- here(run_dir, "results")
dir_create(results_dir)
saveRDS(integrated_LL, file = here(results_dir, "integrated_LL.rds"))

# ---- Run profile likelihood ----
if (.Platform$OS.type == "unix") {
  plan(multicore, workers = I(2))
} else {
  plan(multisession, workers = I(2))
}

profile_LL <- get_profile_LL(config, X_design, model_df)

plan(sequential)

saveRDS(profile_LL, file = here(results_dir, "profile_LL.rds"))

report_objects <- get_report_objects(run_dir)
save_list_objects(report_objects, results_dir)
message("✓ Experiment results saved to /", sub(".*(/?experiments/.*)", "\\1", results_dir))
