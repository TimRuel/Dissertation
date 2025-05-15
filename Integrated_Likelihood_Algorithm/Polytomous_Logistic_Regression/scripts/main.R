# scripts/main.R

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(fs)
  library(yaml)
  library(Rcpp)
  library(nloptr)
  library(PolytomousUtils)
  # sourceCpp("../../polytomous_utils.cpp")
  miceadds::source.all("scripts/helpers")
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript main.R <experiment_id> [iteration_id]")

experiment_id <- args[1]
iteration_id <- if (length(args) >= 2) args[2] else NULL

# Step 1: Ensure experiment config exists
config_path <- file.path("config", paste0(experiment_id, ".yml"))
if (!file.exists(config_path)) {
  message("Creating experiment config...")
  system2("Rscript", c("scripts/make_experiment_config.R", experiment_id))
}

# Step 2: Ensure data and resolved config exist
resolved_config_path <- file.path("results", experiment_id, "resolved_config.yml")
if (!file.exists(resolved_config_path)) {
  message("Generating data and resolved config...")
  system2("Rscript", c("scripts/generate_data.R", experiment_id))
}

# Step 3: Load resolved config
config <- yaml::read_yaml(resolved_config_path)

# Step 4: Determine output mode and directory
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
if (is.null(iteration_id)) {
  mode <- "individual_runs"
  run_id <- paste0("run_", timestamp)
} else {
  mode <- "simulations"
  run_id <- paste0("iter_", sprintf("%03d", as.integer(iteration_id)))
}

output_dir <- file.path("results", experiment_id, mode, run_id)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Step 5: Save config snapshot
write_yaml(config, file.path(output_dir, "config_snapshot.yml"))

# Step 6: Run model and save output
source("scripts/run_model.R")  # defines run_model(config)
result <- run_model(config)
saveRDS(result, file = file.path(output_dir, "output.rds"))

message("✓ Run completed: ", output_dir)
