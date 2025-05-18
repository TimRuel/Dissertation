#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(yaml)
  library(fs)
  library(Rcpp)
  library(nloptr)
  library(PolytomousUtils)
})

proj_subdir <- here("Integrated_Likelihood_Algorithm", "Polytomous_Logistic_Regression")
proj_path <- function(...) here(proj_subdir, ...)

# Load helpers
miceadds::source.all(proj_path("scripts", "helpers"), print.source = FALSE)

# Wrapper to run other scripts
run_script <- function(script_rel_path, args = character()) {
  script_path <- proj_path(script_rel_path)
  system2("Rscript", c(script_path, args))
}

# Parse command line args
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript main.R <experiment_id> [iteration_id] [requested_cores]")

experiment_id <- args[1]
iteration_id <- args[2] %||% NULL
requested_cores <- args[3] %||% NULL

# Step 1: Ensure config exists
config_path <- proj_path("config", paste0(experiment_id, ".yml"))
if (!file.exists(config_path)) {
  message("Creating experiment config...")
  run_script("scripts/make_experiment_config.R", experiment_id)
}

# Step 2: Ensure data and resolved config exist
resolved_config_path <- proj_path("results", experiment_id, "resolved_config.yml")
if (!file.exists(resolved_config_path)) {
  message("Generating data and resolved config...")
  run_script("scripts/generate_data.R", experiment_id)
}

# Step 3: Load resolved config
config <- read_yaml(resolved_config_path)

# Step 4: Determine run mode and output paths
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S", tz = "UTC")
if (is.null(iteration_id)) {
  config$optimization_specs$run_type <- "individual"
  mode <- "individual_runs"
  run_id <- paste0("run_", timestamp)
} else {
  config$optimization_specs$run_type <- "simulation"
  mode <- "simulations"
  run_id <- paste0("iter_", sprintf("%03d", as.integer(iteration_id)))
}

# Core settings
core_info <- get_core_config(requested_cores)
config$optimization_specs$max_cores <- core_info$max_cores
config$optimization_specs$num_workers <- core_info$num_workers

# Output directory
output_dir <- proj_path("results", experiment_id, mode, run_id)
dir_create(output_dir, recurse = TRUE)

# Add run metadata to config
config$run_id <- run_id
config$iteration_id <- iteration_id
config$output_dir <- output_dir

# Step 5: Save config snapshot
config_snapshot_path <- file.path(output_dir, "config_snapshot.yml")
write_yaml(config, config_snapshot_path)

# Step 6: Run experiment
if (mode == "individual_runs") message("Running experiment...")
start_time <- Sys.time()
run_script("scripts/run_experiment.R", config_snapshot_path)
end_time <- Sys.time()
elapsed_time <- round(as.numeric(difftime(end_time, start_time, units = "secs")), 2)

# Step 7: Gather metadata
git_hash <- tryCatch({
  system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
}, error = function(e) NA_character_)

metadata <- list(
  experiment_id = experiment_id,
  run_id = run_id,
  mode = mode,
  iteration_id = iteration_id,
  slurm_array_id = Sys.getenv("SLURM_ARRAY_TASK_ID", unset = NA),
  timestamp_start = format(start_time, "%Y-%m-%d %H:%M:%S"),
  timestamp_end = format(end_time, "%Y-%m-%d %H:%M:%S"),
  elapsed_seconds = elapsed_time,
  git_commit = git_hash
)

# Step 8: Save metadata & append to run log
save_run_metadata(metadata, output_dir)

log_path <- proj_path("results", experiment_id, "run_log.csv")
append_run_log(metadata, log_path)

message("✓ Run completed: ", output_dir)
