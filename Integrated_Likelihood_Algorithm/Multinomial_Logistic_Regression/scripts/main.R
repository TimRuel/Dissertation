# scripts/main.R

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(doFuture)
  library(here)
  library(quarto)
  library(yaml)
  library(fs)
  library(Rcpp)
  library(nloptr)
  library(PolytomousUtils)
})

i_am("Integrated_Likelihood_Algorithm/Polytomous_Logistic_Regression/scripts/main.R")

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
if (length(args) < 2) stop("Usage: Rscript main.R <experiment_id> [requested_cores] [sim_id] [run_id] [--force] [--non-interactive]")

experiment_id <- args[1]
requested_cores <- if (length(args) >= 2 && !grepl("^--", args[2])) as.integer(args[2]) else NULL
sim_id <- if (length(args) >= 3 && !grepl("^--", args[3])) args[3] else NULL
run_id <- if (length(args) >= 4 && !grepl("^--", args[4])) args[4] else NULL

force <- "--force" %in% args
non_interactive <- "--non-interactive" %in% args

# Step 1: Ensure config exists
config_path <- proj_path("config", "exps", paste0(experiment_id, ".yml"))
if (!file.exists(config_path)) {
  message("Creating experiment config...")
  run_script("scripts/make_experiment_config.R", experiment_id)
} else if (is.null(run_id)) {
  if (!should_proceed_or_abort(
    paste0("Experiment config already exists for '", experiment_id, "'. Proceed with it?"),
    force = force, non_interactive = non_interactive
  )) {
    stop("Aborting due to existing config.")
  }
}

# Step 2: Generate true parameters if missing
true_params_dir <- proj_path("experiments", experiment_id, "true_params")
dir_create(true_params_dir)

if (length(dir_ls(true_params_dir, fail = FALSE)) == 0) {
  message("[INFO] No true parameters found — generating...")
  run_script("scripts/generate_true_params.R", experiment_id)
} else {
  message("[INFO] True parameters already exist — skipping generation.")
}

# Step 3: Generate run_id and path to the run directory, creating if necessary 
if (is.null(sim_id)) {
  run_id <- paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S", tz = "UTC"))
  run_dir <- proj_path("experiments", experiment_id, "individual_runs", run_id)
  dir_create(run_dir)
  data_args <- c(experiment_id, run_id)
} else {
  run_dir <- proj_path("experiments", experiment_id, "simulations", sim_id, run_id)
  data_args <- c(experiment_id, sim_id, run_id)
}

# Step 4: Generate data
run_script("scripts/generate_data.R", data_args)

# Step 5: Add optimization config to config snapshot
config_snapshot_path <- here(run_dir, "config_snapshot.yml")
config_snapshot <- read_yaml(config_snapshot_path)
opt_config_path <- proj_path("config", "opt_config.yml")
opt_config <- read_yaml(opt_config_path)
config_snapshot$optimization_specs <- c(config_snapshot$optimization_specs, opt_config)

# Step 6: Core settings
core_info <- get_core_config(requested_cores)
config_snapshot$optimization_specs$IL$max_cores <- core_info$max_cores
config_snapshot$optimization_specs$IL$num_workers <- core_info$num_workers

# Step 7: Save config snapshot
write_strict_yaml(config_snapshot, config_snapshot_path)

# Step 8: Run experiment
message("Running experiment...")
start_time <- Sys.time()
run_script("scripts/run_experiment.R", run_dir)
end_time <- Sys.time()

# Step 9: Metadata
elapsed_time <- round(as.numeric(difftime(end_time, start_time, units = "mins")), 2)

git_hash <- tryCatch(
  system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
  error = function(e) NA_character_
)

metadata <- list(
  experiment_id   = experiment_id,
  sim_id          = sim_id,
  run_id          = run_id,
  slurm_array_id  = Sys.getenv("SLURM_ARRAY_TASK_ID", unset = NA),
  timestamp_start = start_time,
  timestamp_end   = end_time,
  elapsed_seconds = elapsed_time,
  git_commit      = git_hash
)

# Step 10: Save metadata
log_dir <- here(run_dir, "logs")
save_run_metadata(metadata, log_dir)
message("✓ Run completed")

