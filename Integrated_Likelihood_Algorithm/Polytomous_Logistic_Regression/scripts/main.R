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

# Prompt helper
should_proceed_or_abort <- function(prompt, force = FALSE, non_interactive = FALSE) {
  if (force) return(TRUE)
  if (non_interactive) return(FALSE)
  
  # Prompt always, even in non-interactive terminal sessions
  cat(prompt, " [y/N]: ")
  response <- tolower(trimws(readLines("stdin", n = 1)))
  response %in% c("y", "yes")
}

# Parse command line args
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript main.R <experiment_id> [requested_cores] [iteration_id] [--force] [--non-interactive]")

experiment_id <- args[1]
requested_cores <- if (length(args) >= 2 && !grepl("^--", args[2])) as.integer(args[2]) else NULL
iteration_id <- if (length(args) >= 3 && !grepl("^--", args[3])) args[3] else NULL

force <- "--force" %in% args
non_interactive <- "--non-interactive" %in% args

# Step 1: Ensure config exists
config_path <- proj_path("config", paste0(experiment_id, ".yml"))
if (!file.exists(config_path)) {
  message("Creating experiment config...")
  run_script("scripts/make_experiment_config.R", experiment_id)
} else if (is.null(iteration_id)) {
  if (!should_proceed_or_abort(
    paste0("Experiment config already exists for '", experiment_id, "'. Proceed with it?"),
    force = force, non_interactive = non_interactive
  )) {
    stop("Aborting due to existing config.")
  }
}

# Step 2: Ensure data and resolved config exist
resolved_config_path <- proj_path("results", experiment_id, "resolved_config.yml")
if (!file.exists(resolved_config_path)) {
  message("Generating data and resolved config...")
  run_script("scripts/generate_data.R", experiment_id)
} else if (is.null(iteration_id)) {
  if (!should_proceed_or_abort(
    paste0("Data already exists for '", experiment_id, "'. Proceed with it?"),
    force = force, non_interactive = non_interactive
  )) {
    stop("Aborting due to existing data.")
  }
}

# Step 3: Load resolved config
config <- read_yaml(resolved_config_path)

# Step 4: Determine run mode and output paths
start_time <- Sys.time()
timestamp <- format(start_time, "%Y%m%d_%H%M%S", tz = "UTC")

if (is.null(iteration_id)) {
  config$optimization_specs$run_type <- "individual"
  mode <- "individual_runs"
  run_id <- paste0("run_", timestamp)
  config$run_id <- run_id
} else {
  config$optimization_specs$run_type <- "simulation"
  mode <- "simulations"
  run_id <- paste0("iter_", sprintf("%03d", as.integer(iteration_id)))
  config$iteration_id <- run_id
}

# Step 5: Core settings
core_info <- get_core_config(requested_cores)
config$optimization_specs$IL$max_cores <- core_info$max_cores
config$optimization_specs$IL$num_workers <- core_info$num_workers

# Step 6: Output directory and snapshot
output_dir <- proj_path("results", experiment_id, mode, run_id)
config$output_dir <- output_dir
dir_create(output_dir, recurse = TRUE)

config_snapshot_path <- file.path(output_dir, "config_snapshot.yml")
write_yaml(config, config_snapshot_path)

# Step 7: Run experiment
if (mode == "individual_runs") message("Running experiment...")
run_script("scripts/run_experiment.R", config_snapshot_path)
end_time <- Sys.time()
elapsed_time <- round(as.numeric(difftime(end_time, start_time, units = "secs")), 2)

# Step 8: Metadata
git_hash <- tryCatch(
  system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
  error = function(e) NA_character_
)

metadata <- list(
  experiment_id   = experiment_id,
  run_id          = run_id,
  mode            = mode,
  iteration_id    = iteration_id %||% NA_character_,
  slurm_array_id  = Sys.getenv("SLURM_ARRAY_TASK_ID", unset = NA),
  timestamp_start = start_time,  # POSIXct
  timestamp_end   = end_time,    # POSIXct
  elapsed_seconds = elapsed_time,
  git_commit      = git_hash
)

# Step 9: Save metadata & append to run log
save_run_metadata(metadata, output_dir)

log_path <- proj_path("results", experiment_id, "run_log.csv")
append_run_log(metadata, log_path)

message("✓ Run completed: ", output_dir)

# Left off needing to decide which config changes are big enough to warrant 
# new experiment and which can be kept within same
# First inclination - anything not to do with true param generation can be 
# considered "small" enough to remain within same experiment, e.g. step size 
# num std errors, etc.
# Optimization specs in other words
# Maybe break that into separate config template that gets appended to snapshot
# when doing an individual run or a simulation
# Goal for individual run is to have final config file saved in run directory
# Goal for simulations is to have separate snapshot for each iteration saved in
# iteration directory, along with data and plots. Main (only?) difference between
# iteration snapshots should be seed value
# Theoretically, every iteration should be reproducible given config snapshot





