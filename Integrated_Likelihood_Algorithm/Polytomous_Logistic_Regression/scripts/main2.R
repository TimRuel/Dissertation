#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(here)
  library(fs)
})

proj_subdir <- here("Integrated_Likelihood_Algorithm", "Polytomous_Logistic_Regression")
proj_path <- function(...) here(proj_subdir, ...)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript main2.R <experiment_id> <sim_id> <iteration_id> <requested_cores>")
}

experiment_id <- args[1]
sim_id <- args[2]
iteration_id <- args[3]
requested_cores <- as.integer(args[4])

exp_dir <- proj_path("experiments", experiment_id)
sim_dir <- here(exp_dir, "simulations", sim_id)
iter_dir <- here(sim_dir, iteration_id)
log_dir <- here(iter_dir, "logs")
log_path <- here(log_dir, paste0(iteration_id, "_log.out"))

message("▶️ Starting main2.R")
message("Experiment ID: ", experiment_id)
message("Simulation ID: ", sim_id)
message("Iteration ID: ", iteration_id)
message("Requested cores: ", requested_cores)
message("📁 Experiment directory: ", exp_dir)
message("📁 Simulation directory: ", sim_dir)
message("📁 Iteration directory: ", iter_dir)
message("📁 Log directory: ", log_dir)
message("📁 Log file: ", log_path)
message("✓ Iteration completed")
