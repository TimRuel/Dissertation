#!/usr/bin/env Rscript

# ---- CONFIGURATION ----
experiment_id <- "experiment_A"
n_iterations <- 10  # Set this to match the number you want locally
parallel <- TRUE    # Set to FALSE to run serially

# ---- Load libraries ----
if (parallel) {
  suppressPackageStartupMessages(library(parallel))
}

# ---- Define one run ----
run_one <- function(iteration_id) {
  message(sprintf("Running iteration %d", iteration_id))
  system2("Rscript", c("scripts/main.R", experiment_id, iteration_id))
}

# ---- Run loop ----
if (parallel) {
  n_cores <- detectCores() - 1
  message(sprintf("Running in parallel using %d cores", n_cores))
  mclapply(0:(n_iterations - 1), run_one, mc.cores = n_cores)
} else {
  for (i in 0:(n_iterations - 1)) {
    run_one(i)
  }
}
