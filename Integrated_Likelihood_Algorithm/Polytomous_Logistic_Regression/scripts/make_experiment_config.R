# scripts/make_experiment_config.R

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(here)
  library(yaml)
  library(fs)
})

proj_subdir <- here("Integrated_Likelihood_Algorithm", "Polytomous_Logistic_Regression")
proj_path <- function(...) here(proj_subdir, ...)

args <- commandArgs(trailingOnly = TRUE)
experiment_id <- if (length(args) > 0) args[1] else stop("[ERROR] Provide experiment ID")

template_path <- proj_path("config", "template_experiment.yml")
dest_path <- proj_path("config", "exps", paste0(experiment_id, ".yml"))

experiment_config <- read_yaml(template_path)
experiment_config$experiment_id <- experiment_id

write_yaml(experiment_config, dest_path)
message("[INFO] Created and saved config to /", sub(".*(/?config/.*)", "\\1", dest_path))
