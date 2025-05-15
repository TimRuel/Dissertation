# scripts/make_experiment_config.R

#!/usr/bin/env Rscript

if (interactive()) {
  message("Loading packages for interactive debugging...")
  library(yaml)
  library(fs)
}

args <- commandArgs(trailingOnly = TRUE)
experiment_id <- if (length(args) > 0) args[1] else stop("Provide experiment ID")

template_path <- "config/template_experiment.yml"
dest_path <- file.path("config", paste0(experiment_id, ".yml"))

if (file_exists(dest_path)) {
  stop("Experiment config already exists: ", dest_path)
}

experiment_config <- read_yaml(template_path)

# (Optional: modify experiment_config dynamically if needed here)

write_yaml(experiment_config, dest_path)
message("Created config: ", dest_path)
