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

# ---- Run Experiment ----
list2env(config$model_specs, environment())
list2env(config$optimization_specs, environment())

X1_levels <- config$X1_levels
formula <- as.formula(formula)
ml_model <- fit_multinomial_logistic_model(model_df, formula)
Y_design <- get_Y_design(model_df)
Beta_MLE <- get_Beta_MLE(ml_model)
threshold <- get_threshold(Beta_MLE, X_design, Y_design, IL$threshold_offset)
h <- get_X1_level_of_interest(X1_levels)
X_h_design <- get_X_h_design(X_design, X1_levels)
psi_hat <- get_psi_hat_from_model(ml_model, X1_levels)
n_h <- nrow(X_h_design)
psi_endpoints <- get_psi_endpoints(psi_hat, Beta_MLE, X_h_design, IL$num_std_errors, J, n_h)
psi_grid <- get_psi_grid(psi_endpoints, IL$step_size, J)

# ---- Run integrated likelihood ----

if (.Platform$OS.type == "unix") {
  plan(multicore, workers = I(IL$num_workers))
} else {
  plan(multisession, workers = I(IL$num_workers))
}

integrated_LL <- get_integrated_LL(
  X_design = X_design, 
  Y_design = Y_design, 
  X_h_design = X_h_design, 
  Jm1 = J - 1, 
  p = p, 
  n = n, 
  psi_grid = psi_grid, 
  psi_hat = psi_hat, 
  threshold = threshold, 
  init_guess_sd = IL$init_guess_sd, 
  num_workers = IL$num_workers, 
  chunk_size = IL$chunk_size
)

plan(sequential)

IL_plot <- get_LL_plot(integrated_LL$log_L_bar_df)
IL_branches_plot <- get_branches_plot(integrated_LL$branches_matrix)

# ---- Run profile likelihood ----

if (.Platform$OS.type == "unix") {
  plan(multicore, workers = I(2))
} else {
  plan(multisession, workers = I(2))
}

profile_LL <- get_profile_LL(
  step_size = PL$step_size, 
  alpha = PL$alpha,
  psi_hat = psi_hat,
  Beta_MLE = Beta_MLE, 
  X_design = X_design,
  Y_design = Y_design,
  X_h_design = X_h_design, 
  Jm1 = J - 1,
  p = p,
  n = n)

plan(sequential)

PL_plot <- get_LL_plot(profile_LL)

# ---- Store Results ----
results <- list(
  integrated_LL = integrated_LL,
  profile_LL = profile_LL,
  plots = list(IL_plot = IL_plot,
               IL_branches_plot = IL_branches_plot,
               PL_plot = PL_plot)
)

# ---- Save Results ----
results_dir <- here(run_dir, "results")
dir_create(results_dir)
saveRDS(results$integrated_LL, file = here(results_dir, "integrated_LL.rds"))
saveRDS(results$profile_LL, file = here(results_dir, "profile_LL.rds"))
message("✓ Experiment results saved to /", sub(".*(/?experiments/.*)", "\\1", results_dir))

plots_dir <- here(run_dir, "plots")
save_list_plots(results$plots, plots_dir)
message("✓ Log-likelihood plots saved to /", sub(".*(/?experiments/.*)", "\\1", plots_dir))


