# scripts/compare_results.R

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(geomtextpath)
  library(ggnewscale)
  library(ggrepel)
  library(viridis)
  library(here)
  library(yaml)
  library(fs)
  library(zeallot)
  library(kableExtra)
  library(stringr)
  library(splines)
})

proj_subdir <- here("Integrated_Likelihood_Algorithm", "Polytomous_Logistic_Regression")
proj_path <- function(...) here(proj_subdir, ...)
miceadds::source.all(proj_path("scripts", "helpers"), print.source = FALSE)

# Parse Arguments ---------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
run_dir <- if (length(args) > 0) args[1] else stop("[ERROR] A run directory was not provided.")
if (!file.exists(run_dir)) stop("[ERROR] Run directory does not exist at /", sub(".*(/?experiments/.*)", "\\1", run_dir))

# Load Config -------------------------------------------------------------
config_snapshot_path <- here(run_dir, "config_snapshot.yml")
config <- read_yaml(config_snapshot_path)

# Load True Parameters ----------------------------------------------------
true_params_dir <- run_dir |> 
  dirname() |> 
  dirname() |> 
  here("true_params")

H_0 <- readRDS(here(true_params_dir, "H_0.rds"))
h <- get_X1_level_of_interest(config$X1_levels)
psi_0 <- H_0 |>
  filter(X1 == h) |>
  pull(entropy)

# Load Results ------------------------------------------------------------
results_dir <- here(run_dir, "results")
integrated_LL_path <- here(results_dir, "integrated_LL.rds")
profile_LL_path <- here(results_dir, "profile_LL.rds")
integrated_LL <- readRDS(integrated_LL_path)
profile_LL <- readRDS(profile_LL_path)

LL_df <- integrated_LL$log_L_bar_df |> 
  merge(profile_LL, all = TRUE)

# Compare Pseudolikelihoods -----------------------------------------------
LL_df_long <- get_LL_df_long(LL_df)

spline_models <- get_spline_models(LL_df_long)

MLE_data <- get_MLE_data(spline_models)

pseudolikelihoods <- get_pseudolikelihoods(spline_models, MLE_data)

alpha <- config$optimization_specs$PL$alpha

J <- config$model_specs$J

conf_ints <- get_confidence_intervals(
  pseudolikelihoods = pseudolikelihoods, 
  LL_df_long = LL_df_long,
  MLE_data = MLE_data, 
  alpha = alpha, 
  J = J
  )

LL_comparison_table <- get_LL_comparison_table(
  MLE_data = MLE_data, 
  conf_ints = conf_ints, 
  psi_0 = psi_0,
  alpha = alpha
)

stat_fns %<-% get_stat_fns(pseudolikelihoods, LL_df)

LL_comparison_plot <- get_LL_comparison_plot(
  stat_fns = stat_fns, 
  LL_df_long = LL_df_long,
  MLE_data = MLE_data,
  alpha = alpha
  )

# Save Analysis -----------------------------------------------------------
plots_dir <- here(run_dir, "plots")
invisible(
  suppressMessages(
    suppressWarnings({
      save_kable(LL_comparison_table, file = here(plots_dir, "LL_comparison_table.png"), zoom = 3)
      save_list_plots(list(LL_comparison_plot = LL_comparison_plot), plots_dir)
      }
      )
    )
  )
message("✓ Log-likelihood comparison table and plot saved to /", sub(".*(/?experiments/.*)", "\\1", plots_dir))
