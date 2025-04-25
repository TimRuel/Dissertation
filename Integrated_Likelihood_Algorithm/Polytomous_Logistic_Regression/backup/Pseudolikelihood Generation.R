library(future)
library(doFuture)
library(zeallot)
library(purrr)
library(tidyverse)
library(insight)
library(invgamma)
library(stringr)
library(progressr)
library(tictoc)
library(rstudioapi)

setwd(dirname(getActiveDocumentContext()$path))

population_directory <- selectDirectory(caption = "Select population directory")

setwd(population_directory)

source("data.R")

set.seed(seed)

################################################################################
################################## PARAMETERS ################################## 
################################################################################

h <- "A"

step_size <- 0.01

num_std_errors <- 3.5

################################################################################
########################## INTEGRATED LIKELIHOOD - VANILLA MC ##################
################################################################################

threshold <- ceiling(abs(log_likelihood(Beta_MLE, X_design, model.matrix(~ Y)[,-1]))) + 20

init_guess_sd <- 5

# num_workers <- Sys.getenv("SLURM_NPROCS") |>
#   as.numeric()

# num_workers <- availableCores() |>
#   as.numeric()

num_workers <- parallel::detectCores() |>
  as.integer()

chunk_size <- 1

num_branches <- num_workers * chunk_size

IL_maxtime <- 10

tic()

log_integrated_likelihood <- get_log_integrated_likelihood(data,
                                                           formula,
                                                           h,
                                                           step_size,
                                                           num_std_errors,
                                                           init_guess_sd,
                                                           threshold,
                                                           num_workers,
                                                           chunk_size,
                                                           IL_maxtime)

toc()

log_integrated_likelihood_filepath <- glue::glue("IL_objects/R={num_branches}_J={J}_h={h}_stepsize={step_size}_numse={num_std_errors}.Rda")

saveRDS(log_integrated_likelihood, log_integrated_likelihood_filepath)

# log_integrated_likelihood <- readRDS("IL_objects/R=260_h=b_stepsize=0.02.Rda")

################################################################################
############################## PROFILE LIKELIHOOD ############################## 
################################################################################

PL_maxtime <- 10

tic()

log_profile_likelihood <- get_log_profile_likelihood(data,
                                                     formula,
                                                     h,
                                                     step_size,
                                                     num_std_errors,
                                                     init_guess_sd,
                                                     PL_maxtime)

toc()

log_profile_likelihood_filepath <- glue::glue("PL_objects/R={num_branches}_J={J}_h={h}_stepsize={step_size}_numse={num_std_errors}.Rda")

saveRDS(log_profile_likelihood, log_profile_likelihood_filepath)

# log_profile_likelihood <- readRDS("PL_objects/R=260_h=b_stepsize=0.02.Rda")

################################################################################
################################### STORAGE #################################### 
################################################################################

log_integrated_likelihood <- readRDS("Integrated_Likelihood_Algorithm/Polytomous Logistic Regression/Simulations/Marginal Entropy/IL_objects/h=A/Sim1.Rda")

get_log_L_bar(log_integrated_likelihood)

log_likelihood_vals <- log_integrated_likelihood$log_L_bar$df |> 
  merge(log_profile_likelihood, all = TRUE)

# log_likelihood_vals <- log_L_bar_df |> 
#   merge(log_profile_likelihood, all = TRUE)
# 
# log_likelihood_vals_file_path <- glue::glue("log_likelihood_vals/test.Rda")

log_likelihood_vals_file_path <- glue::glue("log_likelihood_vals/R={num_branches}_J={J}_h={h}_stepsize={step_size}_numse={num_std_errors}.Rda")

saveRDS(log_likelihood_vals, log_likelihood_vals_file_path)



for (df in log_integrated_likelihood$log_L_tilde_df) plot(df)








