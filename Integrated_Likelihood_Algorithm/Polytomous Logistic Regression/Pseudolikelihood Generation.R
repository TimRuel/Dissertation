library(future)
library(doFuture)
library(zeallot)
library(purrr)
library(plyr)
library(tidyverse)
library(insight)
library(invgamma)
library(stringr)
library(progressr)
library(tictoc)
library(rstudioapi)

handlers(global = TRUE)
handlers("progress")

setwd(dirname(getActiveDocumentContext()$path))

population_directory <- selectDirectory(caption = "Select population directory")

setwd(population_directory)

source("data.R")

set.seed(seed)

################################################################################
################################## PARAMETERS ################################## 
################################################################################

h <- 1

X_h <- data.frame(X = factor(h))

step_size <- 0.01

alpha <- 0.03

lambda <- 1e-4

max_retries <- 10

################################################################################
########################## INTEGRATED LIKELIHOOD - VANILLA MC ##################
################################################################################

threshold <- 40

init_guess_sd <- 5

# num_workers <- Sys.getenv("SLURM_NPROCS") |>
#   as.numeric()

# num_workers <- availableCores() |>
#   as.numeric()

num_workers <- parallel::detectCores() |>
  as.integer()

chunk_size <- 3

num_branches <- num_workers * chunk_size

tic()

branch_specs <- generate_branch_specs(data,
                                      X_h,
                                      init_guess_sd,
                                      alpha,
                                      num_workers,
                                      chunk_size,
                                      lambda,
                                      max_retries)

branch_specs_filepath <- glue::glue("branch_specs/R={num_branches}_h={h}_alpha={alpha}.Rda")

saveRDS(branch_specs, branch_specs_filepath)

toc()

tic()

log_integrated_likelihood <- get_log_integrated_likelihood(branch_specs,
                                                           data,
                                                           X_h,
                                                           alpha,
                                                           step_size,
                                                           init_guess_sd,
                                                           chunk_size,
                                                           lambda,
                                                           max_retries)

toc()

log_integrated_likelihood_filepath <- glue::glue("IL_obj_R={num_branches}_h={h}_stepsize={step_size}.Rda")

saveRDS(log_integrated_likelihood, log_integrated_likelihood_filepath)

################################################################################
############################## PROFILE LIKELIHOOD ############################## 
################################################################################

tic()

log_profile_likelihood <- get_log_profile_likelihood(data,
                                                     X_h,
                                                     step_size,
                                                     alpha,
                                                     init_guess_sd,
                                                     lambda,
                                                     max_retries)

toc()

log_profile_likelihood_filepath <- glue::glue("PL_obj_R={num_branches}_h={h}_stepsize={step_size}.Rda")

saveRDS(log_profile_likelihood, log_profile_likelihood_filepath)

################################################################################
################################### STORAGE #################################### 
################################################################################

log_likelihood_vals <- log_integrated_likelihood$log_L_bar$df |> 
  merge(log_profile_likelihood, all = TRUE)

# log_likelihood_vals <- values_df |> 
#   merge(log_profile_likelihood$values_df, all = TRUE)

log_likelihood_vals_file_path <- glue::glue("log_likelihood_vals_R={num_branches}_h={h}_stepsize={step_size}.Rda")

saveRDS(log_likelihood_vals, log_likelihood_vals_file_path)
