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

h <- "B"

step_size <- 0.05

alpha <- 0.03

max_retries <- 10

max_eval <- 1e4

maxtime <- 1.5

################################################################################
########################## INTEGRATED LIKELIHOOD - VANILLA MC ##################
################################################################################

quantiles <- c(0.25, 0.5)

# init_guess_sd <- 5 # population C

init_guess_sd <- 0.5 # population D

# num_workers <- Sys.getenv("SLURM_NPROCS") |>
#   as.numeric()

# num_workers <- availableCores() |>
#   as.numeric()

num_workers <- parallel::detectCores() |>
  as.integer()

# num_workers <- 12

chunk_size <- 5

num_branches <- num_workers * chunk_size

tic()

branch_specs <- generate_branch_specs(data,
                                      h,
                                      init_guess_sd,
                                      alpha,
                                      num_workers,
                                      chunk_size,
                                      max_retries,
                                      max_eval,
                                      maxtime)

toc()

branch_specs_filepath <- glue::glue("branch_specs/R={num_branches}_J={J}_h={h}_alpha={alpha}.Rda")

saveRDS(branch_specs, branch_specs_filepath)

# branch_specs <- readRDS("branch_specs/R=60_h=A_alpha=0.03.Rda")

tic()

log_integrated_likelihood <- get_log_integrated_likelihood(branch_specs |> sample(num_branches),
                                                           data,
                                                           h,
                                                           alpha,
                                                           step_size,
                                                           quantiles,
                                                           init_guess_sd,
                                                           num_workers,
                                                           chunk_size,
                                                           max_retries,
                                                           max_eval,
                                                           maxtime)

toc()

log_integrated_likelihood_filepath <- glue::glue("IL_objects/R={num_branches}_J={J}_h={h}_stepsize={step_size}.Rda")

saveRDS(log_integrated_likelihood, log_integrated_likelihood_filepath)

# log_integrated_likelihood <- readRDS("IL_objects/R=260_h=b_stepsize=0.02.Rda")

################################################################################
############################## PROFILE LIKELIHOOD ############################## 
################################################################################

tic()

log_profile_likelihood <- get_log_profile_likelihood(data,
                                                     h,
                                                     step_size,
                                                     alpha,
                                                     init_guess_sd,
                                                     max_retries,
                                                     max_eval,
                                                     maxtime)

toc()

log_profile_likelihood_filepath <- glue::glue("PL_objects/R={num_branches}_J={J}_h={h}_stepsize={step_size}.Rda")

saveRDS(log_profile_likelihood, log_profile_likelihood_filepath)

# log_profile_likelihood <- readRDS("PL_objects/R=260_h=b_stepsize=0.02.Rda")

################################################################################
################################### STORAGE #################################### 
################################################################################

log_likelihood_vals <- log_integrated_likelihood$log_L_bar$df |> 
  merge(log_profile_likelihood, all = TRUE)

# log_likelihood_vals <- log_L_bar_df |> 
#   merge(log_profile_likelihood, all = TRUE)
# 
# log_likelihood_vals_file_path <- glue::glue("log_likelihood_vals/test.Rda")

log_likelihood_vals_file_path <- glue::glue("log_likelihood_vals/R={num_branches}_J={J}_h={h}_stepsize={step_size}.Rda")

saveRDS(log_likelihood_vals, log_likelihood_vals_file_path)
