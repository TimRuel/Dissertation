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
handlers("cli")

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

step_size <- 0.05

alpha <- 0.03

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

chunk_size <- 5

num_branches <- num_workers * chunk_size

tic()

log_integrated_likelihood <- get_log_integrated_likelihood(data,
                                                           X_h,
                                                           step_size,
                                                           threshold,
                                                           alpha,
                                                           prop,
                                                           init_guess_sd,
                                                           num_workers,
                                                           chunk_size)

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
                                                     alpha)

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

# log_integrated_likelihood$log_L_bar$df |> 
#   make_plot()
