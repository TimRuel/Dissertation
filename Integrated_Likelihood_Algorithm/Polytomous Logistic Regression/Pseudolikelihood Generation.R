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

# source("utils.R")

population_directory <- selectDirectory(caption = "Select population directory")

setwd(population_directory)

source("data.R")

set.seed(seed)

################################################################################
################################## PARAMETERS ################################## 
################################################################################

h <- 2

X_h <- data.frame(X = factor(h))

step_size <- 0.1

alpha <- 0.05

################################################################################
########################## INTEGRATED LIKELIHOOD - VANILLA MC ##################
################################################################################

threshold <- 50

init_guess_sd <- 20

prop <- 0.7

# num_workers <- Sys.getenv("SLURM_NPROCS") |>
#   as.numeric()

# num_workers <- availableCores() |>
#   as.numeric()

num_workers <- parallel::detectCores() |>
  as.numeric()

chunk_size <- 1

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

################################################################################
############################## PROFILE LIKELIHOOD ############################## 
################################################################################

tic()

log_profile_likelihood <- get_log_profile_likelihood(data,
                                                     X_h,
                                                     step_size,
                                                     alpha)

toc()

################################################################################
################################### STORAGE #################################### 
################################################################################

log_likelihood_vals <- log_integrated_likelihood$log_L_bar$df |> 
  merge(log_profile_likelihood, all = TRUE)

# log_likelihood_vals <- values_df |> 
#   merge(log_profile_likelihood$values_df, all = TRUE)

log_likelihood_vals_file_path <- glue::glue("seed={seed}_J={J}_p={p}_m={m}_R={num_branches}_h={h}_stepsize={step_size}_alpha={alpha}.Rda")

saveRDS(log_likelihood_vals, log_likelihood_vals_file_path)

# log_integrated_likelihood$log_L_bar$df |> 
#   make_plot()
