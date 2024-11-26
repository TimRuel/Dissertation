library(future)
library(doFuture)
library(zeallot)
library(purrr)
library(plyr)
library(tidyverse)
library(stringr)
library(progressr)
library(tictoc)

handlers(global = TRUE)
handlers("cli")

# num_cores <- Sys.getenv("SLURM_NPROCS") |>
#   as.numeric()

num_cores <- availableCores() |>
  as.numeric()

# num_cores <- parallel::detectCores() |>
#   as.numeric()

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("utils.R")

population_directory <- choose_directory()

setwd(population_directory)

source("data.R")

set.seed(seed)

################################################################################
############################## PROFILE LIKELIHOOD ############################## 
################################################################################

step_size <- 0.01

num_std_errors <- 10

x_h <- 20.95197

profile_log_likelihood_vals <- get_profile_log_likelihood(x, 
                                                          y, 
                                                          x_h,
                                                          step_size, 
                                                          num_std_errors)

################################################################################
################################### STORAGE #################################### 
################################################################################

psi_grid <- get_psi_grid(x, y, x_h, step_size, num_std_errors, split = FALSE)[-1]

log_likelihood_vals <- data.frame(psi = psi_grid,
                                  Profile = profile_log_likelihood_vals)

log_likelihood_vals_file_path <- glue::glue("log_likelihood_vals_seed={seed}_xh={as.character(x_h)}_stepsize={step_size}_numse={num_std_errors}.Rda")

saveRDS(log_likelihood_vals, log_likelihood_vals_file_path)


plot(psi_grid, profile_log_likelihood_vals)



