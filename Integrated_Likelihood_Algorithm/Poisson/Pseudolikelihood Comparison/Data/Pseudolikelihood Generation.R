library(future)
library(doFuture)
library(zeallot)
library(purrr)
library(dplyr)
library(stringr)
library(progressr)
library(tictoc)

handlers(global = TRUE)
handlers("cli")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("Data Generation.R")

num_cores <- Sys.getenv("SLURM_NPROCS") |>
  as.numeric()
# num_cores <- availableCores() |>
#   as.numeric()
# num_cores <- parallel::detectCores() |>
#   as.numeric()

step_size <- 0.01

num_std_errors <- 3

R <- 250

tol <- 0.0001

################################################################################
############################ INTEGRATED LIKELIHOOD ############################# 
################################################################################

plan(sequential)

set.seed(seed)

alpha <- 1/2

beta <- 1/2

u_params <- list(alpha = alpha, beta = beta)

chunk_size <- ceiling(R / num_cores)

plan(multisession, workers = num_cores)

integrated_log_likelihood_vals <- get_integrated_log_likelihood_vals(data,
                                                                     weights,
                                                                     step_size, 
                                                                     num_std_errors,
                                                                     u_params, 
                                                                     R, 
                                                                     tol, 
                                                                     chunk_size)
  
################################################################################
######################## MODIFIED INTEGRATED LIKELIHOOD ########################
################################################################################

plan(sequential)

set.seed(seed)

alpha <- data + 1/2

beta <- 1 + 1/2

u_params <- list(alpha = alpha, beta = beta)

# Q_name <- "euclidean_distance"
Q_name <- "neg_log_likelihood"

Q <- Q_name |>
  get()

chunk_size <- ceiling(R / num_cores)

plan(multisession, workers = num_cores)

mod_integrated_log_likelihood_vals <- get_mod_integrated_log_likelihood_vals(data,
                                                                             weights,
                                                                             Q, 
                                                                             step_size, 
                                                                             num_std_errors, 
                                                                             u_params, 
                                                                             R, 
                                                                             tol, 
                                                                             chunk_size)

################################################################################
############################## PROFILE LIKELIHOOD ############################## 
################################################################################

plan(sequential)

psi_grid_list <- get_psi_grid(data, weights, step_size, num_std_errors, split = TRUE)

profile_log_likelihood_vals <- get_profile_log_likelihood(data, weights, psi_grid_list)

################################################################################
################################### STORAGE #################################### 
################################################################################

psi_grid <- get_psi_grid(data, weights, step_size, num_std_errors, split = FALSE)

log_likelihood_vals <- data.frame(psi = psi_grid,
                                  Integrated = integrated_log_likelihood_vals,
                                  Mod_Integrated = mod_integrated_log_likelihood_vals,
                                  Profile = profile_log_likelihood_vals)

Q_name <- Q_name |> 
  strsplit("_") |> 
  pluck(1) |> 
  substr(1, 1) |> 
  paste(collapse = "")

log_likelihood_vals_file_path <- population_params_file_path |> 
  str_extract("Population\\s[A-Za-z0-9]+") |> 
  paste0("Pseudolikelihoods/", . = _, "/log_likelihood_vals_") |> 
  glue::glue("R={R}_Q={Q_name}__seed={seed}_stepsize={step_size}_numse={num_std_errors}.Rda")

saveRDS(log_likelihood_vals, log_likelihood_vals_file_path)

