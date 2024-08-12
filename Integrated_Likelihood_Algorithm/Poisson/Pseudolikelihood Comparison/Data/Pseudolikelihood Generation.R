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

# num_cores <- Sys.getenv("SLURM_NPROCS") |>
#   as.numeric()
# num_cores <- availableCores() |>
#   as.numeric()
num_cores <- parallel::detectCores() |>
  as.numeric()

step_size <- 0.1

num_std_errors <- 3

R <- 104

# tol <- 0.001

################################################################################
############################ INTEGRATED LIKELIHOOD ############################# 
################################################################################

plan(sequential)

set.seed(seed)

dist <- runif

dist_params <- list(min = 0, max = 100)

# chunk_size <- ceiling(R / num_cores)
chunk_size <- floor(R / num_cores)

plan(multisession, workers = num_cores)

tic()

integrated_log_likelihood_vals <- get_integrated_log_likelihood_vals(data,
                                                                     weights,
                                                                     step_size, 
                                                                     num_std_errors,
                                                                     dist,
                                                                     dist_params, 
                                                                     R, 
                                                                     chunk_size)

toc()
  
################################################################################
######################## MODIFIED INTEGRATED LIKELIHOOD ########################
################################################################################

plan(sequential)

set.seed(seed)

alpha_prior <- 1/2

beta_prior <- 1/2

alpha_posterior <- data |> 
  map_dbl(sum) |>
  (`+`)(alpha_prior)

beta_posterior <- data |> 
  map_dbl(length) |>
  (`+`)(beta_prior)

dist <- rgamma

dist_params <- list(shape = alpha_posterior, rate = beta_posterior)

# Q_name <- "euclidean_distance"
# Q_name <- "neg_log_likelihood"
# 
# Q <- Q_name |>
#   get()

# chunk_size <- ceiling(R / num_cores)
chunk_size <- floor(R / num_cores)

plan(multisession, workers = num_cores)

tic()

mod_integrated_log_likelihood_vals <- get_mod_integrated_log_likelihood_vals(data,
                                                                             weights,
                                                                             step_size, 
                                                                             num_std_errors, 
                                                                             dist,
                                                                             dist_params, 
                                                                             R, 
                                                                             chunk_size)

toc()

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

# Q_name <- Q_name |> 
#   strsplit("_") |> 
#   pluck(1) |> 
#   substr(1, 1) |> 
#   paste(collapse = "")

log_likelihood_vals_file_path <- population_params_file_path |> 
  sub(".*Pseudolikelihoods/(.*)/parameters\\.R", "\\1", x = _) |> 
  paste0("Pseudolikelihoods/", . = _, "/log_likelihood_vals_") |> 
  glue::glue("R={R}_seed={seed}_stepsize={step_size}_numse={num_std_errors}.Rda")

saveRDS(log_likelihood_vals, log_likelihood_vals_file_path)

