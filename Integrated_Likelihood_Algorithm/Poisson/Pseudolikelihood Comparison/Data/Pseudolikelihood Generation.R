library(future)
library(doFuture)
library(zeallot)
library(purrr)
library(dplyr)
library(plyr)
library(stringr)
library(progressr)
library(tictoc)

handlers(global = TRUE)
handlers("cli")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("Data Generation.R")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# num_cores <- Sys.getenv("SLURM_NPROCS") |>
#   as.numeric()

# num_cores <- availableCores() |>
#   as.numeric()

num_cores <- parallel::detectCores() |>
  as.numeric()

step_size <- 0.01

num_std_errors <- 3

psi_grid <- get_psi_grid(data, weights, step_size, num_std_errors, split = FALSE)

R <- 250

################################################################################
############################ INTEGRATED LIKELIHOOD ############################# 
################################################################################

alpha_prior <- 0.5

beta_prior <- 0.5

alpha_posterior <- alpha_prior + data

beta_posterior <- beta_prior + rep(1, length(data))

rng <- rgamma

density <- dgamma

nominal_rng_params <- list(n = n, shape = alpha_posterior, rate = beta_posterior)

nominal_density_params <- list(shape = alpha_posterior, rate = beta_posterior)

importance_rng_rate_vec <- 2 

importance_rng_shape_vec <- importance_rng_rate_vec / beta_posterior * (alpha_posterior - 1) + 1

importance_rng_params <- list(n = n, shape = importance_rng_shape_vec, rate = importance_rng_rate_vec)

importance_density_params <- list(shape = importance_rng_shape_vec, rate = importance_rng_rate_vec)

lambda_method <- "accumulate"

chunk_size <- 5

init_lambda <- 0

################################################################################
########################## INTEGRATED LIKELIHOOD - VANILLA MC ##################
################################################################################

method = "vanilla_MC"

MC_params <- list(method = method, 
                  nominal = list(rng = rng, 
                                 rng_params = nominal_rng_params,
                                 density = density, 
                                 density_params = nominal_density_params))

tic()

plan(multisession, workers = I(50))

integrated_likelihood_vanilla_MC <- get_integrated_likelihood(data, 
                                                              weights, 
                                                              psi_grid, 
                                                              R,
                                                              MC_params,
                                                              lambda_method,
                                                              init_lambda,
                                                              chunk_size)

toc()

################################################################################
##################### INTEGRATED LIKELIHOOD - IMPORTANCE SAMPLING ##############
################################################################################

method = "basic_IS"

MC_params <- list(method = method, 
                  importance = list(rng = rng, 
                                    rng_params = importance_rng_params,
                                    density = density, 
                                    density_params = importance_density_params), 
                  nominal = list(rng = rng, 
                                 rng_params = nominal_rng_params,
                                 density = density, 
                                 density_params = nominal_density_params))

tic()

plan(multisession, workers = I(50))

integrated_likelihood_basic_IS <- get_integrated_likelihood(data, 
                                                            weights, 
                                                            step_size, 
                                                            num_std_errors, 
                                                            R,
                                                            MC_params,
                                                            lambda_method,
                                                            init_lambda,
                                                            chunk_size)

toc()

basic_IS_diagnostics <- get_importance_sampling_diagnostics(integrated_likelihood_basic_IS)

################################################################################
########## INTEGRATED LIKELIHOOD - SELF-NORMALIZED IMPORTANCE SAMPLING #########
################################################################################

method = "self_norm_IS"

MC_params <- list(method = method, 
                  importance = list(rng = rng, 
                                    rng_params = importance_rng_params,
                                    density = density, 
                                    density_params = importance_density_params), 
                  nominal = list(rng = rng, 
                                 rng_params = nominal_rng_params,
                                 density = density, 
                                 density_params = nominal_density_params))

tic()

plan(multisession, workers = I(50))

integrated_likelihood_self_norm_IS <- get_integrated_likelihood(data, 
                                                                weights, 
                                                                step_size, 
                                                                num_std_errors, 
                                                                R,
                                                                MC_params,
                                                                lambda_method,
                                                                init_lambda,
                                                                chunk_size)

toc()

self_norm_IS_diagnostics <- get_importance_sampling_diagnostics(integrated_likelihood_self_norm_IS)

################################################################################
############# INTEGRATED LIKELIHOOD - REGRESSION IMPORTANCE SAMPLING ###########
################################################################################

method = "regression_IS"

MC_params <- list(method = method, 
                  importance = list(rng = rng, 
                                    rng_params = importance_rng_params,
                                    density = density, 
                                    density_params = importance_density_params), 
                  nominal = list(rng = rng, 
                                 rng_params = nominal_rng_params,
                                 density = density, 
                                 density_params = nominal_density_params))

tic()

plan(multisession, workers = I(50))

integrated_likelihood_regression_IS <- get_integrated_likelihood(data, 
                                                                 weights, 
                                                                 step_size, 
                                                                 num_std_errors, 
                                                                 R,
                                                                 MC_params,
                                                                 lambda_method,
                                                                 init_lambda,
                                                                 chunk_size)

toc()

regression_IS_diagnostics <- get_importance_sampling_diagnostics(integrated_likelihood_regression_IS)

################################################################################
############################## PROFILE LIKELIHOOD ############################## 
################################################################################

tic()

plan(multisession, workers = I(2))

profile_log_likelihood_vals <- get_profile_log_likelihood(data, 
                                                          weights, 
                                                          step_size, 
                                                          num_std_errors)

toc()

################################################################################
################################### STORAGE #################################### 
################################################################################

psi_grid <- get_psi_grid(data, weights, step_size, num_std_errors, split = FALSE)

log_likelihood_vals <- data.frame(psi = psi_grid,
                                  Integrated = integrated_likelihood_vanilla_MC$l_bar,
                                  Profile = profile_log_likelihood_vals)

# log_likelihood_vals_file_path <- population_params_file_path |> 
#   sub(".*Pseudolikelihoods\\\\(.*)\\parameters\\.R", "\\1", x = _) |> 
#   paste0("Pseudolikelihoods\\\\", . = _, "\\log_likelihood_vals_") |> 
#   glue::glue("R={R}_seed={seed}_stepsize={step_size}_numse={num_std_errors}.Rda")

log_likelihood_vals_file_path <- population_params_file_path |>
  sub(".*Pseudolikelihoods/(.*)/parameters.R", "\\1", x = _) |>
  paste0("Pseudolikelihoods/", . = _, "/log_likelihood_vals_") |>
  glue::glue("R={R}_seed={seed}_stepsize={step_size}_numse={num_std_errors}.Rda")

saveRDS(log_likelihood_vals, log_likelihood_vals_file_path)

