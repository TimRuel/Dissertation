library(future)
library(doFuture)
library(zeallot)
library(purrr)
library(plyr)
library(tidyverse)
library(invgamma)
library(stringr)
library(progressr)
library(tictoc)
library(rstudioapi)

handlers(global = TRUE)
handlers("cli")

# num_cores <- Sys.getenv("SLURM_NPROCS") |>
#   as.numeric()

# num_cores <- availableCores() |>
#   as.numeric()

num_cores <- parallel::detectCores() |>
  as.numeric()

setwd(dirname(getActiveDocumentContext()$path))

source("utils.R")

population_directory <- selectDirectory(caption = "Select population directory")

setwd(population_directory)

source("data.R")

set.seed(seed)

################################################################################
################################## PARAMETERS ################################## 
################################################################################

beta0_prior <- c(0, 1)
beta0_rng <- rnorm
beta0_rng_params <- list(n = 1, mean = beta0_prior[1], sd = beta0_prior[2])
beta0_density <- dnorm
beta0_density_params <- list(mean = beta0_prior[1], sd = beta0_prior[2])
beta0_dist_list <- list(rng = beta0_rng, 
                        rng_params = beta0_rng_params, 
                        density = beta0_density, 
                        density_params = beta0_density_params)

beta1_prior <- c(0, 1)
beta1_rng <- rnorm
beta1_rng_params <- list(n = 1, mean = beta1_prior[1], sd = beta1_prior[2])
beta1_density <- dnorm
beta1_density_params <- list(mean = beta1_prior[1], sd = beta1_prior[2])
beta1_dist_list <- list(rng = beta1_rng, 
                       rng_params = beta1_rng_params, 
                       density = beta1_density, 
                       density_params = beta1_density_params)

# nominal_rng_params <- list(n = n, shape = alpha_posterior, rate = beta_posterior)
# 
# nominal_density_params <- list(shape = alpha_posterior, rate = beta_posterior)

# importance_rng_rate_vec <- 2 
# 
# importance_rng_shape_vec <- importance_rng_rate_vec / beta_posterior * (alpha_posterior - 1) + 1
# 
# importance_rng_params <- list(n = n, shape = importance_rng_shape_vec, rate = importance_rng_rate_vec)
# 
# importance_density_params <- list(shape = importance_rng_shape_vec, rate = importance_rng_rate_vec)

beta_hat_method <- "accumulate"

chunk_size <- 5

step_size <- 0.01

num_std_errors <- 4

R <- 250

x_h <- 0.89805755

################################################################################
########################## INTEGRATED LIKELIHOOD - VANILLA MC ##################
################################################################################

psi_grid_list <- get_psi_grid(step_size, x, y, x_h, split = TRUE)

method = "vanilla_MC"

init_guess <- c(0, 0)

MC_params <- list(method = method, 
                  nominal = list(beta0_dist_list,
                                 beta1_dist_list))

tic()

plan(multisession, workers = I(50))

log_integrated_likelihood_vanilla_MC <- get_log_integrated_likelihood(x,
                                                                      y,
                                                                      x_h, 
                                                                      psi_grid_list, 
                                                                      R,
                                                                      MC_params,
                                                                      beta_hat_method,
                                                                      init_guess,
                                                                      chunk_size)

toc()

################################################################################
############################## PROFILE LIKELIHOOD ############################## 
################################################################################

init_guess <- c(0, 0)

profile_log_likelihood_vals <- get_profile_log_likelihood(x, 
                                                          y, 
                                                          x_h,
                                                          step_size, 
                                                          init_guess)

################################################################################
################################### STORAGE #################################### 
################################################################################

psi_grid <- get_psi_grid(step_size)

log_likelihood_vals <- data.frame(psi = psi_grid,
                                  Integrated = log_integrated_likelihood_vanilla_MC$log_L_bar$estimate,
                                  Profile = profile_log_likelihood_vals)

log_likelihood_vals_file_path <- glue::glue("log_likelihood_vals_seed={seed}_R={R}_xh={as.character(x_h)}_stepsize={step_size}_numse={num_std_errors}.Rda")

saveRDS(log_likelihood_vals, log_likelihood_vals_file_path)

plot(psi_grid, profile_log_likelihood_vals)
# 
# plot(psi_grid, log_integrated_likelihood_vanilla_MC$log_L_bar$estimate)








