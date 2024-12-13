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

alpha_prior <- c(0, 1)
alpha_rng <- rnorm
alpha_rng_params <- list(n = 1, mean = alpha_prior[1], sd = alpha_prior[2])
alpha_density <- dnorm
alpha_density_params <- list(mean = alpha_prior[1], sd = alpha_prior[2])
alpha_dist_list <- list(rng = alpha_rng, 
                        rng_params = alpha_rng_params, 
                        density = alpha_density, 
                        density_params = alpha_density_params)

beta_prior <- c(0, 1)
beta_rng <- rnorm
beta_rng_params <- list(n = 1, mean = beta_prior[1], sd = beta_prior[2])
beta_density <- dnorm
beta_density_params <- list(mean = beta_prior[1], sd = beta_prior[2])
beta_dist_list <- list(rng = beta_rng, 
                       rng_params = beta_rng_params, 
                       density = beta_density, 
                       density_params = beta_density_params)

sigma_squared_prior <- c(0.1, 0.1)
sigma_squared_rng <- rinvgamma
sigma_squared_rng_params <- list(n = 1, shape = sigma_squared_prior[1], rate = sigma_squared_prior[2])
sigma_squared_density <- dinvgamma
sigma_squared_density_params <- list(shape = sigma_squared_prior[1], rate = sigma_squared_prior[2])
sigma_squared_dist_list <- list(rng = sigma_squared_rng, 
                                rng_params = sigma_squared_rng_params, 
                                density = sigma_squared_density, 
                                density_params = sigma_squared_density_params)

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

theta_hat_method <- "accumulate"

chunk_size <- 5

step_size <- 0.01

num_std_errors <- 4

R <- 250

x_h <- 20

alpha_hat <- get_alpha_hat(x, y)

beta_hat <- get_beta_hat(x, y)

sigma_squared_hat <- get_sigma_squared_hat(x, y)

theta_MLE <- c(alpha_hat, beta_hat, sigma_squared_hat)

# init_guess <- c(rnorm(2), rinvgamma(1, 0.01, 0.01))

init_guess <- c(0, 0, 1)

# init_guess <- theta_MLE

################################################################################
########################## INTEGRATED LIKELIHOOD - VANILLA MC ##################
################################################################################

psi_grid_list <- get_psi_grid(x, y, x_h, step_size, num_std_errors, split = TRUE)

method = "vanilla_MC"

MC_params <- list(method = method, 
                  nominal = list(alpha_dist_list,
                                 beta_dist_list,
                                 sigma_squared_dist_list))

tic()

plan(multisession, workers = I(50))

log_integrated_likelihood_vanilla_MC <- get_log_integrated_likelihood(x,
                                                                      y,
                                                                      x_h, 
                                                                      psi_grid_list, 
                                                                      R,
                                                                      MC_params,
                                                                      theta_hat_method,
                                                                      init_guess,
                                                                      chunk_size)

toc()

################################################################################
############################## PROFILE LIKELIHOOD ############################## 
################################################################################

profile_log_likelihood_vals <- get_profile_log_likelihood(x, 
                                                          y, 
                                                          x_h,
                                                          step_size, 
                                                          num_std_errors)

################################################################################
################################### STORAGE #################################### 
################################################################################

psi_grid <- get_psi_grid(x, y, x_h, step_size, num_std_errors, split = FALSE)

log_likelihood_vals <- data.frame(psi = psi_grid,
                                  Integrated = log_integrated_likelihood_vanilla_MC$log_L_bar$estimate,
                                  Profile = profile_log_likelihood_vals)

log_likelihood_vals_file_path <- glue::glue("log_likelihood_vals_seed={seed}_R={R}_xh={as.character(x_h)}_stepsize={step_size}_numse={num_std_errors}.Rda")

saveRDS(log_likelihood_vals, log_likelihood_vals_file_path)

# plot(psi_grid, profile_log_likelihood_vals)
# 
# plot(psi_grid, log_integrated_likelihood_vanilla_MC$log_L_bar$estimate)








