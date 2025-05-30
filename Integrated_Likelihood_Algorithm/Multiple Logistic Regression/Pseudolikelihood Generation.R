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

setwd(dirname(getActiveDocumentContext()$path))

source("utils.R")

population_directory <- selectDirectory(caption = "Select population directory")

setwd(population_directory)

source("data.R")

set.seed(seed)

################################################################################
################################## PARAMETERS ################################## 
################################################################################

Beta_prior <- c(0, 1)
Beta_rng <- rnorm
Beta_rng_params <- list(n = p, mean = Beta_prior[1], sd = Beta_prior[2])
Beta_density <- dnorm
Beta_density_params <- list(mean = Beta_prior[1], sd = Beta_prior[2])
Beta_dist_list <- list(rng = Beta_rng, 
                       rng_params = Beta_rng_params, 
                       density = Beta_density, 
                       density_params = Beta_density_params)

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

model <- get_logistic_regression_model(data)

h <- 50L

X_h <- data |> 
  select(-Y) |> 
  slice(h) |>
  unname() |> 
  as.matrix() |> 
  t()

plogis(get_Y_hat(model, X_h))

R <- 250

N <- 10^7

step_size <- 0.01

psi_grid <- get_psi_grid(step_size)

################################################################################
########################## INTEGRATED LIKELIHOOD - VANILLA MC ##################
################################################################################

# num_workers <- Sys.getenv("SLURM_NPROCS") |>
#   as.numeric()

# num_workers <- availableCores() |>
#   as.numeric()

# num_workers <- parallel::detectCores() |>
#   as.numeric()

num_workers <- 50

chunk_size <- 5

method <- "vanilla_MC"

MC_params <- list(method = method, 
                  nominal = Beta_dist_list)

tic()

plan(multisession, workers = I(num_workers))

log_integrated_likelihood_vanilla_MC <- get_log_integrated_likelihood(data,
                                                                      X_h, 
                                                                      psi_grid, 
                                                                      R,
                                                                      N,
                                                                      MC_params,
                                                                      chunk_size)

toc()

################################################################################
############################## PROFILE LIKELIHOOD ############################## 
################################################################################

profile_log_likelihood_vals <- get_profile_log_likelihood(data, X_h, step_size)

################################################################################
################################### STORAGE #################################### 
################################################################################

log_likelihood_vals <- data.frame(psi = psi_grid,
                                  Integrated = log_integrated_likelihood_vanilla_MC$log_L_bar$estimate,
                                  Profile = profile_log_likelihood_vals)

log_likelihood_vals_file_path <- glue::glue("log_likelihood_vals_seed={seed}_R={R}_h={h}_stepsize={step_size}.Rda")

saveRDS(log_likelihood_vals, log_likelihood_vals_file_path)

# plot(psi_grid, profile_log_likelihood_vals)
# 
# plot(psi_grid, log_integrated_likelihood_vanilla_MC$log_L_bar$estimate)








