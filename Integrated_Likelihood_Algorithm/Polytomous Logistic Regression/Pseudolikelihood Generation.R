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

Beta_prior <- list(mu = rep(0, J - 1),
                   Sigma = diag(J - 1))
Beta_rng <- MASS::mvrnorm
Beta_rng_params <- list(n = p, 
                        mu = Beta_prior$mu, 
                        Sigma = Beta_prior$Sigma)
# Beta_density <- dnorm
# Beta_density_params <- list(mean = Beta_prior[1], sd = Beta_prior[2])
Beta_dist_list <- list(rng = Beta_rng, 
                       rng_params = Beta_rng_params)
                       # density = Beta_density, 
                       # density_params = Beta_density_params)

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

model <- get_multinomial_logistic_model(data)

h <- 50L

X_h <- data |> 
  select(-Y) |> 
  slice(h) |>
  unname() |> 
  as.matrix() |> 
  (\(mat) cbind(1, mat))()

R <- 100

step_size <- 0.05

psi_grid <- get_psi_grid(step_size, model)

init_guess <- rep(1, length(coef(model)))

################################################################################
########################## INTEGRATED LIKELIHOOD - VANILLA MC ##################
################################################################################

# num_workers <- Sys.getenv("SLURM_NPROCS") |>
#   as.numeric()

# num_workers <- availableCores() |>
#   as.numeric()

# num_workers <- parallel::detectCores() |>
#   as.numeric()

num_workers <- 12

chunk_size <- 4

method <- "vanilla_MC"

MC_params <- list(method = method, 
                  nominal = Beta_dist_list)

tic()

plan(multisession, workers = I(num_workers))

log_integrated_likelihood_vanilla_MC <- get_log_integrated_likelihood(data,
                                                                      X_h, 
                                                                      psi_grid, 
                                                                      R,
                                                                      MC_params,
                                                                      init_guess,
                                                                      chunk_size)

toc()

################################################################################
############################## PROFILE LIKELIHOOD ############################## 
################################################################################

log_profile_likelihood_vals <- get_log_profile_likelihood(data,
                                                          X_h, 
                                                          psi_grid, 
                                                          init_guess)

################################################################################
################################### STORAGE #################################### 
################################################################################

log_likelihood_vals <- data.frame(psi = psi_grid,
                                  Integrated = log_integrated_likelihood_vanilla_MC$log_L_bar$estimate,
                                  Profile = log_profile_likelihood_vals)

log_likelihood_vals_file_path <- glue::glue("log_likelihood_vals_seed={seed}_R={R}_h={h}_stepsize={step_size}.Rda")

saveRDS(log_likelihood_vals, log_likelihood_vals_file_path)

plot(psi_grid, log_profile_likelihood_vals)

plot(psi_grid, log_integrated_likelihood_vanilla_MC$log_L_bar$estimate)
# 
# Rcpp::sourceCpp("accumulate_rcpp.cpp")
# 
# Beta_MLE <- get_Beta_MLE(model)
# 
# U_list <- get_U_list(MC_params, R)
# 
# omega_hat_list <- get_omega_hat_list(U_list, Beta_MLE, X_h)
# 
# accumulate_Beta_hats(psi_grid, omega_hat, X, X_h, init_guess)
# 
# Y_one_hot <- model.matrix( ~ factor(Y))[,-1]
# 
# get_log_L_tilde(psi_grid, omega_hat, X, Y_one_hot, X_h, init_guess)
# 
# registerDoFuture()
# plan(multisession)
# 
# parallel::clusterEvalQ(plan(), Rcpp::sourceCpp("accumulate_rcpp.cpp"))
# 
# get_log_L_tilde_mat(psi_grid, omega_hat_list, X, Y_one_hot, X_h, init_guess, chunk_size)

