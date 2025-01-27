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

X_level <- 2

X_h <- data.frame(X = factor(X_level))

step_size <- 0.01

# threshold <- 80

threshold <- 170

alpha <- 0.03

noise_sd <- 0.1

prop <- 0.9

verbose <- FALSE

################################################################################
########################## INTEGRATED LIKELIHOOD - VANILLA MC ##################
################################################################################

# num_workers <- Sys.getenv("SLURM_NPROCS") |>
#   as.numeric()

# num_workers <- availableCores() |>
#   as.numeric()

num_workers <- parallel::detectCores() |>
  as.numeric()

# num_workers <- 12

R <- num_workers * 2

Beta_MLE <- get_Beta_MLE(model)

init_guess <- c(Beta_MLE)

# num_workers <- 20

chunk_size <- 1

method <- "vanilla_MC"

MC_params <- list(method = method, 
                  nominal = Beta_dist_list)

tic()

plan(multisession, workers = I(num_workers))

log_integrated_likelihood_vanilla_MC <- get_log_integrated_likelihood(data,
                                                                      X_h, 
                                                                      R,
                                                                      MC_params,
                                                                      step_size, 
                                                                      threshold,
                                                                      alpha,
                                                                      prop,
                                                                      chunk_size,
                                                                      noise_sd,
                                                                      verbose)

toc()

################################################################################
############################## PROFILE LIKELIHOOD ############################## 
################################################################################



tic()

log_profile_likelihood <- get_log_profile_likelihood(data,
                                                     X_h,
                                                     step_size,
                                                     alpha,
                                                     noise_sd,
                                                     verbose)

toc()

################################################################################
################################### STORAGE #################################### 
################################################################################

log_likelihood_vals <- log_integrated_likelihood_vanilla_MC$values_df |> 
  merge(log_profile_likelihood$values_df, all = TRUE)

log_likelihood_vals_file_path <- glue::glue("log_likelihood_vals_seed={seed}_R={R}_h={X_level}_stepsize={step_size}.Rda")

saveRDS(log_likelihood_vals, log_likelihood_vals_file_path)

# matplot(log_integrated_likelihood_vanilla_MC$log_L_tilde_mat |> colnames() |> as.numeric(),
#         log_integrated_likelihood_vanilla_MC$log_L_tilde_mat |> t(),
#         type )

# plot(log_profile_likelihood$values_df)
# psi_hat <- get_psi_hat(model, X_h)
# abline(v = psi_hat, col = "green")
# 
# plot(log_integrated_likelihood_vanilla_MC$values_df)
# 
# plot(log_integrated_likelihood_vanilla_MC$psi_grid, log_integrated_likelihood_vanilla_MC$log_L_bar$estimate)
# abline(v = psi_hat, col = "green")
# abline(v = theta_0[[1]] |> get_entropy(), col = "red")
# abline(v = log_integrated_likelihood_vanilla_MC$psi_grid[which.max(log_integrated_likelihood_vanilla_MC$log_L_bar$estimate)], col = "blue")
# 
# for (i in nrow(log_integrated_likelihood_vanilla_MC$log_L_tilde_mat):1) {
# 
#   plot(log_integrated_likelihood_vanilla_MC$psi_grid, log_integrated_likelihood_vanilla_MC$log_L_tilde_mat[i,])
#   title(main = i)
#   abline(v = psi_hat, col = "green")
#   abline(v = theta_0[[1]] |> get_entropy(), col = "red")
#   abline(v = log_integrated_likelihood_vanilla_MC$psi_grid[which.max(log_integrated_likelihood_vanilla_MC$log_L_tilde_mat[i,])], col = "blue")
# }
# 
# for (i in nrow(log_L_tilde_mat):1) {
# 
#   plot(psi_grid, log_L_tilde_mat[i,])
#   title(main = i)
#   abline(v = psi_hat, col = "green")
#   abline(v = theta_0[[1]] |> get_entropy(), col = "red")
#   abline(v = psi_grid[which.max(log_L_tilde_mat[i,])], col = "blue")
# }
# 
# for (i in 2:1) {
# 
#   plot(psi_grid, log_L_tilde_mat[i,])
#   title(main = i)
#   abline(v = psi_hat, col = "green")
#   abline(v = theta_0[[1]] |> get_entropy(), col = "red")
#   abline(v = psi_grid[which.max(log_L_tilde_mat[i,])], col = "blue")
# }
# 

