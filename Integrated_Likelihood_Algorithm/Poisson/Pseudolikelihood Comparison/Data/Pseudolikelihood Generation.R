library(future)
library(zeallot)
library(purrr)
library(dplyr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../../utils.R")

set.seed(38497283)

# Define dimension of parameter
n <- 18

# Define true values of full model parameter
theta_0 <- n |> 
  runif() |> 
  round(2)

theta_0 <- theta_0 / sum(theta_0) * 10

# Define observed data from each population
data <- rpois(n, theta_0)

# Define weights for PoI function
weights <- m |> 
  runif() |>
  round(2)

weights <- weights / mean(weights)

theta_MLE <- data |> 
  map_dbl(mean)

psi_MLE <- weighted_sum(theta_MLE, weights)

step_size <- 0.01

num_std_errors <- 2

psi_grid_list <- data |> 
  get_psi_grid(weights, step_size, num_std_errors, split = TRUE)

R <- 250

tol <- 0.001

################################################################################
############################ INTEGRATED LIKELIHOOD ############################# 
################################################################################

alpha <- rep(1, m)

beta <- rep(0.01, m)

u_params <- list(alpha, beta)

plan(multisession, workers = future::availableCores())

poisson_weighted_sum_values_IL <- neg_log_likelihood |> 
  get_omega_hat_list(psi_MLE, weights, u_params, R, tol) |> 
  pluck("omega_hat") |> 
  get_poisson_weighted_sum_values_IL(data, weights, psi_grid_list)
  
################################################################################
######################## MODIFIED INTEGRATED LIKELIHOOD ########################
################################################################################

plan(sequential)

alpha <- data |> 
  map2_dbl(1/2, sum)

beta <- n

u_params <- list(alpha, beta)

c(u_list, omega_hat_list) %<-% get_omega_hat_list(neg_log_likelihood, psi_MLE, weights, u_params, R, tol)

l <- u_list |> 
  map_dbl(\(u) log_likelihood(u, data)) 

plan(multisession, workers = future::availableCores())

poisson_weighted_sum_values_mod_IL <- omega_hat_list |> 
  get_poisson_weighted_sum_values_modified_IL(data, weights, psi_grid_list, l)

################################################################################
############################## PROFILE LIKELIHOOD ############################## 
################################################################################

plan(sequential)

poisson_weighted_sum_values_PL <- data |> 
  get_poisson_weighted_sum_values_PL(weights, psi_grid_list)

################################################################################
################################### STORAGE #################################### 
################################################################################

psi_grid <- data |> 
  get_psi_grid(weights, step_size, num_std_errors, split = FALSE)

log_likelihood_vals <- data.frame(psi = psi_grid,
                                  Mod_Integrated = poisson_weighted_sum_values_mod_IL,
                                  Integrated = poisson_weighted_sum_values_IL,
                                  Profile = poisson_weighted_sum_values_PL) 

log_likelihood_vals_file_path <- "log_likelihood_vals_2.Rda"

saveRDS(log_likelihood_vals, paste0("Pseudolikelihoods/", log_likelihood_vals_file_path))

