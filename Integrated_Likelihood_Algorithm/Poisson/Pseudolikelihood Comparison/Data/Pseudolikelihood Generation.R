library(future)
library(doFuture)
library(zeallot)
library(purrr)
library(dplyr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../../utils.R")

set.seed(7835)

# Define dimension of parameter
n <- 18

# Define true values of full model parameter
theta_0 <- n |> 
  runif() |> 
  round(2)

theta_0 <- theta_0 / sum(theta_0) * 10

# Define observed data from each population
data <- rpois(n, theta_0) |> 
  as.numeric()

# Define weights for PoI function
weights <- n |> 
  runif() |>
  round(2)

weights <- weights / mean(weights)

psi_MLE <- weighted_sum(data, weights)

step_size <- 0.5

num_std_errors <- 3

psi_grid_list <- data |> 
  get_psi_grid(weights, step_size, num_std_errors, split = TRUE)

R <- 10

tol <- 0.0001

################################################################################
############################ INTEGRATED LIKELIHOOD ############################# 
################################################################################

plan(multisession, workers = future::availableCores())

alpha <- 1/2

beta <- 1/2

omega_hat_list <- foreach(
  
  i = 1:R, 
  .combine = "list", 
  .multicombine = TRUE, 
  .maxcombine = R
  
  ) %dofuture% {
    
    neg_log_likelihood |> 
      get_omega_hat(psi_MLE, weights, alpha, beta, tol)
    }

poisson_weighted_sum_values_IL <- foreach(
  
  omega_hat = omega_hat_list,
  .combine = "rbind",
  .multicombine = TRUE,
  .maxcombine = R,
  .options.future = list(seed = TRUE)
  
  ) %dofuture% {
    
    omega_hat |> 
      get_poisson_weighted_sum_values_IL.aux(data, weights, psi_grid_list)
    } |> 
  matrixStats::colLogSumExps() |> 
  (`-`)(log(length(omega_hat_list)))
  
################################################################################
######################## MODIFIED INTEGRATED LIKELIHOOD ########################
################################################################################

plan(multisession, workers = future::availableCores())

alpha <- data + 1/2

beta <- 1 + 1/2

c(u_list, omega_hat_list) %<-% transpose(
  
  foreach(
    
    i = 1:R, 
    .combine = "list", 
    .multicombine = TRUE,
    .maxcombine = R
    
    ) %dofuture% {
      
      euclidean_distance |> 
        get_omega_hat(psi_MLE, weights, alpha, beta, tol, return_u = TRUE)
      }
  )

l <- u_list |> 
  map_dbl(\(u) log_likelihood(u, data)) 

poisson_weighted_sum_values_mod_IL <- foreach(
  
  omega_hat = omega_hat_list,
  .combine = "rbind",
  .multicombine = TRUE,
  .maxcombine = R
  
  ) %dofuture% {
  
  omega_hat |> 
    get_poisson_weighted_sum_values_IL.aux(data, weights, psi_grid_list)
    } |> 
  (`-`)(l) |>
  matrixStats::colLogSumExps() |> 
  (`-`)(log(length(omega_hat_list)))

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

log_likelihood_vals_file_path <- "log_likelihood_vals_6.Rda"

saveRDS(log_likelihood_vals, paste0("Pseudolikelihoods/", log_likelihood_vals_file_path))

