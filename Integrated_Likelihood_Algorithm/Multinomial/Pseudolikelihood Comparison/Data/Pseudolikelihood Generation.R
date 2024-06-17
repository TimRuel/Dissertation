library(future)
library(doFuture)
library(zeallot)
library(purrr)
library(dplyr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../../utils.R")

population <- "Desert Rodents"
# population <- "Birds in Balrath Woods"
# population <- "Birds in Killarney Woodlands"

switch(population,
       
       "Desert Rodents" = {
         
         data <- c(1, 1, 2, 4, 7, 10)
         },
       
       "Birds in Balrath Woods" = {
         
         data <- c(1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 6, 8)
         },
       
       "Birds in Killarney Woodlands" = {
         
         data <- c(1, 3, 4, 6, 7, 10, 14, 30)
         }
       )  

set.seed(38497283)

n <- sum(data)

m <- length(data)

theta_MLE <- data / n

psi_MLE <- entropy(theta_MLE)

step_size <- 0.01

num_std_errors <- 3

psi_grid_list <- data |> 
  get_psi_grid(step_size, num_std_errors, split = TRUE)

R <- 250

tol <- 0.0001

################################################################################
############################ INTEGRATED LIKELIHOOD ############################# 
################################################################################

plan(multisession, workers = availableCores(method = "system"))

alpha <- 1/2

u_params <- rep(alpha, m)

num_chunks <- 10

stime <- system.time({

integrated_log_likelihood_vals <- foreach(
  
  i = 1:R,
  .combine = "rbind",
  .multicombine = TRUE,
  .maxcombine = R,
  .options.future = list(seed = TRUE,
                         chunk.size = num_chunks)
  
) %dofuture% {
  
  neg_log_likelihood |>
    get_omega_hat(psi_MLE, u_params, tol, return_u = FALSE) |>
    get_integrated_log_likelihood(data, psi_grid_list)
} |> 
  matrixStats::colLogSumExps() |>
  (`-`)(log(R))

})

stime

################################################################################
######################## MODIFIED INTEGRATED LIKELIHOOD ########################
################################################################################

plan(multisession, workers = availableCores(method = "system"))

alpha <- 1/2

u_params <- data + alpha

Q_name <- "euclidean_distance"
# Q_name <- "neg_log_likelihood"

Q <- Q_name |>
  get()

mod_integrated_log_likelihood_vals <- foreach(
  
  i = 1:R,
  .combine = "rbind",
  .multicombine = TRUE,
  .maxcombine = R,
  .options.future = list(seed = TRUE)
  
) %dofuture% {
  
  c(u, omega_hat) %<-%  get_omega_hat(Q, psi_MLE, u_params, tol, return_u = TRUE)
  
  log_like_u <- log_likelihood(u, data)
  
  omega_hat |> 
    get_integrated_log_likelihood(data, psi_grid_list) |> 
    (`-`)(log_like_u)
  } |> 
  matrixStats::colLogSumExps() |> 
  (`-`)(log(R))

################################################################################
############################## PROFILE LIKELIHOOD ############################## 
################################################################################

plan(sequential)

profile_log_likelihood_vals <- data |> 
  get_profile_log_likelihood(psi_grid_list)

################################################################################
################################### STORAGE #################################### 
################################################################################

psi_grid <- data |> 
  get_psi_grid(step_size, num_std_errors, split = FALSE)

log_likelihood_vals <- data.frame(psi = psi_grid,
                                  Mod_Integrated = mod_integrated_log_likelihood_vals,
                                  Integrated = integrated_log_likelihood_vals,
                                  Profile = profile_log_likelihood_vals) 

log_likelihood_vals_file_path <- population |> 
  tolower() |> 
  chartr(" ", "_", x = _) |> 
  glue::glue("_R={R}_step_size={step_size}.Rda")

saveRDS(log_likelihood_vals, paste0("Pseudolikelihoods/", log_likelihood_vals_file_path))

