library(future)
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

plan(multisession, workers = 63)

multinomial_entropy_values_IL <- neg_log_likelihood |> 
    get_omega_hat_list(psi_MLE, rep(1, length(data)), R, tol) |> 
    pluck("omega_hat") |> 
    get_multinomial_entropy_values_IL(psi_grid_list, data)
  
################################################################################
######################## MODIFIED INTEGRATED LIKELIHOOD ########################
################################################################################

plan(sequential)

alpha <- data + 1

euclidean_distance <- function(omega, u) dist(matrix(c(omega, u), 
                                                     nrow = 2, 
                                                     byrow = TRUE),
                                              method = "euclid")[1]

c(u_list, omega_hat_list) %<-% get_omega_hat_list(euclidean_distance, psi_MLE, alpha, R, tol)

L <- u_list |> 
  map_dbl(\(u) likelihood(u, data)) |> 
  unlist() |> 
  as.numeric()

plan(multisession, workers = 63)

multinomial_entropy_values_mod_IL <- omega_hat_list |> 
    get_multinomial_entropy_values_modified_IL(psi_grid_list, L, data)

################################################################################
############################## PROFILE LIKELIHOOD ############################## 
################################################################################

plan(sequential)

multinomial_entropy_values_PL <- data |> 
    get_multinomial_entropy_values_PL(psi_grid_list)

################################################################################
################################### STORAGE #################################### 
################################################################################

psi_grid <- data |> 
  get_psi_grid(step_size, num_std_errors, split = FALSE)

log_likelihood_vals <- data.frame(psi = psi_grid,
                                  Mod_Integrated = multinomial_entropy_values_mod_IL,
                                  Integrated = multinomial_entropy_values_IL,
                                  Profile = multinomial_entropy_values_PL) 

log_likelihood_vals_file_path <- population |> 
  tolower() |> 
  chartr(" ", "_", x = _) |> 
  glue::glue("_R={R}_step_size={step_size}.Rda")

saveRDS(log_likelihood_vals, paste0("Pseudolikelihoods/", log_likelihood_vals_file_path))

