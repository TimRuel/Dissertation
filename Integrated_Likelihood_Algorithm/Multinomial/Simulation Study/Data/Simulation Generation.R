library(dplyr)
library(tidyr)
library(furrr)
library(purrr)
library(dipsaus)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../../utils.R")

population <- "Desert Rodents"
# population <- "Birds in Balrath Woods"
# population <- "Birds in Killarney Woodlands"

switch(population,     
       
       "Desert Rodents" = {

         data <- c(1, 1, 2, 4, 7, 10)
         
         omega_hat_lists_file_path <- "desert_rodents_omega_hat_lists.Rda"
         
         multinomial_entropy_sims_file_path <- "desert_rodents_sims.Rda"
         },
       
       "Birds in Balrath Woods" = {

         data <- c(1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 6, 8)
         
         omega_hat_lists_file_path <- "birds_in_balrath_woods_omega_hat_lists.Rda"
         
         multinomial_entropy_sims_file_path <- "birds_in_balrath_woods_sims.Rda"
         },
       
       "Birds in Killarney Woodlands" = {
         
         data <- c(1, 3, 4, 6, 7, 10, 14, 30)

         omega_hat_lists_file_path <- "birds_in_killarney_woodlands_omega_hat_lists.Rda"
         
         multinomial_entropy_sims_file_path <- "birds_in_killarney_woodlands_sims.Rda"
         }
       )  

set.seed(38497283)

num_sims <- 1000

n <- sum(data)

data_sims <- num_sims |> 
  rmultinom(n, data) |> 
  data.frame() |> 
  as.list()

m <- length(data)

num_std_errors <- 3

step_size <- 0.01

psi_grid_list <- data_sims |> 
  map(\(x) get_psi_grid(x, step_size, num_std_errors))

################################################################################
############################ INTEGRATED LIKELIHOOD ############################# 
################################################################################

R <- 250

tol <- 0.0001

plan(multisession, workers = availableCores())

omega_hat_lists <- data_sims |>
  future_map(\(x) {
    
    neg_log_likelihood <- function(theta, x) -sum(x * log(theta), na.rm = TRUE)
    
    psi_MLE <- entropy(x / sum(x))
    
    prior <- rep(1, length(x))
    
    omega_hat_list <- neg_log_likelihood |> 
      get_omega_hat_list(psi_MLE, prior, R, tol)
    
    return(omega_hat_list)
  },
  .progress = TRUE)

saveRDS(omega_hat_lists, paste0("Omega_hat Preallocations/", omega_hat_lists_file_path))

omega_hat_lists <- readRDS(paste0("Omega_hat Preallocations/", omega_hat_lists_file_path))

plan(multisession, workers = availableCores())

multinomial_entropy_sims_IL <- omega_hat_lists |>
  map2(data_sims,
       \(x, y) get_multinomial_entropy_values_IL(x, y, psi_grid),
       .progress = TRUE)

################################################################################
############################## PROFILE LIKELIHOOD ############################## 
################################################################################

plan(multisession, workers = availableCores())

multinomial_entropy_sims_PL <- data_sims |>
    future_map(\(x) get_multinomial_entropy_values_PL(x, psi_grid),
               .progress = TRUE)

################################################################################
################################### STORAGE #################################### 
################################################################################

multinomial_entropy_sims <- list("Integrated" = multinomial_entropy_sims_IL, 
                                 "Profile" = multinomial_entropy_sims_PL)

saveRDS(multinomial_entropy_sims, paste0("Simulations/", multinomial_entropy_sims_file_path))
