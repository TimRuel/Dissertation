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
  as.list() |> 
  map(as.numeric)

################################################################################
############################ INTEGRATED LIKELIHOOD ############################# 
################################################################################

# R <- 250
# 
# tol <- 0.0001
# 
# plan(multisession, workers = availableCores())
# 
# omega_hat_lists <- data_sims |>
#   future_map(\(x) {
#     
#     neg_log_likelihood <- function(theta, x) -sum(x * log(theta), na.rm = TRUE)
#     
#     psi_MLE <- entropy(x / sum(x))
#     
#     prior <- rep(1, length(x))
#     
#     omega_hat_list <- neg_log_likelihood |> 
#       get_omega_hat_list(psi_MLE, prior, R, tol) |> 
#       pluck("omega_hat") 
#     
#     return(omega_hat_list)
#   },
#   .progress = TRUE)

# seed = 38497283
# saveRDS(omega_hat_lists, paste0("Omega_hat Preallocations/", omega_hat_lists_file_path))

omega_hat_lists <- readRDS(paste0("Omega_hat Preallocations/", omega_hat_lists_file_path))

num_std_errors <- 3

step_size <- 0.01

start <- 1

end <- 1000

plan(multisession, workers = availableCores())

multinomial_entropy_sims_IL <- omega_hat_lists[start:end] |>
  map2(data_sims[start:end],
       \(omega_hat_list, data) get_multinomial_entropy_values_IL(omega_hat_list, data, step_size, num_std_errors),
       .progress = TRUE)

################################################################################
######################## MODIFIED INTEGRATED LIKELIHOOD ########################
################################################################################

# R <- 250
# 
# tol <- 0.0001
# 
# euclidean_distance <- function(u, omega) dist(matrix(c(u, omega), 
#                                                      nrow = 2, 
#                                                      byrow = TRUE),
#                                               method = "euclid")[1]
# 
# test <- data_sims[1:2] |>
#   future_map(\(x) {
# 
#     psi_MLE <- entropy(x / sum(x))
# 
#     alpha <- data + 1
# 
#     out <- get_omega_hat_list(euclidean_distance, psi_MLE, alpha, R, tol) 
#       
#     return(out)
#   },
#   .progress = TRUE)
# 
# c(u_lists, omega_hat_mod_IL_lists) %<-% transpose(test)
# 
# L <- u_list |> 
#   map_dbl(\(u) likelihood(u, data)) |> 
#   unlist() |> 
#   as.numeric()
# 
# plan(multisession, workers = availableCores())
# 
# multinomial_entropy_values_mod_IL <- omega_hat_list |> 
#   get_multinomial_entropy_values_modified_IL(L, data, step_size, num_std_errors)

################################################################################
############################## PROFILE LIKELIHOOD ############################## 
################################################################################

plan(multisession, workers = availableCores())

multinomial_entropy_sims_PL <- data_sims |>
    future_map(\(x) get_multinomial_entropy_values_PL(x, step_size, num_std_errors),
               .progress = TRUE)

################################################################################
################################### STORAGE #################################### 
################################################################################

multinomial_entropy_sims <- list("Integrated" = multinomial_entropy_sims_IL, 
                                 "Profile" = multinomial_entropy_sims_PL)

saveRDS(multinomial_entropy_sims, paste0("Simulations/", multinomial_entropy_sims_file_path))
