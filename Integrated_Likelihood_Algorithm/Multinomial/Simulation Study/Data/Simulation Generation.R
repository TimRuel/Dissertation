library(dplyr)
library(tidyr)
library(furrr)
library(purrr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../../utils.R")

# population <- "Desert Rodents"
population <- "Birds in Balrath Woods"
# population <- "Birds in Killarney Woodlands"

switch(population,     
       
       "Desert Rodents" = {
         
         seed <- 1996
         
         data <- c(1, 1, 2, 4, 7, 10)
         
         step_size <- 0.01
         
         omega_hat_lists_file_path <- "desert_rodents_omega_hat_lists.Rda"
         
         multinomial_entropy_sims_file_path <- "desert_rodents_sims.Rda"
         },
       
       "Birds in Balrath Woods" = {
         
         seed <- 39673
         
         data <- c(1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 6, 8)
         
         step_size <- 0.01
         
         omega_hat_lists_file_path <- "birds_in_balrath_woods_omega_hat_lists.Rda"
         
         multinomial_entropy_sims_file_path <- "birds_in_balrath_woods_sims.Rda"
         },
       
       "Birds in Killarney Woodlands" = {
         
         seed <- 1996
         
         data <- c(1, 3, 4, 6, 7, 10, 14, 30)
         
         step_size <- 0.01
         
         omega_hat_lists_file_path <- "birds_in_killarney_woodlands_omega_hat_lists.Rda"
         
         multinomial_entropy_sims_file_path <- "birds_in_killarney_woodlands_sims.Rda"
         }
       )  

set.seed(seed)

psi_grid <- data |> 
  length() |> 
  log() |> 
  plyr::round_any(step_size, floor) |> 
  seq(0, to = _, step_size)

n_sims <- 1000

data_sims <- n_sims |> 
  rmultinom(sum(data), data) |> 
  data.frame() |> 
  as.list()

################################################################################
############################ INTEGRATED LIKELIHOOD ############################# 
################################################################################

# R <- 250
# 
# u_list <- n_sims |>
#   replicate({
#     R |>
#       LaplacesDemon::rdirichlet(rep(1, length(data)))|>
#       t() |>
#       data.frame() |>
#       as.list()
#     },
#     simplify = FALSE)
# 
# plan(multisession, workers = availableCores())
# 
# omega_hat_lists <- u_list |>
#   future_map2(data_sims,
#               \(x, y) x |>
#                 map(\(z) z |>
#                       get_omega_hat(PoI_fn(y / sum(y)))),
#               .progress = TRUE)
# 
# saveRDS(omega_hat_lists, paste0("Omega_hat Preallocations/", omega_hat_lists_file_path))

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
