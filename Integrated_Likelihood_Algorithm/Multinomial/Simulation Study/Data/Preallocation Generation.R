library(dplyr)
library(tidyr)
library(stringr)
library(furrr)
library(purrr)
library(dipsaus)

plan(sequential)

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

n <- sum(data)

################################################################################
############################ INTEGRATED LIKELIHOOD ############################# 
################################################################################

plan(sequential)

seed <- 38497283

set.seed(seed)

num_sims <- 2

data_sims <- num_sims |> 
  rmultinom(n, data) |> 
  data.frame() |> 
  as.list() |> 
  map(as.numeric)

R <- 250

tol <- 0.0001

plan(multisession, workers = availableCores())

IL_preallocations <- data_sims |>
  future_map(\(x) {
    
    psi_MLE <- entropy(x / sum(x))
    
    prior <- rep(1, length(x))
    
    omega_hat_list_IL <- neg_log_likelihood |>
      get_omega_hat_list(psi_MLE, prior, R, tol) |>
      pluck("omega_hat")
    
    return(omega_hat_list_IL)
  },
  .progress = TRUE)

IL_preallocations_file_path <- population |> 
  tolower() |> 
  str_replace_all(" ", "_") |> 
  glue::glue("_IL_preallocations_seed={seed}_numsims={num_sims}_R={R}_tol={tol}.Rda") |> 
  paste0("Preallocations/Integrated Likelihood/", ... = _)

saveRDS(IL_preallocations, IL_preallocations_file_path)

################################################################################
######################## MODIFIED INTEGRATED LIKELIHOOD ########################
################################################################################

seed <- 38497283

set.seed(seed)

num_sims <- 2

data_sims <- num_sims |> 
  rmultinom(n, data) |> 
  data.frame() |> 
  as.list() |> 
  map(as.numeric)

R <- 250

tol <- 0.0001

euclidean_distance <- function(u, omega) dist(matrix(c(u, omega),
                                                     nrow = 2,
                                                     byrow = TRUE),
                                              method = "euclid")[1]

mod_IL_preallocations <- data_sims |>
  future_map(\(x) {
    
      psi_MLE <- entropy(x / sum(x))
                                      
      alpha <- data + 1
      
      out <- get_omega_hat_list(euclidean_distance, psi_MLE, alpha, R, tol)
      
      return(out)
      },
      .progress = TRUE) |> 
  transpose()

mod_IL_preallocations_file_path <- population |> 
  tolower() |> 
  str_replace_all(" ", "_") |> 
  glue::glue("_mod_IL_preallocations_seed={seed}_numsims={num_sims}_R={R}_tol={tol}.Rda") |> 
  paste0("Preallocations/Modified Integrated Likelihood/", ... = _)

saveRDS(mod_IL_preallocations, mod_IL_preallocations_file_path)
