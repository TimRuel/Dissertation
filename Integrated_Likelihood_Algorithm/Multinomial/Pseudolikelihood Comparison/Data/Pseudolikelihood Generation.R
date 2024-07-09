library(future)
library(doFuture)
library(zeallot)
library(purrr)
library(dplyr)
library(progressr)
library(tictoc)

handlers(global = TRUE)
handlers("cli")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../../utils.R")

num_cores <- Sys.getenv("SLURM_NPROCS") |> 
  as.numeric()
# num_cores <- availableCores() |> 
#   as.numeric()

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

m <- length(data)

step_size <- 0.001

num_std_errors <- 4

R <- 250

tol <- 0.0001

################################################################################
############################ INTEGRATED LIKELIHOOD ############################# 
################################################################################

plan(multisession, workers = num_cores)

alpha <- 1/2

u_params <- rep(alpha, m)

# num_chunks <- ceiling(R / num_cores)

num_chunks <- 15

integrated_log_likelihood_vals <- data |> 
  get_integrated_log_likelihood_vals(step_size, num_std_errors, u_params, R, tol, num_chunks)

################################################################################
######################## MODIFIED INTEGRATED LIKELIHOOD ########################
################################################################################

alpha <- 1/2

u_params <- data + alpha

# Q_name <- "euclidean_distance"
Q_name <- "neg_log_likelihood"

Q <- Q_name |>
  get()

num_chunks <- ceiling(R / num_cores)

# num_chunks <- 15

plan(multisession, workers = num_cores)

tic()

mod_integrated_log_likelihood_vals <- get_mod_integrated_log_likelihood_vals(data,
                                                                             Q, 
                                                                             step_size, 
                                                                             num_std_errors, 
                                                                             u_params, 
                                                                             R, 
                                                                             tol, 
                                                                             num_chunks)

toc()

################################################################################
############################## PROFILE LIKELIHOOD ############################## 
################################################################################

plan(sequential)

psi_grid_list <- data |> 
  get_psi_grid(step_size, num_std_errors, split = TRUE)

profile_log_likelihood_vals <- data |> 
  get_profile_log_likelihood(psi_grid_list)

################################################################################
################################### STORAGE #################################### 
################################################################################

psi_grid <- data |> 
  get_psi_grid(step_size, num_std_errors, split = FALSE)

log_likelihood_vals <- data.frame(psi = psi_grid,
                                  Integrated = integrated_log_likelihood_vals,
                                  Mod_Integrated = mod_integrated_log_likelihood_vals,
                                  Profile = profile_log_likelihood_vals) 

log_likelihood_vals_file_path <- population |> 
  tolower() |> 
  chartr(" ", "_", x = _) |> 
  glue::glue("_R={R}_step_size={step_size}.Rda")

saveRDS(log_likelihood_vals, paste0("Pseudolikelihoods/", log_likelihood_vals_file_path))

