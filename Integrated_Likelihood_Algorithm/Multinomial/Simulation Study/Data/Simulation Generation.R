library(dplyr)
library(tidyr)
library(stringr)
library(doFuture)
library(purrr)
library(zeallot)
library(pushoverr)
library(progressr)
library(tictoc)

handlers(global = TRUE)
handlers("cli")

plan(sequential)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../../utils.R")

num_cores <- Sys.getenv("SLURM_NPROCS") |> 
  as.numeric()# num_cores <- availableCores() |> 
#   as.numeric()

IL_preallocations_file_path <- file.choose()
IL_preallocations <- readRDS(IL_preallocations_file_path)

mod_IL_preallocations_file_path <- file.choose()
mod_IL_preallocations <- readRDS(mod_IL_preallocations_file_path)

population <- IL_preallocations_file_path |>  
  str_remove("^.*/") |> 
  str_remove("_IL.*$") |> 
  str_replace_all("_", " ") |> 
  tools::toTitleCase()

# population <- IL_preallocations_file_path |>  
#   str_remove("^.*\\\\") |> 
#   str_remove("_IL.*$") |> 
#   str_replace_all("_", " ") |> 
#   tools::toTitleCase()

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

m <- length(data)

################################################################################
############################ INTEGRATED LIKELIHOOD ############################# 
################################################################################

plan(sequential)

seed <- IL_preallocations_file_path |>  
  str_remove("^.*seed=") |> 
  str_extract("\\d+") |> 
  as.numeric()

set.seed(seed)

num_sims <- IL_preallocations_file_path |>
  str_remove("^.*numsims=") |>
  str_extract("\\d+") |>
  as.numeric()

data_sims <- num_sims |> 
  rmultinom(n, data) |> 
  data.frame() |> 
  as.list() |> 
  map(as.numeric)

R <- IL_preallocations_file_path |>  
  str_remove("^.*R=") |> 
  str_extract("\\d+") |> 
  as.numeric()

step_size <- 0.01

num_std_errors <- 4

num_chunks <- ceiling(R * num_sims / num_cores)

plan(multisession, workers = num_cores)

integrated_log_likelihood_sims <- get_integrated_log_likelihood_sims(data_sims, omega_hat_lists_IL, R, step_size, num_std_errors, num_chunks)
  
integrated_log_likelihood_sims_file_path <- "seed={seed}_numsims={num_sims}_R={R}_tol={tol}_numse={num_std_errors}_stepsize={step_size}.Rda" |> 
  glue::glue() |> 
  paste0("Simulations/",
         population,
         "/Integrated Likelihood/IL_sims_",
         ... = _)

saveRDS(integrated_log_likelihood_sims, integrated_log_likelihood_sims_file_path)

pushover("Integrated Likelihood Sims Done!")

################################################################################
######################## MODIFIED INTEGRATED LIKELIHOOD ########################
################################################################################

plan(sequential)

seed <- 38497283

set.seed(seed)

num_sims <- 1000

data_sims <- num_sims |> 
  rmultinom(n, data) |> 
  data.frame() |> 
  as.list() |> 
  map(as.numeric)

alpha <- 1/2

R <- 250

tol <- 0.0001

num_std_errors <- 4

step_size <- 0.01

num_chunks <- ceiling(R / num_cores)

# Q_name <- "euclidean_distance"
Q_name <- "neg_log_likelihood"

Q <- Q_name |>
  get()

plan(multisession, workers = num_cores)

tic()

mod_integrated_log_likelihood_sims <-
  
  foreach(
    data = data_sims,
    .combine = "list",
    .multicombine = TRUE,
    .maxcombine = num_sims,
    .options.future = list(seed = TRUE,
                           chunk.size = num_chunks)
    
  ) %:%
  
  foreach(
    i = 1:R,
    .combine = "rbind",
    .multicombine = TRUE,
    .maxcombine = R,
    .options.future = list(seed = TRUE)
    
  ) %dofuture% {
    
    if (i == R) pushover("Still running")
    
    psi_grid_list <-
      get_psi_grid(data, step_size, num_std_errors, split = TRUE)
    
    psi_MLE <- entropy(data / sum(data))
    
    u_params <- data + alpha
    
    c(u, omega_hat) %<-% get_omega_hat(Q, psi_MLE, u_params, tol, return_u = TRUE)
    
    log_like_u <- log_likelihood(u, data)
    
    omega_hat |>
      get_integrated_log_likelihood(data, psi_grid_list) |>
      (`-`)(log_like_u)
  } |>
  map(\(x) x |>
        matrixStats::colLogSumExps() |>
        (`-`)(log(R)))

toc()

Q_name <- Q_name |> 
  strsplit("_") |> 
  pluck(1) |> 
  substr(1, 1) |> 
  paste(collapse = "")

mod_integrated_log_likelihood_sims_file_path <- "seed={seed}_Q={Q_name}_numsims={num_sims}_R={R}_tol={tol}_numse={num_std_errors}_stepsize={step_size}.Rda" |> 
  glue::glue() |> 
  paste0("Simulations/",
         population,
         "/Modified Integrated Likelihood/mod_IL_sims_",
         ... = _)

saveRDS(mod_integrated_log_likelihood_sims, mod_integrated_log_likelihood_sims_file_path)

pushover("Modified Integrated Likelihood Sims Done!")

################################################################################
############################## PROFILE LIKELIHOOD ############################## 
################################################################################

plan(sequential)

handlers(global = TRUE)
handlers("progress")

seed <- 38497283

set.seed(seed)

num_sims <- 1000

data_sims <- num_sims |>
  rmultinom(n, data) |>
  data.frame() |>
  as.list() |>
  map(as.numeric)

num_std_errors <- 4

step_size <- 0.01

num_chunks <- ceiling(R / num_cores)

gen_PL_sims <- function(sims) {
  
  p <- progressor(along = sims)
  
  foreach(
    
    data = sims,
    .combine = "list",
    .multicombine = TRUE,
    .maxcombine = num_sims,
    .options.future = list(seed = TRUE,
                           chunk.size = num_chunks)
    
  ) %dofuture% {
    
    psi_grid_list <- get_psi_grid(data, step_size, num_std_errors, split = TRUE)
    p()
    data |>
      get_profile_log_likelihood(psi_grid_list)
  }
}

plan(multisession, workers = num_cores)

profile_log_likelihood_sims <- gen_PL_sims(data_sims)

profile_log_likelihood_sims_file_path <- "seed={seed}_numsims={num_sims}_numse={num_std_errors}_stepsize={step_size}.Rda" |> 
  glue::glue() |> 
  paste0("Simulations/",
         population,
         "/Profile Likelihood/PL_sims_",
         ... = _)

saveRDS(profile_log_likelihood_sims, profile_log_likelihood_sims_file_path)

pushover("Profile Likelihood Sims Done!")
