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
  as.numeric()
# num_cores <- availableCores() |> 
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

R <- IL_preallocations_file_path |>  
  str_remove("^.*R=") |> 
  str_extract("\\d+") |> 
  as.numeric()

step_size <- 0.01

num_std_errors <- 5

num_chunks <- ceiling(R * num_sims / num_cores)

plan(multisession, workers = num_cores)

integrated_log_likelihood_sims <- get_integrated_log_likelihood_sims(IL_preallocations, step_size, num_std_errors, num_chunks)
  
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

population <- mod_IL_preallocations_file_path |>  
  str_remove("^.*/") |> 
  str_remove("_mod_IL.*$") |> 
  str_replace_all("_", " ") |> 
  tools::toTitleCase()

plan(sequential)

seed <- mod_IL_preallocations_file_path |>  
  str_remove("^.*seed=") |> 
  str_extract("\\d+") |> 
  as.numeric()

set.seed(seed)

num_sims <- mod_IL_preallocations_file_path |>
  str_remove("^.*numsims=") |>
  str_extract("\\d+") |>
  as.numeric()

R <- mod_IL_preallocations_file_path |>  
  str_remove("^.*R=") |> 
  str_extract("\\d+") |> 
  as.numeric()

num_std_errors <- 5

step_size <- 0.01

Q_name <- sub(".*Q=(.*?)_seed.*", "\\1", mod_IL_preallocations_file_path)

num_chunks <- ceiling(R * num_sims/ num_cores)

plan(multisession, workers = num_cores)

mod_integrated_log_likelihood_sims <- get_mod_integrated_log_likelihood_sims(mod_IL_preallocations, step_size, num_std_errors, num_chunks)

mod_integrated_log_likelihood_sims_file_path <- "seed={seed}_Q={Q_name}_numsims={num_sims}_R={R}_numse={num_std_errors}_stepsize={step_size}.Rda" |> 
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



filepaths = list.files("Simulations/Birds in Balrath Woods/Integrated Likelihood",
                       pattern = ".*batch=.*") |>
  sort()

sims_list <- list()

for (path in filepaths) {

  temp <- readRDS(paste0("Simulations/Birds in Balrath Woods/Integrated Likelihood/", path))

  sims_list <- c(sims_list, temp)
}

saveRDS(sims_list, "Simulations/Birds in Balrath Woods/Integrated Likelihood/IL_sims_seed=38497283_numsims=1000_R=250_numse=5_stepsize=0.01.Rda")

filepaths = list.files("Simulations/Birds in Balrath Woods/Modified Integrated Likelihood",
                       pattern = ".*batch=.*") |>
  sort()

sims_list <- list()

for (path in filepaths) {
  
  temp <- readRDS(paste0("Simulations/Birds in Balrath Woods/Modified Integrated Likelihood/", path))
  
  sims_list <- c(sims_list, temp)
}

saveRDS(sims_list, "Simulations/Birds in Balrath Woods/Modified Integrated Likelihood/mod_IL_sims_seed=38497283_Q=nll_numsims=1000_R=250_numse=5_stepsize=0.01.Rda")






