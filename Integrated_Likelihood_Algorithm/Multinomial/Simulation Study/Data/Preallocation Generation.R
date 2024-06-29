library(doFuture)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(pushoverr)
library(progressr)

handlers(global = TRUE)
handlers("cli")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../../utils.R")

num_cores <- Sys.getenv("SLURM_NPROCS") |>
  as.numeric()
# num_cores <- availableCores() |>
#   as.numeric()

# population <- "Desert Rodents"
# population <- "Birds in Balrath Woods"
# population <- "Birds in Killarney Woodlands"

populations <- c("Birds in Killarney Woodlands", "Birds in Balrath Woods")

for (population in populations) {
  
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
  
  if (population == "Birds in Killarney Woodlands") {
    
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
    
    u_params <- rep(alpha, m)
    
    R <- 250
    
    tol <- 0.0001
    
    num_chunks <- ceiling(R / num_cores)
    
    plan(multisession, workers = num_cores)
    
    IL_preallocations <- get_IL_preallocations(data_sims, u_params, R, tol) 
    
    IL_preallocations_file_path <- population |> 
      tolower() |> 
      str_replace_all(" ", "_") |> 
      glue::glue("_IL_preallocations_seed={seed}_numsims={num_sims}_R={R}_tol={tol}.Rda") |> 
      paste0("Preallocations/Integrated Likelihood/", ... = _)
    
    saveRDS(IL_preallocations, IL_preallocations_file_path)
    
    pushover(paste(population, "Integrated Likelihood Preallocation Finished"))
  }
  
  ################################################################################
  ######################## MODIFIED INTEGRATED LIKELIHOOD ########################
  ################################################################################
  
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
  
  Q_name <- "euclidean_distance"
  # Q_name <- "neg_log_likelihood"
  
  Q <- Q_name |>
    get()
  
  num_chunks <- ceiling(R / num_cores)
  
  plan(multisession, workers = num_cores)
  
  mod_IL_preallocations <- get_mod_IL_preallocations(data_sims, Q, alpha, R, tol) 
  
  Q_name <- Q_name |> 
    strsplit("_") |> 
    pluck(1) |> 
    substr(1, 1) |> 
    paste(collapse = "")
  
  mod_IL_preallocations_file_path <- population |> 
    tolower() |> 
    str_replace_all(" ", "_") |> 
    glue::glue("_mod_IL_preallocations_Q={Q_name}_seed={seed}_numsims={num_sims}_R={R}_tol={tol}.Rda") |> 
    paste0("Preallocations/Modified Integrated Likelihood/", ... = _)
  
  saveRDS(mod_IL_preallocations, mod_IL_preallocations_file_path)
  
  pushover(paste(population, "Modified Integrated Likelihood Preallocation Finished"))
}

