library(dplyr)
library(tidyr)
library(stringr)
library(furrr)
library(purrr)
library(zeallot)
library(pushoverr)

plan(sequential)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../../utils.R")

IL_preallocations_file_path <- file.choose()
omega_hat_lists_IL <- readRDS(IL_preallocations_file_path)

mod_IL_preallocations_file_path <- file.choose()
c(u_lists_mod_IL, omega_hat_lists_mod_IL) %<-% readRDS(mod_IL_preallocations_file_path)

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

num_std_errors <- 3

step_size <- 0.01

start <- 1

end <- num_sims

chunk <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 

omega_hat_IL_batches <- chunk(omega_hat_lists_IL, 10)

data_batches <- chunk(data_sims, 10)

plan(multisession, workers = availableCores())

for (batch in 3:3) {
  
  multinomial_entropy_sims_IL <- omega_hat_IL_batches[[batch]] |>
    map2(data_batches[[batch]],
         \(omega_hat_list, data) get_multinomial_entropy_values_IL(omega_hat_list, data, step_size, num_std_errors),
         .progress = TRUE)
  
  multinomial_entropy_sims_IL_file_path <- IL_preallocations_file_path |> 
    str_remove("^.*preallocations_") |> 
    str_remove(".Rda") |> 
    paste0(glue::glue("_numse={num_std_errors}_stepsize={step_size}_batch={batch}.Rda")) |> 
    paste0("Simulations/", 
           population, 
           "/Integrated Likelihood/", 
           population |> 
             tolower() |> 
             str_replace_all(" ", "_"),
           "_IL_sims_",
           ... = _)
  
  saveRDS(multinomial_entropy_sims_IL, multinomial_entropy_sims_IL_file_path)
  
  pushover(paste0("IL Batch ", batch, " done!"))
}

# multinomial_entropy_sims_IL <- omega_hat_lists_IL[start:end] |>
#   map2(data_sims[start:end],
#        \(omega_hat_list, data) get_multinomial_entropy_values_IL(omega_hat_list, data, step_size, num_std_errors),
#        .progress = TRUE)
# 
# multinomial_entropy_sims_IL_file_path <- IL_preallocations_file_path |> 
#   str_remove("^.*preallocations_") |> 
#   str_remove(".Rda") |> 
#   paste0(glue::glue("_numse={num_std_errors}_stepsize={step_size}_batch={batch}.Rda")) |> 
#   paste0("Simulations/", 
#          population, 
#          "/Integrated Likelihood/", 
#          population |> 
#            tolower() |> 
#            str_replace_all(" ", "_"),
#          "_IL_sims_",
#          ... = _)
#   
# saveRDS(multinomial_entropy_sims_IL, multinomial_entropy_sims_IL_file_path)

################################################################################
######################## MODIFIED INTEGRATED LIKELIHOOD ########################
################################################################################

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

data_sims <- num_sims |> 
  rmultinom(n, data) |> 
  data.frame() |> 
  as.list() |> 
  map(as.numeric)

L_lists <- u_lists_mod_IL |> 
  map(\(u_list) u_list |> 
        map_dbl(\(u) likelihood(u, data)) |> 
  unlist() |> 
  as.numeric())

num_std_errors <- 3

step_size <- 0.01

start <- 1

end <- num_sims

omega_hat_mod_IL_batches <- chunk(omega_hat_lists_mod_IL, 10)

L_batches <- chunk(L_lists, 10)

data_batches <- chunk(data_sims, 10)

plan(multisession, workers = availableCores())

for (batch in 3:3) {
  
  multinomial_entropy_sims_mod_IL <- list(omega_hat_mod_IL_batches[[batch]],
                                          L_batches[[batch]],
                                          data_batches[[batch]]) |>
    pmap(\(omega_hat_list, L, data) get_multinomial_entropy_values_modified_IL(omega_hat_list, L, data, step_size, num_std_errors),
         .progress = TRUE)
  
  multinomial_entropy_sims_mod_IL_file_path <- mod_IL_preallocations_file_path |> 
    str_remove("^.*preallocations_") |> 
    str_remove(".Rda") |> 
    paste0(glue::glue("_numse={num_std_errors}_stepsize={step_size}_batch={batch}.Rda")) |> 
    paste0("Simulations/", 
           population, 
           "/Modified Integrated Likelihood/", 
           population |> 
             tolower() |> 
             str_replace_all(" ", "_"),
           "_mod_IL_sims_",
           ... = _)
  
  saveRDS(multinomial_entropy_sims_mod_IL, multinomial_entropy_sims_mod_IL_file_path)
  
  pushover(paste0("Mod IL Batch ", batch, " done!"))
}

# multinomial_entropy_sims_mod_IL <- list(omega_hat_lists_mod_IL[start:end],
#                                         L_lists[start:end],
#                                         data_sims[start:end]) |>
#   pmap(\(omega_hat_list, L, data) get_multinomial_entropy_values_modified_IL(omega_hat_list, L, data, step_size, num_std_errors),
#        .progress = TRUE)
# 
# multinomial_entropy_sims_mod_IL_file_path <- mod_IL_preallocations_file_path |> 
#   str_remove("^.*preallocations_") |> 
#   str_remove(".Rda") |> 
#   paste0(glue::glue("_numse={num_std_errors}_stepsize={step_size}.Rda")) |> 
#   paste0("Simulations/", 
#          population, 
#          "/Modified Integrated Likelihood/", 
#          population |> 
#            tolower() |> 
#            str_replace_all(" ", "_"),
#          "_mod_IL_sims_",
#          ... = _)
# 
# saveRDS(multinomial_entropy_sims_mod_IL, multinomial_entropy_sims_mod_IL_file_path)

################################################################################
############################## PROFILE LIKELIHOOD ############################## 
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

num_std_errors <- 3

step_size <- 0.01

start <- 1

end <- num_sims

data_batches <- chunk(data_sims, 10)

plan(multisession, workers = availableCores())

for (batch in 1:3) {

  multinomial_entropy_sims_PL <- data_batches[[batch]] |>
    future_map(\(data) get_multinomial_entropy_values_PL(data, step_size, num_std_errors),
               .progress = TRUE)

  multinomial_entropy_sims_PL_file_path <- population |>
    tolower() |>
    str_replace_all(" ", "_") |>
    glue::glue("_PL_sims_seed={seed}_numsims={num_sims}_numse={num_std_errors}_stepsize={step_size}_batch={batch}.Rda") |>
    paste0("Simulations/", population, "/Profile Likelihood/", ... = _)

  saveRDS(multinomial_entropy_sims_PL, multinomial_entropy_sims_PL_file_path)

  pushover(paste0("Profile Batch ", batch, " done!"))
}

# multinomial_entropy_sims_PL <- data_sims[start:end] |>
#     future_map(\(data) get_multinomial_entropy_values_PL(data, step_size, num_std_errors),
#                .progress = TRUE)
# 
# multinomial_entropy_sims_PL_file_path <- population |> 
#   tolower() |> 
#   str_replace_all(" ", "_") |> 
#   glue::glue("_PL_sims_seed={seed}_numsims={num_sims}_R={R}_tol={tol}.Rda") |> 
#   paste0("Simulations/", population, "/Profile Likelihood/", ... = _)
# 
# saveRDS(multinomial_entropy_sims_PL, multinomial_entropy_sims_PL_file_path)
