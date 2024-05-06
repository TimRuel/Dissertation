setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../../utils.R")

population <- "Desert Rodents"
# population <- "Birds in Balrath Woods"
# population <- "Birds in Killarney Woodlands"

switch(population,
       
       "Desert Rodents" = {
         
         seed <- 1996
         
         data <- c(1, 1, 2, 4, 7, 10)
         
         step_size <- 0.01
         
         log_likelihood_vals_file_path <- "desert_rodents_R=250_step_size=0.01.Rda"
         },
       
       "Birds in Balrath Woods" = {
         
         seed <- 39673
         
         data <- c(1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 6, 8)
         
         step_size <- 0.05
         
         log_likelihood_vals_file_path <- "birds_in_balrath_woods_R=250_step_size=0.05.Rda"
         },
       
       "Birds in Killarney Woodlands" = {
         
         seed <- 1996
         
         data <- c(1, 3, 4, 6, 7, 10, 14, 30)
         
         step_size <- 0.01
         
         log_likelihood_vals_file_path <- "birds_in_killarney_woodlands_R=250_step_size=0.01.Rda"
         }
       )  

set.seed(seed)

theta_MLE <- data / sum(data)

psi_MLE <- entropy(theta_MLE)

psi_grid <- data |> 
  length() |> 
  log() |> 
  plyr::round_any(step_size, floor) |> 
  seq(0, to = _, step_size)

R <- 250

################################################################################
############################ INTEGRATED LIKELIHOOD ############################# 
################################################################################

future::plan(multisession, workers = availableCores())

u_list_IL <- LaplacesDemon::rdirichlet(R, rep(1, length(data))) |> 
  t() |> 
  data.frame() |> 
  as.list()

omega_hat_list_IL <- u_list_IL |>
  purrr::map(get_omega_hat, psi_MLE, log_likelihood)
  
multinomial_entropy_values_IL <- omega_hat_list_IL |> 
  get_multinomial_entropy_values_IL(data, psi_grid)

################################################################################
######################## MODIFIED INTEGRATED LIKELIHOOD ########################
################################################################################

future::plan(multisession, workers = availableCores())

alpha <- data + 1

u_list_mod_IL <- LaplacesDemon::rdirichlet(R, alpha) |>
  t() |>
  data.frame() |>
  as.list()

omega_hat_list_mod_IL <- u_list_mod_IL |>
  purrr::map(get_omega_hat, psi_MLE, distance)

multinomial_entropy_values_modified_IL <- omega_hat_list_mod_IL |> 
    get_multinomial_entropy_values_modified_IL(u_list2, data, psi_grid)

################################################################################
############################## PROFILE LIKELIHOOD ############################## 
################################################################################

future::plan(sequential)

multinomial_entropy_values_PL <- data |> 
    get_multinomial_entropy_values_PL(psi_grid)

################################################################################
################################### STORAGE #################################### 
################################################################################

log_likelihood_vals <- data.frame(psi = psi_grid,
                                  Mod_Integrated = multinomial_entropy_values_modified_IL,
                                  Integrated = multinomial_entropy_values_IL,
                                  Profile = multinomial_entropy_values_PL) |> 
  tidyr::pivot_longer(cols = c("Mod_Integrated", "Integrated", "Profile"),
                      names_to = "Pseudolikelihood",
                      values_to = "loglikelihood")

saveRDS(log_likelihood_vals, paste0("Pseudolikelihoods/", log_likelihood_vals_file_path))
