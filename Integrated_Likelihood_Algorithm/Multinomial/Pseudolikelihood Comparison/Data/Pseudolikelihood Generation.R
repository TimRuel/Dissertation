library(future)
library(zeallot)
library(purrr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../../utils.R")

population <- "Desert Rodents"
# population <- "Birds in Balrath Woods"
# population <- "Birds in Killarney Woodlands"

switch(population,
       
       "Desert Rodents" = {
         
         seed <- 7835
         
         data <- c(1, 1, 2, 4, 7, 10)
         
         step_size <- 0.01
         
         log_likelihood_vals_file_path <- "desert_rodents_R=250_step_size=0.01.Rda"
         },
       
       "Birds in Balrath Woods" = {
         
         seed <- 7835
         
         data <- c(1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 6, 8)
         
         step_size <- 0.01
         
         log_likelihood_vals_file_path <- "birds_in_balrath_woods_R=250_step_size=0.01.Rda"
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

neg_log_likelihood <- function(theta, data) -sum(data * log(theta), na.rm = TRUE)

plan(multisession, workers = availableCores())

multinomial_entropy_values_IL <- neg_log_likelihood |> 
  get_omega_hat_list(psi_MLE, rep(1, length(data)), R, tol) |> 
  pluck("omega_hat") |> 
  get_multinomial_entropy_values_IL(data, psi_grid)
  
  LaplacesDemon::rdirichlet(R, rep(1, length(data))) |> 
  t() |> 
  data.frame() |> 
  as.list() |> 
  map(\(u) make_objective_fn(u, neg_log_likelihood)) |>
  map(\(objective_fn) get_omega_hat(objective_fn, psi_MLE, init_guess)) |> 
  get_multinomial_entropy_values_IL(data, psi_grid)
  
################################################################################
######################## MODIFIED INTEGRATED LIKELIHOOD ########################
################################################################################

alpha <- data + 1

euclidean_distance <- function(u, omega) dist(matrix(c(u, omega), 
                                                     nrow = 2, 
                                                     byrow = TRUE),
                                              method = "euclid")[1]

c(u_list, omega_hat_list) %<-% get_omega_hat_list(euclidean_distance, psi_MLE, alpha, R, tol)

L <- u_list |> 
  map_dbl(\(u) likelihood(u, data)) |> 
  unlist() |> 
  as.numeric()

plan(multisession, workers = availableCores())

multinomial_entropy_values_mod_IL <- omega_hat_list |> 
  get_multinomial_entropy_values_modified_IL(L, data, psi_grid)

################################################################################
############################## PROFILE LIKELIHOOD ############################## 
################################################################################

plan(sequential)

multinomial_entropy_values_PL <- data |> 
    get_multinomial_entropy_values_PL(psi_grid)

################################################################################
################################### STORAGE #################################### 
################################################################################

log_likelihood_vals <- data.frame(psi = psi_grid,
                                  Mod_Integrated = multinomial_entropy_values_mod_IL,
                                  Integrated = multinomial_entropy_values_IL,
                                  Profile = multinomial_entropy_values_PL) 


saveRDS(log_likelihood_vals, paste0("Pseudolikelihoods/", log_likelihood_vals_file_path))

data.frame(psi = psi_grid,
           Integrated = multinomial_entropy_values_IL) |> 
  ggplot() +
  geom_point(aes(x = psi, y = Integrated))
