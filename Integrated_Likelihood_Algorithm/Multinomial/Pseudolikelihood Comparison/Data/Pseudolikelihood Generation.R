library(future)
library(zeallot)
library(purrr)
library(dipsaus)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../../utils.R")

# population <- "Desert Rodents"
# population <- "Birds in Balrath Woods"
population <- "Birds in Killarney Woodlands"

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

n <- sum(data)

m <- length(data)

theta_MLE <- data / n

psi_MLE <- entropy(theta_MLE)

sigma <- theta_MLE*diag(m) - matrix(theta_MLE) %*% theta_MLE

psi_MLE_SE <- sqrt(sum(matrix(1 + log(theta_MLE)) %*% (1 + log(theta_MLE)) * sigma) / n)

num_std_errors <- 2.5
MoE <- num_std_errors * psi_MLE_SE

step_size <- 0.01

psi_grid <- psi_MLE %+-% MoE |> 
  rev() |> 
  (\(x) c(max(0, x[1]), min(log(m), x[2])))() |> 
  plyr::round_any(step_size, floor) |> 
  (\(x) seq(x[1], x[2], step_size))() 

R <- 500

tol <- 0.0001

################################################################################
############################ INTEGRATED LIKELIHOOD ############################# 
################################################################################

neg_log_likelihood <- function(theta, data) -sum(data * log(theta), na.rm = TRUE)

plan(multisession, workers = availableCores())

multinomial_entropy_values_IL <- neg_log_likelihood |> 
  get_omega_hat_list(psi_MLE, rep(1, length(data)), R, tol) |> 
  pluck("omega_hat") |> 
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

log_likelihood_vals_file_path <- population |> 
  tolower() |> 
  chartr(" ", "_", x = _) |> 
  glue::glue("_R={R}_step_size={step_size}.Rda")

saveRDS(log_likelihood_vals, paste0("Pseudolikelihoods/", log_likelihood_vals_file_path))

data.frame(psi = psi_grid,
           Integrated = multinomial_entropy_values_PL) |> 
  ggplot() +
  geom_point(aes(x = psi, y = Integrated))
