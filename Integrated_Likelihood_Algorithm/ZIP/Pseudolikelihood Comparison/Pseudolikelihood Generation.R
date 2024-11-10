# library(future)
# library(doFuture)
# library(zeallot)
# library(purrr)
# library(dplyr)
# library(plyr)
# library(stringr)
# library(progressr)
# library(tictoc)
library(tidyverse)
library(purrr)

# handlers(global = TRUE)
# handlers("cli")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# num_cores <- Sys.getenv("SLURM_NPROCS") |>
#   as.numeric()

# num_cores <- availableCores() |>
#   as.numeric()

# num_cores <- parallel::detectCores() |>
#   as.numeric()

################################################################################
#################################### DATA ###################################### 
################################################################################

print("Choose your population parameters R script.")
population_params_file_path <- file.choose() 

source(population_params_file_path)
source("../utils.R")

# set.seed(seed)

U <- rbinom(n, 1, 1 - pi_0)
V <- rpois(n, mu_0)
Y <- U * V

################################################################################
############################ INTEGRATED LIKELIHOOD ############################# 
################################################################################

step_size <- 0.01

num_std_errors <- 4

mu_grid <- get_mu_grid(Y, step_size, num_std_errors)

mu_hat <- get_mu_hat(Y)

R <- 3

alpha <- 2

beta <- 2

phi_vals <- rbeta(R, alpha, beta)

x <- phi_vals |> 
  map(\(phi) mu_grid |> 
        map_dbl(\(mu) get_rho(mu, phi, mu_hat) |> 
                  (\(rho) likelihood(mu, rho, Y))()
        )
  )

L_bar <- phi_vals |> 
  map(\(phi) mu_grid |> 
        map_dbl(\(mu) get_rho(mu, phi, mu_hat) |> 
                  (\(rho) likelihood(mu, rho, Y))()
        )
  ) |> 
  unlist() |> 
  matrix(ncol = length(mu_grid),
         byrow = TRUE) |> 
  colMeans()

l_bar <- log(L_bar)

Y_bar <- mean(Y)
pi_bar <- mean(Y == 0)

l_p <- function(mu) n * (Y_bar * log(mu) - (1 - pi_bar) * (log(1 - exp(-mu)) + mu) + pi_bar * log(pi_bar) + (1 - pi_bar) * log(1 - pi_bar))

l_p_vals <- map_dbl(mu_grid, l_p)
mu_hat_p <- mu_vals[which.max(l_p_vals)]
plot(mu_grid, l_p_vals - max(l_p_vals), type = "l", col = "blue")
lines(mu_grid, l_bar - max(l_bar), col = "orange")
# abline(v = mu_hat, col = "green")
# abline(v = mu_0, col = "red")

x <- phi_vals |> 
  map(\(phi) mu_grid |> 
        map_dbl(\(mu) get_rho(mu, phi, mu_hat) |> 
                  (\(rho) likelihood(mu, rho, Y))()
        )
  )
                  
        
log(x[[2]]) - max(log(x[[2]]))

x1 <- (phi_vals[[1]] + (1 - phi_vals[[1]])*exp(-mu_hat))^(n * pi_bar) * (1 - phi_vals[[1]])^(n * (1 - pi_bar))
x2 <- (phi_vals[[2]] + (1 - phi_vals[[2]])*exp(-mu_hat))^(n * pi_bar) * (1 - phi_vals[[2]])^(n * (1 - pi_bar))
  

################################################################################
############################## PROFILE LIKELIHOOD ############################## 
################################################################################

tic()

plan(multisession, workers = I(2))

profile_log_likelihood_vals <- get_profile_log_likelihood(data, 
                                                          weights, 
                                                          step_size, 
                                                          num_std_errors)

toc()

################################################################################
################################### STORAGE #################################### 
################################################################################

psi_grid <- get_psi_grid(data, weights, step_size, num_std_errors, split = FALSE)

log_likelihood_vals <- data.frame(psi = psi_grid,
                                  Integrated = integrated_likelihood_vanilla_MC$l_bar,
                                  Profile = profile_log_likelihood_vals)

log_likelihood_vals_file_path <- population_params_file_path |> 
  sub(".*Pseudolikelihoods/(.*)/parameters\\.R", "\\1", x = _) |> 
  paste0("Pseudolikelihoods/", . = _, "/log_likelihood_vals_") |> 
  glue::glue("R={R}_seed={seed}_stepsize={step_size}_numse={num_std_errors}.Rda")

saveRDS(log_likelihood_vals, log_likelihood_vals_file_path)

