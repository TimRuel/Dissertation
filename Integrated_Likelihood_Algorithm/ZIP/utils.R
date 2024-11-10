################################################################################
#################################### GENERAL ###################################
################################################################################
library(tidyverse)

log_likelihood <- function(mu, rho, y) {
  
  y_bar <- mean(y)
  
  rho_bar <- mean(y == 0)
  
  n <- length(y)
  
  n * (rho_bar * log(rho + (1 - rho) * exp(-mu)) + (1 - rho_bar) * (log(1 - rho) - mu) + y_bar * log(mu))
}

likelihood <- function(mu, rho, y) {
  
  y_bar <- mean(y)
  
  rho_bar <- mean(y == 0)
  
  n <- length(y)
  
  (((rho + (1 - rho) * exp(-mu))^rho_bar) * (((1 - rho) * exp(-mu))^(1 - rho_bar)) * (mu^y_bar))^n
}

neg_log_likelihood <- function(mu, rho, y) -log_likelihood(mu, rho, y)

get_rho <- function(mu, phi, mu_hat) (phi + (1 - phi) * exp(-mu_hat) - exp(-mu)) / (1 - exp(-mu))

get_mu_hat <- function(y) {
  
  y_bar <- mean(y)
  
  rho_bar <- mean(y == 0)
  
  gamma <- y_bar / (1 - rho_bar)
  
  pracma::lambertWp(-gamma * exp(-gamma)) + gamma
}

get_rho_hat <- function(y) {
  
  y_bar <- mean(y)
  
  mu_hat <- get_mu_hat(y)
  
  y_bar / mu_hat
}

get_mu_hat_se <- function(y) {
  
  n <- length(y)
  
  mu_hat <- get_mu_hat(y)
  
  rho_hat <- get_rho_hat(y)
  
  sqrt(mu_hat * (1 - exp(mu_hat)) / (n * rho_hat * (1 + mu_hat - exp(mu_hat))))
}

get_mu_grid <- function(y, step_size, num_std_errors) {
  
  mu_hat <- get_mu_hat(y)
  
  mu_hat_se <- get_mu_hat_se(y)
  
  MoE <- num_std_errors * mu_hat_se
  
  mu_grid <- (mu_hat + MoE*c(-1, 1)) |> 
    (\(x) c(max(0, x[1]), x[2]))() |> 
    plyr::round_any(step_size, floor) |> 
    (\(x) seq(x[1], x[2], step_size))()
  
  return(mu_grid)
}
