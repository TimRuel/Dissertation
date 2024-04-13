library(nloptr)
library(LaplacesDemon)
library(magrittr)
library(purrr)

likelihood <- function(theta, data) prod(theta^data)

PoI_fn <- function(x) -sum(x*log(x), na.rm = TRUE)

get_omega_hat <- function(u, psi_MLE) {
  
  omega_hat <- auglag(x0 = u,
                      fn = function(omega) -sum(u*log(omega)),
                      heq = function(omega) c(sum(omega) - 1, PoI_fn(omega) - psi_MLE),
                      lower = rep(0, length(u)))$par
  
  return(omega_hat)
}

get_theta_hat <- function(psi, omega_hat) {
  
  theta_hat <- auglag(x0 = omega_hat,
                      fn = function(theta) -sum(omega_hat*log(theta)),
                      heq = function(theta) c(sum(theta) - 1, PoI_fn(theta) - psi),
                      lower = rep(0, length(omega_hat)))$par
  
  return(theta_hat)
}

get_multinomial_entropy_IL_values <- function(data, psi_grid) {
  
  m <- length(data)
  
  theta_MLE <- data / sum(data)
  
  psi_MLE <- PoI_fn(theta_MLE)
   
  u <- rdirichlet(1, rep(1, m))
  
  omega_hat <- get_omega_hat(u, psi_MLE)
  
  L <- psi_grid |> 
    map(get_theta_hat, omega_hat) |> 
    sapply(likelihood, data)
  
  return(L)
}