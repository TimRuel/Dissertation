################################################################################
#################################### GENERAL ###################################
################################################################################
library(tidyverse)

choose_directory = function(caption = "Select population directory") {
  if (exists("choose.dir")) {
    choose.dir(caption = caption) 
  } else {
    tcltk::tk_choose.dir(caption = caption)
  }
}

log_likelihood <- function(alpha, beta, sigma_squared, x, y) {
  
  n <- length(x)
  
  SS <- sum((y - alpha - beta * x)^2)
  
  -n / 2 * log(sigma_squared) - SS / (2 * sigma_squared)
}

likelihood <- function(alpha, beta, sigma_squared, x, y) exp(log_likelihood(alpha, beta, sigma_squared, x, y))

neg_log_likelihood <- function(alpha, beta, sigma_squared, x, y) -log_likelihood(alpha, beta, sigma_squared, x, y)

g <- function(alpha, beta, x_h) alpha + beta * x_h

get_alpha_hat <- function(x, y) {
  
  beta_hat <- get_beta_hat(x, y)
  
  x_bar <- mean(x)
  
  y_bar <- mean(y)
  
  y_bar - beta_hat * x_bar
}

get_beta_hat <- function(x, y) {
  
  x_bar <- mean(x)
  
  y_bar <- mean(y)
  
  sum((x - x_bar) * (y - y_bar)) / sum((x - x_bar)^2)
}

get_sigma_squared_hat <- function(x, y) {
  
  alpha_hat <- get_alpha_hat(x, y)
  
  beta_hat <- get_beta_hat(x, y)
  
  mean((y - alpha_hat - beta_hat *x)^2)
}

get_psi_hat <- function(x, y, x_h) {
  
  alpha_hat <- get_alpha_hat(x, y)
  
  beta_hat <- get_beta_hat(x, y)
  
  g(alpha_hat, beta_hat, x_h)
}

get_psi_hat_se <- function(x, y, x_h) {
  
  sigma_squared_hat <- get_sigma_squared_hat(x, y)
  
  x_bar <- mean(x)
  
  sqrt(sigma_squared_hat * mean((x - x_h)^2) / sum((x - x_bar)^2))
}

get_psi_grid <- function(x, y, x_h, step_size, num_std_errors, split = FALSE) {
  
  psi_hat <- get_psi_hat(x, y, x_h)
  
  psi_hat_se <- get_psi_hat_se(x, y, x_h)
  
  MoE <- num_std_errors * psi_hat_se
  
  psi_grid <- (psi_hat + MoE * c(-1, 1)) |> 
    plyr::round_any(step_size, floor) |> 
    (\(x) seq(x[1] - step_size, x[2] + step_size, step_size))()
  
  if (split) {
    
    psi_grid_list <- psi_grid |> 
      split(factor(psi_grid > psi_hat)) |> 
      purrr::modify_in(1, rev) |> 
      unname()
    
    return(psi_grid_list)
  }
  
  return(psi_grid)
}

get_theta_hat <- function(init_guess, psi, omega_hat, x, x_h) {
  
  n <- length(x)

  f <- function(theta) {
    
    n * log(theta[3]) + sum(omega_hat[3] + ((omega_hat[1] - theta[1]) + (omega_hat[2] - theta[2]) * x)^2) / theta[3]
  }
  
  f.gr <- function(theta) nloptr::nl.grad(theta, f)
  fcon <- function(theta) g(theta[1], theta[2], x_h) - psi
  fcon.jac <- function(theta) nloptr::nl.jacobian(theta, fcon)

  theta_hat <- nloptr::auglag(x0 = init_guess,
                              fn = f,
                              gr = f.gr,
                              heq = fcon,
                              heqjac = fcon.jac,
                              lower = c(-Inf, -Inf, 1e-10),
                              localsolver = "LBFGS")$par
  
  if (identical(theta_hat, init_guess)) {
    
    init_guess <- init_guess + c(rnorm(2, sd = 1), init_guess[3])
    
    get_theta_hat(init_guess, psi, omega_hat, x, x_h)
  }

  return(theta_hat)
}

################################################################################
############################ INTEGRATED LIKELIHOOD ############################# 
################################################################################



################################################################################
############################## PROFILE LIKELIHOOD ############################## 
################################################################################

get_profile_log_likelihood <- function(x, y, x_h, step_size, num_std_errors) {
  
  psi_grid <- get_psi_grid(x, y, x_h, step_size, num_std_errors, split = FALSE)
  
  alpha_hat <- get_alpha_hat(x, y)
  
  beta_hat <- get_beta_hat(x, y)
  
  sigma_squared_hat <- get_sigma_squared_hat(x, y)
  
  theta_MLE <- c(alpha_hat, beta_hat, sigma_squared_hat)
  
  psi_grid |> purrr::accumulate(
    \(acc, nxt) {
      get_theta_hat(acc, nxt, theta_MLE, x, x_h)
    },
    .init = c(0, 0, 1)
  ) |> 
    magrittr::extract(-1) |> 
    magrittr::extract(-1) |> 
    purrr::map_dbl(\(theta) log_likelihood(theta[1], theta[2], theta[3], x, y))
}
