################################################################################
#################################### GENERAL ###################################
################################################################################
library(tidyverse)

log_likelihood <- function(theta, data) {
  
  data |> 
    dpois(theta, log = TRUE) |> 
    sum(na.rm = TRUE)
}

likelihood <- function(theta, data) {
  
  data |> 
    dpois(theta) |> 
    prod()
}

neg_log_likelihood <- function(theta, data) -log_likelihood(theta, data)

dot_product <- function(x, y) sum(x * y, na.rm = TRUE)

distance <- function(a, b) sum((a - b)^2)

get_dist_list <- function(MC_params, dist_type = "importance") {
  
  switch(MC_params$method,
         vanilla_MC = MC_params$nominal,
         self_norm_IS = ,
         regression_IS = ,
         basic_IS = switch(dist_type, 
                           importance = MC_params$importance,
                           nominal = MC_params$nominal,
                           stop("Invalid distribution type")),
         stop("Invalid Monte Carlo method")
  )
}

dist_sampler <- function(MC_params, R, simplify = FALSE) {
  
  dist_list <- get_dist_list(MC_params)
  
  rng <- dist_list$rng
  
  rng_params <- dist_list$rng_params
  
  rng |> 
    do.call(args = rng_params) |> 
    replicate(R, expr = _, simplify = simplify)
}

get_importance_weights <- function(u_list, MC_params) {
  
  MC_params$method |> 
    
    switch(
      
      vanilla_MC = 1,
      
      self_norm_IS = ,
      
      regression_IS = ,
      
      basic_IS = {
        
        w <- u_list |>
          map_dbl(\(u) {
            
            p <- MC_params$nominal$density |>
              do.call(args = c(x = list(u), MC_params$nominal$density_params)) |>
              prod()
            
            q <- MC_params$importance$density |>
              do.call(args = c(x = list(u), MC_params$importance$density_params)) |>
              prod()
            
            p / q
          }
          )
        
        if (MC_params$method == "self_norm_IS") w <- w / mean(w)
        
        w
      }
    )
}

Q <- function(u, psi_hat, weights) u / dot_product(u, weights) * psi_hat

get_omega_hat_list <- function(u_list, psi_hat, weights) purrr::map(u_list, \(u) Q(u, psi_hat, weights))

get_theta_hat <- function(lambda, omega_hat, weights) omega_hat / (1 + lambda * weights)

get_lambda <- function(omega_hat, psi, weights) {
  
  psi_dist <- function(lambda) {
    
    lambda |> 
      get_theta_hat(omega_hat, weights) |> 
      dot_product(weights) |> 
      distance(psi)
  }
  
  psi_dist.gr <- function(lambda) nloptr::nl.grad(lambda, psi_dist)
  
  nloptr::auglag(x0 = 0,
                 fn = psi_dist,
                 gr = psi_dist.gr,
                 lower = -min(1 / weights),
                 localsolver = "LBFGS")$par
}

get_L_tilde <- function(omega_hat, psi, weights) {
  
  omega_hat |> 
    get_lambda(psi, weights) |> 
    get_theta_hat(omega_hat, weights) |> 
    likelihood(data)
}

get_l_bar <- function(psi, weights, omega_hat_list, importance_weights) {
  
  L_tilde <- omega_hat_list |> 
    map_dbl(\(omega_hat) get_L_tilde(omega_hat, psi, weights))
  
  l_bar <- log(mean(L_tilde * importance_weights))
  
  return(l_bar)
}

method = "vanilla_MC"

MC_params <- list(method = method, 
                  nominal = list(rng = rng, 
                                 rng_params = nominal_rng_params,
                                 density = density, 
                                 density_params = nominal_density_params))

psi_hat <- dot_product(data, weights)

R <- 250

u_list <- dist_sampler(MC_params, R)

omega_hat_list <- get_omega_hat_list(u_list, psi_hat, weights)

psi_grid <- get_psi_grid(data, weights, step_size, num_std_errors, split = FALSE)

importance_weights <- get_importance_weights(u_list, MC_params)

maximum <- get_l_bar(psi_hat, weights, omega_hat_list, importance_weights)

g <- function(psi) {
  
  l_bar <- get_l_bar(psi, weights, omega_hat_list, importance_weights)
  
  return(l_bar - maximum + crit)
}

tic()

uniroot(g, interval = c(0, psi_hat))$root

toc()

tic()

uniroot(g, interval = c(psi_hat, psi_hat + 3*psi_hat_se))$root

toc()



