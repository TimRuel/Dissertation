################################################################################
#################################### GENERAL ###################################
################################################################################
library(tidyverse)

choose_directory = function(caption = "Select population directory") {
  if (exists("choose.dir")) {
    choose.dir(caption = caption) 
  } else {
    rstudioapi::selectDirectory(caption = caption)
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

get_dist_lists <- function(MC_params, dist_type = "importance") {
  
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

dist_sampler <- function(dist_list) {
  
  rng <- dist_list$rng
  
  rng_params <- dist_list$rng_params
  
  rng |> 
    do.call(args = rng_params) 
}

get_u_list <- function(MC_params, R, simplify = FALSE) {
  
  MC_params |> 
    get_dist_lists() |> 
    purrr::map_dbl(dist_sampler) |> 
    replicate(R, expr = _, simplify = simplify)
}

get_log_importance_weights <- function(u_list, MC_params) {
  
  MC_params$method |> 
    
    switch(
      
      vanilla_MC = 0,
      
      self_norm_IS = ,
      
      regression_IS = ,
      
      basic_IS = {
        
        log_w <- u_list |>
          map_dbl(\(u) {
            
            log_p <- MC_params$nominal$density |>
              do.call(args = c(x = list(u), log = TRUE, MC_params$nominal$density_params)) |>
              sum()
            
            log_q <- MC_params$importance$density |>
              do.call(args = c(x = list(u), log = TRUE, MC_params$importance$density_params)) |>
              sum()
            
            log_p - log_q
          }
          )
        
        if (MC_params$method == "self_norm_IS") log_w <- log_w - matrixStats::logSumExp(log_w) + log(length(u_list))
        
        log_w
      }
    )
}

Q <- function(u, psi_hat, x_h) {
  
  C <- g(u[1], u[2], x_h)
  
  c(u[1] * psi_hat / C, u[2] * psi_hat / C, u[3])
}

get_omega_hat_list <- function(u_list, psi_hat, x_h) purrr::map(u_list, \(u) Q(u, psi_hat, x_h))

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

accumulate_theta_hats <- function(psi_grid_list, omega_hat, x, x_h, init_guess) {
  
  psi_grid_list |>
    purrr::map(\(psi_grid) {
      psi_grid |> 
        purrr::accumulate(
          \(acc, nxt) get_theta_hat(acc, nxt, omega_hat, x, x_h),
          .init = init_guess) |>
        magrittr::extract(-1)
      }
    ) |> 
    purrr::modify_in(1, rev) |> 
    unlist(recursive = FALSE)
}

map_theta_hats <- function(psi_grid_list, omega_hat, x, x_h, init_guess) {
  
  psi_grid_list |>
    purrr::map(\(psi_grid) {
      psi_grid |> 
        purrr::map_dbl(\(psi) get_theta_hat(init_guess, psi, omega_hat, x, x_h))
      }
    ) |> 
    purrr::modify_in(1, rev)
}

################################################################################
############################ INTEGRATED LIKELIHOOD ############################# 
################################################################################

get_log_L_tilde <- function(psi_grid_list, omega_hat, x, y, x_h, init_guess, theta_hat_method = "accumulate") {
  
  theta_hat_method |> 
    
    switch(
      
      accumulate = accumulate_theta_hats(psi_grid_list, omega_hat, x, x_h, init_guess),
      
      map = map_theta_hats(psi_grid_list, omega_hat, x, x_h, init_guess)
      
    ) |> 
    purrr::map_dbl(\(theta_hat) log_likelihood(theta_hat[1], theta_hat[2], theta_hat[3], x, y))
}

get_log_L_tilde_mat <- function(psi_grid_list, omega_hat_list, chunk_size, theta_hat_method, init_guess) {
  
  p <- progressr::progressor(along = omega_hat_list)
  
  foreach(
    
    omega_hat = omega_hat_list,
    .combine = "rbind",
    .multicombine = TRUE,
    .maxcombine = length(omega_hat_list),
    .options.future = list(seed = TRUE,
                           chunk.size = chunk_size)
    
  ) %dofuture% {
    
    p()
    
    get_log_L_tilde(psi_grid_list, omega_hat, x, y, x_h, init_guess, theta_hat_method)
  }
}

get_log_L_bar <- function(log_L_tilde_mat, log_w, MC_params) {
  
  MC_params$method |> 
    
    switch(
      
      vanilla_MC = ,
      
      self_norm_IS = ,
      
      basic_IS = {
        
        log_weighted_vals <- log_L_tilde_mat |> 
          sweep(1, log_w, '+')
      }
      
      # regression_IS = {
      #   
      #   w_bar <- mean(w)
      #   
      #   w_centered <- w - w_bar
      #   
      #   denom <- sum(w_centered^2)
      #   
      #   num <- l_tilde_mat |> 
      #     sweep(1, log_w + w_centered, '*') |> 
      #     colSums()
      #   
      #   beta_hat <- num / denom
      #   
      #   weighted_vals <- sweep(L_tilde_mat, 1, w, '*') - outer(w - 1, beta_hat, '*')
      # }
    )
  
  estimate <- matrixStats::colLogSumExps(log_weighted_vals) - log(nrow(log_L_tilde_mat))
  
  var_estimate <- matrixStats::colVars(log_weighted_vals)
  
  var_exp_estimate <- matrixStats::colVars(exp(log_weighted_vals))
  
  return(list(log_weighted_vals = log_weighted_vals,
              estimate = estimate, 
              var_estimate = var_estimate,
              var_exp_estimate = var_exp_estimate))
}

get_log_integrated_likelihood <- function(x,
                                          y,
                                          x_h, 
                                          psi_grid_list, 
                                          R,
                                          MC_params,
                                          theta_hat_method,
                                          init_guess,
                                          chunk_size) {
  
  alpha_hat <- get_alpha_hat(x, y)
  
  beta_hat <- get_beta_hat(x, y)
  
  psi_hat <- g(alpha_hat, beta_hat, x_h)
  
  u_list <- get_u_list(MC_params, R)
  
  omega_hat_list <- get_omega_hat_list(u_list, psi_hat, x_h)
  
  log_L_tilde_mat <- get_log_L_tilde_mat(psi_grid_list, omega_hat_list, chunk_size, theta_hat_method, init_guess)
  
  log_importance_weights <- get_log_importance_weights(u_list, MC_params)
  
  log_L_bar <- get_log_L_bar(log_L_tilde_mat, log_importance_weights, MC_params)
  
  return(list(log_L_bar = log_L_bar,
              log_L_tilde_mat = log_L_tilde_mat,
              log_importance_weights = log_importance_weights,
              omega_hat_list = omega_hat_list,
              u_list = u_list, 
              psi_hat = psi_hat, 
              MC_params = MC_params,
              x = x,
              y = y,
              x_h =x_h,
              psi_grid_list = psi_grid_list))
}

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
    purrr::map_dbl(\(theta) log_likelihood(theta[1], theta[2], theta[3], x, y))
}
