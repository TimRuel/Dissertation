################################################################################
#################################### GENERAL ###################################
################################################################################
library(tidyverse)
library(pipeR)

log_likelihood <- function(Beta, sigma_squared, X, Y) {
  
  n <- nrow(Y)
  
  SS <- sum((Y - X %*% Beta)^2)
  
  -n / 2 * log(sigma_squared) - SS / (2 * sigma_squared)
}

likelihood <- function(Beta, sigma_squared, X, Y) exp(log_likelihood(Beta, sigma_squared, X, Y))

neg_log_likelihood <- function(Beta, sigma_squared, X, Y) -log_likelihood(Beta, sigma_squared, X, Y)

get_Beta_MLE <- function(X, Y) solve(t(X) %*% X) %*% t(X) %*% Y

get_Y_hat <- function(X_h, b) {
  
  t(X_h) %*% b |> 
    c()
}

get_MSE <- function(X, Y, b) {
  
  e <- Y - X %*% b
  
  SSE <- t(e) %*% e |> 
    c()
  
  n <- nrow(X)
  
  p <- ncol(X)
  
  SSE / (n - p)
}

get_SE_Y_hat <- function(X, Y, X_h, b) {
  
  MSE <- get_MSE(X, Y, b)
  
  sqrt(MSE * t(X_h) %*% solve(t(X) %*% X) %*% X_h)[1]
}

get_CI_Y_hat <- function(X, Y, X_h, b, alpha) {
  
  Y_hat <- get_Y_hat(X_h, b)
  
  SE_Y_hat <- get_SE_Y_hat(X, Y, X_h, b)
  
  n <- nrow(X)
  
  p <- ncol(X)
  
  critical_value <- qt(p = 1 - alpha / 2, df = n - p)
  
  MoE <- critical_value * SE_Y_hat
  
  Y_hat + c(-1, 1) * MoE
}

get_psi_hat <- function(x, y, x_h) {
  
  beta_MLE <- get_beta_MLE(x, y)
  
  get_logistic_mean_response(beta_MLE, x_h)
}

get_psi_grid <- function(step_size, x = NULL, y = NULL, x_h = NULL, split = FALSE) {
  
  psi_grid <- seq(0 + step_size, 1 - step_size, step_size)
  
  if (split) {
    
    psi_hat <- get_psi_hat(x, y, x_h)
    
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

get_omega_hat <- function(u, beta_MLE, x_h) u * get_linear_predictor(beta_MLE, x_h) / get_linear_predictor(u, x_h)

get_omega_hat_list <- function(u_list, beta_MLE, x_h) purrr::map(u_list, \(u) get_omega_hat(u, beta_MLE, x_h))

get_beta_hat <- function(init_guess, psi, omega_hat, x, x_h) {
  
  logistic_mean_response <- get_logistic_mean_response(omega_hat, x)
  
  f <- function(beta) {
    
    linear_predictor <- get_linear_predictor(beta, x)
    
    -sum(logistic_mean_response * linear_predictor - log(1 + exp(linear_predictor)))
  }
  
  f.gr <- function(beta) nloptr::nl.grad(beta, f)
  fcon <- function(beta) get_logistic_mean_response(beta, x_h) - psi
  fcon.jac <- function(beta) nloptr::nl.jacobian(beta, fcon)
  
  out <- nloptr::auglag(x0 = init_guess,
                        fn = f,
                        gr = f.gr,
                        heq = fcon,
                        heqjac = fcon.jac,
                        localsolver = "LBFGS")
  
  beta_hat <- out$par
  
  return(beta_hat)
}

accumulate_beta_hats <- function(psi_grid_list, omega_hat, x, x_h, init_guess) {
  
  psi_grid_list |>
    purrr::map(\(psi_grid) {
      psi_grid |> 
        purrr::accumulate(
          \(acc, nxt) get_beta_hat(acc, nxt, omega_hat, x, x_h),
          .init = init_guess) |>
        magrittr::extract(-1)
    }
    ) |> 
    purrr::modify_in(1, rev) |> 
    unlist(recursive = FALSE)
}

map_beta_hats <- function(psi_grid_list, omega_hat, x, x_h, init_guess) {
  
  psi_grid_list |>
    purrr::map(\(psi_grid) {
      psi_grid |> 
        purrr::map_dbl(\(psi) get_beta_hat(init_guess, psi, omega_hat, x, x_h))
    }
    ) |> 
    purrr::modify_in(1, rev)
}

################################################################################
############################ INTEGRATED LIKELIHOOD ############################# 
################################################################################

get_log_L_tilde <- function(psi_grid_list, omega_hat, x, y, x_h, init_guess, beta_hat_method = "accumulate") {
  
  beta_hat_method |> 
    
    switch(
      
      accumulate = accumulate_beta_hats(psi_grid_list, omega_hat, x, x_h, init_guess),
      
      map = map_beta_hats(psi_grid_list, omega_hat, x, x_h, init_guess)
      
    ) |> 
    purrr::map_dbl(\(beta_hat) log_likelihood(beta_hat, x, y))
}

get_log_L_tilde_mat <- function(psi_grid_list, omega_hat_list, chunk_size, beta_hat_method, init_guess) {
  
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
    
    get_log_L_tilde(psi_grid_list, omega_hat, x, y, x_h, init_guess, beta_hat_method)
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
                                          beta_hat_method,
                                          init_guess,
                                          chunk_size) {
  
  beta_MLE <- get_beta_MLE(x, y)
  
  u_list <- get_u_list(MC_params, R)
  
  omega_hat_list <- get_omega_hat_list(u_list, beta_MLE, x_h)
  
  log_L_tilde_mat <- get_log_L_tilde_mat(psi_grid_list, omega_hat_list, chunk_size, beta_hat_method, init_guess)
  
  log_importance_weights <- get_log_importance_weights(u_list, MC_params)
  
  log_L_bar <- get_log_L_bar(log_L_tilde_mat, log_importance_weights, MC_params)
  
  return(list(log_L_bar = log_L_bar,
              log_L_tilde_mat = log_L_tilde_mat,
              log_importance_weights = log_importance_weights,
              omega_hat_list = omega_hat_list,
              u_list = u_list, 
              beta_MLE = beta_MLE, 
              MC_params = MC_params,
              x = x,
              y = y,
              x_h = x_h,
              psi_grid_list = psi_grid_list))
}

################################################################################
############################## PROFILE LIKELIHOOD ############################## 
################################################################################

get_profile_log_likelihood <- function(x, y, x_h, step_size, init_guess) {
  
  psi_grid <- get_psi_grid(step_size)
  
  beta_MLE <- get_beta_MLE(x, y)
  
  psi_grid |> 
    purrr::accumulate(
      \(acc, nxt) get_beta_hat(acc, nxt, beta_MLE, x, x_h),
      .init = init_guess
    ) |> 
    magrittr::extract(-1) |> 
    purrr::map_dbl(\(beta) log_likelihood(beta, x, y))
}






