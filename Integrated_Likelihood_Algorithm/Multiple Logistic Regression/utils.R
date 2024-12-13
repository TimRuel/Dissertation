################################################################################
#################################### GENERAL ###################################
################################################################################
library(tidyverse)
library(pipeR)

get_logistic_regression_model <- function(data) glm(Y ~ ., data = data, family = binomial(link = logit))

log_likelihood <- function(model) logLik(model)[1]

likelihood <- function(model) exp(log_likelihood(model))

neg_log_likelihood <- function(model) -log_likelihood(model)

get_Beta_MLE <- function(model) {
  
  model |> 
    coef() |> 
    unname() |> 
    matrix()
}

get_Y_hat <- function(model, X_h) {
  
  X_h |> 
    t() |> 
    data.frame() |> 
    predict(model, newdata = _, type = "link") |> 
    unname()
}

get_SE_Y_hat <- function(model, X_h) {
  
  X_h |> 
    t() |> 
    data.frame() |> 
    predict(model, newdata = _, type = "link", se.fit = TRUE) %>>%
    (se.fit)
}

get_CI_Y_hat <- function(model, X_h, alpha) {
  
  Y_hat <- get_Y_hat(model, X_h)
  
  SE_Y_hat <- get_SE_Y_hat(model, X_h)
  
  critical_value <- qnorm(1 - alpha / 2)
  
  MoE <- critical_value * SE_Y_hat
  
  Y_hat + c(-1, 1) * MoE
}

get_psi_hat <- function(model, X_h) {
  
  X_h |> 
    t() |> 
    data.frame() |> 
    predict(model, newdata = _, type = "response") |> 
    unname()
}

get_CI_psi_hat <- function(model, X_h, alpha){
  
  model |> 
    get_CI_Y_hat(X_h, alpha) |>  
    plogis()
}

get_psi_grid <- function(step_size, psi_hat = NULL) {
  
  psi_grid <- seq(0 + step_size, 1 - step_size, step_size)
  
  if (!is.null(psi_hat)) {
    
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

get_omega_hat <- function(U, model, X_h) U * get_Y_hat(model, X_h) / (t(matrix(c(1, X_h))) %*% U)[1]

get_omega_hat_list <- function(U_list, model, X_h) purrr::map(U_list, \(U) get_omega_hat(U, model, X_h))

get_Beta_hat <- function(psi, omega_hat, model, X_h, N) {
  
  w <- model.matrix(model) %*% omega_hat |> 
    plogis()
  
  model_data <- model |> 
    model.frame()
  
  X <- model_data |> 
    select(-Y)
  
  Y <- model_data |> 
    select(Y)
  
  N <- 10^3
  
  data <- X |> 
    sweep(2, X_h) |> 
    tidyr::uncount(weights = N) |> 
    mutate(Y = as.integer(N * w) |> 
             map(\(num_successes) rep(c(1,0), times = c(num_successes, N - num_successes))) |> 
             unlist())
  
  intercept <- log(psi / (1 - psi))
  
  fit <- glm(Y ~ . - 1, data = data, family = binomial(link = logit), offset = rep(intercept, length(Y)))
  
  beta_hat <- fit |> 
    coef() |> 
    unname() |> 
    matrix()
  
  c(intercept - t(beta_hat) %*% X_h, beta_hat) |> 
    matrix()
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






