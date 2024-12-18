################################################################################
#################################### GENERAL ###################################
################################################################################
library(tidyverse)
library(pipeR)

adj_softmax <- function(x) exp(x) / (1 + sum(exp(x)))

get_multinomial_logistic_model <- function(data) nnet::multinom(Y ~ ., data = data)

get_log_likelihood <- function(b, X, Y_one_hot) {
  
  b <- matrix(b,
              nrow = ncol(X),
              ncol = ncol(Y_one_hot),
              byrow = FALSE)
  
  Y_hat <- X %*% b
  
  sum(rowSums(Y_one_hot * Y_hat) - log(1 + rowSums(exp(Y_hat))))
}

get_likelihood <- function(b, X, Y_one_hot) exp(get_log_likelihood(b, X, Y_one_hot))

get_Beta_MLE <- function(model) {
  
  model |> 
    coef() |> 
    unname() |> 
    matrix(ncol = length(model$lev) - 1,
           byrow = TRUE) 
}

get_entropy <- function(p) -sum(p * log(p), na.rm = TRUE)

get_psi_hat <- function(model = NULL, b = NULL, X_h) {
  
  if (!xor(is.null(model), is.null(b))) {
    stop("You must provide either 'model' or 'b', but not both.")
  }
  
  if (!is.null(model)) {
    
    X_h[,-1] |> 
      t() |> 
      data.frame() |> 
      predict(model, newdata = _, type = "probs") |> 
      get_entropy()
    
  } else {
    
    X_h %*% cbind(0, omega_hat) |> 
      LDATS::softmax() |> 
      entropy()
  }
}

get_psi_grid <- function(step_size, model, X_h = NULL) {
  
  J <- length(model$lev)
  
  psi_grid <- seq(0, log(J), step_size)
  
  if (!is.null(X_h)) {
    
    psi_hat <- get_psi_hat(model = model, X_h = X_h)
    
    psi_grid_list <- psi_grid |> 
      split(factor(psi_grid > psi_hat)) |> 
      purrr::modify_in(1, rev) |> 
      unname()
    
    return(psi_grid_list)
  }
  
  return(psi_grid)
}

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

dist_sampler <- function(dist_list) {
  
  rng <- dist_list$rng
  
  rng_params <- dist_list$rng_params
  
  rng |> 
    do.call(args = rng_params) 
}

get_U_list <- function(MC_params, R, simplify = FALSE) {
  
  MC_params |> 
    get_dist_list() |> 
    dist_sampler() |> 
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

get_omega_hat <- function(U, b, X_h) {
  
  psi_hat_logits <- X_h %*% cbind(0, b)
  
  U_logits <- X_h %*% U
  
  sweep(U, 2, psi_hat_logits[-1] / U_logits, "*")
}

get_omega_hat_list <- function(U_list, b, X_h) purrr::map(U_list, \(U) get_omega_hat(U, b, X_h))

get_Beta_hat <- function(init_guess, psi, omega_hat, X, X_h) {

  probs <- X %*% omega_hat |> 
    adj_softmax()

  f <- function(Beta) {
    
    Beta <- matrix(Beta,
                   nrow = nrow(omega_hat),
                   ncol = ncol(omega_hat),
                   byrow = FALSE)
    
    Y_hat <- X %*% Beta

    -sum(rowSums(probs * Y_hat) - log(1 + rowSums(exp(Y_hat))))
  }

  f.gr <- function(Beta) nloptr::nl.grad(Beta, f)

  fcon <- function(Beta) {
    
    Beta <- matrix(Beta,
                   nrow = nrow(omega_hat),
                   ncol = ncol(omega_hat),
                   byrow = FALSE)
    
    entropy <- (X_h %*% cbind(0, Beta)) |> 
      LDATS::softmax() |> 
      get_entropy()
      
    return(entropy - psi)
  }
  
  fcon.jac <- function(Beta) nloptr::nl.jacobian(Beta, fcon)

  nloptr::auglag(x0 = init_guess,
                 fn = f,
                 gr = f.gr,
                 heq = fcon,
                 heqjac = fcon.jac,
                 localsolver = "LBFGS")$par
}

accumulate_Beta_hats <- function(psi_grid, omega_hat, X, X_h, init_guess) {

  psi_grid |>
    purrr::accumulate(
      \(acc, nxt) get_Beta_hat(acc, nxt, omega_hat, X, X_h),
      .init = init_guess) |>
    magrittr::extract(-1)
}

################################################################################
############################ INTEGRATED LIKELIHOOD ############################# 
################################################################################

get_log_L_tilde <- function(psi_grid, omega_hat, X, Y_one_hot, X_h, init_guess) {
  
  psi_grid |> 
    accumulate_Beta_hats(omega_hat, X, X_h, init_guess) |> 
    purrr::map_dbl(\(Beta_hat) get_log_likelihood(b = Beta_hat, X = X, Y_one_hot = Y_one_hot))
}

get_log_L_tilde_mat <- function(psi_grid, omega_hat_list, X, Y_one_hot, X_h, init_guess, chunk_size) {
  
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
    
    get_log_L_tilde(psi_grid, omega_hat, X, Y_one_hot, X_h, init_guess)
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

get_log_integrated_likelihood <- function(data,
                                          X_h, 
                                          psi_grid, 
                                          R,
                                          MC_params,
                                          init_guess,
                                          chunk_size) {
  
  X <- data |> 
    select(-Y) |> 
    as.matrix() |> 
    unname() |> 
    (\(mat) cbind(1, mat))()
  
  Y <- data |> 
    pull(Y)
  
  Y_one_hot <- model.matrix( ~ Y)[,-1]
  
  model <- get_multinomial_logistic_model(data)
  
  Beta_MLE <- get_Beta_MLE(model)
  
  U_list <- get_U_list(MC_params, R)
  
  omega_hat_list <- get_omega_hat_list(U_list, Beta_MLE, X_h)
  
  log_L_tilde_mat <- get_log_L_tilde_mat(psi_grid, omega_hat_list, X, Y_one_hot, X_h, init_guess, chunk_size)
  
  log_importance_weights <- get_log_importance_weights(U_list, MC_params)
  
  log_L_bar <- get_log_L_bar(log_L_tilde_mat, log_importance_weights, MC_params)
  
  return(list(log_L_bar = log_L_bar,
              log_L_tilde_mat = log_L_tilde_mat,
              log_importance_weights = log_importance_weights,
              omega_hat_list = omega_hat_list,
              U_list = U_list, 
              MC_params = MC_params,
              model = model, 
              X_h = X_h,
              psi_grid = psi_grid))
}

################################################################################
############################## PROFILE LIKELIHOOD ############################## 
################################################################################

get_log_profile_likelihood <- function(data,
                                       X_h, 
                                       psi_grid, 
                                       init_guess) {
  
  X <- data |> 
    select(-Y) |> 
    as.matrix() |> 
    unname() |> 
    (\(mat) cbind(1, mat))()
  
  Y <- data |> 
    pull(Y)
  
  Y_one_hot <- model.matrix( ~ Y)[,-1]
  
  model <- get_multinomial_logistic_model(data)
  
  Beta_MLE <- get_Beta_MLE(model)
  
  psi_grid |> 
    accumulate_Beta_hats(Beta_MLE, X, X_h, init_guess) |> 
    purrr::map_dbl(\(Beta_hat) get_log_likelihood(b = Beta_hat, X = X, Y_one_hot = Y_one_hot))
}






