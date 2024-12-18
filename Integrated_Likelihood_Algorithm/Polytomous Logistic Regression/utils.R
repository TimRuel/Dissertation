################################################################################
#################################### GENERAL ###################################
################################################################################
library(tidyverse)
library(pipeR)

adj_softmax <- function(x) exp(x) / (1 + sum(exp(x)))

get_multinomial_logistic_model <- function(data) nnet::multinom(Y ~ ., data = data)

get_log_likelihood <- function(model = NULL, b = NULL, data = NULL) {
  
  if (is.null(model) && (is.null(b) || is.null(data))) {
    stop("If 'model' is not supplied, both 'b' and 'data' must be supplied.")
  }
  
  if (!is.null(model)) {
    
    logLik(model)[1]
  } else {
    
    X <- data |> 
      select(-Y) |> 
      as.matrix() |> 
      unname() |> 
      (\(z) cbind(1, z))()
    
    Y_hat <- X %*% b
    
    Y <- data |> 
      pull(Y)
    
    Y_one_hot <- model.matrix(~ Y)[,-1]
    
    sum(rowSums(Y_one_hot * Y_hat) - log(1 + rowSums(exp(Y_hat))))
  }
}

get_likelihood <- function(model = NULL, b = NULL, data = NULL) exp(get_log_likelihood(model, b, data))

get_neg_log_likelihood <- function(model = NULL, b = NULL, data = NULL) -get_log_likelihood(model, b, data)

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
    
    X_h |> 
      t() |> 
      data.frame() |> 
      predict(model, newdata = _, type = "probs") |> 
      get_entropy()
  } else {
    
    cbind(1, t(X_h)) %*% cbind(0, omega_hat) |> 
      LDATS::softmax() |> 
      entropy()
  }
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
    matrix() |> 
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

get_omega_hat <- function(U, model, X_h) {
  
  b <- get_Beta_MLE(model)
  
  psi_hat_logits <- X_h %*% cbind(0, b)
  
  U_logits <- X_h %*% U
  
  sweep(U, 2, psi_hat_logits[-1] / U_logits, "*")
}

get_omega_hat_list <- function(U_list, model, X_h) purrr::map(U_list, \(U) get_omega_hat(U, model, X_h))

get_Beta_hat <- function(psi, omega_hat, X, X_h, init_guess) {

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

  out <- nloptr::auglag(x0 = init_guess,
                        fn = f,
                        gr = f.gr,
                        heq = fcon,
                        heqjac = fcon.jac,
                        localsolver = "LBFGS")
  
  matrix(out$par,
         nrow = nrow(omega_hat),
         ncol = ncol(omega_hat),
         byrow = FALSE)
}

# accumulate_beta_hats <- function(psi_grid_list, omega_hat, x, x_h, init_guess) {
#   
#   psi_grid_list |>
#     purrr::map(\(psi_grid) {
#       psi_grid |> 
#         purrr::accumulate(
#           \(acc, nxt) get_beta_hat(acc, nxt, omega_hat, x, x_h),
#           .init = init_guess) |>
#         magrittr::extract(-1)
#     }
#     ) |> 
#     purrr::modify_in(1, rev) |> 
#     unlist(recursive = FALSE)
# }

# map_beta_hats <- function(psi_grid_list, omega_hat, x, x_h, init_guess) {
#   
#   psi_grid_list |>
#     purrr::map(\(psi_grid) {
#       psi_grid |> 
#         purrr::map_dbl(\(psi) get_beta_hat(init_guess, psi, omega_hat, x, x_h))
#     }
#     ) |> 
#     purrr::modify_in(1, rev)
# }

################################################################################
############################ INTEGRATED LIKELIHOOD ############################# 
################################################################################

get_log_L_tilde <- function(psi_grid, omega_hat, model, X_h, N) {
  
  data <- model |> 
    model.frame()
  
  psi_grid |> 
    purrr::map_dbl(\(psi) {
      
      psi |> 
        get_Beta_hat(omega_hat, model, X_h, N) |> 
        get_log_likelihood(b = _, data = data)
    }
    )
}

get_log_L_tilde_mat <- function(psi_grid, omega_hat_list, model, X_h, N, chunk_size) {
  
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
    
    get_log_L_tilde(psi_grid, omega_hat, model, X_h, N)
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
                                          N,
                                          MC_params,
                                          chunk_size) {
  
  model <- get_logistic_regression_model(data)
  
  U_list <- get_U_list(MC_params, R)
  
  omega_hat_list <- get_omega_hat_list(U_list, model, X_h)
  
  log_L_tilde_mat <- get_log_L_tilde_mat(psi_grid, omega_hat_list, model, X_h, N, chunk_size)
  
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

get_profile_log_likelihood <- function(data, X_h, step_size) {
  
  psi_grid <- get_psi_grid(step_size)
  
  model <- get_logistic_regression_model(data)
  
  Beta_MLE <- get_Beta_MLE(model)
  
  psi_grid |> 
    get_log_L_tilde(Beta_MLE, model, X_h, N)
}






