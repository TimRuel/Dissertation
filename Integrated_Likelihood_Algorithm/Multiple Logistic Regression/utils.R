################################################################################
#################################### GENERAL ###################################
################################################################################
library(tidyverse)
library(pipeR)

get_logistic_regression_model <- function(data) glm(Y ~ ., data = data, family = binomial(link = logit))

get_log_likelihood <- function(model = NULL, b = NULL, data = NULL) {
  
  if (is.null(model) && (is.null(b) || is.null(data))) {
    stop("If 'model' is not supplied, both 'b' and 'data' must be supplied.")
  }
  
  log_likelihood <- if (!is.null(model)) {
    
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
    
    sum(Y * Y_hat - log(1 + exp(Y_hat)))
  }
  
  log_likelihood
}

get_likelihood <- function(model = NULL, b = NULL, data = NULL) exp(get_log_likelihood(model, b, data))

get_neg_log_likelihood <- function(model = NULL, b = NULL, data = NULL) -get_log_likelihood(model, b, data)

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

get_omega_hat <- function(U, model, X_h) U * get_Y_hat(model, X_h) / (t(matrix(c(1, X_h))) %*% U)[1]

get_omega_hat_list <- function(U_list, model, X_h) purrr::map(U_list, \(U) get_omega_hat(U, model, X_h))

# get_Beta_hat1 <- function(psi, omega_hat, model, X_h) {
#   
#   X <- model.matrix(model)
#   
#   w <- plogis(X %*% omega_hat)
#   
#   f <- function(Beta) {
#     
#     linear_predictor <- X %*% matrix(Beta)
#     
#     -sum(w * linear_predictor - log(1 + exp(linear_predictor)))
#   }
#   
#   f.gr <- function(Beta) nloptr::nl.grad(Beta, f)
#   
#   fcon <- function(Beta) t(X_h) %*% matrix(Beta) - log(psi / (1 - psi))
#   fcon.jac <- function(Beta) nloptr::nl.jacobian(Beta, fcon)
#   
#   out <- nloptr::auglag(x0 = rep(0, ncol(X)),
#                         fn = f,
#                         gr = f.gr,
#                         heq = fcon,
#                         heqjac = fcon.jac,
#                         localsolver = "LBFGS")
#   
#   Beta_hat <- out$par
#   
#   Beta_hat
# }

get_Beta_hat <- function(psi, omega_hat, model, X_h, N) {
  
  model_data <- model |> 
    model.frame()
  
  X <- model_data |> 
    select(-Y)
  
  Y <- model_data |> 
    select(Y)
  
  frac <- model.matrix(model) %*% omega_hat |> 
    plogis()
  
  data <- X |> 
    sweep(2, X_h) |> 
    mutate(successes = round(N * frac),
           failures = N - round(N * frac),
           response = cbind(successes, failures))
    
  intercept <- log(psi / (1 - psi))
  
  fit <- glm(response ~ . - successes - failures - 1, data = data, family = binomial(link = logit), offset = rep(intercept, nrow(Y)))
  
  Beta_hat <- fit |> 
    coef() |> 
    unname() |> 
    matrix()
  
  c(intercept - t(Beta_hat) %*% X_h, Beta_hat) |> 
    matrix()
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






