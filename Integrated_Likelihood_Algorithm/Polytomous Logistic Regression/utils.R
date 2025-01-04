################################################################################
#################################### GENERAL ###################################
################################################################################
library(tidyverse)
library(pipeR)
library(accumulate)
# Rcpp::sourceCpp("accumulate_rcpp.cpp")

adj_softmax <- function(x) exp(x) / (1 + sum(exp(x)))

get_multinomial_logistic_model <- function(data) nnet::multinom(Y ~ ., data = data)

get_log_likelihood <- function(Beta, X_one_hot, Y_one_hot) {
  
  Beta <- matrix(Beta,
                 nrow = ncol(X_one_hot),
                 ncol = ncol(Y_one_hot),
                 byrow = FALSE)
  
  Y_hat <- X_one_hot %*% Beta
  
  sum(rowSums(Y_one_hot * Y_hat) - log(1 + rowSums(exp(Y_hat))))
}

get_likelihood <- function(Beta, X_one_hot, Y_one_hot) exp(get_log_likelihood(Beta, X_one_hot, Y_one_hot))

get_Beta_MLE <- function(model) {
  
  model |> 
    coef() |> 
    unname() |> 
    matrix(ncol = length(model$lev) - 1,
           byrow = TRUE) 
}

get_entropy <- function(p) -sum(p * log(p), na.rm = TRUE)

get_probability_vector <- function(k, target_entropy_range, epsilon = 1e-2, max_iter = 1000) {
  
  # Validate inputs
  if (length(target_entropy_range) != 2 || target_entropy_range[1] > target_entropy_range[2]) {
    stop("target_entropy_range must be a valid range: c(min, max)")
  }
  if (k < 2) {
    stop("k must be at least 2.")
  }
  
  adjust_probabilities <- function(prob, epsilon) {
    while (any(prob < epsilon)) {
      # Identify indices of probabilities below epsilon
      below_epsilon <- which(prob < epsilon)
      deficit <- sum(epsilon - prob[below_epsilon])  # Total deficit to distribute
      
      # Set all values below epsilon to epsilon
      prob[below_epsilon] <- epsilon
      
      # Find the index of the largest probability
      largest_index <- which.max(prob)
      
      # Deduct the deficit from the largest component
      prob[largest_index] <- prob[largest_index] - deficit
      
      # If the largest component becomes too small, redistribute proportionally
      if (prob[largest_index] < epsilon) {
        excess <- prob[largest_index] - epsilon
        prob[largest_index] <- epsilon
        prob[-largest_index] <- prob[-largest_index] + (excess * prob[-largest_index] / sum(prob[-largest_index]))
      }
    }
    return(prob)
  }
  
  # Initialize variables
  lower_bound <- target_entropy_range[1]
  upper_bound <- target_entropy_range[2]
  
  for (i in 1:max_iter) {
    # Generate a random probability vector
    p <- runif(k)
    p <- p / sum(p)  # Normalize to make it a valid probability vector
    
    # Adjust probabilities to ensure all components >= epsilon
    p <- adjust_probabilities(p, epsilon)
    
    # Calculate entropy
    entropy <- get_entropy(p)
    
    # Check if entropy is within the target range
    if (entropy >= lower_bound && entropy <= upper_bound) {
      return(p)
    }
  }
  
  # If no vector found within max_iter
  stop("Failed to generate a probability vector within the target entropy range.")
}

get_psi_hat <- function(model, X_h) {
  
  model |> 
    predict(newdata = X_h, type = "probs") |> 
    get_entropy()
}

get_psi_grid <- function(step_size, num_std_errors, model, X_h, split = FALSE) {
  
  m <- model |> 
    model.frame() |> 
    select(starts_with("X")) |> 
    table() |> 
    (\(table) table[X_h$X])()
  
  J <- length(model$lev)
  
  psi_hat <- get_psi_hat(model, X_h)
  
  probs <- model |> 
    predict(newdata = X_h, type = "probs")
  
  sigma <- probs*diag(J) - matrix(probs) %*% probs
  
  psi_hat_SE <- sqrt(sum(matrix(1 + log(probs)) %*% (1 + log(probs)) * sigma, na.rm = TRUE) / m)
  
  MoE <- num_std_errors * psi_hat_SE
  
  psi_grid <- (psi_hat + MoE*c(-1, 1)) |> 
    (\(x) c(max(0, x[1]), min(log(J), x[2])))() |> 
    plyr::round_any(step_size, floor) |> 
    (\(x) seq(x[1], x[2] + step_size, step_size))()
  
  if (split) {
    
    psi_grid_list <- psi_grid |> 
      split(factor(psi_grid > psi_hat)) |> 
      purrr::modify_in(1, \(x) c(x, tail(x, 1) + step_size / 2) |> rev()) |> 
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

get_log_importance_weights <- function(U_list, MC_params) {
  
  MC_params$method |> 
    
    switch(
      
      vanilla_MC = 0,
      
      self_norm_IS = ,
      
      regression_IS = ,
      
      basic_IS = {
        
        log_w <- U_list |>
          map_dbl(\(U) {
            
            log_p <- MC_params$nominal$density |>
              do.call(args = c(x = list(U), log = TRUE, MC_params$nominal$density_params)) |>
              sum()
            
            log_q <- MC_params$importance$density |>
              do.call(args = c(x = list(U), log = TRUE, MC_params$importance$density_params)) |>
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

get_Beta_hat <- function(init_guess, psi, omega_hat, X_one_hot, X_h_one_hot) {

  probs <- X_one_hot %*% omega_hat |> 
    adj_softmax()

  f <- function(Beta) {
    
    Beta <- matrix(Beta,
                   nrow = nrow(omega_hat),
                   ncol = ncol(omega_hat),
                   byrow = FALSE)
    
    Y_hat <- X_one_hot %*% Beta

    -sum(rowSums(probs * Y_hat) - log(1 + rowSums(exp(Y_hat))))
  }

  f.gr <- function(Beta) nloptr::nl.grad(Beta, f)

  fcon <- function(Beta) {
    
    Beta <- matrix(Beta,
                   nrow = nrow(omega_hat),
                   ncol = ncol(omega_hat),
                   byrow = FALSE)
    
    entropy <- X_h_one_hot %*% cbind(0, Beta) |> 
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

accumulate_Beta_hats <- function(psi_grid, omega_hat, X_one_hot, X_h_one_hot, init_guess) {

  psi_grid |>
    purrr::accumulate(
      \(acc, nxt) get_Beta_hat(acc, nxt, omega_hat, X_one_hot, X_h_one_hot),
      .init = init_guess) |>
    magrittr::extract(-1)
}

# accumulate_Beta_hats <- function(psi_grid, omega_hat, X, X_h, init_guess) {
# 
#   psi_grid |>
#     accumulate_rcpp(
#       \(acc, nxt) get_Beta_hat(acc, nxt, omega_hat, X, X_h),
#       init = init_guess)
# }

################################################################################
############################ INTEGRATED LIKELIHOOD ############################# 
################################################################################

get_log_L_tilde_mat <- function(psi_grid_list, 
                                omega_hat_list, 
                                X_one_hot, 
                                X_h_one_hot, 
                                Y_one_hot, 
                                init_guess, 
                                chunk_size) {
  
  p <- progressr::progressor(steps = 3 * length(omega_hat_list))
  
  foreach(
    
    omega_hat = omega_hat_list,
    .combine = "rbind",
    .multicombine = TRUE,
    .maxcombine = length(omega_hat_list),
    .options.future = list(seed = TRUE,
                           chunk.size = chunk_size,
                           packages = c("nloptr", "accumulate"))
    
  ) %dofuture% {
    
    p()

    Beta_hats_1 <- psi_grid_list[[1]] |>
      accumulate_Beta_hats(omega_hat, X_one_hot, X_h_one_hot, init_guess)

    p()

    Beta_hats_2 <- psi_grid_list[[2]] |>
      accumulate_Beta_hats(omega_hat, X_one_hot, X_h_one_hot, Beta_hats_1[[1]])

    p()

    c(rev(Beta_hats_1[-1]), Beta_hats_2) |>
      purrr::map_dbl(\(Beta_hat) get_log_likelihood(Beta_hat, X_one_hot, Y_one_hot))
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
                                          psi_grid_list, 
                                          R,
                                          MC_params,
                                          init_guess,
                                          chunk_size) {
  
  model <- get_multinomial_logistic_model(data)
  
  X_one_hot <- model.matrix(model)
  
  X_h_one_hot <- X_one_hot[X_h$X,]
  
  Y_one_hot <- data |> 
    pull(Y) |> 
    (\(Y) model.matrix(~ Y)[,-1])()
  
  Beta_MLE <- get_Beta_MLE(model)
  
  U_list <- get_U_list(MC_params, R)
  
  omega_hat_list <- get_omega_hat_list(U_list, Beta_MLE, X_h_one_hot)
  
  log_L_tilde_mat <- get_log_L_tilde_mat(psi_grid_list, omega_hat_list, X_one_hot, X_h_one_hot, Y_one_hot, init_guess, chunk_size)
  
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
              psi_grid_list = psi_grid_list))
}

################################################################################
############################## PROFILE LIKELIHOOD ############################## 
################################################################################

get_log_profile_likelihood <- function(data,
                                       X_h, 
                                       psi_grid_list) {
  
  model <- get_multinomial_logistic_model(data)
  
  X_one_hot <- model.matrix(model)
  
  X_h_one_hot <- X_one_hot[X_h$X,]
  
  Y_one_hot <- data |> 
    pull(Y) |> 
    (\(Y) model.matrix(~ Y)[,-1])()
  
  Beta_MLE <- get_Beta_MLE(model)
  
  Beta_hats_1 <- psi_grid_list[[1]] |> 
    accumulate_Beta_hats(Beta_MLE, X_one_hot, X_h_one_hot, c(Beta_MLE))
  
  Beta_hats_2 <- psi_grid_list[[2]] |> 
    accumulate_Beta_hats(Beta_MLE, X_one_hot, X_h_one_hot, Beta_hats_1[[1]])
  
  c(rev(Beta_hats_1[-1]), Beta_hats_2) |> 
    purrr::map_dbl(\(Beta_hat) get_log_likelihood(Beta_hat, X_one_hot, Y_one_hot))
}