################################################################################
############################### DATA GENERATION ################################
################################################################################
library(tidyverse)
library(pipeR)
library(iterate)
# Rcpp::sourceCpp("../../iterate_while.cpp")

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

get_theta_0 <- function(J, p, epsilon, max_iter) {
  
  seq(0, log(J), length.out = p + 1) |> 
    (\(x) mapply(c, x[-length(x)], x[-1], SIMPLIFY = FALSE))() |> 
    map(\(x) get_probability_vector(k = J, target_entropy_range = x, epsilon = epsilon, max_iter = max_iter))
}

get_Y <- function(theta_0, m) {
  
  J <- length(theta_0[[1]])
  
  theta_0 |> 
    purrr::map(\(prob) sample(1:J, size = m, prob = prob, replace = TRUE)) |> 
    unlist() |> 
    factor(levels = 1:J)
}

get_X <- function(p, m, contrast) {
  
  X <- 1:p |> 
    rep(each = m) |> 
    factor()
  
  contrasts(X) <- contrast
  
  return(X)
}

get_omega_hat_data <- function(omega_hat, X_one_hot, m) {
  
  omega_hat_probs <- X_one_hot %*% cbind(0, omega_hat) |>
    apply(1, LDATS::softmax)
  
  omega_hat_Y <- omega_hat_probs |> 
    unique() |> 
    unname() |> 
    apply(1, \(x) x, simplify = FALSE) |> 
    get_Y(1e6)
  
  omega_hat_data <- data.frame(X = X,
                               Y = omega_hat_Y)
}

################################################################################
#################################### GENERAL ###################################
################################################################################

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

get_psi_hat <- function(model, X_h) {
  
  model |> 
    predict(newdata = X_h, type = "probs") |> 
    get_entropy()
}

get_psi_grid <- function(step_size, num_std_errors, model, X_h, split_at = NULL) {
  
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
    (\(x) seq(x[1] + step_size, x[2] + step_size, step_size))()
  
  if (!is.null(split_at)) {
    
    psi_grid_list <- psi_grid |> 
      split(factor(psi_grid > split_at)) |> 
      purrr::modify_in(1, \(x) x |> rev()) |> 
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

get_U <- function(MC_params) {
  
  MC_params |> 
    get_dist_list() |> 
    dist_sampler()
}

get_U_list <- function(MC_params, R, simplify = FALSE) {
  
  MC_params |> 
    get_U() |> 
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

# get_omega_hat <- function(U, b, X_h) {
#   
#   psi_hat_logits <- X_h %*% cbind(0, b)
#   
#   U_logits <- X_h %*% U
#   
#   sweep(U, 2, psi_hat_logits[-1] / U_logits, "*")
# }
# 
# get_omega_hat_list <- function(U_list, b, X_h) purrr::map(U_list, \(U) get_omega_hat(U, b, X_h))

get_omega_hat <- function(MC_params, X_one_hot, X_h_one_hot, Y_one_hot, psi_hat, threshold) {
  
  U <- get_U(MC_params)
  
  f <- function(Beta) {
    
    Beta <- matrix(Beta,
                   nrow = nrow(U),
                   ncol = ncol(U),
                   byrow = FALSE)
    
    entropy <- X_h_one_hot %*% cbind(0, Beta) |>
      LDATS::softmax() |>
      get_entropy()
    
    return(abs(entropy - psi_hat))
  }
  
  f.gr <- function(Beta) nloptr::nl.grad(Beta, f)
  
  fcon <- function(Beta) -get_log_likelihood(Beta, X_one_hot, Y_one_hot) - threshold
  
  fcon.jac <- function(Beta) nloptr::nl.jacobian(Beta, fcon)
  
  omega_hat <- nloptr::auglag(x0 = U,
                              fn = f,
                              gr = f.gr,
                              # heq = fcon1,
                              # heqjac = fcon1.jac,
                              hin = fcon,
                              hinjac = fcon.jac,
                              localsolver = "LBFGS",
                              deprecatedBehavior = FALSE)$par
  
  if (f(omega_hat) <= 0.01 && fcon(omega_hat) <= 0) {
    
    omega_hat <- omega_hat |> 
      matrix(nrow = nrow(U),
             ncol = ncol(U),
             byrow = FALSE)
    
    return(omega_hat)
  }
  
  else {
    
    get_omega_hat(MC_params, X_one_hot, X_h_one_hot, Y_one_hot, psi_hat, threshold)
  }
}

get_omega_hat_list <- function(MC_params, R, X_one_hot, X_h_one_hot, Y_one_hot, psi_hat, threshold, simplify = FALSE) {
  
  MC_params |> 
    get_omega_hat(X_one_hot, X_h_one_hot, Y_one_hot, psi_hat, threshold) |> 
    replicate(R, expr = _, simplify = simplify)
}

get_Beta_hat <- function(init_guess, 
                         psi, 
                         omega_hat, 
                         X_one_hot, 
                         X_h_one_hot) {


  probs <- X_one_hot %*% omega_hat |>
    apply(1, adj_softmax) |>
    t()

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
                 localsolver = "LBFGS",
                 deprecatedBehavior = FALSE)$par
}

make_branch_fn <- function(omega_hat, X_one_hot, X_h_one_hot, Y_one_hot) {
  
  function(psi) {
    
    omega_hat |> 
      c() |> 
      get_Beta_hat(psi, omega_hat, X_one_hot, X_h_one_hot) |> 
      get_log_likelihood(X_one_hot, Y_one_hot)
  }
}

################################################################################
############################ INTEGRATED LIKELIHOOD ############################# 
################################################################################

get_Beta_hat_matrices <- function(omega_hat_list,
                                  X_one_hot,
                                  X_h,
                                  X_h_one_hot,
                                  Y_one_hot,
                                  step_size, 
                                  model,
                                  crit, 
                                  chunk_size) {

  p <- progressr::progressor(steps = 2 * length(omega_hat_list))

  foreach(

    omega_hat = omega_hat_list,
    .combine = "list",
    .multicombine = TRUE,
    .maxcombine = length(omega_hat_list),
    .errorhandling = "remove",
    .options.future = list(seed = TRUE,
                           chunk.size = chunk_size,
                           packages = c("nloptr", "iterate"))

  ) %dofuture% {

    tryCatch({

      p()
      
      branch_fn <- make_branch_fn(omega_hat, X_one_hot, X_h_one_hot, Y_one_hot)
      
      ga_result <- GA::ga(
        type = "real-valued",
        fitness = branch_fn,
        maxiter = 5,
        run = 2,
        lower = 0,
        upper = log(J))
      
      p()
      
      psi_max <- ga_result@solution[1,] |> 
        unname()
      
      branch_max <- ga_result@fitnessValue
      
      Beta_hat_max <- get_Beta_hat(c(omega_hat), psi_max, omega_hat, X_one_hot, X_h_one_hot)
      
      Beta_hat_matrix <- Beta_hat_max |> 
        matrix(nrow = 1, dimnames = list(psi_max |> plyr::round_any(step_size / 100, floor)))
      
      current_psi_val <- psi_max |> 
        plyr::round_any(step_size, ceiling)
      
      current_branch_val <- branch_max
      
      stopping_branch_val <- branch_max - crit
      
      Beta_hat_matrix <- iterate_while(psi_max = psi_max,
                                       branch_max = branch_max,
                                       initial_guess = Beta_hat_max,
                                       get_Beta_hat = get_Beta_hat,
                                       get_log_likelihood = get_log_likelihood,
                                       omega_hat = omega_hat,
                                       X_one_hot = X_one_hot,
                                       X_h_one_hot = X_h_one_hot,
                                       Y_one_hot = Y_one_hot,
                                       crit = crit,
                                       step_size = step_size)
      
      return(Beta_hat_matrix)
    },
    error = function(e) NULL
    )
  }
}

get_log_L_tilde_mat <- function(Beta_hat_matrices, X_one_hot, Y_one_hot) {
  
  log_likelihood_list <- Beta_hat_matrices |> 
    purrr::compact() |> 
    lapply(\(Beta_hat_matrix) (Beta_hat_matrix |> 
                                 apply(1, \(Beta_hat) (Beta_hat |> 
                                                         get_log_likelihood(X_one_hot, Y_one_hot))))) 
  
  all_labels <- unique(unlist(lapply(log_likelihood_list, names)))
  
  aligned_vectors <- lapply(log_likelihood_list, function(vec) {
    full_vec <- setNames(rep(NA, length(all_labels)), all_labels)  # Create a full vector with NAs
    full_vec[names(vec)] <- vec  # Fill in existing values
    return(full_vec)
  })
  
  do.call(rbind, aligned_vectors)
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
  
  estimate <- matrixStats::colLogSumExps(log_weighted_vals, na.rm = TRUE)
  
  var_estimate <- matrixStats::colVars(log_weighted_vals, na.rm = TRUE)
  
  var_exp_estimate <- matrixStats::colVars(exp(log_weighted_vals), na.rm = TRUE)
  
  return(list(log_weighted_vals = log_weighted_vals,
              estimate = estimate, 
              var_estimate = var_estimate,
              var_exp_estimate = var_exp_estimate))
}

get_log_integrated_likelihood <- function(data,
                                          X_h, 
                                          R,
                                          MC_params,
                                          step_size, 
                                          threshold,
                                          alpha,
                                          chunk_size) {
  
  model <- get_multinomial_logistic_model(data)
  
  m <- model |> 
    model.frame() |> 
    select(starts_with("X")) |> 
    table() |> 
    (\(table) table[X_h$X])()
  
  X_one_hot <- model.matrix(model)
  
  X_level <- X_h |> 
    pull(X) |> 
    as.character() |> 
    as.numeric()
  
  X_h_one_hot <- X_one_hot[1 + m*(X_level - 1),]
  
  Y_one_hot <- data |> 
    pull(Y) |> 
    (\(Y) model.matrix(~ Y)[,-1])()
  
  psi_hat <- get_psi_hat(model, X_h)
  
  omega_hat_list <- get_omega_hat_list(MC_params, R, X_one_hot, X_h_one_hot, Y_one_hot, psi_hat, threshold)
  
  crit <- qchisq(1 - alpha, 1) / 2
  
  Beta_hat_matrices <- get_Beta_hat_matrices(omega_hat_list,
                                             X_one_hot,
                                             X_h,
                                             X_h_one_hot,
                                             Y_one_hot,
                                             step_size, 
                                             model,
                                             crit, 
                                             chunk_size)
  
  log_L_tilde_mat <- get_log_L_tilde_mat(Beta_hat_matrices, X_one_hot, Y_one_hot)
  
  # log_importance_weights <- get_log_importance_weights(U_list, MC_params)
  
  log_importance_weights <- 0
  
  log_L_bar <- get_log_L_bar(log_L_tilde_mat, log_importance_weights, MC_params)
  
  psi_grid <- log_L_tilde_mat |> 
    colnames() |> 
    as.numeric()
  
  return(list(log_L_bar = log_L_bar,
              Beta_hat_matrices = Beta_hat_matrices,
              log_L_tilde_mat = log_L_tilde_mat,
              log_importance_weights = log_importance_weights,
              omega_hat_list = omega_hat_list,
              # U_list = U_list, 
              MC_params = MC_params,
              model = model, 
              X_h = X_h,
              R = nrow(log_L_tilde_mat),
              psi_grid = psi_grid))
}

################################################################################
############################## PROFILE LIKELIHOOD ############################## 
################################################################################

get_log_profile_likelihood <- function(data,
                                       X_h, 
                                       psi_grid_list) {
  
  model <- get_multinomial_logistic_model(data)
  
  X_one_hot <- model.matrix(model)
  
  X_level <- X_h |> 
    pull(X) |> 
    as.character() |> 
    as.numeric()
  
  X_h_one_hot <- X_one_hot[1 + m*(X_level - 1),]
  
  Y_one_hot <- data |> 
    pull(Y) |> 
    (\(Y) model.matrix(~ Y)[,-1])()
  
  Beta_MLE <- get_Beta_MLE(model)
  
  custom_args = list(omega_hat = Beta_MLE, X_one_hot = X_one_hot, X_h_one_hot = X_h_one_hot, use_GA = FALSE)
  
  init_guess <- c(Beta_MLE)
  
  Beta_hats_1 <- iterate_over(psi_grid_list[[1]], init_guess, get_Beta_hat, custom_args, fill_bottom_to_top  = TRUE)
  
  init_guess_2 <- tail(Beta_hats_1, 1)
  
  Beta_hats_2 <- iterate_over(psi_grid_list[[2]], init_guess_2, get_Beta_hat, custom_args, fill_bottom_to_top  = FALSE)
  
  Beta_hats_1 |> 
    rbind(Beta_hats_2) |> 
    apply(1, \(Beta_hat) get_log_likelihood(Beta_hat, X_one_hot, Y_one_hot))
}