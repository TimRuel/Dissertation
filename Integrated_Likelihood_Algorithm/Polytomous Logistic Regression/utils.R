################################################################################
############################### DATA GENERATION ################################
################################################################################
library(tidyverse)
library(pipeR)
library(Rcpp)
library(nloptr)
library(PolytomousUtils)
# sourceCpp("../../polytomous_utils.cpp")

softmax <- function(x) exp(x) / sum(exp(x))

softmax_adj <- function(x) exp(x) / (1 + sum(exp(x)))

entropy <- function(p) -sum(p * log(p), na.rm = TRUE)

# PoI_fn <- function(Beta, X_h_one_hot) {
#   
#   X_h_one_hot %*% cbind(0, Beta) |>
#     softmax() |>
#     entropy()
# }
# 
# PoI_fn2 <- function(Beta, X_h_one_hot) {
#   
#   X_h %*% cbind(0, Beta) |> 
#     apply(1, softmax) |> 
#     t() |> 
#     colMeans() |> 
#     entropy()
# }

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
    p <- k |> 
      runif(min = -10, max = 10) |> 
      softmax() 
    
    # Adjust probabilities to ensure all components >= epsilon
    p <- adjust_probabilities(p, epsilon)
    
    # Calculate entropy
    entropy <- entropy(p)
    
    # Check if entropy is within the target range
    if (entropy >= lower_bound && entropy <= upper_bound) {
      return(sample(p))
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

get_X <- function(c, m, contrast) {
  
  X <- 1:c |> 
    rep(each = m) |> 
    factor()
  
  contrasts(X) <- contrast
  
  return(X)
}

################################################################################
#################################### GENERAL ###################################
################################################################################

fit_multinomial_logistic_model <- function(data, formula) {
  
  # formula <- substitute(formula_expr)
  nnet::multinom(formula, data = data, trace = FALSE)
}

log_likelihood <- function(Beta, X_one_hot, Y_one_hot) {
  
  p <- ncol(X_one_hot)
  
  Jm1 <- ncol(Y_one_hot)
  
  Beta <- matrix(Beta,
                 nrow = p,
                 ncol = Jm1,
                 byrow = FALSE)
  
  Y_hat <- X_one_hot %*% Beta
  
  sum(rowSums(Y_one_hot * Y_hat) - matrixStats::rowLogSumExps(cbind(0, Y_hat)))
}

likelihood <- function(Beta, X_one_hot, Y_one_hot) exp(log_likelihood(Beta, X_one_hot, Y_one_hot))

get_Beta_MLE <- function(model) {
  
  J <- length(model$lev)
  
  model |> 
    coef() |> 
    unname() |> 
    matrix(ncol = J - 1,
           byrow = TRUE) 
}

get_psi_hat <- function(model, data, h) {
  
  model |> 
    predict(data, type = "probs") |>
    as.data.frame() |>
    mutate(X1_level = rep(X1_levels, times = m)) |>
    aggregate(. ~ X1_level, data = _, FUN = mean) |>
    (\(df) setNames(data.frame(t(df[,-1])), df[,1]))() |> 
    apply(2, entropy) |> 
    (\(vec) vec[h])()
}

get_omega_hat <- function(omega_hat_eq_con_fn, omega_hat_ineq_con_fn, Jm1, p, init_guess_sd) {
  
  init_guess <- rnorm(p * Jm1, sd = init_guess_sd)
  
  omega_hat <- nloptr::auglag(x0 = init_guess,
                              fn = function(Beta) 0,
                              heq = omega_hat_eq_con_fn,
                              hin = omega_hat_ineq_con_fn,
                              localsolver = "LBFGS",
                              deprecatedBehavior = FALSE)$par
  
  if (abs(omega_hat_eq_con_fn(omega_hat)) <= 0.1 && omega_hat_ineq_con_fn(omega_hat) <= 0) {
    
    omega_hat <- matrix(omega_hat, nrow = p, ncol = Jm1, byrow = FALSE)
    
    return(omega_hat)
  }
  
  else {
    
    return(get_omega_hat(omega_hat_eq_con_fn, omega_hat_ineq_con_fn, Jm1, p, init_guess_sd))
  }
}

safe_auglag <- purrr::possibly(nloptr::auglag, otherwise = NULL)

get_Beta_hat <- function(X_one_hot, 
                         X_h_one_hot,
                         omega_hat,
                         psi,
                         init_guess,
                         Jm1,
                         p,
                         prev_Beta_hat = NULL,
                         lambda,
                         max_retries,
                         max_eval) {
  
  # Smooth initial guess using a moving average of previous solutions
  if (!is.null(prev_Beta_hat)) {
    phi <- max(0.1, min(0.9, 1 - 1 / max_retries))
    init_guess <- phi * prev_Beta_hat + (1 - phi) * init_guess
  }
  
  # Constraint feasibility check (avoid infeasible starting points)
  if (sum(init_guess^2) > 1000) {
    init_guess <- init_guess * 0.95  # Simple heuristic to move towards feasibility
  }
  
  obj_fn <- function(Beta) Beta_hat_obj_fn_rcpp(Beta, X_one_hot, omega_hat, lambda, Jm1, p, n)
  con_fn <- function(Beta) Beta_hat_con_fn_rcpp(Beta, X_h_one_hot, psi, Jm1, p)
    
  # Optimization attempt with retries
  for (attempt in 1:max_retries) {
    
    result <- safe_auglag(
      x0 = init_guess,
      fn = obj_fn,
      heq = con_fn,
      localsolver = "SLSQP",
      localtol = 1e-8,
      deprecatedBehavior = FALSE,
      control = list(on.error = "ignore",
                     maxeval = max_eval)
    )
    
    if (!is.null(result$par) && result$convergence == 4) {
      
      Beta_hat <- matrix(result$par, nrow = p, ncol = Jm1, byrow = FALSE)
      return(list(Beta_hat = Beta_hat, prev_Beta_hat = Beta_hat))
    }
    
    # Perturbation for next attempt
    init_guess <- init_guess + rnorm(length(init_guess), sd = 0.1)
  }
  
  # If retries failed, attempt a fallback with a different solver
  warning("All SLSQP attempts failed. Trying COBYLA...")
  result_fallback <- safe_auglag(
    x0 = init_guess,
    fn = obj_fn,
    heq = con_fn,
    localsolver = "COBYLA",
    localtol = 1e-6,  # Slightly relaxed tolerance
    deprecatedBehavior = FALSE,
    control = list(on.error = "ignore",
                   maxeval = max_eval)
  )
  
  # If fallback fails, return previous Beta_hat
  if (is.null(result_fallback$par)) {
    warning("All optimization attempts failed. Using previous Beta_hat.")
    return(list(Beta_hat = prev_Beta_hat, prev_Beta_hat = prev_Beta_hat))
  }
  
  Beta_hat <- matrix(result_fallback$par, nrow = p, ncol = Jm1, byrow = FALSE)
  return(list(Beta_hat = Beta_hat, prev_Beta_hat = Beta_hat))
}

make_omega_hat_branch_fn <- function(omega_hat, 
                                     X_one_hot, 
                                     Y_one_hot, 
                                     X_h_one_hot,
                                     Jm1,
                                     p,
                                     lambda,
                                     max_retries,
                                     max_eval) {
  
  init_guess <- c(omega_hat)
  
  prev_Beta_hat <- NULL
  
  function(psi) {
    
    Beta_hat_result <- get_Beta_hat(X_one_hot, 
                                    X_h_one_hot,
                                    omega_hat,
                                    psi,
                                    init_guess,
                                    Jm1,
                                    p,
                                    prev_Beta_hat,
                                    lambda,
                                    max_retries,
                                    max_eval)
    
    Beta_hat <- Beta_hat_result$Beta_hat
    prev_Beta_hat <- Beta_hat_result$prev_Beta_hat
    
    return(log_likelihood(Beta_hat, X_one_hot, Y_one_hot))
  }
}

get_omega_hat_branch_fn_max <- function(omega_hat_branch_fn, J) {
  
  opt_result <- optimize(omega_hat_branch_fn, 
                         interval = c(0, log(J)), 
                         maximum = TRUE)
  branch_argmax <- opt_result$maximum
  branch_max <- opt_result$objective
  
  return(list(branch_argmax = branch_argmax, branch_max = branch_max))
}

get_psi_endpoints <- function(omega_hat_branch_fn, 
                              omega_hat_branch_fn_max, 
                              delta, 
                              J) {
  
  branch_argmax <- omega_hat_branch_fn_max$branch_argmax
  
  branch_max <- omega_hat_branch_fn_max$branch_max
  
  target_value <- branch_max - delta
  
  root_func <- function(psi) omega_hat_branch_fn(psi) - target_value
  
  left_intersection <- uniroot(root_func, interval = c(0.01, branch_argmax))$root
  
  right_intersection <- uniroot(root_func, interval = c(branch_argmax, log(J) - 0.01))$root
  
  return(c(left_intersection, right_intersection))
}

safe_get_psi_endpoints <- purrr::possibly(get_psi_endpoints, otherwise = NULL)

get_psi_grid_1 <- function(psi_endpoints, step_size, J) {
  
  psi_endpoints |>
    (\(x) c(max(0, x[1]), min(log(J), x[2])))() |>
    (\(x) c(plyr::round_any(x[1], step_size, floor),
            plyr::round_any(x[2], step_size, ceiling)))() |>
    (\(x) seq(x[1], x[2], step_size))() |>
    (\(x) c(x[-length(x)], min(log(J), tail(x, 1))))() |>
    round(6)
}

get_psi_grid_2 <- function(branch_specs, step_size, J, quantiles) {
  
  branch_specs |> 
    purrr::map(\(branch_spec) branch_spec$psi_endpoints) |> 
    (\(x) do.call(rbind, x))() |> 
    as.data.frame() |> 
    purrr::set_names(c("Lower", V2 = "Upper")) |> 
    dplyr::summarise(Lower = quantile(Lower, quantiles[1]),
                     Upper = quantile(Upper, quantiles[2])) |> 
    c() |> 
    unlist() |> 
    (\(x) c(plyr::round_any(x[1], step_size, floor), 
            plyr::round_any(x[2], step_size, ceiling)))() |> 
    (\(x) seq(x[1], x[2], step_size))() |>
    (\(x) c(x[-length(x)], min(log(J), tail(x, 1))))() |>
    round(6)
}

################################################################################
############################ INTEGRATED LIKELIHOOD ############################# 
################################################################################

make_plot <- function(df, y_var) {
  
  df |> 
    ggplot(aes(x = psi, y = y_var)) + 
    geom_point(color = "cyan", size = 3, alpha = 0.7) + 
    theme_minimal(base_size = 15) +  # Minimal theme with a larger base font size
    theme(
      plot.background = element_rect(fill = "#2E2E2E", color = NA),  # Dark background for the whole plot
      panel.background = element_rect(fill = "#3A3A3F", color = "#1A1A1A", size = 2),  # Lighter panel with a border
      axis.text = element_text(color = "white"),  # White axis labels
      axis.title = element_text(color = "white"),  # White axis titles
      plot.title = element_text(color = "white", size = 18, face = "bold"),  # White title
      plot.caption = element_text(color = "gray", size = 10),  # Gray caption
      panel.grid = element_line(color = "gray30", linetype = "dashed")  # Subtle grid lines
    ) + 
    labs(
      x = "\u03C8",
      y = "Integrated Log-Likelihood Branch"
    )
}

get_branch_spec <- function(omega_hat_eq_con_fn,
                            omega_hat_ineq_con_fn,
                            X_one_hot, 
                            Y_one_hot, 
                            X_h_one_hot,
                            J,
                            p,
                            init_guess_sd,
                            delta,
                            lambda,
                            max_retries,
                            max_eval) {
  
  psi_endpoints <- NULL
  
  while (is.null(psi_endpoints)) {
    
    omega_hat <- get_omega_hat(omega_hat_eq_con_fn, omega_hat_ineq_con_fn, J - 1, p, init_guess_sd)
    
    omega_hat_branch_fn <- make_omega_hat_branch_fn(omega_hat, 
                                                    X_one_hot, 
                                                    Y_one_hot, 
                                                    X_h_one_hot,
                                                    J - 1,
                                                    p,
                                                    lambda,
                                                    max_retries,
                                                    max_eval)
    
    omega_hat_branch_fn_max <- get_omega_hat_branch_fn_max(omega_hat_branch_fn, J)
    
    psi_endpoints <- safe_get_psi_endpoints(omega_hat_branch_fn, 
                                            omega_hat_branch_fn_max, 
                                            delta,
                                            J)
  }
  
  return(list(omega_hat = omega_hat, 
              psi_endpoints = psi_endpoints))
}

generate_branch_specs <- function(data,
                                  h,
                                  init_guess_sd,
                                  alpha,
                                  num_workers,
                                  chunk_size,
                                  lambda,
                                  max_retries,
                                  max_eval) {
  
  model <- fit_multinomial_logistic_model(data, formula)
  
  X_one_hot <- model.matrix(model)
  
  X_h_one_hot <- extract_X_h(X_one_hot, h, drop_zero_cols = FALSE)
  
  Y_one_hot <- data |>
    pull(Y) |>
    (\(Y) model.matrix(~ Y)[,-1])()
  
  J <- ncol(Y_one_hot) + 1
  
  p <- ncol(X_one_hot)
  
  n <- nrow(X_one_hot)
  
  psi_hat <- get_psi_hat(model, data, h)
  
  omega_hat_eq_con_fn <- function(Beta) omega_hat_eq_con_fn_rcpp(Beta, X_h_one_hot, J - 1, p, psi_hat)
  
  omega_hat_ineq_con_fn <- function(Beta) omega_hat_ineq_con_fn_rcpp(Beta, X_one_hot, Y_one_hot, J - 1, p, n, threshold)
  
  delta <- qchisq(1 - alpha, 1) / 2
  
  num_branches <- num_workers * chunk_size
  
  progress_bar <- progressr::progressor(steps = 2 * num_branches)
  
  plan(multisession, workers = I(num_workers))
  
  branch_specs <- foreach(
    
    i = 1:num_branches,
    .combine = "list",
    .multicombine = TRUE,
    .maxcombine = num_branches,
    .errorhandling = "remove",
    .options.future = list(seed = TRUE,
                           chunk.size = chunk_size,
                           packages = c("PolytomousUtils", "nloptr"))
    
  ) %dofuture% {
    
    progress_bar()
    
    branch_spec <- get_branch_spec(omega_hat_eq_con_fn,
                                   omega_hat_ineq_con_fn,
                                   X_one_hot, 
                                   Y_one_hot, 
                                   X_h_one_hot,
                                   J,
                                   p,
                                   init_guess_sd,
                                   delta,
                                   lambda,
                                   max_retries,
                                   max_eval)
    
    progress_bar()
    
    return(branch_spec)
  }
  
  plan(sequential)
  
  return(branch_specs)
}

get_branch <- function(branch_spec,
                       X_one_hot, 
                       Y_one_hot,
                       X_h_one_hot,
                       J, 
                       p,
                       n,
                       psi_hat,
                       delta,
                       step_size,
                       init_guess_sd,
                       lambda,
                       max_retries,
                       max_eval) {
  
  psi_grid <- get_psi_grid_1(branch_spec$psi_endpoints, step_size, J)
  
  # progress_bar <- progressr::progressor(along = psi_grid)
  
  init_guess <- rnorm(p * (J - 1), sd = init_guess_sd)
  
  log_L_tilde_vec <- numeric(length(psi_grid))
  
  prev_Beta_hat <- NULL 
  
  for (i in seq_along(psi_grid)) {
    
    # progress_bar()
    
    psi <- psi_grid[i]
    
    Beta_hat_result <- get_Beta_hat(X_one_hot, 
                                    X_h_one_hot,
                                    branch_spec$omega_hat,
                                    psi,
                                    init_guess,
                                    J - 1,
                                    p,
                                    prev_Beta_hat,
                                    lambda,
                                    max_retries,
                                    max_eval)
    
    Beta_hat <- Beta_hat_result$Beta_hat
    prev_Beta_hat <- Beta_hat_result$prev_Beta_hat 
    
    log_L_tilde_vec[i] <- log_likelihood_rcpp(Beta_hat, X_one_hot, Y_one_hot, J - 1, p, n)
    
    init_guess <- 0.9 * c(Beta_hat) + 0.1 * init_guess
  }
  
  log_L_tilde_df <- data.frame(psi = psi_grid,
                               Integrated = log_L_tilde_vec)
  
  return(list(omega_hat = branch_spec$omega_hat,
              log_L_tilde_df = log_L_tilde_df))
}

generate_branches <- function(branch_specs,
                              X_one_hot, 
                              Y_one_hot,
                              X_h_one_hot,
                              J, 
                              p,
                              n,
                              psi_hat,
                              delta,
                              step_size,
                              quantiles,
                              init_guess_sd,
                              chunk_size,
                              lambda,
                              max_retries,
                              max_eval) {
  
  psi_grid <- get_psi_grid_2(branch_specs, step_size, J, quantiles)
  
  # num_steps <- 2 * length(psi_grid) * length(branch_specs)
  # 
  # progress_bar <- progressr::progressor(steps = num_steps)
  
  result <- foreach(
    
    branch_spec = branch_specs,
    .combine = "list",
    .multicombine = TRUE,
    .maxcombine = num_branches,
    .errorhandling = "remove",
    .options.future = list(seed = TRUE,
                           chunk.size = chunk_size,
                           packages = c("PolytomousUtils", "nloptr"))
    
  ) %dofuture% {
    
    omega_hat <- branch_spec$omega_hat
    
    init_guess <- rnorm(p * (J - 1), sd = init_guess_sd)
    
    log_L_tilde_vec <- numeric(length(psi_grid))
    
    prev_Beta_hat <- NULL 
    
    for (i in seq_along(psi_grid)) {
      
      # progress_bar()
      
      psi <- psi_grid[i]
      
      Beta_hat_result <- get_Beta_hat(X_one_hot, 
                                      X_h_one_hot,
                                      omega_hat,
                                      psi,
                                      init_guess,
                                      J - 1,
                                      p,
                                      prev_Beta_hat,
                                      lambda,
                                      max_retries,
                                      max_eval)
      
      Beta_hat <- Beta_hat_result$Beta_hat
      prev_Beta_hat <- Beta_hat_result$prev_Beta_hat 
      
      log_L_tilde_vec[i] <- log_likelihood_rcpp(Beta_hat, X_one_hot, Y_one_hot, J - 1, p, n)
      
      init_guess <- 0.9 * c(Beta_hat) + 0.1 * init_guess
      
      # progress_bar()
    }
    
    log_L_tilde_df <- data.frame(psi = psi_grid, 
                                 Integrated = log_L_tilde_vec)
    
    return(list(omega_hat = omega_hat,
                log_L_tilde_df = log_L_tilde_df))
  }
  
  omega_hat <- lapply(result, `[[`, 1)  # Extract first elements
  log_L_tilde_df <- lapply(result, `[[`, 2)  # Extract second elements
  
  list(omega_hat = omega_hat,
       log_L_tilde_df = log_L_tilde_df)
}

get_log_L_bar <- function(branches) {
  
  df_list <- branches$log_L_tilde_df
  
  merged_df <- reduce(df_list, full_join, by = "psi")
  
  branches_matrix <- merged_df[, -1] |>
    as.matrix() |>
    t() |>
    unname()
  
  colnames(branches_matrix) <- merged_df$psi
  
  log_R <- branches_matrix |> 
    nrow() |> 
    log()
  
  log_L_bar <- matrixStats::colLogSumExps(branches_matrix, na.rm = TRUE) - log_R
  
  log_L_bar_df <- data.frame(psi = merged_df$psi,
                             Integrated = log_L_bar)
  
  return(list(df = log_L_bar_df,
              branches_matrix = branches_matrix))
}

get_log_integrated_likelihood <- function(branch_specs,
                                          data,
                                          h,
                                          alpha,
                                          step_size,
                                          quantiles,
                                          init_guess_sd,
                                          num_workers,
                                          chunk_size,
                                          lambda,
                                          max_retries,
                                          max_eval) {
  
  model <- fit_multinomial_logistic_model(data, formula)
  
  X_one_hot <- model.matrix(model)
  
  X_h_one_hot <- extract_X_h(X_one_hot, h, drop_zero_cols = FALSE)
  
  Y_one_hot <- data |>
    pull(Y) |>
    (\(Y) model.matrix(~ Y)[,-1])()
  
  psi_hat <- get_psi_hat(model, data, h)
  
  delta <- qchisq(1 - alpha, 1) / 2
  
  J <- ncol(Y_one_hot) + 1
  
  p <- ncol(X_one_hot)
  
  n <- nrow(X_one_hot)
  
  plan(multisession, workers = I(num_workers))
  
  branches <- generate_branches(branch_specs,
                                X_one_hot, 
                                Y_one_hot,
                                X_h_one_hot,
                                J, 
                                p,
                                n,
                                psi_hat,
                                delta,
                                step_size,
                                quantiles,
                                init_guess_sd,
                                chunk_size,
                                lambda,
                                max_retries,
                                max_eval)
  
  plan(sequential)
  
  log_L_bar <- get_log_L_bar(branches)
  
  return(list(branches = branches,
              log_L_bar = log_L_bar))
}

# ################################################################################
# ############################## PROFILE LIKELIHOOD ############################## 
# ################################################################################

make_profile_plot <- function(df) {
  
  df |> 
    ggplot(aes(x = psi, y = Profile)) + 
    geom_point(color = "cyan", size = 3, alpha = 0.7) + 
    theme_minimal(base_size = 15) +  # Minimal theme with a larger base font size
    theme(
      plot.background = element_rect(fill = "#2E2E2E", color = NA),  # Dark background for the whole plot
      panel.background = element_rect(fill = "#3A3A3F", color = "#1A1A1A", size = 2),  # Lighter panel with a border
      axis.text = element_text(color = "white"),  # White axis labels
      axis.title = element_text(color = "white"),  # White axis titles
      plot.title = element_text(color = "white", size = 18, face = "bold"),  # White title
      plot.caption = element_text(color = "gray", size = 10),  # Gray caption
      panel.grid = element_line(color = "gray30", linetype = "dashed")  # Subtle grid lines
    ) + 
    labs(
      x = "\u03C8",
      y = "Profile Log-Likelihood"
    )
}

get_log_profile_likelihood <- function(data,
                                       h,
                                       step_size,
                                       alpha,
                                       init_guess_sd,
                                       lambda,
                                       max_retries,
                                       max_eval) {
  
  model <- fit_multinomial_logistic_model(data, formula)
  
  X_one_hot <- model.matrix(model)
  
  X_h_one_hot <- extract_X_h(X_one_hot, h, drop_zero_cols = FALSE)
  
  Y_one_hot <- data |>
    pull(Y) |>
    (\(Y) model.matrix(~ Y)[,-1])()
  
  Beta_MLE <- get_Beta_MLE(model)
  
  J <- ncol(Y_one_hot) + 1
  
  p <- ncol(X_one_hot)
  
  n <- nrow(X_one_hot)
  
  delta <- qchisq(1 - alpha, 1) / 2
  
  Beta_MLE_branch_fn <- make_omega_hat_branch_fn(Beta_MLE, 
                                                 X_one_hot, 
                                                 Y_one_hot, 
                                                 X_h_one_hot,
                                                 J - 1,
                                                 p,
                                                 lambda,
                                                 max_retries,
                                                 max_eval)
  
  Beta_MLE_branch_fn_max <- get_omega_hat_branch_fn_max(Beta_MLE_branch_fn, J)
  
  psi_endpoints <- safe_get_psi_endpoints(Beta_MLE_branch_fn, 
                                          Beta_MLE_branch_fn_max, 
                                          delta,
                                          J)
  
  psi_grid <- get_psi_grid_1(psi_endpoints, step_size, J)
  
  init_guess <- rnorm(p * (J - 1), sd = init_guess_sd)
  
  log_L_p_vec <- numeric(length(psi_grid))
  
  prev_Beta_hat <- NULL 
  
  for (i in seq_along(psi_grid)) {
    
    psi <- psi_grid[i]
    
    Beta_hat_result <- get_Beta_hat(X_one_hot, 
                                    X_h_one_hot,
                                    Beta_MLE,
                                    psi,
                                    init_guess,
                                    J - 1,
                                    p,
                                    prev_Beta_hat,
                                    lambda,
                                    max_retries,
                                    max_eval)
    
    Beta_hat <- Beta_hat_result$Beta_hat
    prev_Beta_hat <- Beta_hat_result$prev_Beta_hat
    
    log_L_p_vec[i] <- log_likelihood(Beta_hat, X_one_hot, Y_one_hot)
    
    init_guess <- 0.9 * c(Beta_hat) + 0.1 * init_guess
  }
  
  log_L_p_df <- data.frame(psi = psi_grid, 
                           Profile = log_L_p_vec)
  
  print(make_profile_plot(log_L_p_df))
  
  return(log_L_p_df)
}
