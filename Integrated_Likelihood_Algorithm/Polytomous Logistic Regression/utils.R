################################################################################
############################### DATA GENERATION ################################
################################################################################
library(tidyverse)
library(pipeR)
library(iterate)
library(Rcpp)
library(nloptr)
# Rcpp::sourceCpp("../../iterate_while.cpp")
# sourceCpp("../../get_Beta_hat.cpp")

softmax <- function(x) exp(x) / sum(exp(x))

softmax_adj <- function(x) exp(x) / (1 + sum(exp(x)))

entropy <- function(p) -sum(p * log(p), na.rm = TRUE)

PoI_fn <- function(Beta, X_h_one_hot) {
  
  X_h_one_hot %*% cbind(0, Beta) |>
    softmax() |>
    entropy()
}

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
    entropy <- entropy(p)
    
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

################################################################################
#################################### GENERAL ###################################
################################################################################

fit_multinomial_logistic_model <- function(data) nnet::multinom(Y ~ ., data = data)

log_likelihood <- function(Beta, X_one_hot, Y_one_hot) {
  
  p <- ncol(X_one_hot)
  
  J_minus_one <- ncol(Y_one_hot)
  
  Beta <- matrix(Beta,
                 nrow = p,
                 ncol = J_minus_one,
                 byrow = FALSE)
  
  Y_hat <- X_one_hot %*% Beta
  
  sum(rowSums(Y_one_hot * Y_hat) - log(1 + rowSums(exp(Y_hat))))
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

get_psi_hat <- function(model, X_h) {
  
  model |> 
    predict(newdata = X_h, type = "probs") |> 
    entropy()
}

make_omega_hat_obj_fn <- function(X_h_one_hot, J, p, psi_hat) {
  
  function(Beta) {
    
    Beta <- matrix(Beta,
                   nrow = p,
                   ncol = J - 1,
                   byrow = FALSE)
    
    entropy <- PoI_fn(Beta, X_h_one_hot)
    
    return(abs(entropy - psi_hat))
  }
}

make_omega_hat_con_fn <- function(threshold, X_one_hot, Y_one_hot) {
  
  p <- ncol(X_one_hot)
  
  J_minus_one <- ncol(Y_one_hot)
  
  function(Beta) {
    
    Beta <- matrix(Beta,
                   nrow = p,
                   ncol = J_minus_one,
                   byrow = FALSE)
    
    return(-log_likelihood(Beta, X_one_hot, Y_one_hot) - threshold)
  }
}

get_omega_hat <- function(omega_hat_obj_fn, omega_hat_con_fn, J, p, guess_sd) {
  
  init_guess <- rnorm(p * (J - 1), sd = guess_sd)
  
  omega_hat <- nloptr::auglag(x0 = init_guess,
                              fn = omega_hat_obj_fn,
                              hin = omega_hat_con_fn,
                              localsolver = "LBFGS",
                              deprecatedBehavior = FALSE)$par
  
  if (omega_hat_obj_fn(omega_hat) <= 0.1 && omega_hat_con_fn(omega_hat) <= 0) {
    
    omega_hat <- omega_hat |> 
      matrix(nrow = p,
             ncol = J - 1,
             byrow = FALSE)
    
    return(omega_hat)
  }
  
  else {
    
    get_omega_hat(omega_hat_obj_fn, omega_hat_con_fn, J, p, guess_sd)
  }
}

make_Beta_hat_obj_fn <- function(omega_hat, X_one_hot) {
  
  probs <- X_one_hot %*% cbind(0, omega_hat) |>
    (\(mat) mat[,-1])() |> 
    apply(1, softmax_adj) |>
    t()
  
  p <- nrow(omega_hat)
  
  J_minus_one <- ncol(omega_hat)
  
  function(Beta) {
    
    Beta <- matrix(Beta,
                   nrow = p,
                   ncol = J_minus_one,
                   byrow = FALSE)
    
    Y_hat <- X_one_hot %*% Beta
    
    return(-sum(rowSums(probs * Y_hat) - log(1 + rowSums(exp(Y_hat)))))
  }
}

make_Beta_hat_con_fn <- function(psi, X_h_one_hot, J, p) {
  
  function(Beta) {
    
    Beta <- matrix(Beta,
                   nrow = p,
                   ncol = J - 1,
                   byrow = FALSE)
    
    entropy <- PoI_fn(Beta, X_h_one_hot)
    
    return(entropy - psi)
  }
}

get_Beta_hat <- function(psi,
                         Beta_hat_obj_fn,
                         Beta_hat_con_fn,
                         init_guess,
                         J,
                         p,
                         prev_Beta_hat = NULL,
                         lambda,
                         max_retries) {
  
  # Smooth initial guess using a moving average of previous solutions
  if (!is.null(prev_Beta_hat)) {
    init_guess <- 0.9 * prev_Beta_hat + 0.1 * init_guess
  }
  
  # Constraint feasibility check (avoid infeasible starting points)
  if (max(abs(Beta_hat_con_fn(init_guess))) > 1e-4) {
    warning("Initial guess violates constraints. Adjusting...")
    init_guess <- init_guess * 0.95  # Simple heuristic to move towards feasibility
  }
  
  # Regularized objective function to discourage large jumps
  regularized_Beta_hat_obj_fn <- function(Beta) {
    Beta_hat_obj_fn(Beta) + lambda * sum(Beta^2)
  }
  
  # Optimization attempt with retries
  for (attempt in 1:max_retries) {
    result <- tryCatch({
      nloptr::auglag(
        x0 = init_guess,
        fn = regularized_Beta_hat_obj_fn,
        heq = Beta_hat_con_fn,
        localsolver = "SLSQP",
        localtol = 1e-8,
        deprecatedBehavior = FALSE
      )
    }, error = function(e) {
      warning(sprintf("Optimization attempt %d failed: %s", attempt, e$message))
      return(NULL)
    })
    
    # If optimization succeeded, return the result
    if (!is.null(result$par)) {
      Beta_hat <- matrix(result$par, nrow = p, ncol = J - 1, byrow = FALSE)
      return(list(Beta_hat = Beta_hat, prev_Beta_hat = Beta_hat))
    }
    
    # Retry with a perturbed initial guess
    init_guess <- init_guess + rnorm(length(init_guess), sd = 0.1)
  }
  
  # If retries failed, attempt a fallback with a different solver
  warning("All SLSQP attempts failed. Trying COBYLA...")
  result_fallback <- tryCatch({
    nloptr::auglag(
      x0 = init_guess,
      fn = regularized_Beta_hat_obj_fn,
      heq = Beta_hat_con_fn,
      localsolver = "COBYLA",
      localtol = 1e-6,  # Slightly relaxed tolerance
      deprecatedBehavior = FALSE
    )
  }, error = function(e) {
    warning("Fallback optimization with COBYLA failed. Using previous Beta_hat.")
    return(NULL)
  })
  
  # If fallback fails, return previous Beta_hat
  if (is.null(result_fallback$par)) {
    warning("All optimization attempts failed. Using previous Beta_hat.")
    return(list(Beta_hat = prev_Beta_hat, prev_Beta_hat = prev_Beta_hat))
  }
  
  Beta_hat <- matrix(result_fallback$par, nrow = p, ncol = J - 1, byrow = FALSE)
  return(list(Beta_hat = Beta_hat, prev_Beta_hat = Beta_hat))
}

make_omega_hat_branch_fn <- function(omega_hat, 
                                     X_one_hot, 
                                     Y_one_hot, 
                                     X_h_one_hot,
                                     lambda,
                                     max_retries) {
  
  Beta_hat_obj_fn <- make_Beta_hat_obj_fn(omega_hat, X_one_hot)
  
  init_guess <- c(omega_hat)
  
  J <- ncol(omega_hat) + 1
  
  p <- nrow(omega_hat)
  
  prev_Beta_hat <- NULL
  
  function(psi) {
    
    Beta_hat_con_fn <- make_Beta_hat_con_fn(psi, X_h_one_hot, J, p)
    
    Beta_hat_result <- get_Beta_hat(psi,
                                    Beta_hat_obj_fn,
                                    Beta_hat_con_fn,
                                    init_guess,
                                    J,
                                    p,
                                    prev_Beta_hat,
                                    lambda,
                                    max_retries)
    
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

get_psi_endpoints <- function(omega_hat_branch_fn, omega_hat_branch_fn_max, delta, J) {
  
  branch_argmax <- omega_hat_branch_fn_max$branch_argmax
  
  branch_max <- omega_hat_branch_fn_max$branch_max
  
  target_value <- branch_max - delta
  
  root_func <- function(psi) omega_hat_branch_fn(psi) - target_value
  
  # Find left intersection (x < x_max)
  left_intersection <- uniroot(root_func, interval = c(0, branch_argmax))$root
  
  # Find right intersection (x > x_max)
  right_intersection <- uniroot(root_func, interval = c(branch_argmax, log(J)))$root
  
  return(c(left_intersection, right_intersection))
}

get_psi_grid <- function(psi_endpoints, step_size, J) {
  
  psi_endpoints |> 
    (\(x) c(max(0.01, x[1]), min(log(J) - 0.01, x[2])))() |> 
    (\(x) c(plyr::round_any(x[1], step_size, floor), 
            plyr::round_any(x[2], step_size, ceiling)))() |>
    (\(x) seq(x[1], x[2], step_size))() |> 
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

get_branch <- function(X_one_hot, 
                       Y_one_hot,
                       X_h_one_hot,
                       psi_hat,
                       threshold, 
                       alpha,
                       step_size,
                       init_guess_sd,
                       lambda,
                       max_retries) {
  
  J <- ncol(Y_one_hot) + 1
  
  p <- ncol(X_one_hot)
  
  delta <- qchisq(1 - alpha, 1) / 2
  
  psi_endpoints <- NULL
  
  omega_hat_obj_fn <- make_omega_hat_obj_fn(X_h_one_hot, J, p, psi_hat)
  
  omega_hat_con_fn <- make_omega_hat_con_fn(threshold, X_one_hot, Y_one_hot)
  
  while(is.null(psi_endpoints)) {
    
    omega_hat <- get_omega_hat(omega_hat_obj_fn, omega_hat_con_fn, J, p, init_guess_sd)
    
    omega_hat_branch_fn <- make_omega_hat_branch_fn(omega_hat, 
                                                    X_one_hot, 
                                                    Y_one_hot, 
                                                    X_h_one_hot,
                                                    lambda,
                                                    max_retries)
    
    omega_hat_branch_fn_max <- get_omega_hat_branch_fn_max(omega_hat_branch_fn, J)
    
    psi_endpoints <- tryCatch(get_psi_endpoints(omega_hat_branch_fn, 
                                                omega_hat_branch_fn_max, 
                                                delta,
                                                J),
                              error = function(e) NULL)
  }
  
  psi_grid <- get_psi_grid(psi_endpoints, step_size, J)
  
  Beta_hat_obj_fn <- make_Beta_hat_obj_fn(omega_hat, X_one_hot)
  
  init_guess <- rnorm(p * (J - 1), sd = init_guess_sd)
  
  log_L_tilde_df <- data.frame(psi = psi_grid, Integrated = NA)
  
  prev_Beta_hat <- NULL 
  
  for (i in seq_along(psi_grid)) {
    
    psi_val <- psi_grid[i]
    
    Beta_hat_con_fn <- make_Beta_hat_con_fn(psi_val, X_h_one_hot, J, p)
    
    Beta_hat_result <- get_Beta_hat(psi_val,
                                    Beta_hat_obj_fn,
                                    Beta_hat_con_fn,
                                    init_guess,
                                    J,
                                    p,
                                    prev_Beta_hat,
                                    lambda,
                                    max_retries)
    
    Beta_hat <- Beta_hat_result$Beta_hat
    prev_Beta_hat <- Beta_hat_result$prev_Beta_hat 
    
    log_L_tilde_df$Integrated[i] <- log_likelihood(Beta_hat, X_one_hot, Y_one_hot)
    
    init_guess <- 0.9 * c(Beta_hat) + 0.1 * init_guess
  }
  
  return(list(omega_hat = omega_hat,
              log_L_tilde_df = log_L_tilde_df))
}

generate_branches <- function(X_one_hot, 
                              Y_one_hot,
                              X_h_one_hot,
                              psi_hat,
                              threshold, 
                              alpha,
                              step_size,
                              init_guess_sd,
                              num_branches,
                              chunk_size,
                              lambda,
                              max_retries) {
  
  progress_bar <- progressr::progressor(steps = num_branches)
  
  result <- foreach(
    
    i = 1:num_branches,
    .combine = "list",
    .multicombine = TRUE,
    .maxcombine = num_branches,
    .errorhandling = "remove",
    .options.future = list(seed = TRUE,
                           chunk.size = chunk_size)
    
  ) %dofuture% {
    
    progress_bar()  
    
    get_branch(X_one_hot, 
               Y_one_hot,
               X_h_one_hot,
               psi_hat,
               threshold, 
               alpha,
               step_size,
               init_guess_sd,
               lambda,
               max_retries)
  }
  
  omega_hat <- lapply(result, `[[`, 1)  # Extract first elements
  log_L_tilde_df <- lapply(result, `[[`, 2)  # Extract second elements
  
  list(omega_hat = omega_hat,
       log_L_tilde_df = log_L_tilde_df)
}

get_log_L_bar <- function(branches, alpha, J) {
  
  df_list <- branches$log_L_tilde_df
  
  all_psi <- df_list |>
    lapply(function(df) df$psi) |>
    unlist() |>
    unique() |>
    sort()
  
  df_list_filled <- df_list |>
    lapply(function(df) {
      full_join(data.frame(psi = all_psi), df, by = "psi")
    })
  
  merged_df <- reduce(df_list_filled, full_join, by = "psi")
  
  branches_matrix <- merged_df[, -1] |>
    as.matrix() |>
    t() |>
    unname()
  
  colnames(branches_matrix) <- all_psi
  
  delta <- qchisq(1 - alpha, 1) / 2
  
  best_prop <- 1
  
  best_CI_length <- log(J)
  
  for (prop in seq(0, 1, 0.01)) {
    
    psi_vals_to_keep <- colSums(!is.na(branches_matrix)) >= prop * nrow(branches_matrix)
    
    branches_matrix_new <- branches_matrix[, psi_vals_to_keep]
    
    log_R <- branches_matrix_new |> 
      as.data.frame() |> 
      (\(df) !is.na(df))() |> 
      colSums() |> 
      log()
    
    log_L_bar <- matrixStats::colLogSumExps(branches_matrix_new, na.rm = TRUE) - log_R
    
    log_L_bar_df <- data.frame(psi = psi_vals_to_keep |> 
                                 which() |> 
                                 names() |> 
                                 as.numeric(),
                               Integrated = log_L_bar)
    
    spline_fitted_model <- smooth.spline(log_L_bar_df$psi, log_L_bar_df$Integrated)
    
    MLE_data <- optimize(function(psi) predict(spline_fitted_model, psi)$y,
                           lower = log_L_bar_df |>
                             select(psi) |>
                             min(),
                           upper = log_L_bar_df |>
                             select(psi) |>
                             max(),
                           maximum = TRUE) |> 
      data.frame() |> 
      mutate(MLE = as.numeric(maximum),
             Maximum = as.numeric(objective)) |>
      select(MLE, Maximum)
    
    curve <- function(psi) {
      
      lower_psi_val <- spline_fitted_model$x |> min()
      
      upper_psi_val <- spline_fitted_model$x |> max()
      
      if (psi < lower_psi_val) {
        
        return(head(spline_fitted_model$y, 1) - MLE_data$Maximum)
      }
      
      else if (psi > upper_psi_val) {
        
        return(tail(spline_fitted_model$y, 1) - MLE_data$Maximum)
      }
      
      else {
        
        return(predict(spline_fitted_model, psi)$y - MLE_data$Maximum)
      }
    }
    
    lower_bound <- tryCatch(
      
      uniroot(function(psi) curve(psi) + delta,
              interval = c(0, MLE_data$MLE))$root,
      
      error = function(e) return(0)
    )
    
    upper_bound <- tryCatch(
      
      uniroot(function(psi) curve(psi) + delta,
              interval = c(MLE_data$MLE, log(J)))$root,
      
      error = function(e) return(log(J))
    )
    
    CI_length_new <- upper_bound - lower_bound
    
    if (CI_length_new < best_CI_length) {
      
      best_prop <- prop
      best_CI_length = CI_length_new
    }
  }
  
  psi_vals_to_keep <- colSums(!is.na(branches_matrix)) >= best_prop * nrow(branches_matrix)
  
  branches_matrix <- branches_matrix[, psi_vals_to_keep]
  
  log_R <- branches_matrix |> 
    as.data.frame() |> 
    (\(df) !is.na(df))() |> 
    colSums() |> 
    log()
  
  log_L_bar <- matrixStats::colLogSumExps(branches_matrix, na.rm = TRUE) - log_R
  
  log_L_bar_df <- data.frame(psi = psi_vals_to_keep |> 
                               which() |> 
                               names() |> 
                               as.numeric(),
                             Integrated = log_L_bar)
  
  return(list(df = log_L_bar_df,
              branches_matrix = branches_matrix))
}

get_log_integrated_likelihood <- function(data,
                                          X_h,
                                          step_size,
                                          threshold,
                                          alpha,
                                          prop,
                                          init_guess_sd,
                                          num_workers,
                                          chunk_size,
                                          lambda = 1e-4,
                                          max_retries = 3) {
  
  model <- fit_multinomial_logistic_model(data)
  
  m <- model |>
    model.frame() |>
    filter(X == 1) |>
    nrow()
  
  X_one_hot <- model.matrix(model)
  
  h <- X_h |>
    pull(X) |>
    as.character() |>
    as.numeric()
  
  X_h_one_hot <- X_one_hot[1 + m*(h - 1),]
  
  Y_one_hot <- data |>
    pull(Y) |>
    (\(Y) model.matrix(~ Y)[,-1])()
  
  psi_hat <- get_psi_hat(model, X_h)
  
  num_branches <- num_workers * chunk_size
  
  plan(multisession, workers = I(num_workers))
  
  branches <- generate_branches(X_one_hot,
                                Y_one_hot,
                                X_h_one_hot,
                                psi_hat,
                                threshold,
                                alpha,
                                step_size,
                                init_guess_sd,
                                num_branches,
                                chunk_size,
                                lambda,
                                max_retries)
  
  plan(sequential)
  
  log_L_bar <- get_log_L_bar(branches, alpha, J)
  
  return(list(branches = branches,
              log_L_bar = log_L_bar))
}

# ################################################################################
# ############################## PROFILE LIKELIHOOD ############################## 
# ################################################################################
# 
get_log_profile_likelihood <- function(data,
                                       X_h,
                                       step_size,
                                       alpha,
                                       lambda = 1e-4,
                                       max_retries = 3) {
  
  model <- fit_multinomial_logistic_model(data)
  
  m <- model |>
    model.frame() |>
    filter(X == 1) |>
    nrow()
  
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
  
  J <- ncol(Y_one_hot) + 1
  
  p <- ncol(X_one_hot)
  
  delta <- qchisq(1 - alpha, 1) / 2
  
  Beta_MLE_branch_fn <- make_omega_hat_branch_fn(Beta_MLE, 
                                                 X_one_hot, 
                                                 Y_one_hot, 
                                                 X_h_one_hot,
                                                 lambda,
                                                 max_retries)
  
  Beta_MLE_branch_fn_max <- get_omega_hat_branch_fn_max(Beta_MLE_branch_fn, J)
  
  psi_endpoints <- tryCatch(get_psi_endpoints(Beta_MLE_branch_fn, 
                                              Beta_MLE_branch_fn_max, 
                                              delta,
                                              J),
                            error = function(e) NULL)
  
  psi_grid <- get_psi_grid(psi_endpoints, step_size, J)
  
  Beta_hat_obj_fn <- make_Beta_hat_obj_fn(Beta_MLE, X_one_hot)
  
  init_guess <- rnorm(p * (J - 1), sd = 20)
  
  log_L_p_df <- data.frame(psi = psi_grid, Profile = NA)
  
  prev_Beta_hat <- NULL 
  
  log_L_p_df |>
    make_plot("Profile") |>
    plotly::ggplotly() |>
    print()
  
  for (i in seq_along(psi_grid)) {
    
    psi_val <- psi_grid[i]
    
    Beta_hat_con_fn <- make_Beta_hat_con_fn(psi_val, X_h_one_hot, J, p)
    
    Beta_hat_result <- get_Beta_hat(psi_val,
                                    Beta_hat_obj_fn,
                                    Beta_hat_con_fn,
                                    init_guess,
                                    J,
                                    p,
                                    prev_Beta_hat,
                                    lambda,
                                    max_retries)
    
    Beta_hat <- Beta_hat_result$Beta_hat
    prev_Beta_hat <- Beta_hat_result$prev_Beta_hat
    
    log_L_p_df$Profile[i] <- log_likelihood(Beta_hat, X_one_hot, Y_one_hot)
    
    init_guess <- 0.9 * c(Beta_hat) + 0.1 * init_guess
    
    log_L_p_df |>
      make_plot("Profile") |>
      plotly::ggplotly() |>
      print()
    
    Sys.sleep(0.1)
  }
  
  return(log_L_p_df)
}
