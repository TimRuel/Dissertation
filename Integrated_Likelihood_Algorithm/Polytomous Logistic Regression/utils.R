library(tidyverse)
library(pipeR)
library(Rcpp)
library(nloptr)
library(PolytomousUtils)
# sourceCpp("../../polytomous_utils.cpp")

# True Parameter Generation -----------------------------------------------

get_entropy_ranges <- function(num_classes, levels) {
  
  num_ranges <- length(levels)
  
  entropy_ranges <- seq(0.6, log(num_classes) - 0.2, length.out = num_ranges + 1) |> 
    (\(x) mapply(c, x[-length(x)], x[-1], SIMPLIFY = FALSE))() |> 
    purrr::map(\(x) {
      
      midpoint <- mean(x)
      desired_length <- (x[2] - x[1]) * 0.9
      return(midpoint + c(-1, 1) * desired_length / 2)}
    )
  
  names(entropy_ranges) <- rev(levels)
  
  return(rev(entropy_ranges))
}

compute_probabilities <- function(intercept, slope, X2_vals) {
  
  eta <- matrix(intercept, nrow = length(X2_vals), ncol = length(intercept), byrow = TRUE) + outer(X2_vals, slope, `*`)
  exp_eta <- exp(eta)
  
  denom <- 1 + rowSums(exp_eta)
  probs <- cbind(1 / denom, exp_eta / denom)
  
  return(probs)
}

Beta_0_objective_fn <- function(params, X2_vals, H_min, H_max, Beta2) {
  
  Jm1 <- length(params) / 2
  intercept <- params[1:Jm1]
  slope <- params[(Jm1+1):(2*Jm1)] + Beta2
  
  prob_matrix <- compute_probabilities(intercept, slope, X2_vals)
  entropies <- apply(prob_matrix, 1, entropy)
  
  return(-diff(range(entropies)))
}

Beta_0_constraint_fn <- function(params, X2_vals, H_min, H_max, Beta2) {
  
  Jm1 <- length(params) / 2
  intercept <- params[1:Jm1]
  slope <- params[(Jm1+1):(2*Jm1)] + Beta2
  
  prob_matrix <- compute_probabilities(intercept, slope, X2_vals)
  entropies <- apply(prob_matrix, 1, entropy)
  
  return(c(H_min - min(entropies), max(entropies) - H_max))
}

optimize_Beta_0 <- function(J, X2_interval, H_min, H_max, num_vals, Beta2) {
  
  X2_vals <- seq(X2_interval[1], X2_interval[2], length.out = num_vals) 
  
  Jm1 <- J - 1
  init_params <- rnorm(2 * Jm1)
  
  result <- auglag(
    x0 = init_params,
    fn = function(params) Beta_0_objective_fn(params, X2_vals, H_min, H_max, Beta2),
    lower = rep(-10, length(init_params)),
    upper = rep(10, length(init_params)),
    hin = function(params) Beta_0_constraint_fn(params, X2_vals, H_min, H_max, Beta2),
    deprecatedBehavior = FALSE)
  
  params <- result$par
  
  return(list(intercept = params[1:Jm1], slope = params[(Jm1+1):(2*Jm1)]))
}

get_Beta_0 <- function(X1_levels, p, J, X2_intervals, num_vals) {
  
  X1_ref_level <- X1_levels[1]
  X1_main_effect_names <- paste0("X1", X1_levels)
  X1X2_interaction_names <- paste0("X1", setdiff(X1_levels, X1_ref_level)) |> 
    paste0(":X2")
  
  Beta_0 <- matrix(NA,
                   nrow = p,
                   ncol = J - 1)
  rownames(Beta_0) <- c(X1_main_effect_names, "X2", X1X2_interaction_names)
  
  entropy_ranges <- get_entropy_ranges(J, X1_levels)
  
  for (h in X1_levels) {
    
    X2_interval <- X2_intervals[[h]]
    H_range <- entropy_ranges[[h]]
    
    if (h == X1_ref_level) {
      
      Beta2 <- 0
      
      params <- optimize_Beta_0(J, X2_interval, H_range[1], H_range[2], num_vals, Beta2)
      
      Beta_0[paste0("X1", h), ] <- params$intercept
      Beta_0["X2", ] <- params$slope
    } else {
      
      Beta2 <- Beta_0["X2", ]
      
      params <- optimize_Beta_0(J, X2_interval, H_range[1], H_range[2], num_vals, Beta2)
      
      Beta_0[grepl(h, rownames(Beta_0)), ] <- unname(rbind(params$intercept, params$slope))
    }
  }
  
  return(Beta_0)
}

# General -----------------------------------------------------------------

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
PoI_fn2 <- function(Beta, X_h_one_hot) {

  X_h_one_hot %*% cbind(0, Beta) |>
    apply(1, softmax) |>
    t() |>
    colMeans() |>
    entropy()
}

get_X1_levels <- function(X_design) {
  
  is_X1_main_effect <- attr(X_design, "assign") == 1
  X1_main_effects <- colnames(X_design)[is_X1_main_effect]
  X1_levels <- substr(X1_main_effects, nchar(X1_main_effects), nchar(X1_main_effects))
  return(X1_levels)
}

get_X1_ref_level <- function(X_design) {
  
  X1_levels <- get_X1_levels(X_design)
  
  pattern <- paste0("[", paste(X1_levels, collapse = ""), "]")
  
  is_interaction <- attr(X_design, "assign") == 3
  
  interaction_terms <- colnames(X_design)[is_interaction]
  
  X1_interactions <- pattern |> 
    gregexpr(interaction_terms) |> 
    regmatches(interaction_terms, m = _) |> 
    unlist()
  
  X1_ref_level <- setdiff(X1_levels, X1_interactions)
  
  return(X1_ref_level)
}  

extract_X_h <- function(X_design, h, drop_zero_cols = FALSE) {
  
  X1_levels <- get_X1_levels(X_design)
  
  X1_ref_level <- get_X1_ref_level(X_design)
  
  X1_main_effect_names <- paste0("X1", X1_levels)
  
  X1X2_interaction_names <- paste0("X1", setdiff(X1_levels, X1_ref_level)) |> 
    paste0(":X2")
  
  colnames(X_design) <- c(X1_main_effect_names, "X2", X1X2_interaction_names)
  
  rows_to_keep <- X_design[, paste0("X1", h)] == 1
  
  X_h <- X_design[rows_to_keep,]
  
  if (drop_zero_cols) {
    
    cols_to_keep <- grepl(paste0("^X1", h, "$|^X2$|^X1", h, ":X2$"), colnames(X_design))
    
    X_h <- X_h[, cols_to_keep] |> 
      as.matrix(ncol = length(cols_to_keep))
    
    colnames(X_h) <- colnames(X_design)[cols_to_keep]
  }
  
  return(X_h)
}

generate_X2_samples <- function(arg_list) {
  
  X2_samples <- c()
  
  for (h in names(arg_list)) {
    
    args <- arg_list[[h]]

    A <- args$interval[1]
    B <- args$interval[2]
    
    s1 <- args$shape[[1]]
    s2 <- args$shape[[2]]

    beta_samples <- rbeta(args$n_samples, s1, s2)
    
    X2_transformed <- A + (B - A) * beta_samples
    
    X2_samples <- c(X2_samples, X2_transformed)
  }
  
  return(X2_samples)
}

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
                         max_retries,
                         max_eval,
                         maxtime) {
  
  if (!is.null(prev_Beta_hat)) {
    phi <- max(0.1, min(0.9, 1 - 1 / max_retries))
    init_guess <- phi * prev_Beta_hat + (1 - phi) * init_guess
  }
  
  obj_fn <- function(Beta) Beta_hat_obj_fn_rcpp(Beta, X_one_hot, omega_hat, Jm1, p, n)
  con_fn <- function(Beta) Beta_hat_con_fn_rcpp(Beta, X_h_one_hot, psi, Jm1, p)
  
  result <- safe_auglag(
    x0 = init_guess,
    fn = obj_fn,
    heq = con_fn,
    localsolver = "SLSQP",
    localtol = 1e-3,
    deprecatedBehavior = FALSE,
    control = list(on.error = "ignore",
                   # maxeval = max_eval,
                   maxtime = maxtime)
  )
  
  if (!is.null(result$par) && result$convergence > 0) {
    
    Beta_hat <- matrix(result$par, nrow = p, ncol = Jm1, byrow = FALSE)
    return(list(Beta_hat = Beta_hat, prev_Beta_hat = Beta_hat))
  }
  
  return(list(Beta_hat = prev_Beta_hat, prev_Beta_hat = prev_Beta_hat))
}

# get_Beta_hat <- function(X_one_hot, 
#                          X_h_one_hot,
#                          omega_hat,
#                          psi,
#                          init_guess,
#                          Jm1,
#                          p,
#                          prev_Beta_hat = NULL,
#                          max_retries,
#                          max_eval,
#                          maxtime) {
#   
#     if (!is.null(prev_Beta_hat)) {
#     phi <- max(0.1, min(0.9, 1 - 1 / max_retries))
#     init_guess <- phi * prev_Beta_hat + (1 - phi) * init_guess
#   }
#   
#   obj_fn <- function(Beta) Beta_hat_obj_fn_rcpp(Beta, X_one_hot, omega_hat, Jm1, p, n)
#   con_fn <- function(Beta) Beta_hat_con_fn_rcpp(Beta, X_h_one_hot, psi, Jm1, p)
#     
#   for (attempt in 1:max_retries) {
#     
#     result <- safe_auglag(
#       x0 = init_guess,
#       fn = obj_fn,
#       heq = con_fn,
#       localsolver = "SLSQP",
#       localtol = 1e-3,
#       deprecatedBehavior = FALSE,
#       control = list(on.error = "ignore",
#                      maxeval = max_eval,
#                      maxtime = 1.5)
#     )
#     
#     if (!is.null(result$par) && result$convergence > 0) {
#       
#       Beta_hat <- matrix(result$par, nrow = p, ncol = Jm1, byrow = FALSE)
#       return(list(Beta_hat = Beta_hat, prev_Beta_hat = Beta_hat))
#     }
#     
#     init_guess <- init_guess + rnorm(length(init_guess), sd = 0.1)
#   }
#   
#   result_fallback <- safe_auglag(
#     x0 = init_guess,
#     fn = obj_fn,
#     heq = con_fn,
#     localsolver = "COBYLA",
#     localtol = 1e-3,  
#     deprecatedBehavior = FALSE,
#     control = list(on.error = "ignore",
#                    maxeval = max_eval)
#   )
#   
#   if (is.null(result_fallback$par)) {
#     return(list(Beta_hat = prev_Beta_hat, prev_Beta_hat = prev_Beta_hat))
#   }
#   
#   Beta_hat <- matrix(result_fallback$par, nrow = p, ncol = Jm1, byrow = FALSE)
#   return(list(Beta_hat = Beta_hat, prev_Beta_hat = Beta_hat))
# }

make_omega_hat_branch_fn <- function(omega_hat, 
                                     X_one_hot, 
                                     Y_one_hot, 
                                     X_h_one_hot,
                                     Jm1,
                                     p,
                                     max_retries,
                                     max_eval,
                                     maxtime) {
  
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
                                    max_retries,
                                    max_eval,
                                    maxtime)
    
    Beta_hat <- Beta_hat_result$Beta_hat
    prev_Beta_hat <- Beta_hat_result$prev_Beta_hat
    
    return(log_likelihood(Beta_hat, X_one_hot, Y_one_hot))
  }
}

get_omega_hat_branch_fn_max <- function(omega_hat_branch_fn, J) {
  
  opt_result <- optimize(omega_hat_branch_fn, 
                         interval = c(0, log(J)), 
                         maximum = TRUE,
                         tol = 0.1)
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

# Integrated Likelihood ---------------------------------------------------

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
                            max_retries,
                            max_eval,
                            maxtime) {
  
  psi_endpoints <- NULL

  while (is.null(psi_endpoints)) {
    
    omega_hat <- get_omega_hat(omega_hat_eq_con_fn, omega_hat_ineq_con_fn, J - 1, p, init_guess_sd)
    
    omega_hat_branch_fn <- make_omega_hat_branch_fn(omega_hat, 
                                                    X_one_hot, 
                                                    Y_one_hot, 
                                                    X_h_one_hot,
                                                    J - 1,
                                                    p,
                                                    max_retries,
                                                    max_eval,
                                                    maxtime)
    
    omega_hat_branch_fn_max <- get_omega_hat_branch_fn_max(omega_hat_branch_fn, J)
    
  psi_endpoints <- safe_get_psi_endpoints(omega_hat_branch_fn,
                                          omega_hat_branch_fn_max,
                                          delta,
                                          J)
  }
  
  return(list(omega_hat = omega_hat, 
              omega_hat_branch_fn_max = omega_hat_branch_fn_max,
              psi_endpoints = psi_endpoints))
}

generate_branch_specs <- function(data,
                                  h,
                                  init_guess_sd,
                                  alpha,
                                  num_workers,
                                  chunk_size,
                                  max_retries,
                                  max_eval,
                                  maxtime) {
  
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
                                   max_retries,
                                   max_eval,
                                   maxtime)
    
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
                       max_retries,
                       max_eval,
                       maxtime) {
  
  psi_grid <- get_psi_grid_1(branch_spec$psi_endpoints, step_size, J)
  
  init_guess <- rnorm(p * (J - 1), sd = init_guess_sd)
  
  log_L_tilde_vec <- numeric(length(psi_grid))
  
  prev_Beta_hat <- NULL 
  
  for (i in seq_along(psi_grid)) {
    
    psi <- psi_grid[i]
    
    Beta_hat_result <- get_Beta_hat(X_one_hot, 
                                    X_h_one_hot,
                                    branch_spec$omega_hat,
                                    psi,
                                    init_guess,
                                    J - 1,
                                    p,
                                    prev_Beta_hat,
                                    max_retries,
                                    max_eval,
                                    maxtime)
    
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
                              max_retries,
                              max_eval,
                              maxtime) {
  
  psi_grid <- get_psi_grid_2(branch_specs, step_size, J, quantiles)
  
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
      
      psi <- psi_grid[i]
      
      Beta_hat_result <- get_Beta_hat(X_one_hot, 
                                      X_h_one_hot,
                                      omega_hat,
                                      psi,
                                      init_guess,
                                      J - 1,
                                      p,
                                      prev_Beta_hat,
                                      max_retries,
                                      max_eval,
                                      maxtime)
      
      Beta_hat <- Beta_hat_result$Beta_hat
      prev_Beta_hat <- Beta_hat_result$prev_Beta_hat 
      
      log_L_tilde_vec[i] <- log_likelihood_rcpp(Beta_hat, X_one_hot, Y_one_hot, J - 1, p, n)
      
      init_guess <- 0.9 * c(Beta_hat) + 0.1 * init_guess
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
                                          max_retries,
                                          max_eval,
                                          maxtime) {
  
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
                                max_retries,
                                max_eval,
                                maxtime)
  
  plan(sequential)
  
  log_L_bar <- get_log_L_bar(branches)
  
  return(list(branches = branches,
              log_L_bar = log_L_bar))
}

# Profile Likelihood ------------------------------------------------------

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
                                       max_retries,
                                       max_eval,
                                       maxtime) {
  
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
                                                 max_retries,
                                                 max_eval,
                                                 maxtime)
  
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
                                    max_retries,
                                    max_eval,
                                    maxtime)
    
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
