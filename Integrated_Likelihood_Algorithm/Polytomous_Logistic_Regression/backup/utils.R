library(tidyverse)
library(pipeR)
library(Rcpp)
library(nloptr)
library(PolytomousUtils)
# sourceCpp("../../polytomous_utils.cpp")

# True Parameter Generation -----------------------------------------------

get_entropy_ranges <- function(num_classes, levels) {
  
  num_ranges <- length(levels)
  
  entropy_ranges <- seq(0.3, log(num_classes) - 0.2, length.out = num_ranges + 1) |> 
    (\(x) mapply(c, x[-length(x)], x[-1], SIMPLIFY = FALSE))() |> 
    purrr::map(\(x) {
      
      midpoint <- mean(x)
      desired_length <- (x[2] - x[1]) * 0.75
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

PoI_fn2 <- function(Beta, X_h_design) {

  X_h_design %*% cbind(0, Beta) |>
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

get_X_h_design <- function(X_design, h) {
  
  X1_levels <- get_X1_levels(X_design)
  
  X1_ref_level <- get_X1_ref_level(X_design)
  
  X1_main_effect_names <- paste0("X1", X1_levels)
  
  X1X2_interaction_names <- paste0("X1", setdiff(X1_levels, X1_ref_level)) |> 
    paste0(":X2")
  
  colnames(X_design) <- c(X1_main_effect_names, "X2", X1X2_interaction_names)
  
  rows_to_keep <- X_design[, paste0("X1", h)] == 1
  
  X_h_design <- X_design[rows_to_keep,]
  
  return(X_h_design)
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
  
  nnet::multinom(formula, data = data, trace = FALSE)
}

log_likelihood <- function(Beta, X_design, Y_design) {
  
  p <- ncol(X_design)
  
  Jm1 <- ncol(Y_design)
  
  Beta <- matrix(Beta,
                 nrow = p,
                 ncol = Jm1,
                 byrow = FALSE)
  
  Y_hat <- X_design %*% Beta
  
  sum(rowSums(Y_design * Y_hat) - matrixStats::rowLogSumExps(cbind(0, Y_hat)))
}

likelihood <- function(Beta, X_design, Y_design) exp(log_likelihood(Beta, X_design, Y_design))

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

E_log_like <- function(omega, X_design, Beta) {
  
  Beta <- Beta |> 
    matrix(nrow = p,
           ncol = Jm1,
           byrow = FALSE)
  
  logits <- X_design %*% Beta
  
  probs <- X_design %*% omega |> 
    apply(1, softmax_adj) |> 
    t()
  
  sum(rowSums(probs * logits) - log(1 + rowSums(exp(logits))))
}

safe_auglag <- purrr::possibly(nloptr::auglag)

get_Beta_hat <- function(obj_fn,
                         con_fn,
                         init_guess,
                         prev_Beta_hat = NULL,
                         maxtime) {

   if (!is.null(prev_Beta_hat)) {
   phi <- max(0.1, min(0.9, 1 - 1 / 10))
   init_guess <- 0.9 * prev_Beta_hat + 0.1 * init_guess
 }

 for (attempt in 1:10) {

   result <- safe_auglag(
     x0 = init_guess,
     fn = obj_fn,
     heq = con_fn,
     localsolver = "SLSQP",
     localtol = 1e-3,
     deprecatedBehavior = FALSE,
     control = list(on.error = "ignore",
                    maxtime = maxtime)
   )

   if (!is.null(result$par) && result$convergence > 0) {

     Beta_hat <- result$par
     return(list(Beta_hat = Beta_hat, prev_Beta_hat = Beta_hat))
   }

   init_guess <- init_guess + rnorm(length(init_guess), sd = 0.1)
 }

 result_fallback <- safe_auglag(
   x0 = init_guess,
   fn = obj_fn,
   heq = con_fn,
   localsolver = "COBYLA",
   localtol = 1e-3,
   deprecatedBehavior = FALSE,
   control = list(on.error = "ignore",
                  maxtime = maxtime + 10)
 )

 if (is.null(result_fallback$par)) {
   return(list(Beta_hat = prev_Beta_hat, prev_Beta_hat = prev_Beta_hat))
 }

 Beta_hat <- result_fallback$par
 return(list(Beta_hat = Beta_hat, prev_Beta_hat = Beta_hat))
}

# get_Beta_hat <- function(obj_fn,
#                          con_fn,
#                          init_guess,
#                          maxtime) {
#   
#   Beta_hat <- nloptr::auglag(
#     x0 = init_guess,
#     fn = obj_fn,
#     heq = con_fn,
#     localsolver = "SLSQP",
#     deprecatedBehavior = FALSE,
#     control = list(maxtime = maxtime))$par
#   
#   return(Beta_hat)
# }

get_psi_endpoints <- function(psi_max, Beta_max, X_h_design, num_std_errors, J, n_h) {
  
  theta_MLE <- X_h_design %*% cbind(0, Beta_max) |> 
    apply(1, softmax) |> 
    t() |> 
    colMeans()
  
  sigma <- theta_MLE*diag(J) - matrix(theta_MLE) %*% theta_MLE
  
  psi_max_SE <- sqrt(sum(matrix(1 + log(theta_MLE)) %*% (1 + log(theta_MLE)) * sigma, na.rm = TRUE) / n_h)
  
  psi_endpoints <- psi_max + c(-0.75, 1) * num_std_errors * psi_max_SE
}

get_psi_grid <- function(psi_endpoints, step_size, J) {
  
  lower <- 0
  
  upper <- log(J)
  
  left <- max(psi_endpoints[1], lower)  
  right <- min(psi_endpoints[2], upper)
  
  psi_grid <- seq(ceiling(left / step_size) * step_size, 
                  floor(right / step_size) * step_size, 
                  by = step_size)
  
  if (psi_endpoints[1] < lower) psi_grid <- c(lower, psi_grid[psi_grid > lower])
  if (psi_endpoints[2] > upper) psi_grid <- c(psi_grid[psi_grid < upper], upper)
  
  return(psi_grid)
}

make_plot <- function(df) {
  
  df |> 
    ggplot(aes(x = psi, y = .data[[names(df)[2]]])) + 
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
      y = paste(names(df)[[2]], "Log-Likelihood")
    )
}

# Integrated Likelihood ---------------------------------------------------

get_omega_hat <- function(omega_hat_eq_con_fn, omega_hat_ineq_con_fn, Jm1, p, init_guess_sd) {
  
  init_guess <- rnorm(p * Jm1, sd = init_guess_sd)
  
  result <- safe_auglag(x0 = init_guess,
                        fn = function(Beta) 0,
                        heq = omega_hat_eq_con_fn,
                        hin = omega_hat_ineq_con_fn,
                        localsolver = "SLSQP",
                        deprecatedBehavior = FALSE,
                        control = list(on.error = "ignore"))
  
  if (!is.null(result$par) && abs(omega_hat_eq_con_fn(result$par)) <= 0.1 && omega_hat_ineq_con_fn(result$par) <= 0) {
    
    omega_hat <- matrix(result$par, nrow = p, ncol = Jm1, byrow = FALSE)
    
    return(omega_hat)
  }
  
  else {
    
    return(get_omega_hat(omega_hat_eq_con_fn, omega_hat_ineq_con_fn, Jm1, p, init_guess_sd))
  }
}

make_omega_hat_branch_fn <- function(omega_hat, 
                                     X_design, 
                                     Y_design, 
                                     X_h_design,
                                     Jm1,
                                     p,
                                     n,
                                     maxtime) {
  
  obj_fn <- function(Beta) Beta_hat_obj_fn_rcpp(Beta, X_design, omega_hat, Jm1, p, n)
  
  function(psi) {
    
    con_fn <- function(Beta) Beta_hat_con_fn_rcpp(Beta, X_h_design, psi, Jm1, p)
    
    Beta_hat <- get_Beta_hat(obj_fn,
                             con_fn,
                             c(omega_hat),
                             NULL,
                             maxtime)
    
    return(log_likelihood_rcpp(Beta_hat, X_design, Y_design, Jm1, p, n))
  }
}

get_omega_hat_branch_arg_max <- function(omega_hat_branch_fn, J) {
  
  branch_argmax <- optimize(omega_hat_branch_fn, 
                            interval = c(0, log(J)), 
                            maximum = TRUE,
                            tol = 0.1)$maximum
  
  return(branch_argmax)
}

generate_branches <- function(X_design, 
                              Y_design,
                              X_h_design,
                              Jm1, 
                              p,
                              n,
                              psi_grid,
                              psi_hat,
                              threshold,
                              init_guess_sd,
                              num_workers,
                              chunk_size,
                              num_std_errors,
                              maxtime) {
  
  num_branches <- num_workers * chunk_size
  
  omega_hat_eq_con_fn <- function(Beta) omega_hat_eq_con_fn_rcpp(Beta, X_h_design, Jm1, p, psi_hat)
  
  omega_hat_ineq_con_fn <- function(Beta) omega_hat_ineq_con_fn_rcpp(Beta, X_design, Y_design, Jm1, p, n, threshold)
  
  omega_hat <- get_omega_hat(omega_hat_eq_con_fn, omega_hat_ineq_con_fn, Jm1, p, init_guess_sd)
  
  plan(multisession, workers = I(num_workers))
  
  result <- foreach(
    
    i = 1:num_branches,
    .combine = "list",
    .multicombine = TRUE,
    .maxcombine = num_branches,
    .errorhandling = "remove",
    .options.future = list(seed = TRUE,
                           chunk.size = chunk_size,
                           packages = c("PolytomousUtils", "nloptr"))
    
  ) %dofuture% {
    
    obj_fn <- function(Beta) Beta_hat_obj_fn_rcpp(Beta, X_design, omega_hat, Jm1, p, n)
    
    prev_Beta_hat <- NULL
    
    log_L_tilde_vec <- numeric(length(psi_grid))
    
    log_messages <- character(length(psi_grid))

    for (i in seq_along(psi_grid)) {
      
      psi <- psi_grid[i]
      
      con_fn <- function(Beta) Beta_hat_con_fn_rcpp(Beta, X_h_design, psi, Jm1, p)
      
      init_guess <- rnorm(p * Jm1, sd = init_guess_sd)
      
      Beta_hat_result <- get_Beta_hat(
        obj_fn,
        con_fn,
        init_guess,
        prev_Beta_hat,
        maxtime
      )
      
      Beta_hat <- Beta_hat_result$Beta_hat
      prev_Beta_hat <- Beta_hat_result$prev_Beta_hat
      
      # Evaluate constraint satisfaction
      constraint_val <- tryCatch(
        con_fn(Beta_hat),
        error = function(e) rep(NA_real_, Jm1 * p)
      )
      constraint_violation <- sum(abs(constraint_val), na.rm = TRUE)
      
      # Evaluate objective value (just in case)
      obj_val <- tryCatch(
        obj_fn(Beta_hat),
        error = function(e) NA_real_
      )
      
      # Compute likelihood
      logL <- tryCatch(
        log_likelihood_rcpp(Beta_hat, X_design, Y_design, Jm1, p, n),
        error = function(e) NA_real_
      )
      
      log_L_tilde_vec[i] <- logL
      
      # LOGGING: Save diagnostics for debugging
      log_messages[i] <- sprintf(
        "Branch %d | psi: %.4f | logL: %.4f | obj: %.4f | constraint_violation: %.4e",
        i, psi, logL, obj_val, constraint_violation
      )

      # psi <- psi_grid[i]
      # 
      # con_fn <- function(Beta) Beta_hat_con_fn_rcpp(Beta, X_h_design, psi, Jm1, p)
      # 
      # init_guess <- rnorm(p * Jm1, sd = init_guess_sd)
      # 
      # Beta_hat_result <- get_Beta_hat(obj_fn,
      #                                 con_fn,
      #                                 init_guess,
      #                                 prev_Beta_hat,
      #                                 maxtime)
      # 
      # Beta_hat <- Beta_hat_result$Beta_hat
      # 
      # prev_Beta_hat <- Beta_hat_result$prev_Beta_hat
      # 
      # log_L_tilde_vec[i] <- log_likelihood_rcpp(Beta_hat, X_design, Y_design, Jm1, p, n)
    }
    
    writeLines(log_messages, paste0("branch_log_", Sys.getpid(), ".txt"))
    
    log_L_tilde_df <- data.frame(psi = psi_grid, 
                                 Integrated = log_L_tilde_vec)
    
    list(omega_hat = omega_hat,
         log_L_tilde_df = log_L_tilde_df)
  }
  
  plan(sequential)
  
  omega_hat <- lapply(result, `[[`, 1)
  log_L_tilde_df <- lapply(result, `[[`, 2)
  
  list(omega_hat = omega_hat,
       log_L_tilde_df = log_L_tilde_df)
}

inflection_flag <- function(df) {
  
  spline_fit <- smooth.spline(df$psi, df$Integrated, spar = 0.5)
  
  spline_deriv <- predict(spline_fit, df$psi, deriv = 2)
  
  sign_changes <- which(diff(sign(spline_deriv$y)) != 0)
  
  mid_range <- round(length(df$psi) * c(0.25, 0.75))
  !any(sign_changes > mid_range[1] & sign_changes < mid_range[2])
}

outlier_flag <- function(x) {
  
  iqr <- IQR(x)
  q1 <- quantile(x, 0.25)
  q3 <- quantile(x, 0.75)
  between(x, q1 - 1.5 * iqr, q3 + 1.5 * iqr)
}

get_log_L_bar <- function(branches) {
  
  branch_filter1 <- branches$log_L_tilde_df |> 
    map_lgl(inflection_flag)
  
  branch_filter2 <- branches$log_L_tilde_df |> 
    map_dbl(\(df) max(df$Integrated)) |> 
    outlier_flag()
  
  branch_filter <- branch_filter1 & branch_filter2
  
  df_list <- branches$log_L_tilde_df[branch_filter] 
  
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

get_log_integrated_likelihood <- function(data,
                                          formula,
                                          h,
                                          step_size,
                                          num_std_errors,
                                          init_guess_sd,
                                          threshold,
                                          num_workers,
                                          chunk_size,
                                          maxtime) {
  
  model <- fit_multinomial_logistic_model(data, formula)
  
  Y_design <- data |> 
    pull(Y) |> 
    (\(Y) model.matrix(~ Y)[,-1])()
  
  X_design <- model.matrix(model)
  
  X_h_design <- get_X_h_design(X_design, h)
  
  psi_hat <- get_psi_hat(model, data, h)
  
  Beta_MLE <- get_Beta_MLE(model)
  
  Jm1 <- ncol(Y_design)
  
  J <- Jm1 + 1
  
  p <- ncol(X_design)
  
  n <- nrow(X_design)
  
  n_h <- nrow(X_h_design)
  
  psi_endpoints <- get_psi_endpoints(psi_hat, Beta_MLE, X_h_design, num_std_errors, J, n_h)
  
  psi_grid <- get_psi_grid(psi_endpoints, step_size, J)
  
  branches <- generate_branches(X_design, 
                                Y_design,
                                X_h_design,
                                Jm1, 
                                p,
                                n,
                                psi_grid,
                                psi_hat,
                                threshold,
                                init_guess_sd,
                                num_workers,
                                chunk_size,
                                num_std_errors,
                                maxtime)
  
  return(branches)
  
  # log_L_bar <- get_log_L_bar(branches)
  
  # print(make_plot(log_L_bar$df))
  
  # return(list(branches = branches,
  #             log_L_bar = log_L_bar))
}

# Profile Likelihood ------------------------------------------------------

num_std_errors_profile <- 3

get_log_profile_likelihood <- function(data,
                                       formula,
                                       h,
                                       step_size,
                                       num_std_errors,
                                       init_guess_sd,
                                       maxtime) {
  
  model <- fit_multinomial_logistic_model(data, formula)
  
  Y_design <- data |> 
    pull(Y) |> 
    (\(Y) model.matrix(~ Y)[,-1])()
  
  X_design <- model.matrix(model)
  
  X_h_design <- get_X_h_design(X_design, h)
  
  psi_hat <- get_psi_hat(model, data, h)
  
  Beta_MLE <- get_Beta_MLE(model)
  
  Jm1 <- ncol(Y_design)
  
  J <- Jm1 + 1
  
  p <- ncol(X_design)
  
  n <- nrow(X_design)
  
  n_h <- nrow(X_h_design)
  
  psi_endpoints <- get_psi_endpoints(psi_hat, Beta_MLE, X_h_design, num_std_errors, J, n_h)
  
  psi_grid <- get_psi_grid(psi_endpoints, step_size, J)
  
  obj_fn <- function(Beta) Beta_hat_obj_fn_rcpp(Beta, X_design, Beta_MLE, Jm1, p, n)
  
  prev_Beta_hat <- NULL
  
  log_L_p_vec <- numeric(length(psi_grid))
  
  for (i in seq_along(psi_grid)) {
    
    psi <- psi_grid[i]
    
    con_fn <- function(Beta) Beta_hat_con_fn_rcpp(Beta, X_h_design, psi, Jm1, p)
    
    init_guess <- rnorm(p * Jm1, sd = init_guess_sd)
    
    Beta_hat_result <- get_Beta_hat(obj_fn,
                                    con_fn,
                                    init_guess,
                                    prev_Beta_hat,
                                    maxtime)
    
    Beta_hat <- Beta_hat_result$Beta_hat
    
    prev_Beta_hat <- Beta_hat_result$prev_Beta_hat
    
    log_L_p_vec[i] <- log_likelihood_rcpp(Beta_hat, X_design, Y_design, Jm1, p, n)
  }
  
  log_L_p_df <- data.frame(psi = psi_grid, 
                           Profile = log_L_p_vec)
  
  # print(make_plot(log_L_p_df))
  
  return(log_L_p_df)
}
