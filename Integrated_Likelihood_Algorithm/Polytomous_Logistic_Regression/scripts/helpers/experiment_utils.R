# scripts/helpers/experiment_utils.R

# General -----------------------------------------------------------------

# PoI_fn <- function(Beta, X_h_design) {
#
#   X_h_design %*% cbind(0, Beta) |>
#     apply(1, softmax) |>
#     t() |>
#     colMeans() |>
#     entropy()
# }

fit_multinomial_logistic_model <- function(data, formula) {

  nnet::multinom(formula, data, trace = FALSE)
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

  Jm1 <- length(model$lev) - 1

  model |>
    coef() |>
    unname() |>
    matrix(ncol = Jm1,
           byrow = TRUE)
}

get_psi_hat <- function(model, X1_levels) {
  
  X1_level_names <- names(X1_levels)
  
  m <- map_int(X1_levels, \(x) x$m)
  
  h <- get_X1_level_of_interest(X1_levels)

  model |>
    predict(type = "probs") |>
    as.data.frame() |>
    mutate(X1_level = rep(X1_level_names, times = m)) |>
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

get_threshold <- function(Beta_MLE, X_design, Y_design, threshold_offset) {

  ceiling(abs(log_likelihood(Beta_MLE, X_design, Y_design))) + threshold_offset
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
                              maxtime) {

  num_branches <- num_workers * chunk_size

  omega_hat_eq_con_fn <- function(Beta) omega_hat_eq_con_fn_rcpp(Beta, X_h_design, Jm1, p, psi_hat)

  omega_hat_ineq_con_fn <- function(Beta) omega_hat_ineq_con_fn_rcpp(Beta, X_design, Y_design, Jm1, p, n, threshold)

  omega_hat <- get_omega_hat(omega_hat_eq_con_fn, omega_hat_ineq_con_fn, Jm1, p, init_guess_sd)

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

    # log_messages <- character(length(psi_grid))

    for (i in seq_along(psi_grid)) {

      # psi <- psi_grid[i]
      # 
      # con_fn <- function(Beta) Beta_hat_con_fn_rcpp(Beta, X_h_design, psi, Jm1, p)
      # 
      # init_guess <- rnorm(p * Jm1, sd = init_guess_sd)
      # 
      # Beta_hat_result <- get_Beta_hat(
      #   obj_fn,
      #   con_fn,
      #   init_guess,
      #   prev_Beta_hat,
      #   maxtime
      # )
      # 
      # Beta_hat <- Beta_hat_result$Beta_hat
      # prev_Beta_hat <- Beta_hat_result$prev_Beta_hat
      # 
      # # Evaluate constraint satisfaction
      # constraint_val <- tryCatch(
      #   con_fn(Beta_hat),
      #   error = function(e) rep(NA_real_, Jm1 * p)
      # )
      # constraint_violation <- sum(abs(constraint_val), na.rm = TRUE)
      # 
      # # Evaluate objective value (just in case)
      # obj_val <- tryCatch(
      #   obj_fn(Beta_hat),
      #   error = function(e) NA_real_
      # )
      # 
      # # Compute likelihood
      # logL <- tryCatch(
      #   log_likelihood_rcpp(Beta_hat, X_design, Y_design, Jm1, p, n),
      #   error = function(e) NA_real_
      # )
      # 
      # log_L_tilde_vec[i] <- logL
      # 
      # # LOGGING: Save diagnostics for debugging
      # log_messages[i] <- sprintf(
      #   "Branch %d | psi: %.4f | logL: %.4f | obj: %.4f | constraint_violation: %.4e",
      #   i, psi, logL, obj_val, constraint_violation
      # )

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

      log_L_tilde_vec[i] <- log_likelihood_rcpp(Beta_hat, X_design, Y_design, Jm1, p, n)
    }

    # writeLines(log_messages, paste0("branch_log_", Sys.getpid(), ".txt"))

    log_L_tilde_df <- data.frame(psi = psi_grid,
                                 Integrated = log_L_tilde_vec)

    list(omega_hat = omega_hat,
         log_L_tilde_df = log_L_tilde_df)
  }

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

  # branch_filter1 <- branches$log_L_tilde_df |>
  #   map_lgl(inflection_flag)
  # 
  # branch_filter2 <- branches$log_L_tilde_df |>
  #   map_dbl(\(df) max(df$Integrated)) |>
  #   outlier_flag()
  # 
  # branch_filter <- branch_filter1 & branch_filter2
  # 
  # df_list <- branches$log_L_tilde_df[branch_filter]
  
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

get_integrated_LL <- function(X_design,
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
                              maxtime) {

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
                                maxtime)

  log_L_bar <- get_log_L_bar(branches)

  return(list(branches = branches,
              log_L_bar = log_L_bar))
}

# Profile Likelihood ------------------------------------------------------

get_profile_LL <- function(X_design,
                           Y_design,
                           X_h_design,
                           Beta_MLE,
                           Jm1,
                           p,
                           n,
                           psi_grid,
                           init_guess_sd,
                           maxtime) {

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

  return(log_L_p_df)
}

run_experiment <- function(config, X_design, model_df) {

  list2env(config$model_specs, environment())
  list2env(config$optimization_specs, environment())

  X1_levels <- config$X1_levels
  formula <- as.formula(formula)
  ml_model <- fit_multinomial_logistic_model(model_df, formula)
  Y_design <- get_Y_design(model_df)
  Beta_MLE <- get_Beta_MLE(ml_model)
  threshold <- get_threshold(Beta_MLE, X_design, Y_design, IL$threshold_offset)
  h <- get_X1_level_of_interest(X1_levels)
  X_h_design <- get_X_h_design(X_design, X1_levels)
  psi_hat <- get_psi_hat(ml_model, X1_levels)
  n_h <- nrow(X_h_design)
  psi_endpoints_IL <- get_psi_endpoints(psi_hat, Beta_MLE, X_h_design, IL$num_std_errors, J, n_h)
  psi_grid_IL <- get_psi_grid(psi_endpoints_IL, IL$step_size, J)

  # ---- Run integrated likelihood ----
  integrated_LL <- get_integrated_LL(
    X_design = X_design, 
    Y_design = Y_design, 
    X_h_design = X_h_design, 
    Jm1 = J - 1, 
    p = p, 
    n = n, 
    psi_grid = psi_grid_IL, 
    psi_hat = psi_hat, 
    threshold = threshold, 
    init_guess_sd = IL$init_guess_sd, 
    num_workers = IL$num_workers, 
    chunk_size = IL$chunk_size, 
    maxtime = IL$maxtime
  )
  
  IL_plot <- get_LL_plot(integrated_LL$log_L_bar$df)

  # ---- Run profile likelihood ----
  psi_endpoints_PL <- get_psi_endpoints(psi_hat, Beta_MLE, X_h_design, PL$num_std_errors, J, n_h)
  psi_grid_PL <- get_psi_grid(psi_endpoints_PL, PL$step_size, J)
  
  profile_LL <- get_profile_LL(
    X_design = X_design, 
    Y_design = Y_design, 
    X_h_design = X_h_design, 
    Beta_MLE,
    Jm1 = J - 1, 
    p = p, 
    n = n, 
    psi_grid = psi_grid_PL, 
    init_guess_sd = PL$init_guess_sd, 
    maxtime = PL$maxtime
  )
  
  PL_plot <- get_LL_plot(profile_LL)

  # ---- Return results ----
  list(
    integrated_LL = integrated_LL,
    profile_LL = profile_LL,
    plots = list(IL_plot = IL_plot,
                 PL_plot = PL_plot)
  )
}
