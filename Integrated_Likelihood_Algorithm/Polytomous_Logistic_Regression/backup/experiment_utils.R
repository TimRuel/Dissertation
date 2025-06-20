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

# safe_auglag <- purrr::possibly(nloptr::auglag)
# 
# get_Beta_hat <- function(obj_fn,
#                          con_fn,
#                          init_guess,
#                          prev_Beta_hat = NULL,
#                          maxtime) {
# 
#  if (!is.null(prev_Beta_hat)) init_guess <- prev_Beta_hat
# 
#  for (attempt in 1:100) {
# 
#    result <- safe_auglag(
#      x0 = init_guess,
#      fn = obj_fn,
#      heq = con_fn,
#      localsolver = "SLSQP",
#      localtol = 1e-6,
#      deprecatedBehavior = FALSE,
#      control = list(on.error = "ignore",
#                     maxtime = maxtime)
#    )
# 
#    if (!is.null(result$par) && result$convergence == 0) {
# 
#      Beta_hat <- result$par
#      return(list(Beta_hat = Beta_hat, prev_Beta_hat = Beta_hat))
#    }
# 
#    init_guess <- init_guess + rnorm(length(init_guess), sd = 0.1)
#  }
# 
#  result_fallback <- safe_auglag(
#    x0 = init_guess,
#    fn = obj_fn,
#    heq = con_fn,
#    localsolver = "COBYLA",
#    localtol = 1e-6,
#    deprecatedBehavior = FALSE,
#    control = list(on.error = "ignore",
#                   maxtime = -1)
#  )
# 
#  if (is.null(result_fallback$par)) {
#    return(list(Beta_hat = prev_Beta_hat, prev_Beta_hat = prev_Beta_hat))
#  }
# 
#  Beta_hat <- result_fallback$par
#  return(list(Beta_hat = Beta_hat, prev_Beta_hat = Beta_hat))
# }

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
    
    omega_hat <- get_omega_hat(omega_hat_eq_con_fn, omega_hat_ineq_con_fn, Jm1, p, init_guess_sd)

    obj_fn <- function(Beta) Beta_hat_obj_fn_rcpp(Beta, X_design, omega_hat, Jm1, p, n)

    prev_Beta_hat <- NULL

    log_L_tilde_vec <- numeric(length(psi_grid))

    for (j in seq_along(psi_grid)) {

      psi <- psi_grid[j]

      con_fn <- function(Beta) Beta_hat_con_fn_rcpp(Beta, X_h_design, psi, Jm1, p)

      init_guess <- rnorm(p * Jm1, sd = init_guess_sd)

      Beta_hat_result <- get_Beta_hat(obj_fn,
                                      con_fn,
                                      init_guess,
                                      prev_Beta_hat,
                                      maxtime)

      Beta_hat <- Beta_hat_result$Beta_hat

      prev_Beta_hat <- Beta_hat_result$prev_Beta_hat

      log_L_tilde_vec[j] <- log_likelihood_rcpp(Beta_hat, X_design, Y_design, Jm1, p, n)
    }

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


