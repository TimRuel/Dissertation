# scripts/helpers/profile_LL_utils.R

# Miscellaneous -----------------------------------------------------------

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

get_psi_hat_from_model <- function(model, X1_levels) {
  
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
    (\(vec) vec[h])() |> 
    unname()
}

get_psi_hat_from_MLE <- function(Beta_MLE, X_h_design) {
  
  (X_h_design %*% cbind(0, Beta_MLE)) |> 
    apply(1, softmax) |> 
    t() |> 
    colMeans() |> 
    entropy()
}

safe_auglag <- purrr::possibly(nloptr::auglag)

# Psi Grid ----------------------------------------------------------------

get_psi_endpoints <- function(psi_hat, Beta_MLE, X_h_design, num_std_errors, J, n_h) {
  
  theta_MLE <- X_h_design %*% cbind(0, Beta_MLE) |>
    apply(1, softmax) |>
    t() |>
    colMeans()
  
  sigma <- theta_MLE*diag(J) - matrix(theta_MLE) %*% theta_MLE
  
  psi_hat_SE <- sqrt(sum(matrix(1 + log(theta_MLE)) %*% (1 + log(theta_MLE)) * sigma, na.rm = TRUE) / n_h)
  
  psi_endpoints <- psi_hat + c(-1, 1) * num_std_errors * psi_hat_SE
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

# Beta_hat ----------------------------------------------------------------

get_Beta_hat_template <- function(obj_fn, con_fn, init_guess) {
  
  safe_auglag(
    x0 = init_guess,
    fn = obj_fn,
    heq = con_fn,
    localsolver = "SLSQP",
    deprecatedBehavior = FALSE)$par
}

# Branch Parameters ---------------------------------------------------------------

get_threshold <- function(Beta_MLE, X_design, Y_design, threshold_offset) {
  
  ceiling(abs(log_likelihood(Beta_MLE, X_design, Y_design))) + threshold_offset
}

get_omega_hat_branch_mode <- function(omega_hat,
                                      X_design,
                                      Y_design,
                                      X_h_design, 
                                      Jm1,
                                      p,
                                      n) {
  
  Beta_hat_obj_fn <- function(Beta) Beta_hat_obj_fn_rcpp(Beta, X_design, omega_hat, Jm1, p, n)
  
  omega_hat_branch <- function(psi) {
    
    Beta_hat_con_fn <- function(Beta) Beta_hat_con_fn_rcpp(Beta, X_h_design, psi, Jm1, p)
    
    Beta_hat <- get_Beta_hat_template(Beta_hat_obj_fn, Beta_hat_con_fn, c(omega_hat))
    
    return(log_likelihood_rcpp(Beta_hat, X_design, Y_design, Jm1, p, n))
  }
  
  optimize(omega_hat_branch,
           interval = c(0, log(Jm1 + 1)),
           maximum = TRUE,
           tol = 0.1)$maximum
}

get_branch_params <- function(X_design,
                              Y_design,
                              X_h_design,
                              threshold,
                              psi_hat,
                              psi_grid,
                              Jm1, 
                              p, 
                              n,
                              init_guess_sd) {
  
  omega_hat_obj_fn <- function(Beta) 0
  
  omega_hat_eq_con_fn <- function(Beta) omega_hat_eq_con_fn_rcpp(Beta, X_h_design, Jm1, p, psi_hat)
  
  omega_hat_ineq_con_fn <- function(Beta) omega_hat_ineq_con_fn_rcpp(Beta, X_design, Y_design, Jm1, p, n, threshold)
  
  num_discarded <- 0
  
  while (TRUE) {
    
    init_guess <- rnorm(p * Jm1, sd = init_guess_sd)
    
    omega_hat_candidate <- safe_auglag(
      x0 = init_guess,
      fn = omega_hat_obj_fn,
      heq = omega_hat_eq_con_fn,
      hin = omega_hat_ineq_con_fn,
      localsolver = "SLSQP",
      deprecatedBehavior = FALSE,
      control = list(on.error = "ignore")
      )$par
    
    if (!is.null(omega_hat_candidate) && abs(omega_hat_eq_con_fn(omega_hat_candidate)) <= 0.1 && omega_hat_ineq_con_fn(omega_hat_candidate) <= 0) {
      
      omega_hat_candidate <- matrix(omega_hat_candidate, nrow = p, ncol = Jm1, byrow = FALSE)
      
      omega_hat_candidate_branch_mode <- get_omega_hat_branch_mode(
        omega_hat = omega_hat_candidate,
        X_design = X_design,
        Y_design = Y_design,
        X_h_design = X_h_design, 
        Jm1 = Jm1,
        p = p,
        n = n
        )
      
      if (between(omega_hat_candidate_branch_mode, head(psi_grid, 1), tail(psi_grid, 1))) {
        
        Beta_hat_obj_fn <- function(Beta) Beta_hat_obj_fn_rcpp(Beta, X_design, omega_hat_candidate, Jm1, p, n)
        
        Beta_hat_con_fn <- function(Beta) Beta_hat_con_fn_rcpp(Beta, X_h_design, omega_hat_candidate_branch_mode, Jm1, p)
        
        Beta_hat <- get_Beta_hat_template(Beta_hat_obj_fn, Beta_hat_con_fn, c(omega_hat_candidate)) |> 
          matrix(nrow = p, ncol = Jm1, byrow = FALSE)
        
        return(list(omega_hat = omega_hat_candidate,
                    branch_mode = list(Beta_MLE = Beta_hat,
                                       psi = omega_hat_candidate_branch_mode),
                    num_discarded = num_discarded))
      }
    }
    num_discarded <- num_discarded + 1
  }
}

# Integrated Log-Likelihood ---------------------------------------------------

compute_IL_branch <- function(branch_mode,
                              psi_grid,
                              log_likelihood_fn,
                              get_Beta_hat,
                              con_fn_template) {
  
  psi_grid_left <- rev(psi_grid[psi_grid <= branch_mode$psi])
  psi_grid_right <- psi_grid[psi_grid > branch_mode$psi]
  
  log_L_tilde_vec_left <- numeric(length(psi_grid_left))
  log_L_tilde_vec_right <- numeric(length(psi_grid_right))
  
  init_guess <- c(branch_mode$Beta_MLE)
  
  for (i in seq_along(psi_grid_left)) {
    
    psi <- psi_grid_left[i]
    
    Beta_hat_con_fn <- function(Beta) con_fn_template(Beta, psi)
    
    Beta_hat <- get_Beta_hat(Beta_hat_con_fn, init_guess)
    
    log_L_tilde_vec_left[i] <- log_likelihood_fn(Beta_hat)
    
    init_guess <- Beta_hat
  }
  
  init_guess <- c(branch_mode$Beta_MLE)
  
  for (i in seq_along(psi_grid_right)) {
    
    psi <- psi_grid_right[i]
    
    Beta_hat_con_fn <- function(Beta) con_fn_template(Beta, psi)
    
    Beta_hat <- get_Beta_hat(Beta_hat_con_fn, init_guess)
    
    log_L_tilde_vec_right[i] <- log_likelihood_fn(Beta_hat)
    
    init_guess <- Beta_hat
  }
  
  log_L_tilde_vec <- c(rev(log_L_tilde_vec_left), log_L_tilde_vec_right)
  
  data.frame(psi = psi_grid, Integrated = log_L_tilde_vec)
}

get_log_L_bar <- function(IL_branches) {
  
  merged_df <- reduce(IL_branches, full_join, by = "psi")
  
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
                              chunk_size) {
  
  num_branches <- num_workers * chunk_size
  
  log_likelihood_fn <- function(Beta) log_likelihood_rcpp(Beta, X_design, Y_design, Jm1, p, n)
  
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
    
    branch_params <- get_branch_params(
      X_design = X_design,
      Y_design = Y_design,
      X_h_design = X_h_design, 
      threshold = threshold,
      psi_hat = psi_hat,
      psi_grid = psi_grid,
      Jm1 = Jm1, 
      p = p, 
      n = n,
      init_guess_sd = init_guess_sd
      )

    Beta_hat_obj_fn <- function(Beta) Beta_hat_obj_fn_rcpp(Beta, X_design, branch_params$omega_hat, Jm1, p, n)
    get_Beta_hat <- function(con_fn, init_guess) get_Beta_hat_template(Beta_hat_obj_fn, con_fn, init_guess)
    con_fn_template <- function(Beta, psi) Beta_hat_con_fn_rcpp(Beta, X_h_design, psi, Jm1, p)
    
    IL_branch <- compute_IL_branch(
      branch_mode = branch_params$branch_mode,
      psi_grid = psi_grid,
      log_likelihood_fn = log_likelihood_fn,
      get_Beta_hat = get_Beta_hat,
      con_fn_template = con_fn_template
      )
    
    list(IL_branch = IL_branch,
         branch_params = branch_params)
    }
  
  IL_branches <- lapply(result, `[[`, 1)
  branch_params <- lapply(result, `[[`, 2)
  
  log_L_bar <- get_log_L_bar(IL_branches)
  
  list(log_L_bar_df = log_L_bar$df,
       IL_branches = IL_branches,
       branches_matrix = log_L_bar$branches_matrix,
       branch_params = branch_params)
}

# Profile Log-Likelihood --------------------------------------------------

compute_profile_branch <- function(direction,
                                   step_size,
                                   psi_hat,
                                   Beta_MLE,
                                   PLL_max,
                                   stopping_val,
                                   log_likelihood_fn,
                                   get_Beta_hat,
                                   con_fn_template) {
  
  step <- if (direction == "left") -step_size else step_size
  psi <- if (direction == "left")
    floor(psi_hat / step_size) * step_size
  else
    ceiling(psi_hat / step_size) * step_size
  
  log_L_p <- PLL_max
  init_guess <- Beta_MLE
  
  psi_vals <- list()
  log_L_p_vec <- list()
  
  while (log_L_p > stopping_val) {
    
    Beta_hat_con_fn <- function(Beta) con_fn_template(Beta, psi)
    Beta_hat <- get_Beta_hat(Beta_hat_con_fn, init_guess)
    init_guess <- Beta_hat
    
    log_L_p <- log_likelihood_fn(Beta_hat)
    
    psi_vals[[length(psi_vals) + 1]] <- psi
    log_L_p_vec[[length(log_L_p_vec) + 1]] <- log_L_p
    
    psi <- psi + step
  }
  
  psi_vals <- unlist(psi_vals)
  log_L_p_vec <- unlist(log_L_p_vec)
  
  if (direction == "left") {
    psi_vals <- c(rev(psi_vals), psi_hat)
    log_L_p_vec <- c(rev(log_L_p_vec), PLL_max)
  }
  
  list(psi = psi_vals, Profile = log_L_p_vec)
}

get_profile_LL <- function(step_size, 
                           alpha,
                           psi_hat,
                           Beta_MLE, 
                           X_design,
                           Y_design,
                           X_h_design, 
                           Jm1,
                           p,
                           n) {
  
  log_likelihood_fn <- function(Beta) log_likelihood_rcpp(Beta, X_design, Y_design, Jm1, p, n)
  Beta_hat_obj_fn <- function(Beta) Beta_hat_obj_fn_rcpp(Beta, X_design, Beta_MLE, Jm1, p, n)
  get_Beta_hat <- function(con_fn, init_guess) get_Beta_hat_template(Beta_hat_obj_fn, con_fn, init_guess)
  con_fn_template <- function(Beta, psi) Beta_hat_con_fn_rcpp(Beta, X_h_design, psi, Jm1, p)
  
  PLL_max <- log_likelihood_fn(Beta_MLE)
  crit <- qchisq(1 - alpha, df = 1) / 2
  stopping_val <- PLL_max - crit
  
  result <- foreach(
    dir = c("left", "right"),
    .combine = "list",
    .multicombine = TRUE,
    .maxcombine = 2,
    .errorhandling = "remove",
    .options.future = list(
      seed = TRUE,
      chunk.size = 1,
      packages = c("PolytomousUtils", "nloptr")
    )
  ) %dofuture% {
    
    compute_profile_branch(
      direction = dir,
      step_size = step_size,
      psi_hat = psi_hat,
      Beta_MLE = Beta_MLE,
      PLL_max = PLL_max,
      stopping_val = stopping_val,
      log_likelihood_fn = log_likelihood_fn,
      get_Beta_hat = get_Beta_hat,
      con_fn_template = con_fn_template
      )
  }
  
  profile_LL <- do.call(rbind, lapply(result, function(entry) {
    data.frame(psi = entry$psi, Profile = entry$Profile)
  }))
  
  return(profile_LL)
}

# Experiment --------------------------------------------------------------

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
  psi_hat <- get_psi_hat_from_model(ml_model, X1_levels)
  n_h <- nrow(X_h_design)
  psi_endpoints <- get_psi_endpoints(psi_hat, Beta_MLE, X_h_design, IL$num_std_errors, J, n_h)
  psi_grid <- get_psi_grid(psi_endpoints, IL$step_size, J)
  
  # ---- Run integrated likelihood ----
  integrated_LL <- get_integrated_LL(
    X_design = X_design, 
    Y_design = Y_design, 
    X_h_design = X_h_design, 
    Jm1 = J - 1, 
    p = p, 
    n = n, 
    psi_grid = psi_grid, 
    psi_hat = psi_hat, 
    threshold = threshold, 
    init_guess_sd = IL$init_guess_sd, 
    num_workers = IL$num_workers, 
    chunk_size = IL$chunk_size
  )
  
  IL_plot <- get_LL_plot(integrated_LL$log_L_bar_df)
  IL_branches_plot <- get_branches_plot(integrated_LL$branches_matrix)
  
  # ---- Run profile likelihood ----
  
  profile_LL <- get_profile_LL(
    step_size = PL$step_size, 
    alpha = PL$alpha,
    psi_hat = psi_hat,
    Beta_MLE = Beta_MLE, 
    X_design = X_design,
    Y_design = Y_design,
    X_h_design = X_h_design, 
    Jm1 = Jm1,
    p = p,
    n = n)
  
  PL_plot <- get_LL_plot(profile_LL)
  
  # ---- Return results ----
  list(
    integrated_LL = integrated_LL,
    profile_LL = profile_LL,
    plots = list(IL_plot = IL_plot,
                 IL_branches_plot = IL_branches_plot,
                 PL_plot = PL_plot)
  )
}


