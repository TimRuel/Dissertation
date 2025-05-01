get_num_predictors <- function(X1_level_names, J) {
  
  ncol(model.matrix( ~ factor(X1_level_names)[1] * J - 1))
}

get_entropy_ranges <- function(X1_level_names, J, entropy_range_specs) {
  
  list2env(entropy_range_specs, environment())

  num_ranges <- length(X1_level_names)

  lower <- 0 + offset[1]

  upper <- log(J) - offset[2]

  entropy_ranges <- seq(lower, upper, length.out = num_ranges + 1) |>
    (\(x) mapply(c, x[-length(x)], x[-1], SIMPLIFY = FALSE))() |>
    purrr::map(\(x) {
      midpoint <- mean(x)
      desired_length <- (x[2] - x[1]) * padding
      return(midpoint + c(-1, 1) * desired_length / 2)
    }
    )

  names(entropy_ranges) <- rev(X1_level_names)

  return(rev(entropy_ranges))
}

compute_probabilities <- function(intercept, slope, X2_vals) {

  eta <- matrix(intercept, nrow = length(X2_vals), ncol = length(intercept), byrow = TRUE) + outer(X2_vals, slope, `*`)
  exp_eta <- exp(eta)

  denom <- 1 + rowSums(exp_eta)
  probs <- cbind(1 / denom, exp_eta / denom)

  return(probs)
}

Beta_0_objective_fn <- function(params, X2_vals, Beta2) {

  Jm1 <- length(params) / 2
  intercept <- params[1:Jm1]
  slope <- params[(Jm1+1):(2*Jm1)] + Beta2

  prob_matrix <- compute_probabilities(intercept, slope, X2_vals)
  entropies <- apply(prob_matrix, 1, entropy)

  return(-diff(range(entropies)))
}

Beta_0_constraint_fn <- function(params, X2_vals, entropy_range, Beta2) {

  Jm1 <- length(params) / 2
  intercept <- params[1:Jm1]
  slope <- params[(Jm1+1):(2*Jm1)] + Beta2

  prob_matrix <- compute_probabilities(intercept, slope, X2_vals)
  entropies <- apply(prob_matrix, 1, entropy)

  return(c(entropy_range[1] - min(entropies), max(entropies) - entropy_range[2]))
}

optimize_Beta_0 <- function(J, entropy_range, X2_support, num_X2_test_vals, Beta2) {

  X2_vals <- seq(X2_support[1], X2_support[2], length.out = num_X2_test_vals)

  Jm1 <- J - 1
  init_params <- rnorm(2 * Jm1)

  result <- auglag(
    x0 = init_params,
    fn = function(params) Beta_0_objective_fn(params, X2_vals, Beta2),
    lower = rep(-10, length(init_params)),
    upper = rep(10, length(init_params)),
    hin = function(params) Beta_0_constraint_fn(params, X2_vals, entropy_range, Beta2),
    deprecatedBehavior = FALSE)

  params <- result$par

  return(list(intercept = params[1:Jm1], slope = params[(Jm1+1):(2*Jm1)]))
}

get_experiment_parameters <- function(X1_levels, true_param_specs) {
  
  list2env(true_param_specs, environment())

  X1_level_names <- names(X1_levels)

  X1_ref_level <- X1_levels |>
    purrr::map_lgl(\(x) x$ref_level) |>
    which() |>
    names()

  X1_main_effect_names <- paste0("X1", X1_level_names)
  
  X1X2_interaction_names <- paste0("X1", setdiff(X1_level_names, X1_ref_level)) |>
    paste0(":X2")

  Beta_0 <- matrix(NA,
                   nrow = p,
                   ncol = J - 1)
  
  rownames(Beta_0) <- c(X1_main_effect_names, "X2", X1X2_interaction_names)
  
  entropy_ranges <- get_entropy_ranges(X1_level_names, J, entropy_range_specs)

  for (h in X1_level_names) {

    list2env(X1_levels[[h]]$X2, environment())
    
    entropy_range <- entropy_ranges[[h]]
    
    X1_levels[[h]]$X2$entropy_range <- entropy_range

    if (h == X1_ref_level) {

      Beta2 <- 0

      params <- optimize_Beta_0(J, entropy_range, support, num_X2_test_vals, Beta2)

      Beta_0[paste0("X1", h), ] <- params$intercept
      Beta_0["X2", ] <- params$slope
    } else {

      Beta2 <- Beta_0["X2", ]

      params <- optimize_Beta_0(J, entropy_range, support, num_X2_test_vals, Beta2)

      Beta_0[grepl(h, rownames(Beta_0)), ] <- unname(rbind(params$intercept, params$slope))
    }
  }
  
  experiment_parameters <- list(Beta_0 = Beta_0,
                                X1_levels = X1_levels)

  return(experiment_parameters)
}