################################################################################
#################################### GENERAL ###################################
################################################################################

log_likelihood <- function(theta, data) {
  
  sum_data <- data |> 
    purrr::map_dbl(sum)
  
  n <- data |> 
    purrr::map_dbl(length)
  
  # theta <- theta |>
  #   mpfr(precBits = 256)
  
  l <- sum(sum_data * log(theta) - n * theta)
  
  return(l)
}

neg_log_likelihood <- function(theta, data) -log_likelihood(theta, data)

likelihood <- function(theta, data) exp(log_likelihood(theta, data))

E_log_like <- function(theta, omega) sum((log(theta) * omega - theta) * n, na.rm = TRUE)

weighted_sum <- function(theta, weights) sum(theta * weights, na.rm = TRUE)

euclidean_distance <- function(u, omega) dist(matrix(c(u, omega),
                                                     nrow = 2,
                                                     byrow = TRUE),
                                              method = "euclid")[1]

get_psi_grid <- function(data, weights, step_size, num_std_errors, split = FALSE) {
  
  n <- data |> 
    purrr::map_dbl(length)
  
  theta_MLE <- data |> 
    purrr::map_dbl(mean)
  
  psi_MLE <- weighted_sum(theta_MLE, weights)
  
  psi_MLE_SE <- data |>  
    sapply(var) |>  
    (\(v) v * weights^2 / n)() |> 
    sum() |>  
    sqrt()
  
  MoE <- num_std_errors * psi_MLE_SE
  
  psi_grid <- (psi_MLE + MoE*c(-1, 1)) |> 
    (\(x) c(max(0, x[1]), x[2]))() |> 
    plyr::round_any(step_size, floor) |> 
    (\(x) seq(x[1], x[2], step_size))()
  
  if (split) {
    
    psi_grid_list <- psi_grid |> 
      split(factor(psi_grid > psi_MLE)) |> 
      purrr::modify_in(1, rev) |> 
      unname()
    
    return(psi_grid_list)
  }
  
  return(psi_grid)
}

get_omega_hat_list <- function(objective_fn, psi_MLE, weights, u_params, R, tol) {
  
  u_list <- list()
  
  omega_hat_list <- list()
  
  while (length(omega_hat_list) < R) {
    
    u <- u_params |> 
      purrr::pmap_dbl(\(alpha, beta) rgamma(1, alpha, beta))
    
    omega_hat <- nloptr::auglag(x0 = u,
                                fn = function(omega) objective_fn(omega, u),
                                heq = function(omega) weighted_sum(omega, weights) - psi_MLE,
                                lower = rep(0, length(u)),
                                localsolver = "LBFGS")$par
    
    if (abs(weighted_sum(omega_hat, weights) - psi_MLE) < tol) {
      
      u_list <- c(u_list, list(u))
      
      omega_hat_list <- c(omega_hat_list, list(omega_hat))
      
      if (length(omega_hat_list) %% 10 == 0) print(paste0(length(omega_hat_list), " found!"))
    }
  }
  
  return(list("u" = u_list, "omega_hat" = omega_hat_list))
}

get_theta_hat <- function(init_guess, psi, omega_hat, weights) {
  
  fn <- function(theta) -E_log_like(theta, omega_hat)
  gr <- function(theta) nloptr::nl.grad(theta, fn)
  heq <- function(theta) weighted_sum(theta, weights) - psi
  heqjac <- function(theta) nloptr::nl.jacobian(theta, heq)
  
  theta_hat <- nloptr::auglag(x0 = init_guess,
                              fn = fn,
                              gr = gr,
                              heq = heq,
                              heqjac = heqjac,
                              lower = rep(0, length(omega_hat)),
                              localsolver = "LBFGS")$par
  
  return(theta_hat)
}

################################################################################
############################## PROFILE LIKELIHOOD ############################## 
################################################################################

get_poisson_weighted_sum_values_PL <- function(data, weights, psi_grid_list) {
  
  theta_MLE <- data |> 
    purrr::map_dbl(mean)
  
  l_p <- psi_grid_list |> 
    purrr::map(
      \(x) purrr::accumulate(
        x,
        \(acc, nxt) get_theta_hat(acc, nxt, theta_MLE, weights), 
        .init = theta_MLE
      ) |> 
        magrittr::extract(-1) |> 
        sapply(log_likelihood, data)
    ) |>
    purrr::modify_in(1, rev) |>
    unlist()
  
  return(l_p)
}

################################################################################
############################ INTEGRATED LIKELIHOOD ############################# 
################################################################################

get_poisson_weighted_sum_values_IL.aux <- function(omega_hat, data, weights, psi_grid_list) {
  
  l <- psi_grid_list |> 
    purrr::map(
      \(psi_grid) psi_grid |> 
        purrr::accumulate(
          \(acc, nxt) get_theta_hat(acc, nxt, omega_hat, weights), 
          .init = omega_hat
          ) |> 
        magrittr::extract(-1) |> 
        sapply(log_likelihood, data)
      ) |> 
    purrr::modify_in(1, rev) 
  
  return(l)
}

get_poisson_weighted_sum_values_IL <- function(omega_hat_list, data, weights, psi_grid_list) {
  
  l_bar <- omega_hat_list |>
    furrr::future_map(\(x) get_poisson_weighted_sum_values_IL.aux(x, data, weights, psi_grid_list),
                      .progress = TRUE) |> 
    unlist() |> 
    matrix(ncol = psi_grid_list |> 
             unlist() |> 
             as.numeric() |> 
             length(), 
           byrow = TRUE) |> 
    Rmpfr::mpfr(64) |> 
    exp() |> 
    colMeans() |>
    log() |> 
    as.numeric()
  
  return(l_bar)
}

################################################################################
######################## MODIFIED INTEGRATED LIKELIHOOD ########################
################################################################################

get_poisson_weighted_sum_values_modified_IL <- function(omega_hat_list, data, weights, psi_grid_list, l) {
  
  L_tilde <- omega_hat_list |>
    furrr::future_map(\(x) get_poisson_weighted_sum_values_IL.aux(x, data, weights, psi_grid_list),
                      .progress = TRUE) |> 
    unlist() |> 
    matrix(ncol = psi_grid_list |> 
             unlist() |> 
             as.numeric() |> 
             length(),
           byrow = TRUE)
  
  l_bar <- L_tilde |>
    (`-`)(l) |>
    Rmpfr::mpfr(64) |> 
    exp() |> 
    colMeans() |>
    log() |> 
    as.numeric()

  return(l_bar)
}



