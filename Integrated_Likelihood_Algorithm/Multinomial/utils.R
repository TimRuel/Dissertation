################################################################################
#################################### GENERAL ###################################
################################################################################

log_likelihood <- function(theta, data) sum(data * log(theta), na.rm = TRUE)

neg_log_likelihood <- function(theta, data) -log_likelihood(theta, data)

likelihood <- function(theta, data) exp(log_likelihood(theta, data))

entropy <- function(theta) -sum(theta * log(theta), na.rm = TRUE)

euclidean_distance <- function(omega, u) sqrt(sum((omega - u)^2))

get_psi_grid <- function(data, step_size, num_std_errors, split = FALSE) {
  
  n <- sum(data)
  
  m <- length(data)
  
  theta_MLE <- data / n
  
  psi_MLE <- entropy(theta_MLE)
  
  sigma <- theta_MLE*diag(m) - matrix(theta_MLE) %*% theta_MLE
  
  psi_MLE_SE <- sqrt(sum(matrix(1 + log(theta_MLE)) %*% (1 + log(theta_MLE)) * sigma, na.rm = TRUE) / n)
  
  MoE <- num_std_errors * psi_MLE_SE
  
  psi_grid <- (psi_MLE + MoE*c(-1, 1)) |> 
    (\(x) c(max(0, x[1]), min(log(m), x[2])))() |> 
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

get_omega_hat <- function(objective_fn, psi_MLE, u_params, tol, return_u = FALSE) {
  
  u <- LaplacesDemon::rdirichlet(1, u_params)
  
  omega_hat <- nloptr::auglag(x0 = rep(1 / length(u), length(u)),
                              fn = function(omega) objective_fn(omega, u),
                              heq = function(omega) c(sum(omega) - 1, entropy(omega) - psi_MLE),
                              lower = rep(0, length(u)),
                              localsolver = "LBFGS")$par
  
  if (abs(entropy(omega_hat) - psi_MLE) < tol) {
    
    if (return_u) {
      
      return(list(u = u, omega_hat = omega_hat))
    }
    
    return(omega_hat)
  }
  
  else {
    
    get_omega_hat(objective_fn, psi_MLE, u_params, tol, return_u)
  }
}

get_theta_hat <- function(init_guess, psi, omega_hat) {
  
  fn <- function(theta) neg_log_likelihood(theta, omega_hat)
  gr <- function(theta) nloptr::nl.grad(theta, fn)
  heq <- function(theta) c(sum(theta) - 1, entropy(theta) - psi)
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

get_profile_log_likelihood <- function(data, psi_grid_list) {
  
  theta_MLE <- data / sum(data)
  
  l_p <- psi_grid_list |> 
    purrr::map(
      \(x) purrr::accumulate(
        x,
        \(acc, nxt) get_theta_hat(acc, nxt, theta_MLE), 
        .init = theta_MLE
      ) |> 
        magrittr::extract(-1) |> 
        purrr::map_dbl(likelihood, data)
    ) |> 
    purrr::modify_in(1, rev) |> 
    unlist() |> 
    log()
  
  return(l_p)
}

################################################################################
############################ INTEGRATED LIKELIHOOD ############################# 
################################################################################

get_integrated_log_likelihood <- function(omega_hat, data, psi_grid_list) {
  
  l1 <- psi_grid_list |> 
    purrr::map(
      \(psi_grid) psi_grid |> 
        purrr::accumulate(
          \(acc, nxt) get_theta_hat(acc, nxt, omega_hat), 
          .init = omega_hat
          ) |> 
        magrittr::extract(-1) |> 
        purrr::map_dbl(log_likelihood, data)
      ) |> 
    purrr::modify_in(1, rev) |> 
    unlist()
  
  return(l1)
}

get_multinomial_entropy_values_IL <- function(omega_hat_list, psi_grid_list, data) {
  
  l_bar <- omega_hat_list |>
    furrr::future_map(\(x) get_multinomial_entropy_values_IL.aux(x, data, psi_grid_list),
                      .progress = TRUE) |> 
    unlist() |> 
    matrix(ncol = psi_grid_list |> 
             unlist() |> 
             as.numeric() |> 
             length(), 
           byrow = TRUE) |> 
    colMeans() |> 
    log()
  
  return(l_bar)
}

################################################################################
######################## MODIFIED INTEGRATED LIKELIHOOD ########################
################################################################################

get_multinomial_entropy_values_modified_IL <- function(omega_hat_list, psi_grid_list, L, data) {
  
  L_tilde <- omega_hat_list |>
    furrr::future_map(\(x) get_multinomial_entropy_values_IL.aux(x, data, psi_grid_list),
                      .progress = TRUE) |> 
    unlist() |> 
    matrix(ncol = psi_grid_list |> 
             unlist() |> 
             as.numeric() |> 
             length(),
           byrow = TRUE)
  
  l_bar <- L_tilde |>
      (`/`)(L) |>
      colMeans() |>
      log()

  return(l_bar)
}



