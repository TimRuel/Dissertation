################################################################################
#################################### GENERAL ###################################
################################################################################

likelihood <- function(theta, data) prod(theta^data)

entropy <- function(theta) -sum(theta * log(theta), na.rm = TRUE)

make_objective_fn <- function(u, objective) function(omega) objective(u, omega)

get_omega_hat <- function(objective_fn, psi_MLE, init_guess) {
  
  gr <- function(omega) nloptr::nl.grad(omega, objective_fn)
  heq <- function(omega) c(sum(omega) - 1, entropy(omega) - psi_MLE)
  heqjac <- function(omega) nloptr::nl.jacobian(omega, heq)
  
  omega_hat <- nloptr::auglag(x0 = init_guess,
                              fn = objective_fn,
                              gr = gr,
                              heq = heq,
                              heqjac = heqjac,
                              lower = rep(0, length(init_guess)),
                              localtol = 1e-03,
                              localsolver = "LBFGS")$par
  
  return(omega_hat)
}

get_theta_hat <- function(init_guess, psi, omega_hat) {
  
  fn <- function(theta) -sum(omega_hat * log(theta), na.rm = TRUE)
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

get_multinomial_entropy_values_PL <- function(data, psi_grid) {
  
  theta_MLE <- data / sum(data)
  
  psi_MLE <- entropy(theta_MLE)
  
  psi_grid_list <- psi_grid |> 
    split(factor(psi_grid > psi_MLE)) |> 
    purrr::modify_in(1, rev) |> 
    unname()
  
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

get_multinomial_entropy_values_IL.aux <- function(omega_hat, data, psi_grid) {
  
  theta_MLE <- data / sum(data)
  
  psi_MLE <- entropy(theta_MLE)
  
  psi_grid_list <- psi_grid |> 
    split(factor(psi_grid > psi_MLE)) |> 
    purrr::modify_in(1, rev) |> 
    unname()
  
  L <- psi_grid_list |> 
    purrr::map(
      \(x) purrr::accumulate(
        x,
        \(acc, nxt) get_theta_hat(acc, nxt, omega_hat), 
        .init = omega_hat
      ) |> 
        magrittr::extract(-1) |> 
        purrr::map_dbl(likelihood, data)
    ) |> 
    purrr::modify_in(1, rev) 
  
  return(L)
}

get_multinomial_entropy_values_IL <- function(omega_hat_list, data, psi_grid) {
  
  l_bar <- omega_hat_list |>
    furrr::future_map(get_multinomial_entropy_values_IL.aux, 
                      data, 
                      psi_grid,
                      .progress = TRUE) |> 
    unlist() |> 
    matrix(ncol = length(psi_grid), byrow = TRUE) |> 
    colMeans() |> 
    log()
  
  return(l_bar)
}

################################################################################
######################## MODIFIED INTEGRATED LIKELIHOOD ########################
################################################################################

get_multinomial_entropy_values_modified_IL <- function(omega_hat_list, L, data, psi_grid) {
  
  L_tilde <- omega_hat_list |>
    furrr::future_map(get_multinomial_entropy_values_IL.aux, 
                      data, 
                      psi_grid,
                      .progress = TRUE) |> 
    unlist() |> 
    matrix(ncol = length(psi_grid), byrow = TRUE)
  
  l_bar <- L_tilde |> 
    (`/`)(L) |> 
    colMeans() |> 
    log()
  
  return(l_bar)
}



