likelihood <- function(theta, data) prod(theta^data)

PoI_fn <- function(x) -sum(x*log(x), na.rm = TRUE)

get_omega_hat <- function(u, psi_MLE) {
  
  fn <- function(omega) -sum(u*log(omega), na.rm = TRUE)
  gr <- function(omega) nloptr::nl.grad(omega, fn)
  heq <- function(omega) c(sum(omega) - 1, PoI_fn(omega) - psi_MLE)
  heqjac <- function(omega) nloptr::nl.jacobian(omega, heq)
  
  omega_hat <- nloptr::auglag(x0 = u,
                              fn = fn,
                              gr = gr,
                              heq = heq,
                              heqjac = heqjac,
                              lower = rep(0, length(u)),
                              localsolver = "LBFGS")$par
  
  return(omega_hat)
}

get_theta_hat <- function(init_guess, psi, omega_hat) {
  
  fn <- function(theta) -sum(omega_hat*log(theta), na.rm = TRUE)
  gr <- function(theta) nloptr::nl.grad(theta, fn)
  heq <- function(theta) c(sum(theta) - 1, PoI_fn(theta) - psi)
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

get_multinomial_entropy_values_IL.aux <- function(omega_hat, data, lower_psi_grid, upper_psi_grid) {
  
  lower_L <- lower_psi_grid |>
    purrr::accumulate(\(acc, nxt) get_theta_hat(acc, nxt, omega_hat), 
                      .init = omega_hat) |>
    magrittr::extract(-1) |>
    purrr::map(likelihood, data) |> 
    rev()
  
  upper_L <- upper_psi_grid |>
    purrr::accumulate(\(acc, nxt) get_theta_hat(acc, nxt, omega_hat), 
                      .init = omega_hat) |>
    magrittr::extract(-1) |>
    purrr::map(likelihood, data) 
  
  L <- c(lower_L, upper_L)
  
  return(L)
}

get_multinomial_entropy_values_IL <- function(omega_hat_list, data, lower_psi_grid, upper_psi_grid) {
  
  l_bar <- omega_hat_list |>
    furrr::future_map(get_multinomial_entropy_values_IL.aux, 
                      data, 
                      lower_psi_grid, 
                      upper_psi_grid,
                      .progress = TRUE) |> 
    unlist() |> 
    matrix(ncol = length(c(lower_psi_grid, upper_psi_grid)), byrow = TRUE) |> 
    colMeans() |> 
    log()
  
  return(l_bar)
}

get_multinomial_entropy_values_PL <- function(data, psi_grid) {
  
  theta_MLE <- data / sum(data)
  
  psi_MLE <- PoI_fn(theta_MLE)
  
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


