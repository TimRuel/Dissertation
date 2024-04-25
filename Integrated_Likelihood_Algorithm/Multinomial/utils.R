likelihood <- function(theta, data) prod(theta^data)

PoI_fn <- function(x) -sum(x*log(x), na.rm = TRUE)

get_omega_hat <- function(u, psi_MLE) {
  
  fn <- function(omega) -sum(u*log(omega), na.rm = TRUE)
  gr <- function(omega) nl.grad(omega, fn)
  heq <- function(omega) c(sum(omega) - 1, PoI_fn(omega) - psi_MLE)
  heqjac <- function(omega) nl.jacobian(omega, heq)
  
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
  gr <- function(theta) nl.grad(theta, fn)
  heq <- function(theta) c(sum(theta) - 1, PoI_fn(theta) - psi)
  heqjac <- function(theta) nl.jacobian(theta, heq)
  
  theta_hat <- nloptr::auglag(x0 = init_guess,
                              fn = fn,
                              gr = gr,
                              heq = heq,
                              heqjac = heqjac,
                              lower = rep(0, length(omega_hat)),
                              localsolver = "LBFGS")$par
  
  return(theta_hat)
}

get_multinomial_entropy_values_IL.aux <- function(omega_hat, data, psi_grid) {
  
  L <- psi_grid |> 
    purrr::accumulate(\(acc, nxt) get_theta_hat(acc, nxt, omega_hat), .init = omega_hat) |> 
    magrittr::extract(-1) |> 
    purrr::map_dbl(likelihood, data)
  
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

# Have to change this function to accomodate new init_guess argument in get_theta_hat

get_multinomial_entropy_values_PL <- function(data, psi_grid) {
  
  theta_MLE <- data / sum(data)
  
  l_p <- psi_grid |> 
    purrr::map(get_theta_hat, theta_MLE) |> 
    purrr::map_dbl(likelihood, data) |> 
    log()
  
  return(l_p)
}







