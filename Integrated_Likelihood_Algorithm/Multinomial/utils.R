likelihood <- function(theta, data) {
  
  theta |> 
    Rmpfr::mpfr(256) |> 
    (`^`)(data) |> 
    prod()
}

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

get_multinomial_entropy_values_IL.aux <- function(omega_hat, data, psi_grid) {
  
  L <- psi_grid |>
    purrr::accumulate(\(acc, nxt) get_theta_hat(acc, nxt, omega_hat), .init = omega_hat) |>
    magrittr::extract(-1) |>
    purrr::map(likelihood, data)
  
  return(L)
}

get_multinomial_entropy_values_IL <- function(omega_hat_list, data, psi_grid) {
  
  L_bar <- omega_hat_list |>
    furrr::future_map(get_multinomial_entropy_values_IL.aux, 
                      data, 
                      psi_grid, 
                      .progress = TRUE) |> 
    unlist() |> 
    matrix(ncol = length(omega_hat_list)) |> 
    Rmpfr::mpfr2array(dim = c(length(psi_grid), length(omega_hat_list))) |> 
    apply(1, \(x) sum(x) / length(x))
  
  l_bar <- new("mpfr", unlist(L_bar)) |> 
    log()
  
  return(l_bar)
}

get_multinomial_entropy_values_PL <- function(data, theta_MLE, lower_psi_grid, upper_psi_grid) {
  
  lower_L_p <- lower_psi_grid |> 
    purrr::accumulate(\(acc, nxt) get_theta_hat(acc, nxt, theta_MLE), 
                      .init = theta_MLE) |> 
    magrittr::extract(-1) |> 
    purrr::map(likelihood, data) |> 
    rev()
  
  upper_L_p <- upper_psi_grid |> 
    purrr::accumulate(\(acc, nxt) get_theta_hat(acc, nxt, theta_MLE), 
                      .init = theta_MLE) |> 
    magrittr::extract(-1) |> 
    purrr::map(likelihood, data)
  
  L_p <- c(lower_L_p, upper_L_p)
  
  l_p <- new("mpfr", unlist(L_p)) |>
    log()
  
  return(l_p)
}






