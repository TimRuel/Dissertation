likelihood <- function(theta, data) prod(theta^data)

PoI_fn <- function(x) -sum(x*log(x), na.rm = TRUE)

get_omega_hat <- function(u, psi_MLE) {
  
  omega_hat <- nloptr::auglag(x0 = u,
                              fn = function(omega) -sum(u*log(omega), na.rm = TRUE),
                              heq = function(omega) c(sum(omega) - 1, PoI_fn(omega) - psi_MLE),
                              lower = rep(0, length(u)))$par
  
  return(omega_hat)
}

get_theta_hat <- function(psi, omega_hat) {
  
  theta_hat <- nloptr::auglag(x0 = omega_hat,
                              fn = function(theta) -sum(omega_hat*log(theta), na.rm = TRUE),
                              heq = function(theta) c(sum(theta) - 1, PoI_fn(theta) - psi),
                              lower = rep(0, length(omega_hat)))$par
  
  return(theta_hat)
}

get_multinomial_entropy_values_IL.aux <- function(omega_hat, data, psi_grid) {
  
  L <- psi_grid |> 
    purrr::map(get_theta_hat, omega_hat) |> 
    purrr::map_dbl(likelihood, data)
  
  return(L)
}

get_multinomial_entropy_values_IL <- function(omega_hat_list, data, psi_grid) {
  
  l_bar <- omega_hat_list |>
    future_map(get_multinomial_entropy_values_IL.aux, 
               data, 
               psi_grid, 
               .progress = TRUE) |> 
    unlist() |> 
    matrix(ncol = length(psi_grid), byrow = TRUE) |> 
    colMeans() |> 
    log()
  
  return(l_bar)
}

get_multinomial_entropy_values_PL <- function(data, psi_grid) {
  
  theta_MLE <- data / sum(data)
  
  l_p <- psi_grid |> 
    purrr::map(get_theta_hat, theta_MLE) |> 
    purrr::map_dbl(likelihood, data) |> 
    log()
  
  return(l_p)
}







