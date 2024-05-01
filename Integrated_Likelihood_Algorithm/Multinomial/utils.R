################################################################################
#################################### GENERAL ###################################
################################################################################

likelihood <- function(theta, data) prod(theta^data)

log_likelihood <- function(theta, data) sum(data * log(theta), na.rm = TRUE) 

distance <- function(u, omega) dist(matrix(c(u, omega), nrow = 2, byrow = TRUE))[1] 

entropy <- function(theta) -sum(theta * log(theta), na.rm = TRUE)

get_omega_hat <- function(u, psi_MLE, objective) {
  
  fn <- function(omega) objective(u, omega)
  gr <- function(omega) nloptr::nl.grad(omega, fn)
  heq <- function(omega) c(sum(omega) - 1, entropy(omega) - psi_MLE)
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
  
  fn <- function(theta) -log_likelihood(theta, omega_hat)
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

# E_log_like <- function(theta, omega) -sum(omega*log(theta), na.rm = TRUE)
# 
# get_omega_psi_hat <- function()
# 
# get_multinomial_entropy_values_IL.aux <- function(u, get_omega_psi_hat, data, psi_grid) {
#   
#   theta_MLE <- data / sum(data)
#   
#   psi_MLE <- entropy(theta_MLE)
#   
#   psi_grid_list <- psi_grid |> 
#     split(factor(psi_grid > psi_MLE)) |> 
#     purrr::modify_in(1, rev) |> 
#     unname()
#   
#   L <- psi_grid_list |> 
#     purrr::map(
#       \(x) purrr::accumulate(
#         x,
#         \(acc, nxt) get_theta_hat(acc, nxt, omega_hat), 
#         .init = omega_hat
#       ) |> 
#         magrittr::extract(-1) |> 
#         purrr::map_dbl(likelihood, data)
#     ) |> 
#     purrr::modify_in(1, rev) 
#   
#   return(L)
# }



