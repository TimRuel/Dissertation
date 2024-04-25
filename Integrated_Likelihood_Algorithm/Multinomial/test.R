source("utils.R")

plan(multisession, workers = availableCores())

data <- c(1, 1, 2, 4, 7, 10)

R <- 250

step_size <- 0.01

psi_grid <- data |> 
  length() |> 
  log() |> 
  plyr::round_any(step_size, ceiling) |> 
  seq(0, to = _, step_size)

u_list <- LaplacesDemon::rdirichlet(R, rep(1, length(data))) |> 
  t() |> 
  data.frame() |> 
  as.list()

theta_MLE <- data / sum(data)

psi_MLE <- PoI_fn(theta_MLE)

omega_hat_list <- u_list |>
  furrr::future_map(get_omega_hat, psi_MLE, .progress = TRUE)

omega_hat <- omega_hat_list[[1]]

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

psi_grid |> 
  purrr::accumulate(\(acc, nxt) get_theta_hat(acc, nxt, omega_hat), .init = omega_hat) |> 
  magrittr::extract(-1) |> 
  furrr::future_map_dbl(likelihood, data, .progress = TRUE)


get_multinomial_entropy_values_IL.aux <- function(omega_hat, data, psi_grid) {
  
  L <- psi_grid |> 
    purrr::accumulate(\(acc, nxt) get_theta_hat(acc, nxt, omega_hat), .init = omega_hat) |> 
    magrittr::extract(-1) |> 
    purrr::map_dbl(likelihood, data)
      
  return(L)
}

get_multinomial_entropy_values_IL.aux(omega_hat, data, psi_grid)

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

plan(multisession, workers = availableCores())

s.time1 <- system.time({
  
  multinomial_entropy_values_IL <- omega_hat_list |> 
    get_multinomial_entropy_values_IL(data, psi_grid)
  
})

s.time1

plan(list(tweak(multisession, workers = 6)), tweak(multisession, workers = 2))

s.time2 <- system.time({
  
  multinomial_entropy_values_IL <- omega_hat_list |> 
    get_multinomial_entropy_values_IL(data, psi_grid)
  
})

s.time2

get_theta_hat <- function(psi, omega_hat) {
  
  fn <- function(theta) -sum(omega_hat*log(theta), na.rm = TRUE)
  gr <- function(theta) nl.grad(theta, fn)
  heq <- function(theta) c(sum(theta) - 1, PoI_fn(theta) - psi)
  heqjac <- function(theta) nl.jacobian(theta, heq)
  
  theta_hat <- nloptr::auglag(x0 = omega_hat,
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
    purrr::map(get_theta_hat, omega_hat) |> 
    purrr::map_dbl(likelihood, data)
  
  return(L)
}

s.time1 <- system.time({
  
  get_multinomial_entropy_values_IL.aux(omega_hat, data, psi_grid)
  
})

s.time1


