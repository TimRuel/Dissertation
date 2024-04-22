library(nloptr)

get_omega_hat <- function(u, psi_MLE) {
  
  omega_hat <- nloptr::auglag(x0 = u,
                              fn = function(omega) -sum(u*log(omega), na.rm = TRUE),
                              heq = function(omega) c(sum(omega) - 1, PoI_fn(omega) - psi_MLE),
                              lower = rep(0, length(u)),
                              localsolver = "LBFGS")$par
  
  return(omega_hat)
}

get_omega_hat2 <- function(u, psi_MLE) {
  
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

get_theta_hat <- function(psi, omega_hat) {
  
  theta_hat <- nloptr::auglag(x0 = omega_hat,
                              fn = function(theta) -sum(omega_hat*log(theta), na.rm = TRUE),
                              heq = function(theta) c(sum(theta) - 1, PoI_fn(theta) - psi),
                              lower = rep(0, length(omega_hat)),
                              localsolver = "LBFGS")$par
  
  return(theta_hat)
}

get_theta_hat2 <- function(psi, omega_hat) {
  
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

plan(sequential)

data <- c(1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 6, 8)

R <- 300

step_size <- 0.001

psi_grid <- data |> 
  length() |> 
  log() |> 
  plyr::round_any(step_size, ceiling) |> 
  seq(0, to = _, step_size)

u_list <- LaplacesDemon::rdirichlet(R, rep(1, length(data))) |> 
  t() |> 
  data.frame() |> 
  as.list()

PoI_fn <- function(x) -sum(x*log(x), na.rm = TRUE)

theta_MLE <- data / sum(data)

psi_MLE <- PoI_fn(theta_MLE)

omega_hat <- get_omega_hat2(u, psi_MLE)

plan(multisession, workers = availableCores())

s.time1 <- system.time({
  
  theta_hat1 <- psi_grid |> 
    furrr::future_map(get_theta_hat, omega_hat)
  
})

s.time1

s.time2 <- system.time({
  
  theta_hat2 <- psi_grid |> 
    furrr::future_map(get_theta_hat2, omega_hat)
  
})

s.time2







