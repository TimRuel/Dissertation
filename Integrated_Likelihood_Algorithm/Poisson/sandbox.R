
source("Data Generation.R")

dist <- runif

dist_params <- list(min = 0, max = 100)

omega_hat <- get_omega_hat(data, weights, dist, dist_params)

m <- data |> 
  map_dbl(length)

psi_MLE <- data |>
  map_dbl(mean) |>
  weighted_sum(weights)

psi <- 12.5

fun <- function(lambda, m_i, omega_hat, alpha_i, psi) abs(psi - sum(alpha_i * ((m_i * omega_hat) / (m_i + lambda * alpha_i))))

fun1 <- function(x) fun(x, m_i = m, omega_hat = omega_hat, alpha_i = weights, psi = psi)

out <- nlm(fun1, 0)

lambda <- out$estimate

fun1(lambda)

theta_hat <- (m * omega_hat) / (m + lambda * weights)

get_theta <- function(data, weights, dist, dist_params, psi, return_u = FALSE) {
  
  u <- data |> 
    length() |> 
    c(dist_params) |> 
    do.call(dist, args = _)
  
  omega_hat <- u |> 
    (`/`)(sum(u * weights)) |> 
    (`*`)(psi)
  
  if (return_u) return(rbind(u, omega_hat))
  
  return(omega_hat)
}

Q <- function(theta) log_likelihood(theta, omega_hat)

theta1 <- get_theta(data, weights, dist, dist_params, psi)

Q(theta1)

weighted_sum(theta_hat, weights)

theta_hat <- get_theta_hat(m, psi, omega_hat, weights)



