
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


psi_grid |> 
  purrr::map(\(psi) get_theta_hat(psi, m, omega_hat, weights)) |> 
  purrr::map_dbl(\(theta) weighted_sum(theta, weights))

psi <- 1.7

objective <- function(lambda) abs(psi - sum(weights * ((m * omega_hat) / (m + lambda * weights))))

out <- nlm(objective, 0, fscale = 0, iterlim = 1000)

lambda <- out$estimate

lambda <- psi_grid |>
  purrr::accumulate(
    \(acc, nxt) {
      get_lambda(acc, nxt, m, omega_hat, weights)
    },
    .init = 7) |>
  magrittr::extract(-1)

l <- lambda |> 
  purrr::map_dbl(\(lambda) {
    lambda |>
      get_theta_hat(m, omega_hat, weights) |> 
      log_likelihood(data)
  })

plot(psi_grid, l)

R <- 2

omega_hat_list <- R |> 
  replicate(get_omega_hat(data, weights, dist, dist_params), simplify = FALSE)

l <- list()

tic()

for (i in 1:R) {
  
  l[[i]] <- psi_grid |>
    purrr::accumulate(
      \(acc, nxt) {
        lambda <- get_lambda(acc, nxt, m, omega_hat_list[[i]], weights)
        return(lambda)
      },
      .init = 0) |>
    magrittr::extract(-1) |>
    purrr::map_dbl(\(lambda) {
      lambda |>
        get_theta_hat(m, omega_hat_list[[i]], weights) |> 
        log_likelihood(data)
    })
  print(i)
}

toc()

integrated_log_likelihood_vals <- do.call(rbind, l) |>
  matrixStats::colLogSumExps() |>
  (`-`)(log(R))


library(doFuture)
plan(multisession)

y <- foreach(x = 1:4, y = 1:10) %dofuture% {
  z <- x + y
  slow_sqrt(z)
}