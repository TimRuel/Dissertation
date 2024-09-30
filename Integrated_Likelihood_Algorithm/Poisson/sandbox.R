
dist <- runif

dist_params <- list(min = 0, max = 100)

omega_hat <- get_omega_hat(data, weights, dist, dist_params)

lambdas <- psi_grid |> 
  map_dbl(\(psi) get_lambda(50, psi, m, omega_hat, weights))



hist(lambdas)

weighted_sum(omega_hat, weights)

lambdas <- psi_grid |>
  purrr::accumulate(
    \(acc, nxt) get_lambda(acc, nxt, m, omega_hat, weights),
    .init = 50) |>
  magrittr::extract(-1)

plot(psi_grid, lambdas)

omega_hat_list <- data |> 
  get_omega_hat(weights, dist, dist_params, return_u = FALSE) |> 
  replicate(R, expr = _, simplify = FALSE)

plan(multisession, workers = 50)

lambda_mat <- foreach(
  
  omega_hat = omega_hat_list,
  .combine = "rbind",
  .multicombine = TRUE,
  .maxcombine = R,
  .options.future = list(seed = TRUE,
                         chunk.size = chunk_size)
  
) %dofuture% {
  
  psi_grid |>
    purrr::accumulate(
      \(acc, nxt) get_lambda(acc, nxt, m, omega_hat, weights),
      .init = sum(omega_hat)
      ) |> 
    magrittr::extract(-1) 
}

A <- 1520

B <- 1530

plot(psi_grid[A:B], colMeans(lambda_mat)[A:B])

hist(lambda_mat[,1524])

xx <- lambda_mat[,1524] |> 
  unname() |> 
  map_dbl(\(lambda) {
    
    get_theta_hat(lambda, m, omega_hat, weights) |> 
      log_likelihood(data)
    }) 

for (i in 1520:1530) {
  
  hist(code_mat[,i])
}

hist(xx)

x <- lambda_list |> 
  map_dbl(mean)

plot(psi_grid, x)

lambda_list <- lapply(seq_len(ncol(lambda_mat)), function(i) unname(lambda_mat)[,i])

y <- get_lambda(-1.75, 15.35, m, omega_hat_list[[17]], weights)

get_theta_hat(y, m, omega_hat_list[[17]], weights) |> weighted_sum(weights)


test <- function(lambda) get_theta_hat(lambda, m, omega_hat, weights)


curve(test(x), from = -1000, 1 , xlab="x", ylab="y")

# package Rsolnp




