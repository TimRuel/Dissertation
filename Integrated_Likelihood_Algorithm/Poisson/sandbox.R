
dist <- runif

dist_params <- list(min = 0, max = 100)

omega_hat <- get_omega_hat(data, weights, dist, dist_params, return_u = TRUE)

psi_MLE <- data |>
  map_dbl(mean) |>
  weighted_sum(weights)

u <- omega_hat["u",]

u - ((weighted_sum(u, weights) - psi_MLE) / sum(weights^2)) * weights


weighted_sum(a, weights)



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

m <- data |> 
  map_dbl(length)

psi <- 11.5

objective <- function(lambda) {
  
  theta_hat <- get_theta_hat(lambda, m, omega_hat, weights)
  
  S <- weighted_sum(theta_hat, weights)
  
  return(distance(S, psi))
}

min(m / weights)

objective(0.1)

plot(objective(1:1000), type = "l")

h <- Vectorize(objective); curve(h, -1, 1)

curve(h, from = -0.619, to = 1, xlab="x", ylab="y")

ggplot(data.frame(lambda=c(0, 1)), aes(x=lambda)) + 
  stat_function(fun=objective)

lambda <- 0.5

lambda |> 
  get_theta_hat(m, omega_hat, weights) |> 
  weighted_sum(weights) |> 
  distance(psi)

abline(a = 0.6)

objective <- function(theta) neg_log_likelihood(theta, omega_hat)

psi <- 11.5

h <- function(theta) weighted_sum(theta, weights)

init_guess <- runif(n, 0, 5)

tic()

out <- Rsolnp::solnp(pars = init_guess, 
                     fun = objective, 
                     eqfun = h,
                     eqB = psi,
                     LB = rep(0, n))

toc()

weighted_sum(out$pars, weights)

log_likelihood(out$pars, omega_hat)

log_likelihood(runif(n, 0, 0.5), omega_hat)


omega_hat <- get_omega_hat(data, weights, dist, dist_params)

tic()

theta_hat <- psi_grid |>
  purrr::accumulate(
    \(acc, nxt) get_theta_hat(acc, nxt, m, omega_hat, weights),
    .init = omega_hat)

toc()

theta_mle <- data |> 
  map_dbl(mean)

psi_grid |> 
  purrr::accumulate(
    \(acc, nxt) {
      print(nxt)
      get_theta_hat(acc, nxt, m, theta_mle, weights)
      },
    .init = runif(n, 0, 1))


omega_hat <- get_omega_hat(data, weights, dist, dist_params)

get_lambda <- function(init_guess, psi, m, omega_hat, weights) {
  
  objective <- function(lambda) {
    
    lambda |> 
      get_theta_hat(m, omega_hat, weights) |> 
      weighted_sum(weights) |> 
      distance(psi)
  }
  
  out <- nlm(f = objective, 
             p = init_guess,
             fscale = 0,
             iterlim = 100000)
  
  lambda <- out$estimate
  
  print(out$code)
  
  return(lambda)
}

tic()

test <- psi_grid |>
  purrr::accumulate(
    \(acc, nxt) get_lambda(acc, nxt, m, omega_hat, weights),
    .init = 0) |>
  magrittr::extract(-1) |>
  purrr::map(\(lambda) {
    lambda |>
      get_theta_hat(m, omega_hat, weights)
  })

toc()

omega_hat <- get_omega_hat(data, weights, dist, dist_params)

psi <- 11.5

f <- function(lambda) {
  
  theta_hat <- get_theta_hat(lambda, m, omega_hat, weights)
  
  S <- weighted_sum(theta_hat, weights)
  
  return(distance(S, psi))
}

f.gr <- function(lambda) nloptr::nl.grad(lambda, f)

tic()

lambda <- nloptr::auglag(x0 = 0,
                         fn = f,
                         gr = f.gr,
                         lower = -min(m / weights),
                         localsolver = "LBFGS")$par

toc()

theta_hat <- get_theta_hat(lambda, m, omega_hat, weights)

S <- weighted_sum(theta_hat, weights)

log_likelihood(theta_hat, weights)



get_theta_hat1 <- function(init_guess, psi, omega_hat, weights) {

  f <- function(theta) neg_log_likelihood(theta, omega_hat)
  f.gr <- function(theta) nloptr::nl.grad(theta, f)
  fcon <- function(theta) weighted_sum(theta, weights) - psi
  fcon.jac <- function(theta) nloptr::nl.jacobian(theta, fcon)

  theta_hat <- nloptr::auglag(x0 = init_guess,
                              fn = f,
                              gr = f.gr,
                              heq = fcon,
                              heqjac = fcon.jac,
                              lower = rep(1e-10, length(omega_hat)),
                              localsolver = "LBFGS")$par

  return(theta_hat)
}

tic()

theta_hat1 <- get_theta_hat1(data |> map_dbl(mean) + 0.5, psi, omega_hat, weights)

toc()

S <- weighted_sum(theta_hat1, weights)

log_likelihood(theta_hat1, weights)


omega_hat <- get_omega_hat(data, weights, dist, dist_params)

theta_hat_list <- psi_grid |>
  purrr::accumulate(
    \(acc, nxt) get_lambda(acc, nxt, m, omega_hat, weights),
    .init = 0) |>
  magrittr::extract(-1) |>
  purrr::map(\(lambda) {
    lambda |>
      get_theta_hat(m, omega_hat, weights)
  })

theta_hat_list |> map_dbl(\(theta_hat) weighted_sum(theta_hat, weights))

omega_hat_list <- data |> 
  get_omega_hat(weights, dist, dist_params) |> 
  replicate(R, expr = _, simplify = FALSE)

omega_hat_list |> map_dbl(weighted_sum, weights)

plan(sequential)

set.seed(seed)

alpha_prior <- 10

beta_prior <- 10

alpha_posterior <- data |> 
  map_dbl(sum) |>
  (`+`)(alpha_prior)

beta_posterior <- data |> 
  map_dbl(length) |>
  (`+`)(beta_prior)

dist <- rgamma

dist_params <- list(shape = alpha_posterior, rate = beta_posterior)

chunk_size <- 5

tic()

plan(multisession, workers = 50)

mod_integrated_log_likelihood_vals <- get_mod_integrated_log_likelihood_vals(data,
                                                                             weights,
                                                                             step_size, 
                                                                             num_std_errors, 
                                                                             dist,
                                                                             dist_params, 
                                                                             R, 
                                                                             chunk_size)

toc()

plot(psi_grid, mod_integrated_log_likelihood_vals_1$variance)

data <- tibble(psi = psi_grid,
           variance = L_ratio |> apply(2, \(x) var(x) * 249 /250),
           cv = cv,
           index = seq_along(psi_grid)) 
  

data |> 
  # filter(index %% 32 == 0) |> 
  ggplot() +
  geom_point(aes(x = psi, y = variance)) 
  # coord_cartesian(xlim = c(0.25, 0.3))

data |> 
  ggplot() +
  geom_point(aes(x = psi, y = cv))


L_ratio <- mod_integrated_log_likelihood_vals_1$L_ratio

L_ratio |> apply(2, \(x) var(x) * 249 /250)

x <- L_ratio[,2] |> mean()

mean((L_ratio[,2] - x)^2)

var(L_ratio[,2]) * 249/250

mod_integrated_log_likelihood_vals_1$variance

l_diff <- log(L_ratio)

log_I_hat <- l_diff |>
  matrixStats::colLogSumExps() |>
  (`-`)(log(R))

L_ratio |> colMeans()

I_hat <- log_I_hat

plan(sequential)

set.seed(seed)

dist <- runif

dist_params <- list(min = 0, max = 100)

chunk_size <- 5

tic()

plan(multisession, workers = 50)

integrated_log_likelihood_vals <- get_integrated_log_likelihood_vals(data,
                                                                     weights,
                                                                     step_size,
                                                                     num_std_errors,
                                                                     dist,
                                                                     dist_params,
                                                                     R,
                                                                     chunk_size)

toc()

data1 <- tibble(psi = psi_grid,
               variance = integrated_log_likelihood_vals$s_squared,
               cv = integrated_log_likelihood_vals$cv,
               index = seq_along(psi_grid)) 


data1 |> 
  # filter(index %% 32 == 0) |> 
  ggplot() +
  geom_point(aes(x = psi, y = variance)) 
# coord_cartesian(xlim = c(0.25, 0.3))

data1 |> 
  ggplot() +
  geom_point(aes(x = psi, y = cv))

w_star_sum <- mod_integrated_log_likelihood_vals$L_ratio |> 
  sweep(2, colSums(mod_integrated_log_likelihood_vals$L_ratio), FUN = '/') |>
  apply(2, \(x) x^2) |> 
  colSums()

data2 <- tibble(psi = psi_grid,
                variance = mod_integrated_log_likelihood_vals$s_squared,
                cv = mod_integrated_log_likelihood_vals$cv,
                w_star_sum = w_star_sum,
                index = seq_along(psi_grid)) 


data2 |> 
  # filter(index %% 32 == 0) |> 
  ggplot() +
  geom_point(aes(x = psi, y = variance)) 
# coord_cartesian(xlim = c(0.25, 0.3))

data2 |> 
  ggplot() +
  geom_point(aes(x = psi, y = cv))

data2 |> 
  ggplot() +
  geom_point(aes(x = psi, y = w_star_sum))

#############################################################################################################################
#############################################################################################################################
#############################################################################################################################

omega_hat_mat <- get_omega_hat(data, weights, dist, dist_params, return_u = TRUE)
u <- omega_hat_mat["u",]
omega_hat <- omega_hat_mat["omega_hat",]
log_like_u <- log_likelihood(u, data)

init_guess <- 0

psi_grid |>
  purrr::accumulate(
    \(acc, nxt) get_lambda(acc, nxt, m, omega_hat, weights),
    .init = init_guess) |>
  magrittr::extract(-1) |>
  purrr::map_dbl(\(lambda) {
    lambda |>
      get_theta_hat(m, omega_hat, weights) |> 
      log_likelihood(data) |> 
      (`-`)(log_like_u)
  })

omega_hat_list <- data |> 
  get_omega_hat(weights, dist, dist_params, return_u = TRUE) |> 
  replicate(R, expr = _, simplify = FALSE)

likelihood1 <- function(theta, data) exp(log_likelihood(theta, data)) #* 0.2

vlikelihood <- Vectorize(likelihood1, vectorize.args="theta")

beta <- 0.5

alpha <- data[[3]] * beta + 1

f <- function(x, shape, rate) dgamma(x, shape = shape, rate = rate) * gamma(0.3)

vf <- Vectorize(f, vectorize.args="x")

ggplot() +
  stat_function(fun = vlikelihood,
                args = list(data = data[[3]]),
                xlim = c(0, 15),
                color = "red") +
  stat_function(fun = vf, 
                args = list(shape = alpha, rate = beta),
                xlim = c(0, 15), 
                color = "blue") +
  geom_vline(xintercept = 2)


f1 <- function(x, shape, rate) dgamma(x, shape = shape, rate = rate) 

vf1 <- Vectorize(f1, vectorize.args="x")


ggplot() +
  stat_function(fun = vf1,
                args = list(shape = 2, rate = 2),
                xlim = c(0, 30))







