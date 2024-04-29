source("utils.R")

data <- c(1, 1, 2, 4, 7, 10)

theta_MLE <- data / sum(data)

psi_MLE <- PoI_fn(theta_MLE)

step_size <- 0.01

psi_grid <- data |> 
  length() |> 
  log() |> 
  plyr::round_any(step_size, floor) |> 
  seq(0, to = _, step_size)

lower_psi_grid <- psi_grid[psi_grid < psi_MLE] |> 
  rev()

upper_psi_grid <- psi_grid[psi_grid >= psi_MLE]

R <- 250

u_list <- LaplacesDemon::rdirichlet(R, rep(1, length(data))) |>
  t() |>
  data.frame() |>
  as.list()

omega_hat_list <- u_list |>
  furrr::future_map(get_omega_hat, psi_MLE, .progress = TRUE)


get_theta_hat2 <- function(psi, omega_hat) {
  
  theta_hat <- nloptr::auglag(x0 = omega_hat,
                              fn = function(theta) -sum(omega_hat*log(theta), na.rm = TRUE),
                              heq = function(theta) c(sum(theta) - 1, PoI_fn(theta) - psi),
                              lower = rep(0, length(omega_hat)))$par
  
  return(theta_hat)
}


get_multinomial_entropy_values_PL2 <- function(data, psi_grid) {
  
  m <- length(data)
  
  theta_MLE <- data / sum(data)
  
  l_p <- psi_grid |> 
    furrr::future_map(get_theta_hat2, theta_MLE) |> 
    sapply(likelihood, data) |> 
    log() |>
    as.double()
  
  return(l_p)
}

plan(multisession, workers = availableCores())


s.time1 <- system.time({
  test <- get_multinomial_entropy_values_PL2(data, psi_grid)
})

s.time2 <- system.time({
  test2 <- get_multinomial_entropy_values_PL(data, psi_grid)
})


plan(sequential)

stime1 <- system.time({
  
  test1 <- sims[1:10] |>
    purrr::map(\(x) get_multinomial_entropy_values_PL(x, psi_grid),
               .progress = TRUE)
})

stime1


stime2 <- system.time({
  
  test2 <- sims[1:10] |>
    purrr::map(\(x) get_multinomial_entropy_values_PL2(x, psi_grid),
               .progress = TRUE)
})

stime2



