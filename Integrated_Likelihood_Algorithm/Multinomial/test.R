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

likelihood <- function(theta, data) prod(theta^data)

get_multinomial_entropy_values_IL.aux <- function(omega_hat, data, lower_psi_grid, upper_psi_grid) {
  
  lower_L <- lower_psi_grid |>
    purrr::accumulate(\(acc, nxt) get_theta_hat(acc, nxt, omega_hat), 
                      .init = omega_hat) |>
    magrittr::extract(-1) |>
    purrr::map(likelihood, data) |> 
    rev()
  
  upper_L <- upper_psi_grid |>
    purrr::accumulate(\(acc, nxt) get_theta_hat(acc, nxt, omega_hat), 
                      .init = omega_hat) |>
    magrittr::extract(-1) |>
    purrr::map(likelihood, data) 
  
  L <- c(lower_L, upper_L)
  
  return(L)
}

get_multinomial_entropy_values_IL <- function(omega_hat_list, data, lower_psi_grid, upper_psi_grid) {
  
  l_bar <- omega_hat_list |>
    furrr::future_map(get_multinomial_entropy_values_IL.aux, 
                      data, 
                      lower_psi_grid, 
                      upper_psi_grid,
                      .progress = TRUE) |> 
    unlist() |> 
    matrix(ncol = length(c(lower_psi_grid, upper_psi_grid)), byrow = TRUE) |> 
    colMeans() |> 
    log()
  
  return(l_bar)
}

plan(multisession, workers = availableCores())

multinomial_entropy_values_IL <- omega_hat_list |> 
  get_multinomial_entropy_values_IL(data, lower_psi_grid, upper_psi_grid)

data.frame(psi = psi_grid,
           loglikelihood = multinomial_entropy_values_IL |> as.numeric()) |> 
  ggplot() +
  geom_point(aes(x = psi, y = loglikelihood))





