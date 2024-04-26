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

# R <- 250

# u_list <- LaplacesDemon::rdirichlet(R, rep(1, length(data))) |> 
#   t() |> 
#   data.frame() |> 
#   as.list()

# omega_hat_list <- u_list |>
#   furrr::future_map(get_omega_hat, psi_MLE, .progress = TRUE)


get_multinomial_entropy_values_IL.aux <- function(omega_hat, data, lower_psi_grid, upper_psi_grid) {
  
  L <- psi_grid |>
    purrr::accumulate(\(acc, nxt) get_theta_hat(acc, nxt, omega_hat), .init = omega_hat) |>
    magrittr::extract(-1) |>
    purrr::map(likelihood, data)
  
  return(L)
}

multinomial_entropy_values_PL <- data |> 
  get_multinomial_entropy_values_PL(theta_MLE, lower_psi_grid, upper_psi_grid) 

data.frame(psi = psi_grid,
           loglikelihood = multinomial_entropy_values_PL |> as.numeric()) |> 
  ggplot() +
  geom_point(aes(x = psi, y = loglikelihood))





