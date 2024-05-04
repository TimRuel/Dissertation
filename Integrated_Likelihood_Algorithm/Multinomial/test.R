setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")

set.seed(1996)

data <- c(1, 1, 2, 4, 7, 10)

theta_MLE <- data / sum(data)

psi_MLE <- entropy(theta_MLE)

alpha <- data + 1

step_size <- 0.01

psi_grid <- data |> 
  length() |> 
  log() |> 
  plyr::round_any(step_size, floor) |> 
  seq(0, to = _, step_size)

R <- 250

u_list <- LaplacesDemon::rdirichlet(R, alpha) |>
  t() |>
  data.frame() |>
  as.list()

omega_hat_list <- u_list |>
  purrr::map(get_omega_hat, psi_MLE, log_likelihood)

plan(multisession, workers = 12)

L_tilde <- omega_hat_list |>
  furrr::future_map(get_multinomial_entropy_values_IL.aux, 
                    data, 
                    psi_grid,
                    .progress = TRUE) |> 
  unlist() |> 
  matrix(ncol = length(psi_grid), byrow = TRUE) 

L <- u_list |> 
  purrr::map_dbl(likelihood, data) |> 
  unlist() |> 
  as.numeric()

L_tilde |> 
  (`/`)(L) |> 
  colMeans() |> 
  log()

L_tilde / L


