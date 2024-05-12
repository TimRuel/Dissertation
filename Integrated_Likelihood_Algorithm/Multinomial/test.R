setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")

set.seed(1996)

data <- c(1, 1, 2, 4, 7, 10)

rmultinom(1, sum(data), (data + 1) / (6 + sum(data)))

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

u_list_mod_IL <- LaplacesDemon::rdirichlet(R, alpha) |>
  t() |>
  data.frame() |>
  as.list()

plan(multisession, workers = 12)

multinomial_entropy_values_modified_IL <- u_list_mod_IL |> 
  get_multinomial_entropy_values_modified_IL(u_list_mod_IL, data, psi_grid)

objective <- function(u, t) sum(u * log(t), na.rm = TRUE)

omega_hat_list_mod_IL <- u_list_mod_IL |>
  map(get_omega_hat, psi_MLE, objective)

plan(multisession, workers = 12)

multinomial_entropy_values_modified_IL1 <- omega_hat_list_mod_IL |> 
  get_multinomial_entropy_values_modified_IL(u_list_mod_IL, data, psi_grid)

plan(multisession, workers = availableCores())

u_list_IL <- LaplacesDemon::rdirichlet(R, rep(1, length(data))) |> 
  t() |> 
  data.frame() |> 
  as.list()

u <- LaplacesDemon::rdirichlet(1, rep(1, length(data)))

omega_hat_list_IL <- u_list_IL |>
  map(get_omega_hat, psi_MLE, log_likelihood)

omega_hat <- get_omega_hat(u, psi_MLE, log_likelihood)

multinomial_entropy_values_IL <- u |> list() |> rep(R) |> 
  get_multinomial_entropy_values_IL(data, psi_grid)

multinomial_entropy_values_IL1 <- omega_hat_list_IL |> 
  get_multinomial_entropy_values_IL(data, psi_grid)

plan(multisession, workers = 12)

multinomial_entropy_values_modified_IL2 <- u_list_IL |> 
  get_multinomial_entropy_values_modified_IL(u_list_IL, data, psi_grid)

multinomial_entropy_values_IL2 <- u_list_mod_IL |> 
  get_multinomial_entropy_values_IL(data, psi_grid)

data.frame(psi = psi_grid,
           loglikelihood = multinomial_entropy_values_IL2) |>
  ggplot(aes(x = psi, y = loglikelihood)) +
  geom_point()

data.frame(psi = psi_grid,
           loglikelihood = multinomial_entropy_values_IL) |>
  ggplot(aes(x = psi, y = loglikelihood)) +
  geom_point()

