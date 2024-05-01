setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")

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

get_omega_hat(u_list[[1]], psi_MLE, distance)


