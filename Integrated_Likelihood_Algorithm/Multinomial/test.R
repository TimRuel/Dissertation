setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
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




get_multinomial_entropy_values_PL2 <- function(data, psi_grid) {
  
  theta_MLE <- data / sum(data)
  
  psi_MLE <- PoI_fn(theta_MLE)
  
  psi_grid_list <- psi_grid |> 
    split(factor(psi_grid > psi_MLE)) |> 
    purrr::modify_in(1, rev) |> 
    unname()
  
  l_p <- psi_grid_list |> 
    furrr::future_map(
      \(x) purrr::accumulate(
        x,
        \(acc, nxt) get_theta_hat(acc, nxt, theta_MLE), 
        .init = theta_MLE
        ) |> 
        magrittr::extract(-1) |> 
        purrr::map_dbl(likelihood, data),
      .progress = TRUE
      ) |> 
    purrr::modify_in(1, rev) |> 
    unlist() |> 
    log()
  
  return(l_p)
}

plan(sequential)

s.time1 <- system.time({
  test <- get_multinomial_entropy_values_PL(data, psi_grid)
})

s.time1

plan(multisession, workers = 6)

s.time2 <- system.time({
  test2 <- get_multinomial_entropy_values_PL2(data, psi_grid)
})

s.time2


plan(multisession, workers = 12)

stime1 <- system.time({
  
  test1 <- sims[1:20] |>
    furrr::future_map(\(x) get_multinomial_entropy_values_PL(x, psi_grid),
               .progress = TRUE)
})

stime1

plan(multisession, workers = 2)

stime2 <- system.time({
  
  test2 <- sims[1:20] |>
    purrr::map(\(x) get_multinomial_entropy_values_PL2(x, psi_grid),
               .progress = TRUE)
})

stime2


all.equal(test1, test2)
