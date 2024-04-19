library(parallelly)
library(furrr)

plan(list(tweak(multisession, workers = 4)), tweak(multisession, workers = 3))

data <- c(1, 1, 2, 4, 7, 10)

step_size <- 0.01

psi_grid <- data |> 
  length() |> 
  log() |> 
  plyr::round_any(step_size, ceiling) |> 
  seq(0, to = _, step_size)

R <- 250

n_sims <- 1000

sims <- n_sims |> 
  rmultinom(length(data), data) |> 
  data.frame() |> 
  as.list()

u_list <- n_sims |> 
  replicate({
    R |> 
      LaplacesDemon::rdirichlet(rep(1, length(data)))|> 
      t() |> 
      data.frame() |> 
      as.list()
    }, 
    simplify = FALSE)

omega_hat_lists <- u_list |> 
  purrr::map2(sims, 
              \(x, y) x |> 
                purrr::map(\(z) z |> 
                             get_omega_hat(PoI_fn(y / sum(y)))), 
              .progress = TRUE)

# saveRDS(omega_hat_lists, "omega_hat_lists.Rda")
omega_hat_lists <- readRDS("omega_hat_lists.Rda")


stime <- system.time({
  
  result <- omega_hat_lists |>
    future_map2(sims,
                \(x, y) get_multinomial_entropy_values_IL(x, y, psi_grid),
                .progress = TRUE)
})

stime

mods <- result |> 
  data.frame() |> 
  mutate(psi = psi_grid) |> 
  pivot_longer(cols = -psi,
               names_to = "Iteration",
               values_to = "loglikelihood") |> 
  group_by(Iteration) |> 
  group_map(~ smooth.spline(.x$psi, .x$loglikelihood))

MLE_data <- mods |>
  sapply(
    function(mod) {
      optimize(
        function(psi) predict(mod, psi)$y, 
        lower = psi_grid |> head(1), 
        upper = psi_grid |> tail(1), 
        maximum = TRUE
      )}) |> 
  t() |> 
  data.frame() |> 
  rownames_to_column("Iteration") |> 
  dplyr::rename(MLE = maximum,
                Maximum = objective)

curves <- mapply(
  function(mod, maximum) function(psi) predict(mod, psi)$y - maximum,
  mods,
  MLE_data$Maximum)

crit <- qchisq(0.95, 1) / 2

mapply(function(curve, MLE) {
  uniroot(function(psi) curve(psi) + crit,
          interval = c(psi_grid |> head(1), MLE))$root |> round(3)
  },
  curves,
  MLE_data$MLE)

mapply(function(curve, MLE) {
  uniroot(function(psi) curve(psi) + crit,
          interval = c(MLE, psi_grid |> tail(1)))$root |> round(3)
},
curves,
MLE_data$MLE)

i <- 2

uniroot(function(psi) curves[[i]](psi) + crit,
        interval = c(psi_grid |> head(1), MLE_data$MLE[[i]]))$root |> round(3)


ggplot() + 
  geom_function(fun = function(psi) curves[[i]](psi) + crit) +
  scale_x_continuous(limits = c(0, 1.9))

curves[[i]](psi_grid |> head(1)) + crit
curves[[i]](MLE_data$MLE[[i]]) + crit
curves[[i]](psi_grid |> tail(1)) + crit

