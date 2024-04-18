library(LaplacesDemon)
library(future.apply)
library(parallelly)
library(furrr)
library(progressr)

plan(list(tweak(multisession, workers = 4)), tweak(multisession, workers = 3))

#handlers("cli")

data <- c(1, 1, 2, 4, 7, 10)

step_size <- 0.01

psi_grid <- data |> 
  length() |> 
  log() |> 
  round_any(step_size, ceiling) |> 
  seq(0, to = _, step_size)

# R <- 10

n_sims <- 1000

sims <- rmultinom(n_sims, length(data), data) |> 
  data.frame() |> 
  as.list()

test <- sims |>
  future_map(get_multinomial_entropy_values_PL, psi_grid, .progress = TRUE)

get_multinomial_entropy_values_PL(sims[[2]], psi_grid)

a <-  sample(0:6, 6)
a <- a / sum(a)
psi_grid |> 
  purrr::map(get_theta_hat, a) |> 
  sapply(likelihood, data) |> 
  log() |>
  as.double()

get_theta_hat(1.3, a)
  

# test <- sims |> 
#   future_apply(1, \(x) {
#     u <- rdirichlet(R, rep(1, length(x)))
#     get_multinomial_entropy_values_IL(u, x, psi_grid)
#     },
#     future.seed = TRUE)

mods <- test |> 
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

