library(tidyverse)
library(parallelly)
library(furrr)
library(purrr)

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")

set.seed(38498984)

plan(list(tweak(multisession, workers = 4)), tweak(multisession, workers = 3))
# plan(multisession, workers = availableCores())
# plan(sequential)

# Desert Rodents
data <- c(1, 1, 2, 4, 7, 10)

# Birds in Balrath Woods
# data <- c(1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 6, 8)

# Birds in Killarney Woodlands
# data <- c(1, 3, 4, 6, 7, 10, 14, 30)

psi_0 <- PoI_fn(data / sum(data))

step_size <- 0.01

psi_grid <- data |> 
  length() |> 
  log() |> 
  plyr::round_any(step_size, ceiling) |> 
  seq(0, to = _, step_size)

R <- 250

n_sims <- 1000

sims <- n_sims |> 
  rmultinom(sum(data), data) |> 
  data.frame() |> 
  as.list()

# u_list <- n_sims |> 
#   replicate({
#     R |> 
#       LaplacesDemon::rdirichlet(rep(1, length(data)))|> 
#       t() |> 
#       data.frame() |> 
#       as.list()
#     }, 
#     simplify = FALSE)

# omega_hat_lists <- u_list |> 
#   future_map2(sims, 
#               \(x, y) x |> 
#                 purrr::map(\(z) z |> 
#                              get_omega_hat(PoI_fn(y / sum(y)))), 
#               .progress = TRUE)

# seed = 38498984, set at top of script
# saveRDS(omega_hat_lists, "desert_rodents_omega_hat_lists.Rda")
omega_hat_lists <- readRDS("desert_rodents_omega_hat_lists.Rda")

# seed = 38498984, set at top of script
# saveRDS(omega_hat_lists, "birds_in_balrath_woods_omega_hat_lists.Rda")
# omega_hat_lists <- readRDS("birds_in_balrath_woods_omega_hat_lists.Rda")

# seed = 38498984, set at top of script
# saveRDS(omega_hat_lists, "birds_in_killarney_woodlands_omega_hat_lists.Rda")
# omega_hat_lists <- readRDS("birds_in_killarney_woodlands_omega_hat_lists.Rda")

K <- 60

stime <- system.time({
  
  multinomial_entropy_values_IL <- omega_hat_lists[11:K] |>
    future_map2(sims[11:K],
                \(x, y) get_multinomial_entropy_values_IL(x, y, psi_grid),
                .progress = TRUE)
})

stime

# saveRDS(multinomial_entropy_values_IL, "desert_rodents_IL_sims.Rda")
multinomial_entropy_values_IL <- readRDS("desert_rodents_IL_sims.Rda")

# saveRDS(multinomial_entropy_values_IL, "birds_in_balrath_woods_IL_sims.Rda")
# multinomial_entropy_values_IL <- readRDS("birds_in_balrath_woods_IL_sims.Rda")

# saveRDS(multinomial_entropy_values_IL, "birds_in_killarney_woodlands_IL_sims.Rda")
# multinomial_entropy_values_IL <- readRDS("birds_in_killarney_woodlands_IL_sims.Rda")

stime <- system.time({
  
  multinomial_entropy_values_PL <- sims |>
    future_map(\(x) get_multinomial_entropy_values_PL(x, psi_grid),
               .progress = TRUE)
})

stime

# saveRDS(multinomial_entropy_values_PL, "desert_rodents_PL_sims.Rda")
multinomial_entropy_values_PL <- readRDS("desert_rodents_PL_sims.Rda")

# saveRDS(multinomial_entropy_values_PL, "birds_in_balrath_woods_PL_sims.Rda")
# multinomial_entropy_values_PL <- readRDS("birds_in_balrath_woods_PL_sims.Rda")

# saveRDS(multinomial_entropy_values_PL, "birds_in_killarney_woodlands_PL_sims.Rda")
# multinomial_entropy_values_PL <- readRDS("birds_in_killarney_woodlands_PL_sims.Rda")

mods_PL <- multinomial_entropy_values_PL |> 
  data.frame() |> 
  mutate(psi = psi_grid) |> 
  pivot_longer(cols = -psi,
               names_to = "Iteration",
               values_to = "loglikelihood") |> 
  group_by(Iteration) |> 
  filter(is.finite(loglikelihood)) |> 
  dplyr::group_map(~ smooth.spline(.x$psi, .x$loglikelihood))

mods_IL <- multinomial_entropy_values_IL |> 
  data.frame() |> 
  mutate(psi = psi_grid) |> 
  pivot_longer(cols = -psi,
               names_to = "Iteration",
               values_to = "loglikelihood") |> 
  group_by(Iteration) |> 
  filter(is.finite(loglikelihood)) |> 
  dplyr::group_map(~ smooth.spline(.x$psi, .x$loglikelihood))

MLE_data_PL <- mods_PL |>
  sapply(
    function(mod) {
      optimize(
        function(psi) predict(mod, psi)$y, 
        lower = psi_grid |> head(1), 
        upper = psi_grid |> tail(1), 
        maximum = TRUE
      )}) |> 
  t() |> 
  apply(2, as.numeric) |> 
  data.frame() |> 
  dplyr::rename(MLE = maximum,
                Maximum = objective)

MLE_data_IL <- mods_IL |>
  sapply(
    function(mod) {
      optimize(
        function(psi) predict(mod, psi)$y, 
        lower = psi_grid |> head(1), 
        upper = psi_grid |> tail(1), 
        maximum = TRUE
      )}) |> 
  t() |> 
  apply(2, as.numeric) |> 
  data.frame() |> 
  dplyr::rename(MLE = maximum,
                Maximum = objective)

curves_PL <- mapply(
  function(mod, maximum) function(psi) predict(mod, psi)$y - maximum,
  mods_PL,
  MLE_data_PL$Maximum)

curves_IL <- mapply(
  function(mod, maximum) function(psi) predict(mod, psi)$y - maximum,
  mods_IL,
  MLE_data_IL$Maximum)

crit <- qchisq(0.95, 1) / 2

CI_lower_PL <- mapply(function(curve, MLE) {
  tryCatch(
    uniroot(function(psi) curve(psi) + crit,
          interval = c(psi_grid |> head(1), MLE))$root |> round(3),
    error = function(e) return(NA)
    )
  },
  curves_PL,
  MLE_data_PL$MLE)

CI_upper_PL <- mapply(function(curve, MLE) {
  tryCatch(
    uniroot(function(psi) curve(psi) + crit,
          interval = c(MLE, psi_grid |> tail(1)))$root |> round(3),
    error = function(e) return(NA)
    )
  },
  curves_PL,
  MLE_data_PL$MLE)

sim_results_PL <- MLE_data_PL |> 
  mutate(Lower = CI_lower_PL, 
         Upper = CI_upper_PL,
         psi_0 = psi_0,
         Covered = psi_0 |> between(Lower, Upper)) |> 
  summarise(Bias = mean(MLE - psi_0),
            SD = sd(MLE),
            RMSE = sqrt(mean((MLE - psi_0)^2)),
            Coverage = mean(Covered, na.rm = TRUE),
            Length = mean(Upper - Lower, na.rm = TRUE)) |> 
  t() |> 
  as.data.frame() |> 
  round(3) |> 
  setNames("Profile")

# saveRDS(sim_results_PL, "desert_rodents_sim_results_PL.Rda")
# sim_results_PL <- readRDS("desert_rodents_sim_results_PL.Rda")

# saveRDS(sim_results_PL, "birds_in_balrath_woods_sim_results_PL.Rda")
# sim_results_PL <- readRDS("birds_in_balrath_woods_sim_results_PL.Rda")

# saveRDS(sim_results_PL, "birds_in_killarney_woodlands_sim_results_PL.Rda")
# sim_results_PL <- readRDS("birds_in_killarney_woodlands_sim_results_PL.Rda")

sim_results_PL

CI_lower_IL <- mapply(function(curve, MLE) {
  tryCatch(
    uniroot(function(psi) curve(psi) + crit,
            interval = c(psi_grid |> head(1), MLE))$root |> round(3),
    error = function(e) return(NA)
  )
},
curves_IL,
MLE_data_IL$MLE)

CI_upper_IL <- mapply(function(curve, MLE) {
  tryCatch(
    uniroot(function(psi) curve(psi) + crit,
            interval = c(MLE, psi_grid |> tail(1)))$root |> round(3),
    error = function(e) return(NA)
  )
},
curves_IL,
MLE_data_IL$MLE)

sim_results_IL <- MLE_data_IL |> 
  mutate(Lower = CI_lower_IL, 
         Upper = CI_upper_IL,
         psi_0 = psi_0,
         Covered = psi_0 |> between(Lower, Upper)) |> 
  summarise(Bias = mean(MLE - psi_0),
            SD = sd(MLE),
            RMSE = sqrt(mean((MLE - psi_0)^2)),
            Coverage = mean(Covered, na.rm = TRUE),
            Length = mean(Upper - Lower, na.rm = TRUE)) |> 
  t() |> 
  as.data.frame() |> 
  round(3) |> 
  setNames("Integrated")

saveRDS(sim_results_IL, "desert_rodents_sim_results_IL.Rda")
sim_results_IL <- readRDS("desert_rodents_sim_results_IL.Rda")

# saveRDS(sim_results_IL, "birds_in_balrath_woods_sim_results_IL.Rda")
# sim_results_IL <- readRDS("birds_in_balrath_woods_sim_results_IL.Rda")

# saveRDS(sim_results_IL, "birds_in_killarney_woodlands_sim_results_IL.Rda")
# sim_results_IL <- readRDS("birds_in_killarney_woodlands_sim_results_IL.Rda")

sim_results_IL










