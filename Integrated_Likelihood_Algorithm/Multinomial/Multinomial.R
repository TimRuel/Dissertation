library(tidyverse)
library(geomtextpath)
library(viridis)
library(ggnewscale)
library(future)
library(zeallot)
library(parallelly)

plan(multisession, workers = availableCores())

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")

set.seed(38498984)

# Desert Rodents
# data <- c(1, 1, 2, 4, 7, 10)

# Birds in Balrath Woods
# data <- c(1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 6, 8)
 
# Birds in Killarney Woodlands
data <- c(1, 3, 4, 6, 7, 10, 14, 30)

step_size <- 0.01

psi_grid <- data |> 
  length() |> 
  log() |> 
  plyr::round_any(step_size, ceiling) |> 
  seq(0, to = _, step_size)

R <- 250

u_list <- LaplacesDemon::rdirichlet(R, rep(1, length(data))) |> 
  t() |> 
  data.frame() |> 
  as.list()

theta_MLE <- data / sum(data)

psi_MLE <- PoI_fn(theta_MLE)

omega_hat_list <- u_list |>
  purrr::map(get_omega_hat, psi_MLE)

multinomial_entropy_values_IL <- omega_hat_list |> 
  get_multinomial_entropy_values_IL(data, psi_grid)

multinomial_entropy_values_PL <- data |> 
  get_multinomial_entropy_values_PL(psi_grid) 

log_likelihood_vals_tidy <- data.frame(psi = psi_grid,
                                       Integrated = multinomial_entropy_values_IL,
                                       Profile = multinomial_entropy_values_PL) |> 
  pivot_longer(cols = c("Integrated", "Profile"),
               names_to = "Pseudolikelihood",
               values_to = "loglikelihood")

# Desert Rodents
# R = 250, step_size = 0.01, seed = 38498984
# saveRDS(log_likelihood_vals_tidy, "desert_rodents_R=250_step_size=0.01.Rda")
# log_likelihood_vals_tidy <- readRDS("desert_rodents_R=250_step_size=0.01.Rda")

# Birds in Balrath Woods
# R = 250, step_size = 0.01, seed = 38498984
# saveRDS(log_likelihood_vals_tidy, "birds_in_balrath_woods_R=250_step_size=0.01.Rda")
# log_likelihood_vals_tidy <- readRDS("birds_in_balrath_woods_R=250_step_size=0.01.Rda")

# Birds in Killarney Woodlands
# R = 250, step_size = 0.01, seed = 38498984
# saveRDS(log_likelihood_vals_tidy, "birds_in_killarney_woodlands_R=250_step_size=0.01.Rda")
# log_likelihood_vals_tidy <- readRDS("birds_in_killarney_woodlands_R=250_step_size=0.01.Rda")

spline_fitted_models <- log_likelihood_vals_tidy |>
  group_by(Pseudolikelihood) |> 
  group_map(~ smooth.spline(.x$psi, .x$loglikelihood)) |> 
  set_names(c("Integrated", "Profile"))

MLE_data <- spline_fitted_models |>
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
  rownames_to_column("Pseudolikelihood") |> 
  dplyr::rename(MLE = maximum,
                Maximum = objective) |> 
  mutate(MLE_label = c("hat(psi)[IL]", "hat(psi)[P]"))

c(IL_curve, PL_curve) %<-% mapply(
  function(mod, maximum) function(psi) predict(mod, psi)$y - maximum,
  spline_fitted_models,
  MLE_data$Maximum)

ggplot() +
  stat_function(fun = PL_curve,
                geom = "textpath",
                label = "Profile",
                aes(color = "Profile"),
                linewidth = 1,
                hjust = 0.1,
                show.legend = FALSE) +
  stat_function(fun = IL_curve,
                geom = "textpath",
                label = "Integrated",
                aes(color = "Integrated"),
                linewidth = 1,
                hjust = 0.1,
                show.legend = FALSE) +
  geom_hline(yintercept = 0,
             linetype = 5) +
  geom_labelvline(aes(xintercept = as.numeric(MLE),
                      label = MLE_label,
                      color = Pseudolikelihood),
                  data = MLE_data,
                  parse = TRUE,
                  show.legend = FALSE) +
  ylab("Log-Likelihood") +
  scale_x_continuous(expand = c(0, 0),
                     # limits = c(1, 2) # Desert Rodents
                     # limits = c(2, 3) # Birds in Balrath Woods
                     limits = c(1.4, 2.1) # Birds in Killarney Woodlands
                     ) + 
  scale_y_continuous(expand = c(0.1, 0),
                     # limits = c(-3, 0.1) # Desert Rodents
                     # limits = c(-4, 0.1) # Birds in Balrath Woods
                     limits = c(-5.5, 0.1) # Birds in Killarney Woodlands
                     ) +
  scale_color_brewer(palette = "Set1") +
  xlab(expression(psi)) +
  theme_minimal() +
  theme(axis.line = element_line())

crit <- qchisq(0.95, 1) / 2

c(psi_hat_IL, psi_hat_PL) %<-% MLE_data$MLE

CI_lower_PL <- uniroot(function(psi) PL_curve(psi) + crit,
                      interval = c(psi_grid |> head(1), psi_hat_PL))$root |> 
  round(3)

CI_upper_PL <- uniroot(function(psi) PL_curve(psi) + crit,
                      interval = c(psi_hat_PL, psi_grid |> tail(1)))$root |> 
  round(3)

CI_lower_IL <- uniroot(function(psi) IL_curve(psi) + crit,
                       interval = c(psi_grid |> head(1), psi_hat_IL))$root |> 
  round(3)

CI_upper_IL <- uniroot(function(psi) IL_curve(psi) + crit,
                       interval = c(psi_hat_IL, psi_grid |> tail(1)))$root |> 
  round(3)

data.frame(MLE = c(psi_hat_IL, psi_hat_PL) |> round(3),
           CI_95 = c(paste0("(", CI_lower_IL, ", ", CI_upper_IL, ")"),
                     paste0("(", CI_lower_PL, ", ", CI_upper_PL, ")")),
           row.names = c("Integrated", "Profile"))
