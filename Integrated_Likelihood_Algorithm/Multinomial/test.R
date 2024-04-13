library(tidyverse)
library(plyr)
library(Rmpfr)
library(geomtextpath)
library(viridis)
library(ggnewscale)
library(zeallot)
library(pbapply)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("auxiliary functions.R")

#set.seed(1996)

data <- c(1, 1, 2, 4, 7, 10)

step_size <- 0.1

R <- 10

psi_grid <- data |> 
  length() |> 
  log() |> 
  round_any(step_size, ceiling) |> 
  seq(0, to = _, step_size)
 
multinomial_entropy_IL_values <- data |> 
  get_multinomial_entropy_IL_values(psi_grid) |> 
  pbreplicate(R, expr = _, simplify = FALSE) |> 
  do.call(rbind, args = _)

multinomial_entropy_IL_values |> 
  apply(2, mean)

# Initialize vector for holding values of profile likelihood
L_p <- mpfrArray(NA, precBits = 106, dim = N)

# Initialize progress bar for for loop
pb = txtProgressBar(min = 0, max = N, initial = 0, style = 3) 

theta_hat_p <- theta.hat

for (j in N_lower:1) {
  
  # Update progress bar
  setTxtProgressBar(pb, j)
  
  # Find value of theta that minimizes objective function subject to constraints
  theta_hat_p <- auglag(x0 = theta_hat_p,
                        fn = function(theta) -sum(theta.hat*log(theta)),
                        heq = function(theta) c(sum(theta) - 1, g(theta) - psi1[j]),
                        lower = rep(0, m))$par
  
  L_p[j] <- likelihood(theta_hat_p)
}

theta_hat_p <- theta.hat

for (j in (N_lower + 1):N) {
  
  # Update progress bar
  setTxtProgressBar(pb, j)
  
  # Find value of theta that minimizes objective function subject to constraints
  theta_hat_p <- auglag(x0 = theta_hat_p,
                        fn = function(theta) -sum(theta.hat*log(theta)),
                        heq = function(theta) c(sum(theta) - 1, g(theta) - psi1[j]),
                        lower = rep(0, m))$par
  
  L_p[j] <- likelihood(theta_hat_p)
}

log_likelihood_vals <- data.frame(psi = psi1,
                                  Integrated = L |>
                                    apply(2, mean) |>
                                    log() |>
                                    as.double(),
                                  Profile = L_p |> 
                                    log() |> 
                                    as.double())

log_likelihood_vals_tidy <- log_likelihood_vals |> 
  pivot_longer(cols = c("Integrated", "Profile"),
               names_to = "Pseudolikelihood",
               values_to = "loglikelihood")

spline_fitted_models <- log_likelihood_vals_tidy |>
  group_by(Pseudolikelihood) |> 
  group_map(~ smooth.spline(.x$psi, .x$loglikelihood)) |> 
  set_names(c("Integrated", "Profile"))

MLE_data <- spline_fitted_models |>
  sapply(
    function(mod) {
      optimize(
        function(psi) predict(mod, psi)$y, 
        lower = psi1 |> head(1), 
        upper = psi1 |> tail(1), 
        maximum = TRUE
      )}) |> 
  t() |> 
  data.frame() |> 
  rownames_to_column("Pseudolikelihood") |> 
  dplyr::rename(MLE = maximum,
                Maximum = objective) |> 
  mutate(MLE_label = c("hat(psi)[IL]", "hat(psi)[P]"))

c(IL_curve, P_curve) %<-% mapply(
  function(mod, maximum) function(psi) predict(mod, psi)$y - maximum,
  spline_fitted_models,
  MLE_data$Maximum)

ggplot() +
  stat_function(fun = P_curve,
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
  scale_x_continuous(limits = c(1, 2),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(-4, 0.1),
                     expand = c(0.1, 0)) +
  scale_color_brewer(palette = "Set1") +
  xlab(expression(psi)) +
  theme_minimal() +
  theme(axis.line = element_line())

crit <- qchisq(0.95, 1) / 2

c(psi_hat_IL, psi_hat_P) %<-% MLE_data$MLE

CI_lower_P <- uniroot(function(psi) P_curve(psi) + crit,
                      interval = c(psi1 |> head(1), psi_hat_P))$root |> 
  round(3)

CI_upper_P <- uniroot(function(psi) P_curve(psi) + crit,
                      interval = c(psi_hat_P, psi1 |> tail(1)))$root |> 
  round(3)

CI_lower_IL <- uniroot(function(psi) IL_curve(psi) + crit,
                       interval = c(psi1 |> head(1), psi_hat_IL))$root |> 
  round(3)

CI_upper_IL <- uniroot(function(psi) IL_curve(psi) + crit,
                       interval = c(psi_hat_IL, psi1 |> tail(1)))$root |> 
  round(3)

data.frame(MLE = c(psi_hat_IL, psi_hat_P) |> round(3),
           CI_95 = c(paste0("(", CI_lower_IL, ", ", CI_upper_IL, ")"),
                     paste0("(", CI_lower_P, ", ", CI_upper_P, ")")),
           row.names = c("Integrated", "Profile"))
