library(future)
library(zeallot)
library(purrr)
library(dipsaus)
library(tidyverse)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../../utils.R")

# population <- "Desert Rodents"
population <- "Birds in Balrath Woods"
# population <- "Birds in Killarney Woodlands"

switch(population,
       
       "Desert Rodents" = {
         
         data <- c(1, 1, 2, 4, 7, 10)
         
         x_range <- c(1, 2)
         
         y_range <- c(-3, 0.1)
       },
       
       "Birds in Balrath Woods" = {
         
         data <- c(1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 6, 8)
         
         x_range <- c(2, 2.7)
         
         y_range <- c(-4, 0.1)
       },
       
       "Birds in Killarney Woodlands" = {
         
         data <- c(1, 3, 4, 6, 7, 10, 14, 30)
         
         x_range <- c(1.48, 2)
         
         y_range <- c(-5, 0.1)
       }
)  

set.seed(38497283)

theta_MLE <- data / sum(data)

psi_MLE <- entropy(theta_MLE)

n <- sum(data)

m <- length(data)

sigma <- theta_MLE*diag(m) - matrix(theta_MLE) %*% theta_MLE

psi_MLE_SE <- sqrt(sum(matrix(1 + log(theta_MLE)) %*% (1 + log(theta_MLE)) * sigma) / n)

num_std_errors <- 2.1
MoE <- num_std_errors * psi_MLE_SE

step_size <- 0.01

psi_grid <- psi_MLE %+-% MoE |> 
  rev() |> 
  (\(x) c(max(0, x[1]), min(log(m), x[2])))() |> 
  plyr::round_any(step_size, floor) |> 
  (\(x) seq(x[1], x[2], step_size))() 

R <- 250

tol <- 0.0001

alpha <- data + 1

################################################################################
############################## EUCLIDEAN DISTANCE ##############################
################################################################################

euclidean_distance <- function(u, omega) dist(matrix(c(u, omega), 
                                                     nrow = 2, 
                                                     byrow = TRUE),
                                              method = "euclid")[1]

c(u_list, omega_hat_list_euclid) %<-% get_omega_hat_list(euclidean_distance, psi_MLE, alpha, R, tol)

L <- u_list |> 
  map_dbl(\(u) likelihood(u, data)) |> 
  unlist() |> 
  as.numeric()

plan(multisession, workers = availableCores())

multinomial_entropy_values_modified_IL_euclid <- omega_hat_list_euclid |> 
  get_multinomial_entropy_values_modified_IL(L, data, psi_grid)

################################################################################
############################## CANBERRA DISTANCE ##############################
################################################################################

# plan(sequential)
# 
# canberra_distance <- function(u, t) dist(matrix(c(u, t), 
#                                                 nrow = 2, 
#                                                 byrow = TRUE),
#                                          method = "canberra")[1]
# 
# tol <- 0.01
# 
# c(u_list, omega_hat_list_canberra) %<-% get_omega_hat_list(canberra_distance, psi_MLE, alpha, R, tol)
# 
# L <- u_list |> 
#   map_dbl(\(u) likelihood(u, data)) |> 
#   unlist() |> 
#   as.numeric()
# 
# plan(multisession, workers = availableCores())
# 
# multinomial_entropy_values_modified_IL_canberra <- omega_hat_list_canberra |> 
#   get_multinomial_entropy_values_modified_IL(L, data, psi_grid)
# 
# ################################################################################
# ############################## EXPECTED LOGLIKELIHOOD ##########################
# ################################################################################
# 
# plan(sequential)
# 
# E_log_like <- function(u, t) sum(u * log(t), na.rm = TRUE) 
# 
# tol <- 0.01
# 
# c(u_list, omega_hat_list_E_log_like) %<-% get_omega_hat_list(E_log_like, psi_MLE, alpha, R, tol)
# 
# L <- u_list |> 
#   map_dbl(\(u) likelihood(u, data)) |> 
#   unlist() |> 
#   as.numeric()
# 
# plan(multisession, workers = availableCores())
# 
# multinomial_entropy_values_modified_IL_E_log_like <- omega_hat_list_E_log_like |> 
#   get_multinomial_entropy_values_modified_IL(L, data, psi_grid)

################################################################################
############################## u_list ##########################################
################################################################################

plan(sequential)

u_list <- LaplacesDemon::rdirichlet(R, alpha) |> 
  t() |> 
  data.frame() |> 
  as.list() 

L <- u_list |> 
  map_dbl(\(u) likelihood(u, data)) |> 
  unlist() |> 
  as.numeric()

plan(multisession, workers = availableCores())

multinomial_entropy_values_modified_IL_u_list <- u_list |> 
  get_multinomial_entropy_values_modified_IL(L, data, psi_grid)

u_list <- LaplacesDemon::rdirichlet(R, rep(1, length(alpha))) |> 
  t() |> 
  data.frame() |> 
  as.list()

multinomial_entropy_values_IL_u_list <- u_list |> 
  get_multinomial_entropy_values_IL(data, psi_grid)

################################################################################
############################## COMPARISON ##########################
################################################################################

plan(sequential)

log_likelihood_vals <- data.frame(psi = psi_grid,
                                  Euclidean = multinomial_entropy_values_modified_IL_euclid,
                                  mod_u = multinomial_entropy_values_modified_IL_u_list,
                                  u = multinomial_entropy_values_IL_u_list) 

Q_fns <- c("Euclidean", "mod_u", "u")

spline_fitted_models <- log_likelihood_vals |>
  tidyr::pivot_longer(cols = c("Euclidean", "mod_u", "u"),
                      names_to = "Q",
                      values_to = "loglikelihood") |> 
  group_by(Q) |> 
  group_map(~ smooth.spline(.x$psi, .x$loglikelihood)) |> 
  set_names(Q_fns)

MLE_data <- spline_fitted_models |>
  sapply(
    function(mod) {
      optimize(
        function(psi) predict(mod, psi)$y, 
        lower = 0, 
        upper = log(length(data)), 
        maximum = TRUE
      )}) |> 
  t() |> 
  data.frame() |> 
  rownames_to_column("Q") |> 
  dplyr::rename(MLE = maximum,
                Maximum = objective) |> 
  mutate(MLE_label = c("Euclidean", "mod_u", "u"))

pseudo_log_likelihood_curves <- spline_fitted_models |> 
  map2(MLE_data$Maximum,
       function(mod, maximum) function(psi) predict(mod, psi)$y - maximum)

c(stat_fn_euclid, stat_fn_mod_u, stat_fn_u) %<-% map2(
  pseudo_log_likelihood_curves,
  Q_fns,
  function(curve, Q_fn) {
    
    stat_fn <- stat_function(fun = curve,
                             geom = "textpath",
                             label = Q_fn,
                             aes(color = Q_fn),
                             linewidth = 1,
                             hjust = 0.1,
                             show.legend = FALSE,
                             xlim = c(0, log(length(data))))
    return(stat_fn)
  }
)

ggplot() +
  stat_fn_euclid +
  stat_fn_mod_u +
  stat_fn_u +
  geom_hline(yintercept = 0,
             linetype = 5) +
  geom_vline(aes(xintercept = as.numeric(MLE),
                 color = Q),
             data = MLE_data,
             show.legend = FALSE) +
  ggrepel::geom_label_repel(aes(x = as.numeric(MLE),
                                y = y_range[1] + 1.5,
                                label = MLE_label,
                                color = Q),
                            data = MLE_data,
                            direction = "y",
                            parse = TRUE,
                            show.legend = FALSE) +
  ylab("Log-Likelihood") +
  scale_x_continuous(expand = c(0, 0),
                     limits = x_range) + 
  scale_y_continuous(expand = c(0, 0),
                     limits = y_range) +
  scale_color_brewer(palette = "Set1") +
  xlab(expression(psi)) +
  theme_minimal() +
  theme(axis.line = element_line())

crit <- qchisq(0.95, 1) / 2

conf_ints <- pseudo_log_likelihood_curves |> 
  map2(MLE_data$MLE,
       function(curve, MLE) {
         
         lower_bound <- tryCatch(
           
           uniroot(function(psi) curve(psi) + crit,
                   interval = c(0, MLE))$root,
           
           error = function(e) return(0)
         ) |> 
           round(3)
         
         upper_bound <- tryCatch(
           
           uniroot(function(psi) curve(psi) + crit,
                   interval = c(MLE, log(length(data))))$root,
           
           error = function(e) return(log(length(data)))
         ) |> 
           round(3)
         
         return(c(lower_bound, upper_bound))
       }
  )

MLE_data |> 
  select(Q, MLE) |> 
  mutate(MLE = as.numeric(MLE) |> round(3),
         conf_int = map(conf_ints, \(x) paste0("(", x[1], ", ", x[2], ")")),
         length = map(conf_ints, diff)) |> 
  kbl(col.names = c("Q", 
                    "MLE",
                    "95% Confidence Interval",
                    "CI Length"),   
      align = "c",
      caption = population) |> 
  kable_styling(bootstrap_options = c("striped", "hover")) 
