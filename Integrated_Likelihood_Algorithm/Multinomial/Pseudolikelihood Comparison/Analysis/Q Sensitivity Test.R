library(future)
library(purrr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../../utils.R")

population <- "Desert Rodents"
# population <- "Birds in Balrath Woods"
# population <- "Birds in Killarney Woodlands"

switch(population,
       
       "Desert Rodents" = {
         
         seed <- 1996
         
         data <- c(1, 1, 2, 4, 7, 10)
         
         step_size <- 0.01
         
         log_likelihood_vals_file_path <- "desert_rodents_R=250_step_size=0.01.Rda"
         
         x_range <- c(1, 2)
         
         y_range <- c(-3, 0.1)
       },
       
       "Birds in Balrath Woods" = {
         
         seed <- 7835
         
         data <- c(1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 6, 8)
         
         step_size <- 0.01
         
         log_likelihood_vals_file_path <- "birds_in_balrath_woods_R=250_step_size=0.01.Rda"
       },
       
       "Birds in Killarney Woodlands" = {
         
         seed <- 1996
         
         data <- c(1, 3, 4, 6, 7, 10, 14, 30)
         
         step_size <- 0.01
         
         log_likelihood_vals_file_path <- "birds_in_killarney_woodlands_R=250_step_size=0.01.Rda"
       }
)  

set.seed(seed)

theta_MLE <- data / sum(data)

psi_MLE <- entropy(theta_MLE)

psi_grid <- data |> 
  length() |> 
  log() |> 
  plyr::round_any(step_size, floor) |> 
  seq(0, to = _, step_size)

R <- 250

plan(multisession, workers = availableCores())

alpha <- data + 1

u_list_mod_IL <- LaplacesDemon::rdirichlet(R, alpha) |>
  t() |>
  data.frame() |>
  as.list()

distance_euclid <- function(u, t) dist(matrix(c(u, t), 
                                              nrow = 2, 
                                              byrow = TRUE),
                                       method = "euclid")[1] 

fn <- function(omega) distance_euclid(u_list_mod_IL[[100]], omega)
gr <- function(omega) nloptr::nl.grad(omega, fn)
heq <- function(omega) c(sum(omega) - 1, entropy(omega) - psi_MLE)
heqjac <- function(omega) nloptr::nl.jacobian(omega, heq)

omega_hat <- nloptr::auglag(x0 = rep(1, 6) / 6,
                            fn = fn,
                            gr = gr,
                            heq = heq,
                            heqjac = heqjac,
                            lower = rep(0, length(u)),
                            localsolver = "LBFGS")$par














################################################################################
############################## EUCLIDEAN DISTANCE ##############################
################################################################################

distance_euclid <- function(u, t) dist(matrix(c(u, t), 
                                              nrow = 2, 
                                              byrow = TRUE),
                                       method = "euclidean")[1] 

omega_hat_list_mod_IL_euclid <- u_list_mod_IL |>
  map(get_omega_hat, psi_MLE, distance_euclid)

multinomial_entropy_values_modified_IL_euclid <- omega_hat_list_mod_IL_euclid |> 
  get_multinomial_entropy_values_modified_IL(u_list_mod_IL, data, psi_grid)

################################################################################
############################## CANBERRA DISTANCE ##############################
################################################################################

distance_canberra <- function(u, t) dist(matrix(c(u, t), 
                                                nrow = 2, 
                                                byrow = TRUE),
                                         method = "canberra")[1]   

omega_hat_list_mod_IL_canberra <- u_list_mod_IL |>
  map(get_omega_hat, psi_MLE, distance_canberra)

multinomial_entropy_values_modified_IL_canberra <- omega_hat_list_mod_IL_canberra |> 
  get_multinomial_entropy_values_modified_IL(u_list_mod_IL, data, psi_grid)

################################################################################
############################## EXPECTED LOGLIKELIHOOD ##########################
################################################################################

E_log_like <- function(u, t) sum(u * log(t), na.rm = TRUE) 

omega_hat_list_mod_IL_E_log_like <- u_list_mod_IL |>
  map(get_omega_hat, psi_MLE, E_log_like)

multinomial_entropy_values_modified_IL_E_log_like <- omega_hat_list_mod_IL_E_log_like |> 
  get_multinomial_entropy_values_modified_IL(u_list_mod_IL, data, psi_grid)

################################################################################
############################## COMPARISON ##########################
################################################################################

log_likelihood_vals <- data.frame(psi = psi_grid,
                                  Euclidean = multinomial_entropy_values_modified_IL_euclid,
                                  Canberra = multinomial_entropy_values_modified_IL_canberra,
                                  E_log_like = multinomial_entropy_values_modified_IL_E_log_like) 

Q_fns <- c("Euclidean", "Expected Loglikelihood")

spline_fitted_models <- log_likelihood_vals |>
  tidyr::pivot_longer(cols = c("Euclidean", "E_log_like"),
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
  mutate(MLE_label = c("Euclidean", "E_log_like"))

pseudo_log_likelihood_curves <- spline_fitted_models |> 
  map2(MLE_data$Maximum,
       function(mod, maximum) function(psi) predict(mod, psi)$y - maximum)

c(stat_fn_euclid, stat_fn_E_log_like) %<-% map2(
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
  stat_fn_E_log_like +
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
         
         return(paste0("(", lower_bound, ", ", upper_bound, ")"))
       }
  ) |> 
  unlist()

MLE_data |> 
  select(Q, MLE) |> 
  mutate(MLE = as.numeric(MLE) |> round(3),
         conf_int = conf_ints) |> 
  kbl(col.names = c("Q", 
                    "MLE",
                    "95% Confidence Interval"),   
      align = "c",
      caption = population) |> 
  kable_styling(bootstrap_options = c("striped", "hover")) 
