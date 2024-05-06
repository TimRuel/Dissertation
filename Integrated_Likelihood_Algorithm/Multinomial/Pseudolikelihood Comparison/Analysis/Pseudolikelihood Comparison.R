library(tidyverse)
library(geomtextpath)
library(viridis)
library(ggnewscale)
library(zeallot)
library(kableExtra)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

population <- "Desert Rodents"
# population <- "Birds in Balrath Woods"
# population <- "Birds in Killarney Woodlands"

switch(population,
       
       "Desert Rodents" = {
         
         data <- c(1, 1, 2, 4, 7, 10)
         
         step_size <- 0.01
         
         log_likelihood_vals_file_path <- "desert_rodents_R=250_step_size=0.01.Rda"
         
         x_range <- c(1, 2)
         
         y_range <- c(-3, 0.1)
       },
       
       "Birds in Balrath Woods" = {
      
         data <- c(1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 6, 8)
         
         step_size <- 0.05
         
         log_likelihood_vals_file_path <- "birds_in_balrath_woods_R=250_step_size=0.05.Rda"
         
         x_range <- c(2, 2.7)
         
         y_range <- c(-4, 0.1)
       },
       
       "Birds in Killarney Woodlands" = {
         
         data <- c(1, 3, 4, 6, 7, 10, 14, 30)
         
         step_size <- 0.01
         
         log_likelihood_vals_file_path <- "birds_in_killarney_woodlands_R=250_step_size=0.01.Rda"
         
         x_range <- c(1.48, 2)
         
         y_range <- c(-5, 0.1)
       }
)  

max_psi_val <- data |> 
  length() |> 
  log() |> 
  plyr::round_any(step_size, floor)

log_likelihood_vals <- readRDS(paste0("../Data/Pseudolikelihoods/", log_likelihood_vals_file_path))

pseudolikelihoods <- c("Modified Integrated", "Integrated", "Profile")

spline_fitted_models <- log_likelihood_vals |>
  group_by(Pseudolikelihood) |> 
  group_map(~ smooth.spline(.x$psi, .x$loglikelihood)) |> 
  set_names(pseudolikelihoods)

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
  mutate(MLE_label = c("hat(psi)[m-IL]", "hat(psi)[IL]", "hat(psi)[PL]"))

log_likelihood_curves <- mapply(
  function(mod, maximum) function(psi) predict(mod, psi)$y - maximum,
  spline_fitted_models,
  MLE_data$Maximum)

c(stat_fn_mod_IL, stat_fn_IL, stat_fn_PL) %<-% purrr::map2(
  log_likelihood_curves,
  pseudolikelihoods,
  function(curve, pseudolikelihood) {
    
    stat_fn <- stat_function(fun = curve,
                             geom = "textpath",
                             label = pseudolikelihood,
                             aes(color = pseudolikelihood),
                             linewidth = 1,
                             hjust = 0.1,
                             show.legend = FALSE,
                             xlim = c(0, log(length(data))))
    return(stat_fn)
    }
  )

ggplot() +
  stat_fn_mod_IL +
  stat_fn_IL + 
  stat_fn_PL +
  geom_hline(yintercept = 0,
             linetype = 5) +
  geom_vline(aes(xintercept = as.numeric(MLE),
                 color = Pseudolikelihood),
             data = MLE_data,
             show.legend = FALSE) +
  ggrepel::geom_label_repel(aes(x = as.numeric(MLE),
                                y = y_range[1] + 1.5,
                                label = MLE_label,
                                color = Pseudolikelihood),
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

c(psi_hat_mod_IL, psi_hat_IL, psi_hat_PL) %<-% MLE_data$MLE

conf_ints <- purrr::map2(
  log_likelihood_curves,
  MLE_data$MLE,
  function(curve, MLE) {
    
    lower_bound <- tryCatch(
      
      uniroot(function(psi) curve(psi) + crit,
              interval = c(0, MLE))$root,
      
      error = function(e) return(0)
      ) |> 
      round(3)
    
    upper_bound <- tryCatch(
      
      uniroot(function(psi) curve(psi) + crit,
              interval = c(MLE, max_psi_val))$root,
      
      error = function(e) return(log(length(data)))
      ) |> 
      round(3)
    
    return(paste0("(", lower_bound, ", ", upper_bound, ")"))
    }
  ) |> 
  unlist()

MLE_data |> 
  select(Pseudolikelihood, MLE) |> 
  mutate(MLE = as.numeric(MLE) |> round(3),
         conf_int = conf_ints) |> 
  kbl(col.names = c("Pseudolikelihood", 
                    "MLE",
                    "95% Confidence Interval"), 
      align = "c",
      caption = population) |> 
  kable_styling(bootstrap_options = c("striped", "hover")) 



# log_likelihood_vals |>
#   filter(Pseudolikelihood == "Mod_Integrated") |>
#   ggplot(aes(x = psi, y = loglikelihood)) +
#   geom_point()





