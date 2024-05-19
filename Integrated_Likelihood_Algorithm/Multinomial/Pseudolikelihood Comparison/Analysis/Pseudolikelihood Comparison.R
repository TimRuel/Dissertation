library(tidyverse)
library(purrr)
library(geomtextpath)
library(viridis)
library(ggnewscale)
library(zeallot)
library(kableExtra)
library(stringr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

log_likelihood_vals_file_path <- file.choose()

log_likelihood_vals <- readRDS(log_likelihood_vals_file_path)

population <- log_likelihood_vals_file_path |>  
  str_remove("^.*\\\\") |> 
  str_remove("_R.*$") |> 
  str_replace_all("_", " ") |> 
  tools::toTitleCase()

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

pseudolikelihood_names <- c("Modified Integrated", "Integrated", "Profile")

spline_fitted_models <- log_likelihood_vals |>
  tidyr::pivot_longer(cols = c("Mod_Integrated", "Integrated", "Profile"),
                      names_to = "Pseudolikelihood",
                      values_to = "loglikelihood") |> 
  group_by(Pseudolikelihood) |> 
  group_map(~ smooth.spline(.x$psi, .x$loglikelihood)) |> 
  set_names(pseudolikelihood_names)

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
  rownames_to_column("Pseudolikelihood") |> 
  dplyr::rename(MLE = maximum,
                Maximum = objective) |> 
  mutate(MLE_label = c("hat(psi)[m-IL]", "hat(psi)[IL]", "hat(psi)[PL]"))

pseudo_log_likelihood_curves <- spline_fitted_models |> 
  map2(MLE_data$Maximum,
       function(mod, maximum) function(psi) predict(mod, psi)$y - maximum)

c(stat_fn_mod_IL, stat_fn_IL, stat_fn_PL) %<-% map2(
  pseudo_log_likelihood_curves,
  pseudolikelihood_names,
  function(curve, pseudolikelihood_name) {
    
    stat_fn <- stat_function(fun = curve,
                             geom = "textpath",
                             label = pseudolikelihood_name,
                             aes(color = pseudolikelihood_name),
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

log_likelihood_vals |> 
  tidyr::pivot_longer(cols = c("Mod_Integrated", "Integrated", "Profile"),
                      names_to = "Pseudolikelihood",
                      values_to = "loglikelihood") |> 
  ggplot() +
  geom_point(aes(x = psi, y = loglikelihood, color = Pseudolikelihood),
             size = 1) +
  ylab("Log-Likelihood") +
  scale_color_brewer(palette = "Set1") +
  # scale_x_continuous(limits = c(2, 2.7)) +
  # scale_y_continuous(limits = c(-110, -85)) +
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
  select(Pseudolikelihood, MLE) |> 
  mutate(MLE = as.numeric(MLE) |> round(3),
         conf_int = map(conf_ints, \(x) paste0("(", x[1], ", ", x[2], ")")),
         length = map(conf_ints, diff)) |> 
  kbl(col.names = c("Pseudolikelihood", 
                    "MLE",
                    "95% Confidence Interval",
                    "CI Length"),   
      align = "c",
      caption = population) |> 
  kable_styling(bootstrap_options = c("striped", "hover")) 





