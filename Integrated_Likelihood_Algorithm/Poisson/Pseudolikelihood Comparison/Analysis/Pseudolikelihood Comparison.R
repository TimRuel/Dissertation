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
                             xlim = c(psi_grid |> head(1), psi_grid |> tail(1)))
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
                                y = -1,
                                label = MLE_label,
                                color = Pseudolikelihood),
                            data = MLE_data,
                            direction = "y",
                            parse = TRUE,
                            show.legend = FALSE) +
  ylab("Log-Likelihood") +
  # scale_x_continuous(expand = c(0, 0),
  #                    limits = c(70, 78)) +
  # scale_y_continuous(expand = c(0, 0),
  #                    limits = c(-1, 0)) +
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
         
         # upper_bound <- uniroot(function(psi) curve(psi) + crit,
         #                        interval = c(MLE, psi_grid |> tail(1)))$root |> 
         #   round(3)
         
         upper_bound <- tryCatch(

           uniroot(function(psi) curve(psi) + crit,
                   interval = c(MLE, psi_grid |> tail(1)))$root,

           error = function(e) return(psi_grid |> tail(1))
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
      align = "c") |> 
  kable_styling(bootstrap_options = c("striped", "hover")) 

log_likelihood_vals |> 
  ggplot(aes(x = psi, y = Mod_Integrated)) +
  geom_point()


