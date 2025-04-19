library(tidyverse)
library(purrr)
library(geomtextpath)
library(viridis)
library(ggnewscale)
library(zeallot)
library(kableExtra)
library(stringr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("../Data/Data Generation.R")

print("Choose your pseudolikelihood data file.")
log_likelihood_vals_file_path <- file.choose()

pseudolikelihood_names <- c("Integrated", "Profile")

log_likelihood_vals <- readRDS(log_likelihood_vals_file_path) |>
  tidyr::pivot_longer(cols = all_of(pseudolikelihood_names),
                      names_to = "Pseudolikelihood",
                      values_to = "loglikelihood") |>
  mutate(Pseudolikelihood = Pseudolikelihood |>
           as_factor() |>
           fct_inorder()) |> 
  drop_na()

# log_likelihood_vals <- log_likelihood_vals |>
#   tidyr::pivot_longer(cols = all_of(pseudolikelihood_names),
#                       names_to = "Pseudolikelihood",
#                       values_to = "loglikelihood") |>
#   mutate(Pseudolikelihood = Pseudolikelihood |>
#            as_factor() |>
#            fct_inorder())

spline_fitted_models <- log_likelihood_vals |>
  group_by(Pseudolikelihood) |> 
  group_map(~ smooth.spline(.x$psi, .x$loglikelihood)) |> 
  set_names(pseudolikelihood_names)

MLE_data <- spline_fitted_models |>
  map(
    function(mod) {
      optimize(function(psi) predict(mod, psi)$y,
               lower = log_likelihood_vals |>
                 select(psi) |>
                 min(),
               upper = log_likelihood_vals |>
                 select(psi) |>
                 max(),
               maximum = TRUE) |>
        data.frame() |>
        mutate(MLE = as.numeric(maximum),
               Maximum = as.numeric(objective)) |>
        select(MLE, Maximum)
    }) |>
  do.call(rbind, args = _) |>
  mutate(MLE_label = c("hat(psi)[IL]", "hat(psi)[PL]")) |>
  rownames_to_column("Pseudolikelihood")

pseudo_log_likelihood_curves <- spline_fitted_models |> 
  map2(MLE_data$Maximum,
       function(mod, maximum) function(psi) predict(mod, psi)$y - maximum)

base_plot <- ggplot() +   
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "#444444", color = NA),  
    plot.background = element_rect(fill = "#2A2A2A", color = NA),  
    panel.grid.major = element_line(color = "#4C4C4C"),  
    panel.grid.minor = element_line(color = "#333333"),  
    axis.ticks = element_line(color = "white"),
    axis.text = element_text(color = "white"),
    axis.title = element_text(color = "white"),
    strip.text = element_text(color = "white"),
    plot.title = element_text(color = "white", face = "bold"),
    legend.background = element_rect(fill = "#444444"),
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white")
  )

base_plot +
  geom_point(aes(x = psi, y = loglikelihood, color = Pseudolikelihood),
             data = log_likelihood_vals,
             size = 3) +
  ylab("Log-Likelihood") +
  scale_color_brewer(palette = "Set1") +
  xlab(expression(psi))

crit <- qchisq(0.95, 1) / 2

min_psi <- min(log_likelihood_vals$psi)

max_psi <- max(log_likelihood_vals$psi)

conf_ints <- pseudo_log_likelihood_curves |> 
  map2(MLE_data$MLE,
       function(curve, MLE) {
         
         lower_bound <- tryCatch(
           
           uniroot(function(psi) curve(psi) + crit,
                   interval = c(min_psi, MLE))$root,
           
           error = function(e) return(NA)
         ) |> 
           round(3)
         
         upper_bound <- tryCatch(
           
           uniroot(function(psi) curve(psi) + crit,
                   interval = c(MLE, max_psi))$root,
           
           error = function(e) return(NA)
         ) |>
           round(3)
         
         return(c(lower_bound, upper_bound))
       }
  )

MLE_data |> 
  select(Pseudolikelihood, MLE) |> 
  mutate(MLE = MLE |> 
           round(3),
         conf_int = conf_ints |> 
           map(\(x) paste0("(", x[1], ", ", x[2], ")")),
         length = conf_ints |> 
           map(diff) |> 
           as.numeric()) |> 
  arrange(length) |> 
  kbl(col.names = c("Pseudolikelihood", 
                    "MLE",
                    "95% Confidence Interval",
                    "CI Length"),   
      align = "c") |> 
  kable_styling(bootstrap_options = c("striped", "hover")) 

c(stat_fn_IL, stat_fn_PL) %<-% map2(
  pseudo_log_likelihood_curves,
  pseudolikelihood_names,
  function(curve, pseudolikelihood_name) {

    stat_fn <- stat_function(fun = curve,
                             geom = "line",
                             aes(color = pseudolikelihood_name,
                                 linewidth = pseudolikelihood_name),
                             show.legend = FALSE)
    return(stat_fn)
  }
)

x_range <- conf_ints |> 
  unlist() |> 
  (\(x) c(floor(min(x)), ceiling(max(x))))()

y_min <- pseudo_log_likelihood_curves |> 
  map((\(curve) c(curve(x_range[1]), curve(x_range[2])))) |> 
  unlist() |> 
  min() |> 
  round()

label_data <- MLE_data |> 
  add_row(Pseudolikelihood = "Truth",
          MLE = dot_product(theta_0, weights),
          MLE_label = "psi[0]") 

base_plot +
  stat_fn_IL +
  stat_fn_PL +
  geom_hline(yintercept = 0,
             linetype = 5) +
  geom_hline(aes(yintercept = -crit),
             linewidth = 1.5) +
  geom_text(aes(x_range[2] - 0.25 * diff(x_range), -crit, label = "95% CI Line", vjust = 1.3)) +
  geom_vline(aes(xintercept = MLE,
                 color = Pseudolikelihood),
             data = label_data,
             show.legend = FALSE) +
  ggrepel::geom_label_repel(aes(x = MLE,
                                y = y_min / 2,
                                label = MLE_label,
                                color = Pseudolikelihood),
                            data = label_data,
                            direction = "y",
                            parse = TRUE,
                            show.legend = FALSE) +
  ylab("Log-Likelihood") +
  scale_x_continuous(expand = c(0, 0),
                     limits = x_range) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(y_min, 0.1)) +
  scale_color_brewer(palette = "Set1") +
  scale_linewidth_manual(values = c(4, 2.5, 1)) +
  xlab(expression(psi))






