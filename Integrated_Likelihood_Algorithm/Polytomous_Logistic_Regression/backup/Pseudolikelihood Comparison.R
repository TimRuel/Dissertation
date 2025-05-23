library(tidyverse)
library(purrr)
library(geomtextpath)
library(viridis)
library(ggnewscale)
library(zeallot)
library(kableExtra)
library(stringr)
library(rstudioapi)
library(splines)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("utils.R")

print("Select population directory")
population_directory <- selectDirectory(path = getwd())

setwd(population_directory)

print("Choose your pseudolikelihood data file.")
log_likelihood_vals_file_path <- selectFile(path = getwd())

model <- fit_multinomial_logistic_model(data, formula)

step_size <- log_likelihood_vals_file_path |>  
  str_remove("^.*stepsize=") |> 
  str_extract("\\d+\\.\\d+") |> 
  as.numeric()

h <- log_likelihood_vals_file_path |>
  str_extract("(?<=h=).") 

psi_hat <- get_psi_hat(model, data, h)

psi_0 <- H_0 |> 
  filter(X1 == h) |> 
  pull(entropy)

pseudolikelihood_names <- c("Integrated", "Profile")

log_likelihood_vals <- readRDS(log_likelihood_vals_file_path) |>
  tidyr::pivot_longer(cols = all_of(pseudolikelihood_names),
                      names_to = "Pseudolikelihood",
                      values_to = "loglikelihood") |>
  mutate(Pseudolikelihood = Pseudolikelihood |>
           as_factor())

# log_likelihood_vals <- readRDS("Simulations/Marginal Entropy/log_likelihood_dfs/h=A/Sim2.Rda") |>
#   tidyr::pivot_longer(cols = all_of(pseudolikelihood_names),
#                       names_to = "Pseudolikelihood",
#                       values_to = "loglikelihood") |>
#   mutate(Pseudolikelihood = Pseudolikelihood |>
#            as_factor())

# log_likelihood_vals <- log_likelihood_vals |>
#   tidyr::pivot_longer(cols = all_of(pseudolikelihood_names),
#                       names_to = "Pseudolikelihood",
#                       values_to = "loglikelihood") |>
#   mutate(Pseudolikelihood = Pseudolikelihood |>
#            as_factor())

spline_fitted_models <- log_likelihood_vals |>
  drop_na(loglikelihood) |>
  group_by(Pseudolikelihood) |>
  group_map(~ smooth.spline(.x$psi, .x$loglikelihood, spar = 0.3)) |>
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
  rownames_to_column("Source")

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
             size = 3,
             show.legend = FALSE) +
  ggrepel::geom_label_repel(aes(x = psi,
                                y = loglikelihood,
                                label = Pseudolikelihood,
                                color = Pseudolikelihood),
                            data = log_likelihood_vals |> 
                              filter(psi == round(psi_hat, 1)),
                            direction = "y",
                            nudge_y = 2,
                            parse = TRUE,
                            show.legend = FALSE) +
  ylab("Log-Likelihood") +
  scale_color_brewer(palette = "Set1") +
  xlab(expression(psi))

alpha <- 0.05

crit <- qchisq(1 - alpha, 1) / 2

psi_grid <- log_likelihood_vals$psi |> 
  unique()

conf_ints <- pseudo_log_likelihood_curves |> 
  map2(MLE_data$MLE,
       function(curve, MLE) {
         
         lower_bound <- tryCatch(
           
           uniroot(function(psi) curve(psi) + crit,
                   interval = c(psi_grid |> head(1), MLE))$root,
           
           error = function(e) return(0)
         ) |> 
           round(3)
         
         upper_bound <- tryCatch(
           
           uniroot(function(psi) curve(psi) + crit,
                   interval = c(MLE, psi_grid |> tail(1)))$root,
           
           error = function(e) return(log(length(model$lev)))
         ) |>
           round(3)
         
         return(c(lower_bound, upper_bound))
       }
  )

MLE_data |> 
  select(Source, MLE) |> 
  mutate(MLE = MLE |> 
           round(3),
         conf_int = conf_ints |> 
           map(\(x) paste0("(", x[1], ", ", x[2], ")")),
         length = conf_ints |> 
           map(diff) |> 
           as.numeric()) |> 
  arrange(length) |> 
  add_row(Source = "Truth",
          MLE = psi_0) |>
  kbl(col.names = c("Source", 
                    "MLE",
                    "95% Confidence Interval",
                    "CI Length"),   
      align = "c") |> 
  kable_styling(bootstrap_options = c("striped", "hover")) |> 
  row_spec(3, color = "green", bold = TRUE)

c(stat_fn_IL, stat_fn_PL) %<-% map2(
  pseudo_log_likelihood_curves,
  pseudolikelihood_names,
  function(curve, pseudolikelihood_name) {
    
    stat_fn <- stat_function(fun = curve,
                             geom = "line",
                             aes(color = pseudolikelihood_name),
                             linewidth = 1.5,
                             show.legend = FALSE,
                             xlim = c(psi_grid |> head(1), psi_grid |> tail(1)))
    return(stat_fn)
  }
)

psi_range <- c(min(log_likelihood_vals$psi) - step_size, max(conf_ints$psi[2], conf_ints$Integrated[2]) + step_size)

y_min <- pseudo_log_likelihood_curves[[1]](conf_ints$Integrated[1]) - 0.5

label_data <- MLE_data |>
  add_row(Source = "Truth",
          MLE = psi_0,
          MLE_label = "psi[0]")

base_plot +
  stat_fn_IL +
  stat_fn_PL +
  geom_hline(yintercept = 0,
             linetype = 5) +
  geom_hline(aes(yintercept = -crit),
             linewidth = 1.5) +
  geom_text(aes(1.5, -crit, label = "95% CI Line", vjust = 1.3)) +
  geom_vline(aes(xintercept = MLE,
                 color = Source),
             data = label_data,
             show.legend = FALSE) +
  ggrepel::geom_label_repel(aes(x = MLE,
                                y = y_min / 2,
                                label = MLE_label,
                                color = Source),
                            data = label_data,
                            direction = "y",
                            parse = TRUE,
                            show.legend = FALSE) +
  geom_point(aes(x = psi, y = loglikelihood - max(loglikelihood, na.rm = TRUE)),
             data = log_likelihood_vals |> 
               filter(Pseudolikelihood == "Integrated"),
             size = 1) +
  geom_point(aes(x = psi, y = loglikelihood - max(loglikelihood, na.rm = TRUE)),
             data = log_likelihood_vals |> 
               filter(Pseudolikelihood == "Profile"),
             size = 1) +
  ylab("Log-Likelihood") +
  scale_x_continuous(expand = c(0, 0),
                     limits = psi_range) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(y_min, 0.1)) +
  scale_color_brewer(palette = "Set1") +
  xlab(expression(psi))
          





