library(tidyverse)
library(purrr)
library(geomtextpath)
library(viridis)
library(ggnewscale)
library(zeallot)
library(kableExtra)
library(stringr)
library(rstudioapi)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("utils.R")

print("Select population directory")
population_directory <- selectDirectory(path = getwd())

setwd(population_directory)

source("data.R")

print("Choose your pseudolikelihood data file.")
log_likelihood_vals_file_path <- selectFile(path = getwd())

model <- get_multinomial_logistic_model(data)

step_size <- log_likelihood_vals_file_path |>  
  str_remove("^.*stepsize=") |> 
  str_extract("\\d+\\.\\d+") |> 
  as.numeric()

X_level <- log_likelihood_vals_file_path |>
  str_remove("^.*h=") |>
  str_extract("-?\\d+(\\.\\d+)?") |>
  as.numeric()

X_h <- data.frame(X = factor(X_level))

psi_hat <- get_psi_hat(model, X_h)

psi_0 <- theta_0[[X_level]] |> 
  get_entropy()

# pseudolikelihood_names <- c("Integrated", "Mod_Integrated", "Profile")

pseudolikelihood_names <- c("Integrated", "Profile")

# pseudolikelihood_names <- "Profile"

log_likelihood_vals <- readRDS(log_likelihood_vals_file_path) |>
  tidyr::pivot_longer(cols = all_of(pseudolikelihood_names),
                      names_to = "Pseudolikelihood",
                      values_to = "loglikelihood") |>
  mutate(Pseudolikelihood = Pseudolikelihood |>
           as_factor() |>
           fct_inorder())

# log_likelihood_vals <- log_likelihood_vals |>
#   tidyr::pivot_longer(cols = all_of(pseudolikelihood_names),
#                       names_to = "Pseudolikelihood",
#                       values_to = "loglikelihood") |>
#   mutate(Pseudolikelihood = Pseudolikelihood |>
#            as_factor() |>
#            fct_inorder())

spline_fitted_models <- log_likelihood_vals |>
  drop_na(loglikelihood) |> 
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
  # mutate(MLE_label = c("hat(psi)[PL]")) |>
  mutate(MLE_label = c("hat(psi)[IL]", "hat(psi)[PL]")) |>
  # mutate(MLE_label = c("hat(psi)[IL]", "hat(psi)[m-IL]", "hat(psi)[PL]")) |>
  rownames_to_column("Source")

pseudo_log_likelihood_curves <- spline_fitted_models |> 
  map2(MLE_data$Maximum,
       function(mod, maximum) function(psi) predict(mod, psi)$y - maximum)

log_likelihood_vals |> 
  # filter(Pseudolikelihood == "Integrated") |> 
  ggplot() +
  geom_point(aes(x = psi, y = loglikelihood, color = Pseudolikelihood),
             size = 1) +
  ylab("Log-Likelihood") +
  scale_color_brewer(palette = "Set1") +
  xlab(expression(psi)) +
  # scale_x_continuous(expand = c(0, 0),
  #                    # limits = c(31.5, 32.5)) +
  #                    limits = c(15, 16)) +
  # scale_y_continuous(expand = c(0, 0),
  #                    # limits = c(-0.025, 0)) +
  #                    limits = c(0, 1)) +
  theme_minimal() +
  theme(axis.line = element_line())

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

# psi_hat <- get_psi_hat(model, X_h)
# 
# conf_ints$Classical <- get_CI_psi_hat(model, X_h, alpha) |> 
#   round(3)
# 
# MLE_data <- MLE_data |>
#   add_row(Source = "Classical",
#           MLE = psi_hat,
#           MLE_label = "hat(psi)")

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
  add_row(Source = "Classical",
          MLE = psi_hat) |>
  add_row(Source = "Truth",
          MLE = psi_0) |>
  kbl(col.names = c("Source", 
                    "MLE",
                    "95% Confidence Interval",
                    "CI Length"),   
      align = "c") |> 
  kable_styling(bootstrap_options = c("striped", "hover")) |> 
  row_spec(4, color = "green", bold = TRUE)

# stat_fn_PL <- map2(
c(stat_fn_IL, stat_fn_PL) %<-% map2(
  # c(stat_fn_IL, stat_fn_mod_IL, stat_fn_PL) %<-% map2(
  pseudo_log_likelihood_curves,
  pseudolikelihood_names,
  function(curve, pseudolikelihood_name) {
    
    stat_fn <- stat_function(fun = curve,
                             geom = "line",
                             # label = pseudolikelihood_name,
                             aes(color = pseudolikelihood_name),
                             # linewidth = pseudolikelihood_name),
                             linewidth = 1.5,
                             # hjust = 0.1,
                             show.legend = FALSE,
                             xlim = c(psi_grid |> head(1), psi_grid |> tail(1)))
    return(stat_fn)
  }
)

psi_range <- c(min(psi_grid), max(psi_grid))

y_min <- pseudo_log_likelihood_curves |> 
  map((\(curve) c(curve(psi_range[1]), curve(psi_range[2])))) |> 
  unlist() |> 
  min() |> 
  round()

label_data <- MLE_data |>
  add_row(Source = "Classical",
          MLE = psi_hat,
          MLE_label = "hat(psi)") |> 
  add_row(Source = "Truth",
          MLE = psi_0,
          MLE_label = "psi[0]")

ggplot() +
  stat_fn_IL +
  # stat_fn_mod_IL +
  stat_fn_PL +
  geom_hline(yintercept = 0,
             linetype = 5) +
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
  annotate("rect",
           xmin = conf_ints$Integrated[1],
           xmax = min(conf_ints$Integrated[2], psi_range[2]),
           ymin = -Inf,
           ymax = Inf,
           fill = "blue",
           alpha = 0.5) +
  annotate("rect",
           xmin = conf_ints$Profile[1],
           xmax = conf_ints$Profile[2],
           ymin = -Inf,
           ymax = Inf,
           fill = "green",
           alpha = 0.5) +
  # annotate("rect",
  #          xmin = conf_ints$Classical[1],
  #          xmax = conf_ints$Classical[2],
  #          ymin = -Inf,
  #          ymax = Inf,
  #          fill = "red",
  #          alpha = 0.5) +
  ylab("Log-Likelihood") +
  scale_x_continuous(expand = c(0, 0),
                     limits = psi_range) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(y_min, 0.1)) +
  scale_color_brewer(palette = "Set1") +
  # scale_linewidth_manual(values = c(4, 2.5, 1)) +
  xlab(expression(psi)) +
  theme_minimal() +
  theme(axis.line = element_line())






