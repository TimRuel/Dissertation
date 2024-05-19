library(dplyr)
library(purrr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../../utils.R")

multinomial_entropy_sims_file_path <- file.choose()

multinomial_entropy_sims <- readRDS(multinomial_entropy_sims_file_path)

population <- multinomial_entropy_sims |>  
  str_remove("^.*\\\\") |> 
  str_remove("_R.*$") |> 
  str_replace_all("_", " ") |> 
  tools::toTitleCase()

step_size <- log_likelihood_vals_file_path |>  
  str_remove("^.*step_size=") |> 
  str_remove(".Rda") |> 
  as.numeric()

switch(population,
       
       "Desert Rodents" = {
         
         data <- c(1, 1, 2, 4, 7, 10)
         
         step_size <- 0.01
         
         multinomial_entropy_sims_file_path <- "desert_rodents_sims.Rda"
         
         multinomial_entropy_sim_results_file_path <- "desert_rodents_sim_results.Rda"
         },
       
       "Birds in Balrath Woods" = {
         
         data <- c(1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 6, 8)
         
         multinomial_entropy_sim_results_file_path <- "birds_in_balrath_woods_sim_results.Rda"
         },
       
       "Birds in Killarney Woodlands" = {
         
         data <- c(1, 3, 4, 6, 7, 10, 14, 30)
         
         multinomial_entropy_sim_results_file_path <- "birds_in_killarney_woodlands_sim_results.Rda"
         }
       )  

psi_grid <- data |> 
  length() |> 
  log() |> 
  plyr::round_any(step_size, floor) |> 
  seq(0, to = _, step_size)

multinomial_entropy_sims <- readRDS(paste0("../Data/Simulations/", multinomial_entropy_sims_file_path))

spline_fitted_models_list <- multinomial_entropy_sims |> 
  map(
    function(sims) {
      sims |> 
        data.frame() |> 
        mutate(psi = psi_grid) |> 
        tidyr::pivot_longer(cols = -psi,
                            names_to = "Iteration",
                            values_to = "loglikelihood") |>
        group_by(Iteration) |>
        filter(is.finite(loglikelihood)) |>
        group_map(~ smooth.spline(.x$psi, .x$loglikelihood))
      }
    )

MLE_data_list <- spline_fitted_models_list |> 
  map(
    function(spline_fitted_models) {
      spline_fitted_models |> 
        sapply(
          function(mod) {
            optimize(
              function(psi) predict(mod, psi)$y, 
              lower = 0, 
              upper = log(length(data)), 
              maximum = TRUE)
            }
          ) |> 
        t() |> 
        data.frame() |> 
        mutate(MLE = as.numeric(maximum),
               Maximum = as.numeric(objective)) |> 
        select(MLE, Maximum)
      }
    )

pseudo_log_likelihood_curves_list <- spline_fitted_models_list |> 
  map2(MLE_data_list,
       function(spline_fitted_models, MLE_data) {
         spline_fitted_models |>
           map2(MLE_data$Maximum,
                function(mod, maximum) function(psi) predict(mod, psi)$y - maximum)
         }
       )

crit <- qchisq(0.95, 1) / 2

conf_ints_list <- pseudo_log_likelihood_curves_list |> 
  map2(MLE_data_list,
       function(pseudo_log_likelihood_curves, MLE_data) {
         pseudo_log_likelihood_curves |> 
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
       }
  ) |> 
  map(
    function(conf_ints) {
      conf_ints |> 
        do.call(rbind, args = _) |> 
        as.data.frame() |> 
        setNames(c("Lower", "Upper"))
    }
  )

psi_0 <- entropy(data / sum(data))

multinomial_entropy_sim_results <- MLE_data_list |> 
  map2(conf_ints_list,
       function(MLE_data, conf_ints) {
         
         MLE_data |> 
           mutate(Lower = conf_ints$Lower, 
                  Upper = conf_ints$Upper,
                  psi_0 = psi_0,
                  Covered = psi_0 |> between(Lower, Upper)) |> 
           summarise(Bias = mean(MLE - psi_0),
                     SD = sd(MLE),
                     RMSE = sqrt(mean((MLE - psi_0)^2)),
                     Coverage = mean(Covered, na.rm = TRUE),
                     Length = mean(Upper - Lower, na.rm = TRUE))
         }
       ) |> 
  do.call(rbind, args = _) |> 
  t() |> 
  as.data.frame() |> 
  rownames_to_column("Metric")

saveRDS(multinomial_entropy_sim_results, paste0("Results/", multinomial_entropy_sim_results_file_path))












