library(dplyr)
library(purrr)
library(stringr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../../utils.R")

multinomial_entropy_sims_IL_file_path <- file.choose()

multinomial_entropy_sims_IL <- readRDS(multinomial_entropy_sims_IL_file_path)

multinomial_entropy_sims_mod_IL_file_path <- file.choose()

multinomial_entropy_sims_mod_IL <- readRDS(multinomial_entropy_sims_mod_IL_file_path)

multinomial_entropy_sims_PL_file_path <- file.choose()

multinomial_entropy_sims_PL <- readRDS(multinomial_entropy_sims_PL_file_path)

# population <- multinomial_entropy_sims_IL_file_path |>
#   str_remove("^.*/") |>
#   str_remove("_IL.*$") |>
#   str_replace_all("_", " ") |>
#   tools::toTitleCase()

population <- gsub(".*/Simulations/(.*)/Integrated Likelihood/.*", "\\1", multinomial_entropy_sims_IL_file_path)

switch(population,
       
       "Desert Rodents" = {
         
         data <- c(1, 1, 2, 4, 7, 10)
         
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

n <- sum(data)

################################################################################
############################ INTEGRATED LIKELIHOOD ############################# 
################################################################################

seed_IL <- multinomial_entropy_sims_IL_file_path |>  
  str_remove("^.*seed=") |> 
  str_extract("\\d+") |> 
  as.numeric()

set.seed(seed_IL)

num_sims_IL <- multinomial_entropy_sims_IL_file_path |>  
  str_remove("^.*numsims=") |> 
  str_extract("\\d+") |> 
  as.numeric()

data_sims_IL <- num_sims_IL |> 
  rmultinom(n, data) |> 
  data.frame() |> 
  as.list() |> 
  map(as.numeric)

step_size_IL <- multinomial_entropy_sims_IL_file_path |>  
  str_remove("^.*stepsize=") |> 
  str_extract("\\d+\\.\\d+") |> 
  as.numeric()

num_std_errors_IL <- multinomial_entropy_sims_IL_file_path |>  
  str_remove("^.*numse=") |> 
  str_extract("\\d+") |> 
  as.numeric()

psi_grid_list_IL <- data_sims_IL |> 
  map(\(data) get_psi_grid(data, step_size_IL, num_std_errors_IL))

################################################################################
######################## MODIFIED INTEGRATED LIKELIHOOD ########################
################################################################################

seed_mod_IL <- multinomial_entropy_sims_mod_IL_file_path |>  
  str_remove("^.*seed=") |> 
  str_extract("\\d+") |> 
  as.numeric()

set.seed(seed_mod_IL)

num_sims_mod_IL <- multinomial_entropy_sims_mod_IL_file_path |>  
  str_remove("^.*numsims=") |> 
  str_extract("\\d+") |> 
  as.numeric()

data_sims_mod_IL <- num_sims_mod_IL |> 
  rmultinom(n, data) |> 
  data.frame() |> 
  as.list() |> 
  map(as.numeric)

step_size_mod_IL <- multinomial_entropy_sims_mod_IL_file_path |>  
  str_remove("^.*stepsize=") |> 
  str_extract("\\d+\\.\\d+") |> 
  as.numeric()

num_std_errors_mod_IL <- multinomial_entropy_sims_mod_IL_file_path |>  
  str_remove("^.*numse=") |> 
  str_extract("\\d+") |> 
  as.numeric()

psi_grid_list_mod_IL <- data_sims_mod_IL |> 
  map(\(data) get_psi_grid(data, step_size_mod_IL, num_std_errors_mod_IL))

################################################################################
############################## PROFILE LIKELIHOOD ############################## 
################################################################################

seed_PL <- multinomial_entropy_sims_PL_file_path |>  
  str_remove("^.*seed=") |> 
  str_extract("\\d+") |> 
  as.numeric()

set.seed(seed_PL)

num_sims_PL <- multinomial_entropy_sims_PL_file_path |>  
  str_remove("^.*numsims=") |> 
  str_extract("\\d+") |> 
  as.numeric()

data_sims_PL <- num_sims_PL |> 
  rmultinom(n, data) |> 
  data.frame() |> 
  as.list() |> 
  map(as.numeric)

step_size_PL <- multinomial_entropy_sims_PL_file_path |>  
  str_remove("^.*stepsize=") |> 
  str_extract("\\d+\\.\\d+") |> 
  as.numeric()

num_std_errors_PL <- multinomial_entropy_sims_PL_file_path |>  
  str_remove("^.*numse=") |> 
  str_extract("\\d+") |> 
  as.numeric()

psi_grid_list_PL <- data_sims_PL |> 
  map(\(data) get_psi_grid(data, step_size_PL, num_std_errors_PL))

################################################################################
################################## SYNTHESIS ################################### 
################################################################################

multinomial_entropy_sims_lists <- list("Integrated" = multinomial_entropy_sims_IL,
                                       "Mod_Integrated" = multinomial_entropy_sims_mod_IL,
                                       "Profile" = multinomial_entropy_sims_PL)

psi_grid_lists <- list("Integrated" = psi_grid_list_IL,
                       "Mod_Integrated" = psi_grid_list_mod_IL,
                       "Profile" = psi_grid_list_PL)

spline_fitted_models_list <- multinomial_entropy_sims_lists |> 
  map2(psi_grid_lists,
       \(sims_list, psi_grid_list) {
         sims_list |> 
           map2(psi_grid_list,
                \(sims, psi_grid) {
                  sims |> 
                    data.frame() |> 
                    mutate(psi = psi_grid) |> 
                    tidyr::pivot_longer(cols = -psi,
                                        names_to = "Iteration",
                                        values_to = "loglikelihood") |>
                    group_by(Iteration) |>
                    # filter(is.finite(loglikelihood)) |>
                    group_map(~ smooth.spline(.x$psi, .x$loglikelihood))
                  }
           )
         }
       )

MLE_data_list <- spline_fitted_models_list |> 
  map(
    function(spline_fitted_models) {
      spline_fitted_models |> 
        sapply(
          function(mod) {
            optimize(
              function(psi) predict(mod[[1]], psi)$y, 
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
                function(mod, maximum) function(psi) predict(mod[[1]], psi)$y - maximum)
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
                  
                  interval = c(lower_bound, upper_bound)
                  
                  # if (upper_bound < 1) return(NA)
                  
                  return(interval)
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
           drop_na() |> 
           summarise(Bias = mean(MLE - psi_0),
                     SD = sqrt(mean((MLE - mean(MLE))^2)),
                     RMSE = sqrt(mean((MLE - psi_0)^2)),
                     Coverage = mean(Covered, na.rm = TRUE),
                     Length = mean(Upper - Lower, na.rm = TRUE))
         }
       ) |> 
  do.call(rbind, args = _) |> 
  t() |> 
  as.data.frame() |> 
  tibble::rownames_to_column("Metric")

saveRDS(multinomial_entropy_sim_results, paste0("Results/", multinomial_entropy_sim_results_file_path))












