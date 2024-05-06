library(dplyr)
library(purrr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

population <- "Desert Rodents"
# population <- "Birds in Balrath Woods"
# population <- "Birds in Killarney Woodlands"

switch(population,
       
       "Desert Rodents" = {
         
         data <- c(1, 1, 2, 4, 7, 10)
         
         step_size <- 0.01
         
         multinomial_entropy_sims_IL_file_path <- "desert_rodents_IL_sims.Rda"
         
         multinomial_entropy_sims_PL_file_path <- "desert_rodents_PL_sims.Rda"
         },
       
       "Birds in Balrath Woods" = {
         
         data <- c(1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 6, 8)
         
         step_size <- 0.05
         
         multinomial_entropy_sims_IL_file_path <- "birds_in_balrath_woods_IL_sims.Rda"
         
         multinomial_entropy_sims_PL_file_path <- "birds_in_balrath_woods_PL_sims.Rda"
         },
       
       "Birds in Killarney Woodlands" = {
         
         data <- c(1, 3, 4, 6, 7, 10, 14, 30)
         
         step_size <- 0.01
         
         multinomial_entropy_sims_IL_file_path <- "birds_in_killarney_woodlands_IL_sims.Rda"
         
         multinomial_entropy_sims_PL_file_path <- "birds_in_killarney_woodlands_PL_sims.Rda"
         }
       )  

psi_grid <- data |> 
  length() |> 
  log() |> 
  plyr::round_any(step_size, ceiling) |> 
  seq(0, to = _, step_size)

multinomial_entropy_sims_IL <- readRDS(paste0("../Data/Simulations/", multinomial_entropy_sims_IL_file_path))

multinomial_entropy_sims_PL <- readRDS(paste0("../Data/Simulations/", multinomial_entropy_sims_PL_file_path))

multinomial_entropy_sims <- list("Integrated" = multinomial_entropy_sims_IL, "Profile" = multinomial_entropy_sims_PL)

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
              upper = psi_grid |> tail(1), 
              maximum = TRUE)
            }
          ) |> 
        t() |> 
        data.frame() |> 
        dplyr::rename(MLE = maximum,
               Maximum = objective)
      }
    )

curves_PL <- mapply(
  function(mod, maximum) function(psi) predict(mod, psi)$y - maximum,
  mods_PL,
  MLE_data_PL$Maximum)

curves_IL <- mapply(
  function(mod, maximum) function(psi) predict(mod, psi)$y - maximum,
  mods_IL,
  MLE_data_IL$Maximum)

crit <- qchisq(0.95, 1) / 2

CI_lower_PL <- mapply(function(curve, MLE) {
  tryCatch(
    uniroot(function(psi) curve(psi) + crit,
            interval = c(psi_grid |> head(1), MLE))$root |> round(3),
    error = function(e) return(NA)
  )
},
curves_PL,
MLE_data_PL$MLE)

CI_upper_PL <- mapply(function(curve, MLE) {
  tryCatch(
    uniroot(function(psi) curve(psi) + crit,
            interval = c(MLE, psi_grid |> tail(1)))$root |> round(3),
    error = function(e) return(NA)
  )
},
curves_PL,
MLE_data_PL$MLE)

psi_0 <- entropy(data / sum(data))

sim_results_PL <- MLE_data_PL |> 
  mutate(Lower = CI_lower_PL, 
         Upper = CI_upper_PL,
         psi_0 = psi_0,
         Covered = psi_0 |> between(Lower, Upper)) |> 
  summarise(Bias = mean(MLE - psi_0),
            SD = sd(MLE),
            RMSE = sqrt(mean((MLE - psi_0)^2)),
            Coverage = mean(Covered, na.rm = TRUE),
            Length = mean(Upper - Lower, na.rm = TRUE)) |> 
  t() |> 
  as.data.frame() |> 
  round(3) |> 
  setNames("Profile")

# saveRDS(sim_results_PL, "desert_rodents_sim_results_PL.Rda")
# sim_results_PL <- readRDS("desert_rodents_sim_results_PL.Rda")

saveRDS(sim_results_PL, "birds_in_balrath_woods_sim_results_PL.Rda")
# sim_results_PL <- readRDS("birds_in_balrath_woods_sim_results_PL.Rda")

# saveRDS(sim_results_PL, "birds_in_killarney_woodlands_sim_results_PL.Rda")
# sim_results_PL <- readRDS("birds_in_killarney_woodlands_sim_results_PL.Rda")

sim_results_PL

CI_lower_IL <- mapply(function(curve, MLE) {
  tryCatch(
    uniroot(function(psi) curve(psi) + crit,
            interval = c(psi_grid |> head(1), MLE))$root |> round(3),
    error = function(e) return(NA)
  )
},
curves_IL,
MLE_data_IL$MLE)

CI_upper_IL <- mapply(function(curve, MLE) {
  tryCatch(
    uniroot(function(psi) curve(psi) + crit,
            interval = c(MLE, psi_grid |> tail(1)))$root |> round(3),
    error = function(e) return(NA)
  )
},
curves_IL,
MLE_data_IL$MLE)

sim_results_IL <- MLE_data_IL |> 
  mutate(Lower = CI_lower_IL, 
         Upper = CI_upper_IL,
         psi_0 = psi_0,
         Covered = psi_0 |> between(Lower, Upper)) |> 
  summarise(Bias = mean(MLE - psi_0),
            SD = sd(MLE),
            RMSE = sqrt(mean((MLE - psi_0)^2)),
            Coverage = mean(Covered, na.rm = TRUE),
            Length = mean(Upper - Lower, na.rm = TRUE)) |> 
  t() |> 
  as.data.frame() |> 
  round(3) |> 
  setNames("Integrated")

# saveRDS(sim_results_IL, "desert_rodents_sim_results_IL.Rda")
# sim_results_IL <- readRDS("desert_rodents_sim_results_IL.Rda")

# saveRDS(sim_results_IL, "birds_in_balrath_woods_sim_results_IL.Rda")
# sim_results_IL <- readRDS("birds_in_balrath_woods_sim_results_IL.Rda")

saveRDS(sim_results_IL, "birds_in_killarney_woodlands_sim_results_IL.Rda")
sim_results_IL <- readRDS("birds_in_killarney_woodlands_sim_results_IL.Rda")

sim_results_IL










