library(kableExtra)
library(dplyr)
library(purrr)

desert_rodents_sim_results_PL <- readRDS("desert_rodents_sim_results_PL.Rda") |> 
  setNames("Desert Rodents (Ruel)") |> 
  mutate("Desert Rodents (Severini)" = c(-0.115, 0.156, 0.194, 0.851, 0.522))

birds_in_balrath_woods_sim_results_PL <- readRDS("birds_in_balrath_woods_sim_results_PL.Rda") |> 
  setNames("Birds in Balrath Woods (Ruel)") |> 
  mutate("Birds in Balrath Woods (Severini)" = c(-0.197, 0.133, 0.238, 0.598, 0.472))

birds_in_killarney_woodlands_sim_results_PL <- readRDS("birds_in_killarney_woodlands_sim_results_PL.Rda") |> 
  setNames("Birds in Killarney Woodlands (Ruel)") |> 
  mutate("Birds in Killarney Woodlands (Severini)" = formatC(c(-0.0503, 0.0927, 0.105, 0.851, 0.352), 3, format = "fg")) 

sim_results_df <- list(desert_rodents_sim_results_PL, 
                       birds_in_balrath_woods_sim_results_PL, 
                       birds_in_killarney_woodlands_sim_results_PL) |>
  map(\(x) rownames_to_column(x, "Metric")) |> 
  reduce(join, by = "Metric") |> 
  column_to_rownames("Metric")

sim_results_df |> 
  kbl(caption = "Profile Likelihood",
      col.names = rep(c("Ruel", "Severini"), 3), 
      align = "c") |> 
  kable_styling(bootstrap_options = c("striped", "hover"),
                full_width = FALSE) |> 
  add_header_above(c(" " = 1, 
                     "Desert Rodents" = 2, 
                     "Birds in Balrath Woods" = 2, 
                     "Birds in Killarney Woodlands" = 2)) |> 
  column_spec(1, bold = TRUE)




