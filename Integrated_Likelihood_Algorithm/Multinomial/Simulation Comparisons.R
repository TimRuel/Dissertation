library(kableExtra)
library(tidyverse)
library(purrr)
library(gdata)

desert_rodents_sim_results_PL <- readRDS("desert_rodents_sim_results_PL.Rda") |> 
  setNames("Desert Rodents (Profile)") 

desert_rodents_sim_results_IL <- readRDS("desert_rodents_sim_results_IL.Rda") |> 
  setNames("Desert Rodents (Integrated)")

birds_in_balrath_woods_sim_results_PL <- readRDS("birds_in_balrath_woods_sim_results_PL.Rda") |> 
  setNames("Birds in Balrath Woods (Profile)")

birds_in_balrath_woods_sim_results_IL <- data.frame(rep(NA, 5)) |> 
  setNames("Birds in Balrath Woods (Integrated)") |> 
  mutate(Metric = c("Bias", "SD", "RMSE", "Coverage", "Length")) |> 
  column_to_rownames("Metric")

birds_in_killarney_woodlands_sim_results_PL <- readRDS("birds_in_killarney_woodlands_sim_results_PL.Rda") |> 
  setNames("Birds in Killarney Woodlands (Profile)")

birds_in_killarney_woodlands_sim_results_IL <- readRDS("birds_in_killarney_woodlands_sim_results_IL.Rda") |> 
  setNames("Birds in Killarney Woodlands (Integrated)")

ruel_sim_results_df <- list(desert_rodents_sim_results_IL,
                            desert_rodents_sim_results_PL, 
                            birds_in_balrath_woods_sim_results_IL,
                            birds_in_balrath_woods_sim_results_PL, 
                            birds_in_killarney_woodlands_sim_results_IL,
                            birds_in_killarney_woodlands_sim_results_PL) |>
  map(\(x) rownames_to_column(x, "Metric")) |> 
  reduce(join, by = "Metric") |> 
  add_column(Author = "Ruel", 
             .after = "Metric")

severini_sim_results_df <- data.frame(c("Bias", "SD", "RMSE", "Coverage", "Length"),
                                      c(-0.018, 0.151, 0.153, 0.958, 0.547),
                                      c(-0.115, 0.156, 0.194, 0.851, 0.522),
                                      c(-0.058, 0.102, 0.117, 0.930, 0.403), 
                                      c(-0.197, 0.133, 0.238, 0.598, 0.472),
                                      c(0.0497, 0.0902, 0.103, 0.902, 0.342), 
                                      c(-0.0503, 0.0927, 0.105, 0.851, 0.352)) |>
  setNames(c("Metric",
             "Desert Rodents (Integrated)",
             "Desert Rodents (Profile)",
             "Birds in Balrath Woods (Integrated)",
             "Birds in Balrath Woods (Profile)",
             "Birds in Killarney Woodlands (Integrated)",
             "Birds in Killarney Woodlands (Profile)")) |> 
  add_column(Author = "Severini", 
             .after = "Metric")

sim_results_df1 <- ruel_sim_results_df |> 
  interleave(severini_sim_results_df) 

sim_results_df1 |> 
  kbl(col.names = c("Metric", "Author", rep(c("Integrated", "Profile"), 3)), 
      align = "c",
      row.names = FALSE) |> 
  kable_styling(bootstrap_options = c("striped", "hover"),
                full_width = FALSE) |> 
  add_header_above(c(" " = 2, 
                     "Desert Rodents" = 2, 
                     "Birds in Balrath Woods" = 2, 
                     "Birds in Killarney Woodlands" = 2)) |> 
  column_spec(1, bold = TRUE) |> 
  collapse_rows(columns = 1, valign = "top")

sim_results_df2 <- ruel_sim_results_df |> 
  rbind(severini_sim_results_df) 

sim_results_df2 |> 
  select(-Author) |> 
  kbl(col.names = c("Metric", rep(c("Integrated", "Profile"), 3)), 
      align = "c",
      row.names = FALSE) |> 
  kable_styling(bootstrap_options = c("striped", "hover"),
                full_width = FALSE) |> 
  add_header_above(c(" " = 1, 
                     "Desert Rodents" = 2, 
                     "Birds in Balrath Woods" = 2, 
                     "Birds in Killarney Woodlands" = 2)) |> 
  column_spec(1, bold = TRUE) |> 
  pack_rows(index = c("Ruel" = 5, "Severini" = 5),
            label_row_css = "background-color: #666; color: #fff;")



