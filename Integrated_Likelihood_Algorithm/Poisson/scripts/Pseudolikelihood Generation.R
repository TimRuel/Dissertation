setwd("C:/Northwestern/Dissertation/Integrated_Likelihood_Algorithm/Poisson")
source("scripts/utils.R")

population <- "A"

cfg <- load_config(population)

set.seed(cfg$seed)

theta_0 <- get_theta_0(cfg)
Y <- get_Y(theta_0)
weights <- get_weights(cfg)

cfg$theta_0$values <- theta_0
cfg$data$Y <- Y
cfg$weights$values <- weights
cfg$lambda_grid$interval$stop <- min(1 / weights)

yaml::write_yaml(cfg, file.path(cfg$result_dir, "config.yml"))

integrated_LL <- get_integrated_LL(cfg)

profile_LL <- get_profile_LL(cfg)

PLL_df_merged <- merge(integrated_LL$log_L_bar$df, profile_LL$log_L_p_df, all = TRUE) 
  
PLL_df_merged_filepath <- file.path(cfg$result_dir, "PLL_df_merged.rda")

saveRDS(PLL_df_merged, PLL_df_merged_filepath)

