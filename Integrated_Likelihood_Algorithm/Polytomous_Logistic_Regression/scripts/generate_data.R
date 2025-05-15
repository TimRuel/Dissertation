experiment_id <- "experiment_A"
seed <- 7835L
set.seed(seed)

X1_levels <- list(
  
  A = list(
    X2 = list(dist = list(name = "rbeta", 
                          params = c(2, 5)),
              support = c(0, 1)),
    m = 60,
    ref_level = TRUE,
    level_of_interest = TRUE
  ),
  
  B = list(
    X2 = list(dist = list(name = "rbeta", 
                          params = c(3, 3)),
              support = c(0, 1)),
    m = 60,
    ref_level = FALSE,
    level_of_interest = FALSE
  ),
  
  C = list(
    X2 = list(dist = list(name = "rbeta", 
                          params = c(5, 2)),
              support = c(0, 1)),
    m = 60,
    ref_level = FALSE,
    level_of_interest = FALSE
  )
)

ep_specs <- list(J = 6L,
                 entropy_range_specs = list(offset = c(0.3, 0.2), padding = 0.75),
                 reps = list(Beta_0 = 1000L, pY_0 = 1e6L),
                 formula = "~.^2 - 1")

experiment_parameters <- get_experiment_parameters(X1_levels, ep_specs)

X1_levels <- experiment_parameters$X1_levels

data <- get_data(X1_levels, ep_specs$formula, experiment_parameters$true_params$Beta_0)

plots <- get_plots(X1_levels, experiment_parameters$pY_0, data$Y_probs, data$model_df)

experiment_config <- list(
  experiment_id = experiment_id,
  seed = seed,
  X1_levels = X1_levels,
  specs = experiment_parameters$specs,
  true_params_dir = file.path("results", experiment_id, "true_params"),
  data_dir = file.path("results", experiment_id, "data"),
  plots_dir = file.path("results", experiment_id, "plots")
)

save_list_objects(experiment_parameters$true_params, experiment_config$true_params_dir)
save_list_objects(data, experiment_config$data_dir)
save_list_plots(plots, experiment_config$plots_dir)
write_yaml(experiment_config, file.path("config", paste0(experiment_id, ".yml")))
