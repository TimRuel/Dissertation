source("scripts/helpers/config_utils.R")
source("scripts/helpers/truth_utils.R")

experiment_id <- "experiment_A"
seed <- 7835L
set.seed(seed)

X1_levels <- list(
  
  A = list(
    X2 = list(dist = list(name = "rbeta", 
                          params = c(2, 5)),
              support = c(0, 1)),
    ref_level = TRUE,
    level_of_interest = TRUE
  ),
  
  B = list(
    X2 = list(dist = list(name = "rbeta", 
                          params = c(3, 3)),
              support = c(0, 1)),
    ref_level = FALSE,
    level_of_interest = FALSE
  ),
  
  C = list(
    X2 = list(dist = list(name = "rbeta", 
                          params = c(5, 2)),
              support = c(0, 1)),
    ref_level = FALSE,
    level_of_interest = FALSE
  )
)

true_param_specs <- list(J = 6L,
                         entropy_range_specs = list(offset = c(0.3, 0.2), padding = 0.75),
                         num_X2_test_vals = 1000L)
true_param_specs$p <- get_num_predictors(names(X1_levels), true_param_specs$J)

experiment_parameters <- get_experiment_parameters(X1_levels, true_param_specs)

experiment_config <- list(
  experiment_id = experiment_id,
  seed = seed,
  X1_levels = experiment_parameters$X1_levels,
  true_param_specs = true_param_specs,
  Beta_0_path = file.path("results", experiment_id, "Beta_0.rds")
)

dir.create(dirname(experiment_config$Beta_0_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(experiment_parameters$Beta_0, experiment_config$Beta_0_path)
yaml::write_yaml(experiment_config, file.path("config", paste0(experiment_id, ".yml")))