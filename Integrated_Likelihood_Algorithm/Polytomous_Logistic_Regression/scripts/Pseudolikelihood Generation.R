setwd("C:/Northwestern/Dissertation/Integrated_Likelihood_Algorithm/Polytomous_Logistic_Regression")
source("scripts/utils.R")

population <- "A"

cfg <- load_config(population)

set.seed(cfg$seed)

# Population Parameters -----------------------------------------------

entropy_ranges <- get_entropy_ranges(cfg)

cfg$X1$levels <- purrr::imap(cfg$X1$levels, function(level, X1_level) {
  level$X2$specs$entropy_range <- entropy_ranges[[X1_level]]
  level
})

Beta_0 <- get_Beta_0(cfg)

cfg$Beta_0$values <- Beta_0

# Data Parameters ---------------------------------------------------------



X2_arg_list <- list(X2_intervals, X2_shape_params, m) |>
  unlist(recursive = FALSE) |>
  split(X1_levels) |>
  map(unname) |>
  map(~ setNames(., c("interval", "shape", "n_samples")))

# Pseudolikelihood Parameters --------------------------------------------

step_size <- 0.01

num_std_errors <- 3.5

init_guess_sd <- 5

# num_workers <- parallel::detectCores() |>
#   as.integer()

num_workers <- Sys.getenv("MC_CORES") |>
  as.integer()

# num_workers <- 12

chunk_size <- 1

IL_maxtime <- 10

PL_maxtime <- 10

# Data Generation ---------------------------------------------------------

X1 <- get_X1(cfg)

X2 <- get_X2(cfg)

cfg$X1$levels <- purrr::imap(cfg$X1$levels, function(level, X1_level) {
  level$X2$values <- X2[names(X2) == X1_level] |> 
    unname()
  level
})

X_design <- model.matrix(~ X1*X2 - 1)

Y_probs <- get_Y_probs(cfg)

cfg$data$Y$true_probs <- Y_probs

Y <- Y_probs |>
  select(-c(X1, X2)) |>
  apply(1, \(prob) sample(1:J, size = 1, prob = prob)) |>
  unlist() |>
  unname() |>
  factor(levels = 1:J)

Y_design <- model.matrix(~ Y)[,-1]

data <- data.frame(X1 = X1,
                   X2 = X2,
                   Y)

formula <- Y ~ .^2 - 1

model <- fit_multinomial_logistic_model(data, formula)

Beta_MLE <- get_Beta_MLE(model)

threshold <- ceiling(abs(log_likelihood(Beta_MLE, X_design, model.matrix(~ Y)[,-1]))) + 20

data_filepath <- glue::glue("sim_dfs/Sim{i}.Rda")

saveRDS(data, data_filepath)


################################################################################
########################## INTEGRATED LIKELIHOOD - VANILLA MC ##################
################################################################################

threshold <- ceiling(abs(log_likelihood(Beta_MLE, X_design, model.matrix(~ Y)[,-1]))) + 20

init_guess_sd <- 5

# num_workers <- Sys.getenv("SLURM_NPROCS") |>
#   as.numeric()

# num_workers <- availableCores() |>
#   as.numeric()

num_workers <- parallel::detectCores() |>
  as.integer()

chunk_size <- 1

num_branches <- num_workers * chunk_size

IL_maxtime <- 10

tic()

log_integrated_likelihood <- get_log_integrated_likelihood(data,
                                                           formula,
                                                           h,
                                                           step_size,
                                                           num_std_errors,
                                                           init_guess_sd,
                                                           threshold,
                                                           num_workers,
                                                           chunk_size,
                                                           IL_maxtime)

toc()

log_integrated_likelihood_filepath <- glue::glue("IL_objects/R={num_branches}_J={J}_h={h}_stepsize={step_size}_numse={num_std_errors}.Rda")

saveRDS(log_integrated_likelihood, log_integrated_likelihood_filepath)

# log_integrated_likelihood <- readRDS("IL_objects/R=260_h=b_stepsize=0.02.Rda")

################################################################################
############################## PROFILE LIKELIHOOD ############################## 
################################################################################

PL_maxtime <- 10

tic()

log_profile_likelihood <- get_log_profile_likelihood(data,
                                                     formula,
                                                     h,
                                                     step_size,
                                                     num_std_errors,
                                                     init_guess_sd,
                                                     PL_maxtime)

toc()

log_profile_likelihood_filepath <- glue::glue("PL_objects/R={num_branches}_J={J}_h={h}_stepsize={step_size}_numse={num_std_errors}.Rda")

saveRDS(log_profile_likelihood, log_profile_likelihood_filepath)

# log_profile_likelihood <- readRDS("PL_objects/R=260_h=b_stepsize=0.02.Rda")

################################################################################
################################### STORAGE #################################### 
################################################################################

log_integrated_likelihood <- readRDS("Integrated_Likelihood_Algorithm/Polytomous Logistic Regression/Simulations/Marginal Entropy/IL_objects/h=A/Sim1.Rda")

get_log_L_bar(log_integrated_likelihood)

log_likelihood_vals <- log_integrated_likelihood$log_L_bar$df |> 
  merge(log_profile_likelihood, all = TRUE)

# log_likelihood_vals <- log_L_bar_df |> 
#   merge(log_profile_likelihood, all = TRUE)
# 
# log_likelihood_vals_file_path <- glue::glue("log_likelihood_vals/test.Rda")

log_likelihood_vals_file_path <- glue::glue("log_likelihood_vals/R={num_branches}_J={J}_h={h}_stepsize={step_size}_numse={num_std_errors}.Rda")

saveRDS(log_likelihood_vals, log_likelihood_vals_file_path)



for (df in log_integrated_likelihood$log_L_tilde_df) plot(df)








