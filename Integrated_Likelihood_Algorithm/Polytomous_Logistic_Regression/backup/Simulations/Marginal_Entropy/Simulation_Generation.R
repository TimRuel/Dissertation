library(doFuture)

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("/gpfs/home/tbr0780/Dissertation/Integrated_Likelihood_Algorithm/Polytomous Logistic Regression/Simulations/Marginal Entropy")
source("../../utils.R")
seed <- 7835
set.seed(seed)

TASK_ID <- Sys.getenv("TASK_ID")
num_cores <- Sys.getenv("MC_CORES") |> 
  as.integer()

# Population Parameters -----------------------------------------------

J <- 6 # number of levels of response variable

X1_levels <- c("A", "B", "C") # Levels of categorical predictor

X1_ref_level <- X1_levels[1]

p <- ncol(model.matrix( ~ factor(X1_levels)[1] * J - 1)) # number of terms in model after dummy encoding

entropy_ranges <- get_entropy_ranges(J, X1_levels)

X2_intervals <- list("A" = c(0, 1),
                     "B" = c(0, 1),
                     "C" = c(0, 1))

X2_shape_params <- list("A" = c(2, 5),
                        "B" = c(3, 3),
                        "C" = c(5, 2))

num_vals <- 1000

Beta_0 <- get_Beta_0(X1_levels, p, J, X2_intervals, num_vals)

# Data Parameters ---------------------------------------------------------

m <- c(60, 60, 60)  # number of observations at each level of categorical predictor

names(m) <- X1_levels

C <- length(X1_levels) # number of levels of categorical predictor

n <- sum(m) # total number of observations

X1 <- X1_levels |>
  rep(times = m) |>
  factor() |>
  relevel(X1_ref_level)

X2_arg_list <- list(X2_intervals, X2_shape_params, m) |>
  unlist(recursive = FALSE) |>
  split(X1_levels) |>
  map(unname) |>
  map(~ setNames(., c("interval", "shape", "n_samples")))

# Pseudolikelihood Parameters ----------------------------------------

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

X2 <- generate_X2_samples(X2_arg_list)

X_design <- model.matrix(~ X1*X2 - 1)

Y_probs <- X_design %*% cbind(0, Beta_0) |>
  apply(1, softmax) |>
  t() |>
  data.frame() |>
  rename_with( ~ paste0("Y", 1:J)) |>
  mutate(X1 = rep(X1_levels, times = m),
         X2 = X2) |>
  select(X1, X2, everything())

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

# Pseudolikelihood Generation ---------------------------------------------

for (h in X1_levels) {

  # Integrated Likelihood Generation ----------------------------------------

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

  log_integrated_likelihood_filepath <- glue::glue("IL_objects/h={h}/Sim{TASK_ID}.Rda")

  saveRDS(log_integrated_likelihood, log_integrated_likelihood_filepath)

  # Profile Likelihood Generation -------------------------------------------

  log_profile_likelihood <- get_log_profile_likelihood(data,
                                                       formula,
                                                       h,
                                                       step_size,
                                                       num_std_errors,
                                                       init_guess_sd,
                                                       PL_maxtime)

  log_profile_likelihood_filepath <- glue::glue("PL_objects/h={h}/Sim{TASK_ID}.Rda")

  saveRDS(log_profile_likelihood, log_profile_likelihood_filepath)

  # Storage ------------------------------------------------------------

  log_likelihood_df <- log_integrated_likelihood$log_L_bar$df |>
    merge(log_profile_likelihood, all = TRUE)

  log_likelihood_df_filepath <- glue::glue("log_likelihood_dfs/h={h}/Sim{TASK_ID}.Rda")

  saveRDS(log_likelihood_df, log_likelihood_df_filepath)
}

pushoverr::pushover(paste("Sim", TASK_ID, "done!"))

