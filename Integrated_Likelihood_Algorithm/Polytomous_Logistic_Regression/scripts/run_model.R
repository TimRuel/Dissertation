
threshold <- ceiling(abs(log_likelihood(Beta_MLE, X_design, model.matrix(~ Y)[,-1]))) + 20

init_guess_sd <- 5

chunk_size <- 1

num_branches <- num_workers * chunk_size

IL_maxtime <- 10

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

PL_maxtime <- 10

log_profile_likelihood <- get_log_profile_likelihood(data,
                                                     formula,
                                                     h,
                                                     step_size,
                                                     num_std_errors,
                                                     init_guess_sd,
                                                     PL_maxtime)

