library(nnet)

set.seed(123)
model_interaction <- multinom(Y ~ X1 * X2 - 1, data = data, trace = FALSE)
model_main_effects <- multinom(Y ~ X1 + X2 - 1, data = data, trace = FALSE)
model_reduced <- multinom(Y ~ X1 - 1, data = data, trace = FALSE)

n_test <- 1000

X1_test <- X1_levels |> 
  rep(each = n_test) |> 
  factor() |> 
  relevel(X1_ref_level)

X2_arg_list_test <- list(X2_intervals, X2_shape_params, n_test) |> 
  modify_if(is.numeric, ~ setNames(rep(., length(X1_levels)), X1_levels)) |> 
  unlist(recursive = FALSE) |> 
  split(X1_levels) |> 
  map(unname) |> 
  map(~ setNames(., c("interval", "shape", "n_samples")))

X2_test <- generate_X2_samples(X2_arg_list_test)

test_data <- data.frame(X1 = X1_test,
                        X2 = X2_test)

entropy_rmse <- function(model, test_data, true_entropy) {
  probs_pred <- predict(model, newdata = test_data, type = "probs")
  entropy_pred <- probs_pred |> 
    colMeans() |> 
    entropy()
  print(entropy_pred)
  return(entropy_pred - true_entropy)
}

error_df <- data.frame(X1 = X1_levels,
                      interaction = NA,
                      main_effects = NA,
                      reduced = NA)

for (h in X1_levels) {
  
  true_entropy <- H_0 |> 
    filter(X1 == h) |> 
    pull(entropy)
  
  test_data_h <- test_data |> 
    filter(X1 == h)
  
  entropy_errors_interaction <- entropy_rmse(model_interaction, test_data_h, true_entropy)
  entropy_errors_main_effects <- entropy_rmse(model_main_effects, test_data_h, true_entropy)
  entropy_errors_reduced <- entropy_rmse(model_reduced, test_data_h, true_entropy)
  
  error_df[error_df$X1 == h,]$interaction <- entropy_errors_interaction
  error_df[error_df$X1 == h,]$main_effects <- entropy_errors_main_effects
  error_df[error_df$X1 == h,]$reduced <- entropy_errors_reduced
}

error_df

# Perform LRT
lrt_stat <- 2 * (logLik(model_interaction) - logLik(model_main_effects))
p_value_lrt <- pchisq(lrt_stat, df = length(coef(model_interaction)) - length(coef(model_main_effects)), lower.tail = FALSE)

print(p_value_lrt)

lrt_stat <- 2 * (logLik(model_main_effects) - logLik(model_reduced))
p_value_lrt <- pchisq(lrt_stat, df = length(coef(model_main_effects)) - length(coef(model_reduced)), lower.tail = FALSE)

print(p_value_lrt)



