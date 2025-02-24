source("data_utils.R")
source("../../utils.R")
library(dplyr)

seed <- 6991

set.seed(seed)

J <- 6 # number of levels of response variable

m <- c(100, 100, 100)  # number of observations at each level of categorical predictor

C <- length(m) # number of levels of categorical predictor

n <- sum(m) # total number of observations

X1_levels <- LETTERS[1:C]

names(m) <- X1_levels

X1_ref_level <- X1_levels[1]

X1 <- X1_levels |> 
  rep(times = m) |> 
  factor() |> 
  relevel(X1_ref_level)

sigma <- 1

X2 <- get_X2_samples(X1, sigma)

X_design <- model.matrix(~ X1*X2 - 1)

p <- ncol(X_design) # Number of terms in model after dummy encoding

Beta2 <- rnorm(J - 1, sd = 0.3)

n_samples <- 10000

Beta_0 <- get_Beta_0(X_design, Beta2, n_samples, sigma)

true_probs <- X_design %*% cbind(0, Beta_0) |> 
  apply(1, softmax) |> 
  t() |> 
  data.frame() |> 
  rename_with( ~ paste0("Y", 1:J)) |> 
  mutate(X1_level = rep(X1_levels, times = m)) |> 
  select(X1_level, everything())

theta_0 <- true_probs |> 
  group_by(X1_level) |> 
  summarise(across(everything(), \(x) mean(x))) |> 
  data.frame()

H_0 <- theta_0 |> 
  group_by(X1_level) |> 
  rowwise() |> 
  mutate(entropy = entropy(c_across(everything()))) |> 
  select(X1_level, entropy) |> 
  data.frame()

Y <- true_probs |>
  select(-X1_level) |> 
  apply(1, \(prob) sample(1:J, size = 1, prob = prob)) |> 
  unlist() |>
  unname() |> 
  factor(levels = 1:J)

Y_one_hot <- model.matrix(~ Y)[,-1]

data <- data.frame(Y = Y,
                   X1 = X1, 
                   X2 = X2)

formula <- Y ~ .^2 - 1

model <- fit_multinomial_logistic_model(data, formula)

Beta_MLE <- get_Beta_MLE(model)

threshold <- ceiling(abs(log_likelihood(Beta_MLE, X, model.matrix(~ Y)[,-1])) + 20)

predict(model, data, type = "probs") |>
  as.data.frame() |>
  mutate(X1_level = rep(X1_levels, times = m)) |>
  aggregate(. ~ X1_level, data = _, FUN = mean) |>
  select(-X1_level) |>
  apply(1, entropy)

fit <- nnet::multinom(Y ~ X1 - 1, data = data, trace = FALSE)

predict(fit, data, type = "probs") |>
  unique() |>
  apply(1, entropy)

H_0

