source("data_utils.R")
source("../../utils.R")

seed <- 6991

set.seed(seed)

J <- 6 # number of levels of response variable

m <- c(25, 25, 25) # number of observations at each level of categorical predictor

C <- length(m) # number of levels of categorical predictor

n <- sum(m) # total number of observations

X1_levels <- letters[1:C]

X1_ref_level <- X1_levels[1]

X1 <- X1_levels |> 
  rep(times = m) |> 
  factor() |> 
  relevel(X1_ref_level)

# contrasts(X1) <- contr.sum

X2 <- runif(n, 0, 90) |> 
  scale(center = FALSE) |> 
  round(1)

X <- model.matrix(~ X1*X2 - 1)

p <- ncol(X) # Number of terms in model after dummy encoding

X2_coef <- rnorm(J - 1, sd = 0.3)

Beta_0 <- get_Beta_0(X2_coef)

true_probs <- X %*% cbind(0, Beta_0) |> 
  apply(1, softmax) |> 
  t()

theta_0 <- true_probs |> 
  as.data.frame() |> 
  mutate(X1_level = rep(X1_levels, times = m)) |> 
  aggregate(. ~ X1_level, data = _, FUN = mean) |> 
  select(-X1_level) |> 
  apply(1, entropy) |> 
  setNames(X1_levels)

Y <- true_probs |>
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

theta_0

