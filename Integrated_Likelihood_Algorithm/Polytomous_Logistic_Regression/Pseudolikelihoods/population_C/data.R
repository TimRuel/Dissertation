source("data_utils.R")
source("../../utils.R")

seed <- 348257

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

X <- model.matrix(~ X1 - 1)

p <- ncol(X) # Number of terms in model after dummy encoding

Beta_0 <- get_Beta_0()

true_probs <- X %*% cbind(0, Beta_0) |> 
  apply(1, softmax) |> 
  t() |> 
  unique() |> 
  tapply(rep(1:C, times = J), 
         function(i) i) |> 
  setNames(X1_levels)

theta_0 <- true_probs |> 
  purrr::map_dbl(entropy)

Y <- true_probs |> 
  purrr::map2(m, \(prob, size) sample(1:J, size = size, prob = prob, replace = TRUE)) |> 
  unlist() |> 
  unname() |> 
  factor(levels = 1:J)

Y_one_hot <- model.matrix(~ Y)[,-1]

data <- data.frame(X1 = X1,
                   Y = Y)

formula <- Y ~ . - 1

model <- fit_multinomial_logistic_model(data, formula)

Beta_MLE <- get_Beta_MLE(model)

threshold <- ceiling(abs(log_likelihood(Beta_MLE, X, model.matrix(~ Y)[,-1])) + 20)
