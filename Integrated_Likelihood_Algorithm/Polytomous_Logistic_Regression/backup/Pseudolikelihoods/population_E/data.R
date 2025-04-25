source("../../utils.R")

seed <- 348257

set.seed(seed)

J <- 8 # number of unique levels of response variable

# number of observations of response at each combination of level of categorical predictor with each unique value of continuous predictor
m1 <- rbind(
  
  c(10, 10, 10, 10, 10),
  c(10, 10, 10, 10, 10),
  c(10, 10, 10, 10, 10),
  c(10, 10, 10, 10, 10)
  
) 

m2 <- rep(ncol(m1), nrow(m1)) # number of unique values of continuous predictor within each level of categorical predictor

c <- length(m2) # number of levels of categorical predictor

n <- sum(m1)

epsilon <- 0.01

max_iter <- 1e6

theta_0 <- get_theta_0(J, c, epsilon, max_iter)

Y <- m1 |> 
  split(row(m1)) |> 
  purrr::map2(theta_0,
              \(vec, prob) purrr::map(vec, \(num) sample(1:J, size = num, prob = prob, replace = TRUE))
  ) |> 
  unlist() |> 
  unname() |> 
  factor(levels = 1:J)

contrast <- contr.sum

X1 <- get_X(c, m, contrast)

X2 <- m1 |> 
  split(row(m1)) |> 
  purrr::map(\(vec, prob) purrr::map(vec, \(num) sample(1:J, size = num, prob = prob, replace = TRUE))
  ) |> 
  unlist() |> 
  unname() |> 
  factor(levels = 1:J)

data <- data.frame(X1 = X1,
                   X2 = X2,
                   Y = Y)

formula <- Y ~ .^2

model <- fit_multinomial_logistic_model(data, formula)

Beta_MLE <- get_Beta_MLE(model)

p <- nrow(Beta_MLE) # Number of coefficient terms in model

threshold <- 170

