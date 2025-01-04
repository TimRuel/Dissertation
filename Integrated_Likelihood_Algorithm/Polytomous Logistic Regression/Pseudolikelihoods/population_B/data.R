library(dplyr)

source("../../utils.R")

seed <- 212532

set.seed(seed)

J <- 5 # number of levels of response variable

p <- 3 # number of levels of predictor

m <- 50 # number of observations at each level of predictor

theta_0 <- seq(0, log(J), length.out = p + 1) |> 
  (\(x) mapply(c, x[-length(x)], x[-1], SIMPLIFY = FALSE))() |> 
  map(\(x) get_probability_vector(k = J, target_entropy_range = x, max_iter = 100000))

get_entropy(theta_0[[1]])
get_entropy(theta_0[[2]])
get_entropy(theta_0[[3]])

Y <- theta_0 |> 
  map(\(p) sample(1:J, size = m, prob = p, replace = TRUE)) |> 
  unlist() |> 
  factor(levels = 1:J)

X <- 1:p |> 
  rep(each = m) |> 
  factor()

contrasts(X) <- contr.sum

data <- data.frame(X = X,
                   Y = Y)