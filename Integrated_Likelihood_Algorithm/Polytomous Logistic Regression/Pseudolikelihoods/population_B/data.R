library(dplyr)

source("../../utils.R")

seed <- 46990395

set.seed(seed)

J <- 6 # number of levels of response variable

p <- 3 # number of levels of predictor

m <- 25 # number of observations at each level of predictor

theta_0 <- get_theta_0(J, p, 0.01, 1e6)

get_entropy(theta_0[[1]])
get_entropy(theta_0[[2]])
get_entropy(theta_0[[3]])

Y <- get_Y(theta_0, m)

X <- get_X(p, m, contr.sum)

data <- data.frame(X = X,
                   Y = Y)