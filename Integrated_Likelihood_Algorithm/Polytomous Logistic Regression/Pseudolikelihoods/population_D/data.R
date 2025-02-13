source("../../utils.R")

seed <- 348257

set.seed(seed)

J <- 6 # number of levels of response variable

c <- 3 # number of levels of categorical predictor

m <- 50 # number of observations at each level of predictor

n <- sum(m * c)

epsilon <- 0.01

max_iter <- 1e6

theta_0 <- get_theta_0(J, c, epsilon, max_iter)

Y <- get_Y(theta_0, m)

contrast <- contr.sum

X1 <- get_X(c, m, contrast)

# X2 <- rnorm(n, mean = 5*(1:c), sd = 5)

X2 <- rnorm(n, sd = 5)

data <- data.frame(X1 = X1,
                   X2 = X2,
                   Y = Y)

formula <- Y ~ .^2

model <- fit_multinomial_logistic_model(data, formula)

Beta_MLE <- get_Beta_MLE(model)

p <- nrow(Beta_MLE) # Number of coefficient terms in model

threshold <- 170

