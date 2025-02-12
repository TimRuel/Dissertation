source("../../utils.R")

seed <- 348257

set.seed(seed)

J <- 6 # number of levels of response variable

c <- 3 # number of levels of predictor

m <- 10 # number of observations at each level of predictor

n <- sum(m * c)

epsilon <- 0.01

max_iter <- 1e6

theta_0 <- get_theta_0(J, c, epsilon, max_iter)

Y <- get_Y(theta_0, m)

contrast <- contr.sum

X <- get_X(c, m, contrast)

data <- data.frame(X = X,
                   Y = Y)

formula <- Y ~ .

model <- fit_multinomial_logistic_model(data, formula)

Beta_MLE <- get_Beta_MLE(model)

p <- nrow(Beta_MLE) # Number of parameters in model

threshold <- 40

# Add interactions and/or continuous variables

