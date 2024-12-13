seed <- 3409836

set.seed(seed)

beta_0 <- -1

beta_1 <- 1

n <- 15

x <- rnorm(n, 0, 1)

y <- rbinom(n, 1, sigmoid::sigmoid(beta_0 + beta_1 * x))

data <- data.frame(x = x, y = y)