seed <- 205683

set.seed(seed)

alpha_0 <- 1

beta_0 <- -2

sigma_squared_0 <- 10

n <- 200

x <- rnorm(n, 0, 10)

eps <- rnorm(n, 0, sqrt(sigma_squared_0))

y <- alpha_0 + beta_0 * x + eps

