seed <- 398724

set.seed(seed)

alpha_0 <- 3

beta_0 <- 5

sigma_squared_0 <- 1

n <- 100

x <- sample(1:n)

eps <- rnorm(n, 0, sqrt(sigma_squared_0))

y <- alpha_0 + beta_0 * x + eps