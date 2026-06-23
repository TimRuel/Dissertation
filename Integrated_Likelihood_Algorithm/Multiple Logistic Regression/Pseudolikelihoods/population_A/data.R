seed <- 5463498

set.seed(seed)

p <- 10

Beta_0 <- runif(p, -1, 1) |>
  matrix(ncol = 1)

n <- 100

X <- rep(1, n) |>
  c(rnorm(n * (p - 1), 0, 0.1)) |>
  matrix(nrow = n, ncol = p, byrow = FALSE)

Y <- rbinom(n, 1, plogis(X %*% Beta_0))

data <- X[, -1] |>
  data.frame(Y = Y)
