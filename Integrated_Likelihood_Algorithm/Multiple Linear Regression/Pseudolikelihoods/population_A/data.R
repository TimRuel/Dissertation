seed <- 956103

set.seed(seed)

p <- 10

Beta_0 <- runif(p, -5, 5) |> 
  matrix(ncol = 1)

n <- 100

X <- rep(1, n) |> 
  c(rnorm(n * (p - 1), 0, 1)) |> 
  matrix(nrow = n,
         ncol = p, 
         byrow = FALSE)

sigma_squared_0 <- 0.1

eps <- rnorm(n, 0, sqrt(sigma_squared_0)) |> 
  matrix(ncol = 1)

Y <- X %*% Beta_0 + eps

X_h <- X[4,] |> 
  matrix()
