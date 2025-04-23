library(dplyr)

seed <- 54323872

set.seed(seed)

J <- 4

p <- 5

n <- 100

Beta_0 <- runif(p, -1, 1) |> 
  replicate(J - 1, expr = _)

X <- rep(1, n) |> 
  c(rnorm(n * (p - 1), 0, 0.1)) |> 
  matrix(nrow = n,
         ncol = p, 
         byrow = FALSE)

logits <- cbind(0, X %*% Beta_0)

probs <- LDATS::softmax(logits)

Y <- apply(probs, 1, function(p) sample(1:J, size = 1, prob = p))

data <- X[,-1] |> 
  data.frame(Y = factor(Y))
