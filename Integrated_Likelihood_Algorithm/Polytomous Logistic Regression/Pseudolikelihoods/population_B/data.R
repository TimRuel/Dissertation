library(dplyr)

seed <- 54323872

set.seed(seed)

J <- 3 # number of levels of response variable

p <- 3 # number of levels of predictor

m <- 200 # number of observations at each level of predictor

Beta_0 <- runif(p, -1, 1) |> 
  replicate(J - 1, expr = _)

X <- 1:p |> 
  rep(each = m) |> 
  factor() |> 
  relevel(ref = p) |> 
  (\(X) model.matrix(~ X))()

logits <- cbind(0, X %*% Beta_0)

probs <- LDATS::softmax(logits)

Y <- apply(probs, 1, function(p) sample(1:J, size = 1, prob = p))

data <- X[,-1] |> 
  data.frame(Y = factor(Y))
