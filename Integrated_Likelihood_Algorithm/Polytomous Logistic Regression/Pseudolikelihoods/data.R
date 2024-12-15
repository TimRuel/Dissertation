seed <- 54323872

set.seed(seed)

J <- 4

p <- 5

n <- 10000

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

model <- get_multinomial_logistic_model(data)

b <- t(coef(model))
# 
# get_log_likelihood(model)
# 
# get_log_likelihood(b = b, data = data)
# 
# get_Beta_MLE(model)

h <- 500L

X_h <- data |> 
  select(-Y) |> 
  slice(h) |>
  unname() |> 
  as.matrix() |> 
  t()

X_h |> 
  t() |> 
  data.frame() |> 
  predict(model, newdata = _, type = "probs") 

get_psi_hat(model, X_h)

cbind(1, t(X_h)) %*% b |> softmax()

