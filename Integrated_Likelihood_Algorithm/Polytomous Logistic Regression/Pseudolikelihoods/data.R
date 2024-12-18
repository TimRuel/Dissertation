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

model <- get_multinomial_logistic_model(data)

b <- t(coef(model))
# 
# get_log_likelihood(model)
# 
# get_log_likelihood(b = b, data = data)
# 
# get_Beta_MLE(model)

h <- 50L

X_h <- data |> 
  select(-Y) |> 
  slice(h) |>
  unname() |> 
  as.matrix() |> 
  (\(mat) cbind(1, mat))()

U <- MASS::mvrnorm(p, rep(0, J - 1), diag(J - 1))

omega_hat <- get_omega_hat(U, model, X_h)

Beta_MLE <- get_Beta_MLE(model)

Beta_hat <- get_Beta_hat(0.5, omega_hat, X, X_h, rep(1, length(omega_hat)))

X_h %*% cbind(0, Beta_hat) |> 
  LDATS::softmax() |> 
  get_entropy()

X_h |> 
  t() |> 
  data.frame() |> 
  predict(model, newdata = _, type = "probs") 



omega_hat

get_psi_hat(b = omega_hat, X_h = X_h)

get_psi_hat(model = model, X_h = X_h)




