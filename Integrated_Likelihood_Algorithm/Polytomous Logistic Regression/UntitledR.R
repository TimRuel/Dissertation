
num <- 40

psi_grid[num]

omega_hat <- log_integrated_likelihood_vanilla_MC$omega_hat_list[[30]]

probs <- X_one_hot %*% omega_hat |>
  apply(1, adj_softmax) |>
  t()

f <- function(Beta) {
  
  Beta <- matrix(Beta,
                 nrow = nrow(omega_hat),
                 ncol = ncol(omega_hat),
                 byrow = FALSE)
  
  Y_hat <- X_one_hot %*% Beta
  
  sum(rowSums(probs * Y_hat) - log(1 + rowSums(exp(Y_hat))))
}

Beta_hat <- log_integrated_likelihood_vanilla_MC$Beta_hat_matrices[[1]][num,] |> 
  matrix(nrow = nrow(omega_hat),
         ncol = ncol(omega_hat),
         byrow = FALSE)

psi_logits <- X_h_one_hot %*% cbind(0, Beta_hat)

U <- rnorm(15) |> matrix(nrow = 3)

U_logits <- X_h_one_hot %*% U

a <- sweep(U, 2, psi_logits[-1] / U_logits, "*")

X_h_one_hot %*% cbind(0, a) |>
  LDATS::softmax() |>
  get_entropy()

psi_grid[num]

f(a)

f(Beta_hat)

data |> 
  group_by(X, Y) |> 
  summarise(n = n(), .groups = "drop") |> 
  complete(X, Y, fill = list(n = 0)) |> 
  group_by(X) |> 
  mutate(proportion = n / sum(n)) |> 
  ungroup() |> 
  pull(proportion) |> 
  matrix(nrow = 3,
         byrow = TRUE)

unique(X_one_hot) %*% cbind(0, Beta_MLE) |> 
  apply(1, LDATS::softmax) |> 
  t() |> 
  round(3)

i <- 5

sqrt(sum((U_list[[i]] - omega_hat_list[[i]])^2))

X_h_one_hot %*% U_list[[1]]



dim(mat)

a <- init_guess |> 
  unname() |> 
  matrix(nrow = nrow(omega_hat),
         ncol = ncol(omega_hat),
         byrow = FALSE)

X_h_one_hot %*% cbind(0, a) |> 
  LDATS::softmax() |> 
  get_entropy()

X_h_one_hot %*% cbind(0, Beta) |> 
  LDATS::softmax() |> 
  get_entropy()

c <- b |> 
  unname() |> 
  matrix(nrow = nrow(omega_hat),
         ncol = ncol(omega_hat),
         byrow = FALSE)

X_h_one_hot %*% cbind(0, c) |> 
  LDATS::softmax() |> 
  get_entropy()

-f(init_guess)

-f(b)

get_log_likelihood(init_guess, X_one_hot, Y_one_hot)





-f(omega_hat)

-f(omega_hat_list[[2
                   ]])





