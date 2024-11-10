n <- 30
pi_0 <- 0.2
mu_0 <- 5
U <- rbinom(n, 1, 1 - pi_0)
V <- rpois(n, mu_0)
Y <- U * V



Y_bar <- mean(Y)
pi_bar <- mean(Y == 0)

l_p <- function(mu) n * (Y_bar * log(mu) - (1 - pi_bar) * (log(1 - exp(-mu)) + mu))

mu_grid <- get_mu_grid(Y, 0.01, 4)
l_p_vals <- purrr::map_dbl(mu_grid, l_p)
mu_hat_p <- mu_vals[which.max(l_p_vals)]
plot(mu_grid, l_p_vals)
abline(v = mu_hat, col = "green")
abline(v = mu_0, col = "red")


