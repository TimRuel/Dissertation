n <- 10000

mu <- 0.1

pi <- 0.01

u <- rbinom(n, 1, 1 - pi)
v <- rpois(n, mu)

y <- u * v

y_bar <- mean(y)

pi_bar <- mean(y == 0)

pi_bar < exp(-mu)

pi_hat <- (pi_bar - exp(-mu)) / (1 - exp(-mu))

l <- function(p) {
  
  n * pi_bar * log(p + (1-p)*exp(-mu)) + n*(1 - pi_bar) * log(1-p)
  
}

profile <- function(mu) {
  
  pi_hat <- (pi_bar - exp(-mu)) / (1 - exp(-mu))
  
  n * pi_bar * log(pi_hat + (1 - pi_hat) * exp(-mu)) + n * (1 - pi_bar) * (log(1 - pi_hat) - mu) + n * y_bar * log(mu)
}

profile2 <- function(mu) {
  
  -n * (1 - pi_bar) * (mu + log(1 - exp(-mu))) + n * y_bar * log(mu)
  
}


mu_vals <- seq(0.001, 50, 0.001)

profile_vals <- purrr::map_dbl(mu_vals, profile2)

plot(mu_vals, profile_vals, pch = 16)

p_vals <- seq(0, 1, 0.001)

l_vals <- purrr::map_dbl(p_vals, l)

plot(p_vals, l_vals, pch = 16)
points(pi_hat, l(pi_hat), col = "red", pch = 16)
points(pi, l(pi), col = "blue", pch = 16)
abline(v = pi_bar, col = "green")
abline(v = exp(-mu), col = "orange")



l(pi)
l(pi_hat)

p_vals[which.max(l_vals)]


