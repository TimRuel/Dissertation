# Store observations
data <- c(1, 1, 2, 4, 7, 10)

# Define parameter of interest function
entropy <- function(x) {
  
  y <- ifelse(x != 0, x*log(x), 0)
  
  return(-sum(y))
}

# Define function that generates a random variate from the set omega_psihat with a distribution not depending on psi
gen_omega_hat <- function(n, g) {
  
  # Store number of observations
  m <- length(n)
  
  # Draw random variate from the m-dimensional probability simplex
  u <- rexp(m) 
  u <- u / sum(u)
  
  # Define MLE for parameter of interest
  psi_hat <- g(n / sum(n))
  
  # Find value of omega that minimizes objective function subject to constraints
  omega_hat <- nloptr::auglag(x0 = u,
                              fn = function(omega) -sum(u*log(omega)),
                              heq = function(omega) c(sum(omega) - 1, g(omega) - psi_hat),
                              lower = rep(0, m))$par
  
  return(omega_hat)
}

omega_hat <- gen_omega_hat(data, entropy)

# Define function that generates theta_hat from omega_hat

gen_theta_hat <- function(omega_hat, g, psi) {
  
  # Store dimension of omega_hat (same as number of observations)
  m <- length(omega_hat)
  
  # Find value of theta that minimizes objective function subject to constraints
  out <- nloptr::auglag(x0 = omega_hat,
                        fn = function(theta) -sum(omega_hat*log(theta)),
                        heq = function(theta) c(sum(theta) - 1, g(theta) - psi),
                        lower = rep(0, m))
  
  return(out$par)
}

theta_hat <- gen_theta_hat(omega_hat, entropy, 0.5)

# Define likelihood function
likelihood <- function(theta, n) prod(theta^n)

likelihood(theta_hat, data)


