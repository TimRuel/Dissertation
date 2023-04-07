multinom_MC <- function(data, param, reps, step_size) {
  
  n <- data
  
  m <- length(n)
  
  # Define MLE for theta
  theta.hat <- n / sum(n)
  
  # Define MLE for parameter of interest
  psi_hat <- param(theta.hat)
  
  # Define likelihood function
  likelihood <- function(theta, n) prod(theta^n)
  
  omega_hat <- list()
  
  for (j in 1:reps) {
    
    # Draw random variartes from the m-dimensional probability simplex
    u <- rdirichlet(1, rep(1, m)) 
    # u <- rexp(m) 
    # u <- u / sum(u)
    
    # Define first objective function to be minimized
    entropy1 <- function(omega) -sum(u*log(omega))
    
    # Define first constraint function
    constraint1 <- function(omega) c(sum(omega) - 1, param(omega) - psi_hat)
    
    # Initialize starting point for searching for optimal omega value
    omega0 <- rdirichlet(1, rep(1, m)) 
    # omega0 <- sample(1:10, m)
    # omega0 <- omega0 / sum(omega0) 
    
    # Find value of omega that minimizes objective function subject to constraints
    omega_hat[[j]] <- auglag(x0 = omega0,
                             fn = entropy1,
                             heq = constraint1,
                             lower = rep(0, m))$par
  }
  
  # Define values for parameter of interest at which to evaluate the
  # integrated likelihood function 
  psi1 <- seq(1, log(m), step_size)
  
  # Initialize vector for holding values of integrated likelihood function
  L_bar <- c()
  
  L_p <- c()
  
  # Calculate value of integrated likelihood for each parameter of interest value
  for (i in 1:length(psi1)) {
    
    # Initialize vector for holding values of likelihood function
    L <- c()
    
    L1 <- c()
    
    # Initialize starting point for searching for optimal theta value
    theta_hat0 <- sample(1:10, m)
    theta_hat0 <- theta_hat0 / sum(theta_hat0) 
    
    for (j in 1:reps) {
      
      # Define second objective function to be minimized
      entropy2 <- function(theta) -sum(omega_hat[[j]]*log(theta))
      
      # Define second constraint function
      constraint2 <- function(theta) c(sum(theta) - 1, param(theta) - psi1[i])
      
      # Find value of theta that minimizes objective function subject to constraints
      out <- auglag(x0 = theta_hat0,
                    fn = entropy2,
                    heq = constraint2,
                    lower = rep(0, m))
      
      # Store optimal theta value
      theta_hat <- out$par
      
      # Evaluate likelihood function at optimal theta value and store result
      L[j] <- likelihood(theta_hat, n)
      
      # Define third objective function to be minimized
      entropy3 <- function(theta) -sum(theta.hat*log(theta))
      
      # Define second constraint function
      constraint3 <- function(theta) c(sum(theta) - 1, param(theta) - psi1[i])
      
      # Find value of theta that minimizes objective function subject to constraints
      out <- auglag(x0 = theta_hat0,
                    fn = entropy3,
                    heq = constraint3,
                    lower = rep(0, m))
      
      # Store optimal theta value
      theta_hat_p <- out$par
      
      # Evaluate likelihood function at optimal theta value and store result
      L1[j] <- likelihood(theta_hat_p, n)
    }
    
    print(out$convergence)
    
    L_bar[i] <- mean(L)
    
    L_p[i] <- mean(L1)
  }
  
  return(list(integrated = L_bar,
              profile = L_p))
}

g <- function(x) {
  
  y <- ifelse(x != 0, x*log(x), 0)
  
  return(-sum(y))
}

out <- multinom_MC(data = c(1, 1, 2, 4, 7, 10),
                   param = g,
                   reps = 250,
                   step_size = 0.01)



