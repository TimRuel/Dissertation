library(tidyverse)
library(nloptr)
library(LaplacesDemon)

# Define observed data
n <- c(1, 1, 2, 4, 7, 10)

m <- length(n)

# Define likelihood function
likelihood <- function(theta, n) prod(theta^n)

# Define MLE for theta
theta_hat <- n / sum(n)

# Define parameter of interest function (entropy in this case)
g <- function(x) {
  
  y <- ifelse(x != 0, x*log(x), 0)
  
  return(-sum(y))
}

# Define MLE for parameter of interest
psi_hat <- g(theta_hat)

# Define parameters of Dirichlet posterior based on uniform Dirichlet prior and multinomial likelihood
alpha <- n + 1

# Define values for parameter of interest at which to evaluate the integrated likelihood  
psi <- seq(1, log(m), 0.01)

# Define first constraint function
constraint1 <- function(omega) c(sum(omega) - 1, g(omega) - psi_hat)

# Define entropy function to be minimized
entropy <- function(theta) -sum(theta_hat*log(theta))

# Initialize vector for holding values of integrated likelihood
L_bar <- c()

# Initialize vector for holding values of profile likelihood
L_p <- c()

# Number of replications for each value of psi
R <- 250

for (i in 1:length(psi)) {
  
  # Define second constraint function
  constraint2 <- function(theta) c(sum(theta) - 1, g(theta) - psi[i])
  
  # Initialize vector for holding values of likelihood ratio
  L_ratio <- c()
  
  # Initialize vector for holding values of likelihood 
  L <- c()
  
  for (j in 1:R) {
    
    # Random draw from posterior for theta
    u <- rdirichlet(1, alpha)
    
    # Define distance function to be minimized
    distance <- function(omega) dist(matrix(c(u, omega), nrow = 2, byrow = TRUE))[1]
    
    # Find value of omega that minimizes distance function subject to constraints
    Q <- auglag(x0 = u,
                fn = distance,
                heq = constraint1,
                lower = rep(0, m))$par
    
    # Define log-likelihood expectation function to be minimized
    E_log_like <- function(theta) -sum(Q*log(theta))
    
    # Initialize starting point for searching for optimal theta value
    theta0 <- sample(1:10, m)
    theta0 <- theta0 / sum(theta0) 
    
    # Find value of theta that minimizes log-likelihood expectation function subject to constraints
    T_psi <- auglag(x0 = theta0,
                    fn = E_log_like,
                    heq = constraint2,
                    lower = rep(0, m))$par
    
    # Calculate ratio of likelihood at optimal theta to likelihood at initial random draw for theta
    L_ratio[j] <- likelihood(T_psi, n) / likelihood(u, n)
    
    # Find value of theta that minimizes entropy function subject to constraints
    theta_hat_p <- auglag(x0 = theta0,
                          fn = entropy,
                          heq = constraint2,
                          lower = rep(0, m))$par
    
    # Evaluate likelihood function at optimal theta value and store result
    L[j] <- likelihood(theta_hat_p, n)
  }
  
  # Calculate value of integrated likelihood for current value of psi
  L_bar[i] <- mean(L_ratio)
  
  # Calculate value of profile likelihood for current value of psi
  L_p[i] <- mean(L)
}

likelihood_vals <- data.frame(psi = psi, 
                              Integrated = log(L_bar / max(L_bar)),
                              Profile = log(L_p / max(L_p)))

likelihood_vals %>% 
  pivot_longer(cols = c("Integrated", "Profile"),
               names_to = "Pseudolikelihood",
               values_to = "log-likelihood") %>% 
  ggplot() +
  scale_y_continuous(limits = c(-3, 0.1)) +
  geom_smooth(aes(x = psi, y = `log-likelihood`, color = Pseudolikelihood),
              se = FALSE,
              linewidth = 0.9,
              fullrange = TRUE) +
  theme_minimal() +
  theme(legend.position = c(0.2, 0.8),
        legend.background = element_rect())






