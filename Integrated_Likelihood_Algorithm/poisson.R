library(tidyverse)
library(nloptr)
library(LaplacesDemon)

set.seed(1996)

# Define observed x data
lambda_0 = 5
n = 10
x = rpois(n, lambda_0)

# Define observed y data
mu_0 = 10
m = 15
y = rpois(m, mu_0)

# Define likelihood function
likelihood <- function(theta, x, y) {
  exp(-(length(x)*theta[1] + length(y)*theta[2]))*lambda^(sum(x))*mu^(sum(y))
  }

# Define parameter of interest function 
g <- function(theta) theta[1] + 2*theta[2]

# Define MLE for theta
theta_hat <- c(mean(x), mean(y))

# Define MLE for parameter of interest
psi_hat <- g(theta_hat)

# Calculate standard error of MLE
psi_hat_se = sqrt(var(x) / n + 4*var(y) / m)

# Define values for parameter of interest at which to evaluate the integrated likelihood  
psi <- seq(psi_hat - 3*psi_hat_se, psi_hat + 3*psi_hat_se, 0.01)

# Define first constraint function
constraint1 <- function(omega) g(omega) - psi_hat

# Define second distance function to be minimized
distance1 <- function(theta) -sum(theta_hat*log(theta))

# Initialize vector for holding values of integrated likelihood
L_bar <- c()

# Initialize vector for holding values of profile likelihood
L_p <- c()

# Number of replications for each value of psi
R <- 250

for (i in 1:length(psi)) {
  
  # Define second constraint function
  constraint2 <- function(theta) g(theta) - psi[i]
  
  # Initialize vector for holding values of likelihood ratio
  L_ratio <- c()
  
  # Initialize vector for holding values of likelihood 
  L <- c()
  
  for (j in 1:R) {
    
    # Random draw from posterior for theta
    alpha_1 = 1
    beta_1 = 1
    u1 = rgamma(1, sum(x) + alpha_1, n + beta_1)
    
    alpha_2 = 1
    beta_2 = 1
    u2 = rgamma(1, sum(y) + alpha_2, m + beta_2)
    
    u = u1 * u2
    
    # Define first distance function to be minimized
    distance2 <- function(omega) dist(matrix(c(u, omega), nrow = 2, byrow = TRUE))[1]
    
    # Find value of omega that minimizes distance function subject to constraints
    Q <- auglag(x0 = u,
                fn = distance1,
                heq = constraint1,
                lower = 0)$par
    
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

