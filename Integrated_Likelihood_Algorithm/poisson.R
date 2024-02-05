library(tidyverse)
library(nloptr)
library(LaplacesDemon)

set.seed(1996)

# Define observed x data
lambda_0 <- 5
n <- 100
x <- rpois(n, lambda_0)

# Define observed y data
mu_0 <- 100
m <- 15
y <- rpois(m, mu_0)

# Define hyperparameters for u1 (lambda)
alpha_1 <- 1
beta_1 <- 1

# Define hyperparameters for u2 (mu)
alpha_2 <- 1
beta_2 <- 1

# Define likelihood function
likelihood <- function(theta) {
  exp(-(n*theta[1] + m*theta[2]))*theta[1]^(sum(x))*theta[2]^(sum(y))
}

# Define log-likelihood function
log_likelihood <- function(theta) {
  -(n*theta[1] + m*theta[2]) + sum(x)*log(theta[1]) + sum(y)*log(theta[2])
}

# Define parameter of interest function 
g <- function(theta) theta[1] + 2*theta[2]

theta_0 <- c(lambda_0, mu_0)
psi_0 <- g(theta_0)

# Define MLE for theta
theta_hat <- c(mean(x), mean(y))

# Define MLE for parameter of interest
psi_hat <- g(theta_hat)

# Calculate standard error of MLE
psi_hat_se <- sqrt(var(x) / n + 4*var(y) / m)

# Calculate margin of error
n_std_errors = 3
MoE = n_std_errors * psi_hat_se

# Define values for parameter of interest at which to evaluate the integrated likelihood  
psi <- seq(psi_hat - MoE, psi_hat + MoE, 0.1)

# Define log-likelihood expectation function to be minimized
E_log_like <- function(theta, omega) sum((-theta + log(theta)*omega)*c(n, m))

# Initialize vector for holding values of integrated likelihood
L_bar <- c()

# Initialize vector for holding values of profile likelihood
log_L_p <- c()

# Number of replications for each value of psi
R <- 1

# Initialize progress bar for for loop
pb = txtProgressBar(min = 0, max = length(psi), initial = 0, style = 3) 

for (i in 1:length(psi)) {
  
  # Update progress bar
  setTxtProgressBar(pb, i)
  
  # Initialize vector for holding values of likelihood ratio
  L_ratio <- c()
  
  for (j in 1:R) {
    
    # Random draw from posterior for theta
    u1 <- rgamma(1, shape = sum(x) + alpha_1, rate = n + beta_1)
    u2 <- rgamma(1, shape = sum(y) + alpha_2, rate = m + beta_2)
    u <- c(u1, u2)
    
    # Find value of omega that minimizes distance from u subject to constraints
    Q <- auglag(x0 = u,
                fn = function(omega) dist(matrix(c(u, omega), nrow = 2, byrow = TRUE))[1],
                heq = function(omega) g(omega) - psi_hat,
                lower = c(0, 0))$par
    
    # Initialize starting point for searching for optimal theta value
    theta0 <- c(rgamma(1, shape = sum(x), rate = n), rgamma(1, shape = sum(y), rate = m)) 
    
    # Find value of theta that minimizes log-likelihood expectation function subject to constraints
    T_psi <- auglag(x0 = theta0,
                    fn = function(theta) -E_log_like(theta, Q),
                    heq = function(theta) g(theta) - psi[i],
                    lower = c(0, 0))$par
    
    # Calculate ratio of likelihood at optimal theta to likelihood at initial random draw for theta
    L_ratio[j] <- log_likelihood(T_psi) - log_likelihood(u)
  }
  
  # Calculate value of integrated likelihood for current value of psi
  L_bar[i] <- mean(exp(L_ratio))
  
  # Find value of theta for obtaining profile log-likelihood
  theta_hat_p <- auglag(x0 = theta0,
                        fn = function(theta) -E_log_like(theta, theta_hat),
                        heq = function(theta) g(theta) - psi[i],
                        lower = c(0, 0))$par
  
  # Calculate value of profile likelihood for current value of psi
  log_L_p[i] <- log_likelihood(theta_hat_p)
}

likelihood_vals <- data.frame(psi = psi, 
                              Integrated = log(L_bar / max(L_bar)),
                              Profile = log_L_p - max(log_L_p)) %>% 
  pivot_longer(cols = c("Integrated", "Profile"),
               names_to = "Pseudolikelihood",
               values_to = "log-likelihood")

likelihood_vals %>% 
  ggplot() +
  geom_point(aes(x = psi, y = `log-likelihood`, color = Pseudolikelihood),
             size = 0.5, ) +
  geom_smooth(aes(x = psi, y = `log-likelihood`, color = Pseudolikelihood),
              se = FALSE,
              linewidth = 0.9,
              fullrange = TRUE) +
  theme_minimal() +
  theme(legend.position = c(0.5, 0.2),
        legend.background = element_rect())


likelihood_vals %>% 
  group_by(Pseudolikelihood) %>% 
  mutate(variance = sd(`log-likelihood`))
