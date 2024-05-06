library(tidyverse)
library(nloptr)
library(LaplacesDemon)
library(broom)

set.seed(7835)

# Define observed x data
lambda_0 <- 5
n <- 10
x <- rpois(n, lambda_0)

# Define observed y data
mu_0 <- 10
m <- 15
y <- rpois(m, mu_0)

# Define hyperparameters for u1 (lambda)
alpha_1 <- 0
beta_1 <- 0

# Define hyperparameters for u2 (mu)
alpha_2 <- 0
beta_2 <- 0

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
R <- 10

# Initialize progress bar for for loop
pb = txtProgressBar(min = 0, max = length(psi), initial = 0, style = 3) 

for (i in 1:length(psi)) {
  
  # Update progress bar
  setTxtProgressBar(pb, i)
  
  # Initialize vector for holding values of log of likelihood ratio
  log_L_ratio <- c()
  
  # Random draw from posterior for theta
  u1 <- rgamma(1, shape = sum(x) + alpha_1, rate = n + beta_1)
  u2 <- rgamma(1, shape = sum(y) + alpha_2, rate = m + beta_2)
  u <- c(u1, u2)
  
  # Initialize starting point for searching for optimal theta value
  theta0 <- c(rgamma(1, shape = sum(x), rate = n), rgamma(1, shape = sum(y), rate = m)) 
  
  for (j in 1:R) {
    
    # Find value of omega that minimizes distance from u subject to constraints
    Q <- auglag(x0 = u,
                fn = function(omega) dist(matrix(c(u, omega), nrow = 2, byrow = TRUE))[1],
                heq = function(omega) g(omega) - psi_hat,
                lower = c(0, 0))$par
    
    # Find value of theta that minimizes log-likelihood expectation function subject to constraints
    T_psi <- auglag(x0 = theta0,
                    fn = function(theta) -E_log_like(theta, Q),
                    heq = function(theta) g(theta) - psi[i],
                    lower = c(0, 0))$par
    
    # Calculate ratio of likelihood at optimal theta to likelihood at initial random draw for theta
    log_L_ratio[j] <- log_likelihood(T_psi) - log_likelihood(u)
  }
  
  # Estimate value of integrated likelihood for current value of psi
  L_bar[i] <- mean(exp(log_L_ratio))
  
  # Find value of theta for obtaining profile log-likelihood
  theta_hat_p <- auglag(x0 = theta0,
                        fn = function(theta) -E_log_like(theta, theta_hat),
                        heq = function(theta) g(theta) - psi[i],
                        lower = c(0, 0))$par
  
  # Calculate value of profile likelihood for current value of psi
  log_L_p[i] <- log_likelihood(theta_hat_p)
}

likelihood_vals <- data.frame(psi = psi, 
                              Integrated = L_bar,
                              Profile = exp(log_L_p))

# save(likelihood_vals, file = "likelihood_df_1000_iter.Rda")

log_likelihood_vals <- likelihood_vals %>% 
  # filter(Integrated < max(Integrated)) %>%
  mutate(Integrated = log(Integrated / max(Integrated)),
         Profile = log(Profile / max(Profile))) 
  
log_likelihood_vals_tidy <- log_likelihood_vals %>% 
  pivot_longer(cols = c("Integrated", "Profile"),
               names_to = "Pseudolikelihood",
               values_to = "loglikelihood") 

log_likelihood_vals_tidy %>% 
  ggplot() +
  geom_point(aes(x = psi, y = loglikelihood, color = Pseudolikelihood),
             size = 0.5,
             alpha = 0.5) +
  geom_smooth(aes(x = psi, y = loglikelihood, color = Pseudolikelihood),
              linewidth = 0.9,
              se = TRUE,
              fullrange = TRUE) +
  ylab("Log-Likelihood") +
  xlab(expression(psi)) +
  theme_minimal()

fitted_models <- log_likelihood_vals_tidy %>%
  group_by(Pseudolikelihood) %>% 
  nest() %>% 
  mutate(model = map(data, ~ loess(loglikelihood ~ psi, data = .))) %>% 
  select(-data) 

fitted_models %>% 
  mutate(RSE = map_dbl(model, \(x) x$s)) %>% 
  select(-model)





  
  








