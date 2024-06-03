library(tidyverse)
library(nloptr)
library(LaplacesDemon)
library(Rmpfr)
library(geomtextpath)

set.seed(7835)

# Define observed x data
lambda_0 <- 0.5
n <- 30
x <- rpois(n, lambda_0)

# Define observed y data
mu_0 <- 0.3
m <- 20
y <- rpois(m, mu_0)

# Define weights for PoI function
w <- c(0.6, 0.4)

# Define hyperparameters for u1 (lambda)
alpha_1 <- 0
beta_1 <- 0

# Define hyperparameters for u2 (mu)
alpha_2 <- 0
beta_2 <- 0

# Define likelihood function
likelihood <- function(theta) {
  
  theta <- theta %>% 
    mpfr(precBits = 106)
  
  (-theta * c(n, m)) %>% 
    sum() %>% 
    exp() %>% 
    "*"((theta^c(sum(x), sum(y))) %>%
          prod()) 
}

# Define parameter of interest function 
g <- function(theta) sum(theta*w)

theta_0 <- c(lambda_0, mu_0)
psi_0 <- g(theta_0)

# Define MLE for theta
theta_hat <- c(mean(x), mean(y))

# Define MLE for parameter of interest
psi_hat <- g(theta_hat)

# Calculate standard error of MLE
psi_hat_se <- sqrt(sum(c(var(x), var(y)) / c(n, m) * w^2))

# Calculate margin of error
n_std_errors = 1
MoE = n_std_errors * psi_hat_se

# Define values for parameter of interest at which to evaluate the integrated likelihood  
psi <- seq(max(0, psi_hat - MoE), psi_hat + MoE, 0.001)

# Define log-likelihood expectation function to be minimized
E_log_like <- function(theta, omega) sum((-theta + log(theta)*omega)*c(n, m))

# Initialize matrix for holding each iteration's values for likelihood ratio at each value of psi 
L_ratio <- mpfrArray(NA, precBits = 106, dim = c(R, length(psi)))

# Initialize vector for holding values of profile likelihood
L_p <- mpfrArray(NA, precBits = 106, dim = length(psi))

# Number of replications for each value of psi
R <- 10

# Initialize progress bar for for loop
pb = txtProgressBar(min = 0, max = R, initial = 0, style = 3) 

for (i in 1:R) {
  
  # Update progress bar
  setTxtProgressBar(pb, i)
  
  # Random draw from posterior for theta
  u1 <- rgamma(1, shape = sum(x) + alpha_1, rate = n + beta_1)
  u2 <- rgamma(1, shape = sum(y) + alpha_2, rate = m + beta_2)
  u <- c(u1, u2)
  
  # Find value of omega that minimizes distance from u subject to constraints
  Q <- auglag(x0 = u,
              fn = function(omega) dist(matrix(c(u, omega), nrow = 2, byrow = TRUE))[1],
              heq = function(omega) g(omega) - psi_hat,
              lower = c(0, 0))$par
  
  for (j in 1:length(psi)) {
    
    # Initialize starting point for searching for optimal theta value
    theta0 <- c(rgamma(1, shape = sum(x), rate = n), rgamma(1, shape = sum(y), rate = m))
    
    # Find value of theta that minimizes log-likelihood expectation function subject to constraints
    T_psi <- auglag(x0 = theta0,
                    fn = function(theta) -E_log_like(theta, Q),
                    heq = function(theta) g(theta) - psi[j],
                    lower = c(0, 0))$par
    
    # Calculate ratio of likelihood at optimal theta to likelihood at initial random draw for theta
    L_ratio[i, j] <- likelihood(T_psi) / likelihood(u)
    
    if (i == 1) {
      
      # Find value of theta for obtaining profile log-likelihood
      theta_hat_p <- auglag(x0 = theta_hat,
                            fn = function(theta) -E_log_like(theta, theta_hat),
                            heq = function(theta) g(theta) - psi[j],
                            lower = c(0, 0))$par
      
      # Calculate value of profile likelihood for current value of psi
      L_p[j] <- likelihood(theta_hat_p)
    }
  }
}

log_likelihood_vals <- data.frame(psi = psi, 
                                  Integrated = L_ratio %>% 
                                    apply(2, mean) %>% 
                                    log() %>%
                                    as.double(),
                                  Profile = L_p %>% 
                                    log() %>% 
                                    as.double()
                                  )

log_likelihood_vals_tidy <- log_likelihood_vals %>% 
  pivot_longer(cols = c("Integrated", "Profile"),
               names_to = "Pseudolikelihood",
               values_to = "loglikelihood") 

fitted_models <- log_likelihood_vals_tidy %>%
  group_by(Pseudolikelihood) %>% 
  nest() %>% 
  mutate(model = map(data, ~ loess(loglikelihood ~ psi, data = .))) %>% 
  select(-data) 

fitted_models %>% 
  mutate(shift = map(model, predict)) %>% 
  unnest(c(Pseudolikelihood, shift)) %>% 
  group_by(Pseudolikelihood) %>% 
  slice_max(shift) %>% 
  select(-model) %>% 
  ungroup() %>% 
  inner_join(log_likelihood_vals_tidy) %>% 
  mutate(loglikelihood = loglikelihood - shift) %>% 
  select(-shift) %>% 
  ggplot(aes(x = psi, y = loglikelihood, color = Pseudolikelihood, linetype = Pseudolikelihood)) +
  geom_smooth(linewidth = 1.3,
              se = FALSE) +
  geom_hline(aes(yintercept = max(loglikelihood)),
             linetype = 2) +
  geom_textvline(xintercept = psi_0, 
                 label = "psi[0]",
                 parse = TRUE,
                 color = "green4") +
  geom_textvline(aes(xintercept = psi[loglikelihood == max(loglikelihood)]),
                 label = "MPLE\nMILE",
                 color = "red3") +
  ylab("Log-Likelihood") +
  scale_x_continuous(limits = c(0.4, 0.4817),
                     expand = c(0, 0)) +
  scale_y_continuous(expand = c(0.1, 0)) +
  xlab(expression(psi)) +
  theme_minimal() +
  theme(axis.line = element_line())

fitted_models %>% 
  mutate(RSE = map_dbl(model, \(x) x$s)) %>% 
  select(-model)
  










