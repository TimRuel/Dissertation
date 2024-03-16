library(tidyverse)
library(nloptr)
library(LaplacesDemon)
library(DescTools)

n <- c(1, 1, 2, 4, 7, 10)

m <- length(n)

# Define MLE for theta
theta.hat <- n / sum(n)

# Define parameter of interest function
g <- function(x) {
  
  y <- ifelse(x != 0, x*log(x), 0)
  
  return(-sum(y))
}

# Define MLE for parameter of interest
psi_hat <- g(theta.hat)

# Define likelihood function
likelihood <- function(theta, n) prod(theta^n)

R <- 250

omega_hat <- list()

for (j in 1:R) {
  
  # Draw random variartes from the m-dimensional probability simplex
  u <- rdirichlet(1, rep(1, m)) 
  # u <- rexp(m) 
  # u <- u / sum(u)
  
  # Define first objective function to be minimized
  entropy1 <- function(omega) -sum(u*log(omega))
  
  # Define first constraint function
  constraint1 <- function(omega) c(sum(omega) - 1, g(omega) - psi_hat)
  
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
psi1 <- seq(0, log(m), 0.01)

# Initialize vector for holding values of integrated likelihood function
L_bar <- c()

L_p <- c()

# Initialize progress bar for for loop
pb = txtProgressBar(min = 0, max = length(psi1), initial = 0, style = 3) 

# Calculate value of integrated likelihood for each parameter of interest value
for (i in 1:length(psi1)) {
  
  # Update progress bar
  setTxtProgressBar(pb, i)
  
  # Initialize vector for holding values of likelihood function
  L <- c()
  
  L1 <- c()
  
  # Initialize starting point for searching for optimal theta value
  theta_hat0 <- sample(1:10, m)
  theta_hat0 <- theta_hat0 / sum(theta_hat0) 
  
  for (j in 1:R) {
    
    # Define second objective function to be minimized
    entropy2 <- function(theta) -sum(omega_hat[[j]]*log(theta))
    
    # Define second constraint function
    constraint2 <- function(theta) c(sum(theta) - 1, g(theta) - psi1[i])
    
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
    
    # Define third constraint function
    constraint3 <- function(theta) c(sum(theta) - 1, g(theta) - psi1[i])
    
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
  
  L_bar[i] <- mean(L)
  
  L_p[i] <- mean(L1)
}

likelihood_dat <- data.frame(psi = psi1, 
                             L_bar = L_bar / max(L_bar), 
                             L_p = L_p / max(L_p)) %>% 
  mutate(Integrated = log(L_bar),
         Profile = log(L_p))

likelihood_dat %>% 
  pivot_longer(cols = c("Integrated", "Profile"),
               names_to = "Pseudolikelihood",
               values_to = "loglikelihood") %>% 
  ggplot() +
  scale_y_continuous(limits = c(-3.5, 0.3)) +
  scale_x_continuous(limits = c(1.1, 1.8)) +
  #geom_point(aes(x = psi, y = loglikelihood, color = Pseudolikelihood)) +
  geom_smooth(aes(x = psi, y = loglikelihood, color = Pseudolikelihood),
              se = FALSE,
              linewidth = 0.9,
              fullrange = TRUE) +
  theme_minimal() +
  theme(legend.position = c(0.5, 0.3),
        legend.background = element_rect())

log_likelihood_vals <- data.frame(psi = psi1, 
                                  Integrated = log(L_bar),
                                  Profile = log(L_p))

fit_IL <- loess(Integrated ~ psi, data = log_likelihood_vals)
fit_P <- loess(Profile ~ psi, data = log_likelihood_vals)

crit <- qchisq(0.95, 1) / 2

x <- optimize(
  function(psi) predict(fit_IL, psi), 
  lower = psi1 %>% head(1), 
  upper = psi1 %>% tail(1), 
  maximum = TRUE
)$objective

psi_max <- optimize(
  function(psi) predict(fit_IL, psi), 
  lower = psi1 %>% head(1), 
  upper = psi1 %>% tail(1), 
  maximum = TRUE
)$maximum

l <- uniroot(function(psi) predict(fit_IL, psi) - x + crit,
             interval = c(0, psi_max))$root

u <- uniroot(function(psi) predict(fit_IL, psi) - x + crit,
             interval = c(psi_max, psi1 %>% tail(1)))$root

curve <- function(psi) predict(fit_IL, psi) - x + crit

c(l, u)

ggplot() +
  geom_function(fun = curve) +
  scale_x_continuous(limits = c(0,2))

x <- optimize(
  function(psi) predict(fit_P, psi), 
  lower = psi1 %>% head(1), 
  upper = psi1 %>% tail(1), 
  maximum = TRUE
)$objective

psi_max <- optimize(
  function(psi) predict(fit_P, psi), 
  lower = psi1 %>% head(1), 
  upper = psi1 %>% tail(1), 
  maximum = TRUE
)$maximum

l <- uniroot(function(psi) predict(fit_P, psi) - x + crit,
             interval = c(0, psi_max))$root

u <- uniroot(function(psi) predict(fit_P, psi) - x + crit,
             interval = c(psi_max, psi1 %>% tail(1)))$root

curve <- function(psi) predict(fit_P, psi) - x + crit

ggplot() +
  geom_function(fun = curve) +
  scale_x_continuous(limits = c(0,2))

