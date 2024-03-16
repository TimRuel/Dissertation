library(tidyverse)
library(nloptr)
library(LaplacesDemon)
library(Rmpfr)
library(geomtextpath)
library(viridis)
library(ggnewscale)
library(purrr)

set.seed(1996)

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
likelihood <- function(theta) prod(theta^n)

# Define values for parameter of interest at which to evaluate the integrated likelihood  
psi <- seq(0, log(m), 0.01)

# Number of replications for each value of psi
R <- 250

u <- rdirichlet(R, rep(1, m)) 

# Initialize matrix for holding each iteration's values for likelihood ratio at each value of psi 
L <- mpfrArray(NA, precBits = 106, dim = c(R, length(psi)))

# Initialize vector for holding values of profile likelihood
L_p <- mpfrArray(NA, precBits = 106, dim = length(psi))

omega_hat <- list()

# Initialize progress bar for for loop
pb = txtProgressBar(min = 0, max = R, initial = 0, style = 3) 

for (i in 1:R) {
  
  # Update progress bar
  setTxtProgressBar(pb, i)
  
  # Find value of omega that minimizes objective function subject to constraints
  omega_hat[[i]] <- auglag(x0 = rdirichlet(1, rep(1, m)),
                           fn = function(omega) -sum(u[i,]*log(omega)),
                           heq = function(omega) c(sum(omega) - 1, g(omega) - psi_hat),
                           lower = rep(0, m))$par
  
  for (j in 1:length(psi)) {
    
    # Find value of theta that minimizes objective function subject to constraints
    theta_hat <- auglag(x0 = rdirichlet(1, rep(1, m)),
                        fn = function(theta) -sum(omega_hat[[i]]*log(theta)),
                        heq = function(theta) c(sum(theta) - 1, g(theta) - psi[j]),
                        lower = rep(0, m))$par
    
    # Calculate ratio of likelihood at optimal theta to likelihood at initial random draw for theta
    L[i, j] <- likelihood(theta_hat)
    
    if (i == 1) {
      
      # Find value of theta for obtaining profile log-likelihood
      theta_hat_p <- auglag(x0 = rdirichlet(1, rep(1, m)),
                            fn = function(theta) -sum(theta.hat*log(theta)),
                            heq = function(theta) g(theta) - psi[j],
                            lower = rep(0, m))$par
      
      # Calculate value of profile likelihood for current value of psi
      L_p[j] <- likelihood(theta_hat_p)
    }
  }
}

log_likelihood_vals <- data.frame(psi = psi, 
                              Integrated = L %>% 
                                apply(2, mean) %>% 
                                log() %>% 
                                as.double(),
                              Profile = L_p %>% 
                                log() %>% 
                                as.double())

log_likelihood_vals %>% 
  mutate(Integrated = Integrated - max(Integrated),
         Profile = Profile - max(Profile)) %>% 
  pivot_longer(cols = c("Integrated", "Profile"),
               names_to = "Pseudolikelihood",
               values_to = "loglikelihood") %>% 
  ggplot() +
  # scale_y_continuous(limits = c(-3, 0)) +
  geom_smooth(aes(x = psi, y = loglikelihood, color = Pseudolikelihood),
              se = FALSE,
              linewidth = 0.9,
              fullrange = TRUE) +
  theme_minimal() +
  theme(legend.position = c(0.2, 0.8),
        legend.background = element_rect())

fit_IL <- loess(Integrated ~ psi, data = log_likelihood_vals)
fit_P <- loess(Profile ~ psi, data = log_likelihood_vals)

crit <- qchisq(0.95, 1) / 2

x <- optimize(
  function(psi) predict(fit_IL, psi), 
  lower = psi %>% head(1), 
  upper = psi %>% tail(1), 
  maximum = TRUE
)$objective

psi_max <- optimize(
  function(psi) predict(fit_IL, psi), 
  lower = psi %>% head(1), 
  upper = psi %>% tail(1), 
  maximum = TRUE
)$maximum

l <- uniroot(function(psi) predict(fit_IL, psi) - x + crit,
             interval = c(0, psi_max))$root

u <- uniroot(function(psi) predict(fit_IL, psi) - x + crit,
             interval = c(psi_max, psi1 %>% tail(1)))$root

curve <- function(psi) predict(fit_IL, psi) - x + crit

c(l, u)


l <- uniroot(function(psi) predict(fit_P, psi) - log(likelihood(theta.hat)) + crit,
             interval = c(0, psi1))$root

u <- uniroot(function(psi) predict(fit_P, psi) - x + crit,
             interval = c(psi1, psi %>% tail(1))$root
             
c(l, u)


curve <- function(psi) predict(fit_IL, psi) - x + crit

ggplot() +
  geom_function(fun = curve) +
  scale_x_continuous(limits = c(0,2))

log(likelihood(theta.hat))







