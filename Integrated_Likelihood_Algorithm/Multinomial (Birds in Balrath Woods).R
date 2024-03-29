library(tidyverse)
library(nloptr)
library(LaplacesDemon)
library(Rmpfr)
library(geomtextpath)
library(viridis)
library(ggnewscale)
library(purrr)
library(zeallot)

#set.seed(1996)

n <- c(1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 6, 8)

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
psi1 <- seq(0, log(m), 0.001)

# Number of replications for each value of psi
R <- 250

u <- rdirichlet(R, rep(1, m)) 

# Initialize matrix for holding each iteration's values for likelihood ratio at each value of psi 
L <- mpfrArray(NA, precBits = 106, dim = c(R, length(psi1)))

omega_hat <- list()

# Initialize progress bar for for loop
pb = txtProgressBar(min = 0, max = R, initial = 0, style = 3) 

for (i in 1:R) {
  
  # Update progress bar
  setTxtProgressBar(pb, i)
  
  # Find value of omega that minimizes objective function subject to constraints
  omega_hat[[i]] <- auglag(x0 = u[i,],
                           fn = function(omega) -sum(u[i,]*log(omega)),
                           heq = function(omega) c(sum(omega) - 1, g(omega) - psi_hat),
                           lower = rep(0, m))$par
  
  for (j in 1:length(psi1)) {
    
    # Find value of theta that minimizes objective function subject to constraints
    theta_hat <- auglag(x0 = rdirichlet(1, rep(1, m)),
                        fn = function(theta) -sum(omega_hat[[i]]*log(theta)),
                        heq = function(theta) c(sum(theta) - 1, g(theta) - psi1[j]),
                        lower = rep(0, m))$par
    
    # Calculate ratio of likelihood at optimal theta to likelihood at initial random draw for theta
    L[i, j] <- likelihood(theta_hat)
  }
}

# Initialize vector for holding values of profile likelihood
L_p <- mpfrArray(NA, precBits = 106, dim = length(psi1))

# Initialize progress bar for for loop
pb = txtProgressBar(min = 0, max = length(psi1), initial = 0, style = 3) 

for (j in 1:length(psi1)) {
  
  # Update progress bar
  setTxtProgressBar(pb, j)
  
  # Find value of theta that minimizes objective function subject to constraints
  theta_hat_p <- auglag(x0 = theta.hat,
                        fn = function(theta) -sum(theta.hat*log(theta)),
                        heq = function(theta) c(sum(theta) - 1, g(theta) - psi1[j]),
                        lower = rep(0, m))$par
  
  L_p[j] <- likelihood(theta_hat_p)
}

log_likelihood_vals <- data.frame(psi = psi1,
                                  # Integrated = L %>% 
                                  #   apply(2, mean) %>% 
                                  #   log() %>% 
                                  #   as.double(),
                                  Profile = L_p %>% 
                                    log() %>% 
                                    as.double())

# saveRDS(log_likelihood_vals, "balrath_woods_profile_log_likelihood_vals_0.001_step_size.Rda")
# log_likelihood_vals <- readRDS(file = "balrath_woods_profile_log_likelihood_vals_0.001_step_size.Rda")

# fit_P <- loess(Profile ~ psi, data = log_likelihood_vals)

spline_fit_P <- log_likelihood_vals %>% 
  with(smooth.spline(psi, Profile))

log_likelihood_vals %>% 
  #mutate(Profile = Profile - max(Profile)) %>% 
  # mutate(Integrated = Integrated - max(Integrated),
  #        Profile = Profile - max(Profile)) %>%
  # pivot_longer(cols = c("Integrated", "Profile"),
  #              names_to = "Pseudolikelihood",
  #              values_to = "loglikelihood") %>% 
  # filter(Pseudolikelihood == "Profile") %>%
  ggplot() +
  # scale_y_continuous(limits = c(-4, 1)) +
  scale_x_continuous(limits = c(2, 2.7)) +
  # geom_point(aes(x = psi, y = Profile)) +
  # geom_smooth(aes(x = psi, y = Profile),
  #             se = FALSE,
  #             linewidth = 0.9,
  #             fullrange = TRUE) +
  geom_function(fun = function(psi) predict(spline_fit_P, psi)$y) +
  theme_minimal() +
  theme(legend.position = c(0.2, 0.2),
        legend.background = element_rect())

crit <- qchisq(0.95, 1) / 2

c(l_p_maximizer, l_p_maximum) %<-% optimize(
  function(psi) predict(spline_fit_P, psi)$y, 
  lower = psi1 %>% head(1), 
  upper = psi1 %>% tail(1), 
  maximum = TRUE)

l <- uniroot(function(psi) predict(spline_fit_P, psi)$y - l_p_maximum + crit,
             interval = c(psi1 %>% head(1), l_p_maximizer))$root

u <- uniroot(function(psi) predict(spline_fit_P, psi)$y - l_p_maximum + crit,
             interval = c(l_p_maximizer, psi1 %>% tail(1)))$root

print("Profile")
c(l, u)

# fit_IL <- loess(Integrated ~ psi, data = log_likelihood_vals)
# 
# l_bar_psi_hat <- predict(fit_IL, psi_hat)
# 
# l_bar_maximum <- optimize(
#   function(psi) predict(fit_IL, psi),
#   lower = psi1 %>% head(1),
#   upper = psi1 %>% tail(1),
#   maximum = TRUE
# )$objective
# 
# l_bar_maximizer <- optimize(
#   function(psi) predict(fit_IL, psi), 
#   lower = psi1 %>% head(1), 
#   upper = psi1 %>% tail(1), 
#   maximum = TRUE
# )$maximum
# 
# l <- uniroot(function(psi) predict(fit_IL, psi) - l_bar_psi_hat + crit,
#              interval = c(psi1 %>% head(1), l_bar_maximizer))$root
# 
# u <- uniroot(function(psi) predict(fit_IL, psi) - l_bar_psi_hat + crit,
#              interval = c(l_bar_maximizer, psi1 %>% tail(1)))$root
# 
# print("Integrated")
# c(l, u)

# curve <- function(psi) predict(fit_IL, psi) - l_bar_psi_hat + crit
# curve(psi1 %>% head(1))
# curve(l_bar_maximizer)
# curve(psi1 %>% tail(1))
# 
# ggplot() + 
#   geom_function(fun = curve) +
#   scale_x_continuous(limits = c(0, 2)) +
#   scale_y_continuous(limits = c(-50, 10))

