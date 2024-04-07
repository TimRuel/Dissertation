library(tidyverse)
library(plyr)
library(nloptr)
library(LaplacesDemon)
library(Rmpfr)
library(geomtextpath)
library(viridis)
library(ggnewscale)
library(purrr)
library(zeallot)

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
step_size <- 0.05
psi1 <- seq(0, round_any(log(m), step_size, ceiling), step_size)

N <- length(psi1)

psi1_lower <- psi1[psi1 <= psi_hat] %>% rev()

N_lower <- length(psi1_lower)

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
    theta_hat <- auglag(x0 = omega_hat[[i]],
                        fn = function(theta) -sum(omega_hat[[i]]*log(theta)),
                        heq = function(theta) c(sum(theta) - 1, g(theta) - psi1[j]),
                        lower = rep(0, m))$par
    
    # Calculate ratio of likelihood at optimal theta to likelihood at initial random draw for theta
    L[i, j] <- likelihood(theta_hat)
  }
}

# Initialize vector for holding values of profile likelihood
L_p <- mpfrArray(NA, precBits = 106, dim = N)

# Initialize progress bar for for loop
pb = txtProgressBar(min = 0, max = N, initial = 0, style = 3) 

theta_hat_p <- theta.hat

for (j in N_lower:1) {
  
  # Update progress bar
  setTxtProgressBar(pb, j)
  
  # Find value of theta that minimizes objective function subject to constraints
  theta_hat_p <- auglag(x0 = theta_hat_p,
                        fn = function(theta) -sum(theta.hat*log(theta)),
                        heq = function(theta) c(sum(theta) - 1, g(theta) - psi1[j]),
                        lower = rep(0, m))$par
  
  L_p[j] <- likelihood(theta_hat_p)
}

theta_hat_p <- theta.hat

for (j in (N_lower + 1):N) {
  
  # Update progress bar
  setTxtProgressBar(pb, j)
  
  # Find value of theta that minimizes objective function subject to constraints
  theta_hat_p <- auglag(x0 = theta_hat_p,
                        fn = function(theta) -sum(theta.hat*log(theta)),
                        heq = function(theta) c(sum(theta) - 1, g(theta) - psi1[j]),
                        lower = rep(0, m))$par
  
  L_p[j] <- likelihood(theta_hat_p)
}

log_likelihood_vals <- data.frame(psi = psi1,
                                  Integrated = L %>%
                                    apply(2, mean) %>%
                                    log() %>%
                                    as.double(),
                                  Profile = L_p %>% 
                                    log() %>% 
                                    as.double())

log_likelihood_vals_tidy <- log_likelihood_vals %>% 
  pivot_longer(cols = c("Integrated", "Profile"),
               names_to = "Pseudolikelihood",
               values_to = "loglikelihood")

spline_fitted_models <- log_likelihood_vals_tidy %>%
  group_by(Pseudolikelihood) %>% 
  group_map(~ smooth.spline(.x$psi, .x$loglikelihood)) %>% 
  set_names(c("Integrated", "Profile"))

MLE_data <- spline_fitted_models %>%
  sapply(
    function(mod) {
      optimize(
        function(psi) predict(mod, psi)$y, 
        lower = psi1 %>% head(1), 
        upper = psi1 %>% tail(1), 
        maximum = TRUE
      )}) %>% 
  t() %>% 
  data.frame() %>% 
  rownames_to_column("Pseudolikelihood") %>% 
  dplyr::rename(MLE = maximum,
         Maximum = objective) %>% 
  mutate(MLE_label = c("hat(psi)[IL]", "hat(psi)[P]"))

c(IL_curve, P_curve) %<-% mapply(
  function(mod, maximum) function(psi) predict(mod, psi)$y - maximum,
  spline_fitted_models,
  MLE_data$Maximum)

ggplot() +
  stat_function(fun = P_curve,
                geom = "textpath",
                label = "Profile",
                aes(color = "Profile"),
                linewidth = 1,
                hjust = 0.1,
                show.legend = FALSE) +
  stat_function(fun = IL_curve,
                geom = "textpath",
                label = "Integrated",
                aes(color = "Integrated"),
                linewidth = 1,
                hjust = 0.1,
                show.legend = FALSE) +
  geom_hline(yintercept = 0,
             linetype = 5) +
  geom_labelvline(aes(xintercept = as.numeric(MLE),
                      label = MLE_label,
                      color = Pseudolikelihood),
                  data = MLE_data,
                  parse = TRUE,
                  show.legend = FALSE) +
  ylab("Log-Likelihood") +
  scale_x_continuous(limits = c(1, 2),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(-4, 0.1),
                     expand = c(0.1, 0)) +
  scale_color_brewer(palette = "Set1") +
  xlab(expression(psi)) +
  theme_minimal() +
  theme(axis.line = element_line())

crit <- qchisq(0.95, 1) / 2

c(psi_hat_IL, psi_hat_P) %<-% MLE_data$MLE

CI_lower_P <- uniroot(function(psi) P_curve(psi) + crit,
                    interval = c(psi1 %>% head(1), psi_hat_P))$root %>% 
  round(3)

CI_upper_P <- uniroot(function(psi) P_curve(psi) + crit,
             interval = c(psi_hat_P, psi1 %>% tail(1)))$root %>% 
  round(3)

CI_lower_IL <- uniroot(function(psi) IL_curve(psi) + crit,
                      interval = c(psi1 %>% head(1), psi_hat_IL))$root %>% 
  round(3)

CI_upper_IL <- uniroot(function(psi) IL_curve(psi) + crit,
                      interval = c(psi_hat_IL, psi1 %>% tail(1)))$root %>% 
  round(3)

data.frame(MLE = c(psi_hat_IL, psi_hat_P) %>% round(3),
           CI_95 = c(paste0("(", CI_lower_IL, ", ", CI_upper_IL, ")"),
                     paste0("(", CI_lower_P, ", ", CI_upper_P, ")")),
           row.names = c("Integrated", "Profile"))
