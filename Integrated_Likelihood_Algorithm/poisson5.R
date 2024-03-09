library(tidyverse)
library(nloptr)
library(LaplacesDemon)
library(Rmpfr)
library(geomtextpath)
library(viridis)
library(ggnewscale)
library(purrr)

set.seed(1996)

# Define dimension of parameter
d <- 20

# Define true values of full model parameter
theta_0 <- runif(d) %>% round(2)

# Define sample sizes from each population
n <- sample(4:5, d, replace = TRUE)

# Define observed data from each population
x <- mapply(rpois, n, theta_0)

# Define weights for PoI function
w <- runif(d) %>% round(2)

# Define hyperparameters for u's
alpha <- rep(1, d)
beta <- rep(1, d)

# Define likelihood function
likelihood <- function(theta) {
  
  theta <- theta %>% mpfr(precBits = 106)
  
  sum_x <- x %>% sapply(sum)
  
  (-theta * n) %>% 
    sum() %>% 
    exp() %>% 
    "*"((theta^sum_x) %>%
          prod()) 
}

# Define parameter of interest function 
g <- function(theta) sum(theta*w)

# Define parameter of interest
psi_0 <- g(theta_0)

# Define MLE for theta
theta_hat <- x %>% sapply(mean)

# Define MLE for parameter of interest
psi_hat <- g(theta_hat)

# Calculate standard error of MLE
psi_hat_se <- x %>% 
  sapply(var) %>% 
  "*"(w^2 / n) %>% 
  sum() %>% 
  sqrt()

# Calculate margin of error
n_std_errors = 3
MoE = n_std_errors * psi_hat_se

# Define values for parameter of interest at which to evaluate the integrated likelihood  
psi <- seq(max(0, psi_hat - MoE), psi_hat + MoE, 0.001)

# Define log-likelihood expectation function to be minimized
E_log_like <- function(theta, omega) sum((-theta + log(theta)*omega)*n)

# Number of replications for each value of psi
R <- 10

# Initialize matrix for holding each iteration's values for likelihood ratio at each value of psi 
L_ratio <- mpfrArray(NA, precBits = 106, dim = c(R, length(psi)))

# Initialize vector for holding values of profile likelihood
L_p <- mpfrArray(NA, precBits = 106, dim = length(psi))

# Initialize progress bar for for loop
pb = txtProgressBar(min = 0, max = R, initial = 0, style = 3) 

for (i in 1:R) {
  
  # Update progress bar
  setTxtProgressBar(pb, i)
  
  # Random draw from posterior for theta
  u <- x %>% 
    map2_dbl(alpha, sum) %>% 
    mapply(rgamma, n = 1, shape = ., rate = n + beta)
  
  # Find value of omega that minimizes distance from u subject to constraints
  Q <- auglag(x0 = u,
              fn = function(omega) dist(matrix(c(u, omega), nrow = 2, byrow = TRUE))[1],
              heq = function(omega) g(omega) - psi_hat,
              lower = rep(0, d))$par
  
  for (j in 1:length(psi)) {
    
    # Find value of theta that minimizes log-likelihood expectation function subject to constraints
    T_psi <- auglag(x0 = rgamma(d, 1, 1),
                    fn = function(theta) -E_log_like(theta, Q),
                    heq = function(theta) g(theta) - psi[j],
                    lower = rep(0, d))$par
    
    # Calculate ratio of likelihood at optimal theta to likelihood at initial random draw for theta
    L_ratio[i, j] <- likelihood(T_psi) / likelihood(u)
    
    if (i == 1) {
      
      # Find value of theta for obtaining profile log-likelihood
      theta_hat_p <- auglag(x0 = rgamma(d, 1, 1),
                            fn = function(theta) -E_log_like(theta, theta_hat),
                            heq = function(theta) g(theta) - psi[j],
                            lower = rep(0, d))$par
      
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

log_likelihood_vals_tidy_shifted <- fitted_models %>% 
  ungroup() %>% 
  pull(model) %>% 
  sapply(
    function(mod) {
      optimize(
        function(val) predict(mod, val), 
        lower = psi %>% head(1), 
        upper = psi %>% tail(1), 
        maximum = TRUE
      )
      }
    ) %>% 
  t() %>% 
  data.frame() %>% 
  mutate(Pseudolikelihood = c("Integrated", "Profile"),
         MLE = unlist(maximum),
         maximum = unlist(objective)) %>% 
  merge(log_likelihood_vals_tidy) %>% 
  mutate(loglikelihood = loglikelihood - maximum) %>% 
  select(Pseudolikelihood, psi, loglikelihood, MLE, maximum)

MLE_data <- log_likelihood_vals_tidy_shifted %>% 
  group_by(Pseudolikelihood) %>% 
  distinct(MLE, maximum) %>% 
  ungroup() %>% 
  mutate(MLE_label = c("hat(psi)[IL]", "hat(psi)[P]")) %>% 
  add_row(data.frame("MLE" = psi_0, "MLE_label" = "psi[0]"))

ggplot() +
  stat_smooth(aes(x = psi, 
                  y = loglikelihood, 
                  color = Pseudolikelihood, 
                  linetype = Pseudolikelihood,
                  label = Pseudolikelihood),
              data = log_likelihood_vals_tidy_shifted,
              linewidth = 0.7,
              geom = "textpath",
              se = FALSE,
              hjust = 0.1,
              show.legend = FALSE) +
  geom_hline(yintercept = 0,
             linetype = 5,
             show.legend = FALSE) +
  geom_labelvline(aes(xintercept = MLE,
                      label = MLE_label,
                      color = MLE_label,
                      angle = 0),
                  data = MLE_data,
                  parse = TRUE,
                  angle = -90,
                  show.legend = FALSE) +
  ylab("Log-Likelihood") +
  # scale_x_continuous(limits = c(2.5, 4.5),
  #                    expand = c(0, 0)) +
  # scale_y_continuous(expand = c(0.1, 0)) +
  scale_color_brewer(palette = "Set1") +
  xlab(expression(psi)) +
  theme_minimal() +
  theme(axis.line = element_line())

fitted_models %>% 
  mutate(RSE = map_dbl(model, \(x) x$s)) %>% 
  select(-model)







