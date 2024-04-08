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

g <- function(x) {
  
  y <- ifelse(x != 0, x*log(x), 0)
  
  return(-sum(y))
}

get_multinom_entropy_IL <- function(n, step_size, R) {
  
  likelihood <- function(theta) prod(theta^n)
  
  m <- length(n)
  
  theta.hat <- n / sum(n)
  
  psi_hat <- g(theta.hat)
  
  psi1 <- seq(0, round_any(log(m), step_size, ceiling), step_size)
  
  N <- length(psi1)
  
  N_lower <- psi1[psi1 <= psi_hat] %>% length()
  
  u <- rdirichlet(R, rep(1, m)) 
  
  L <- mpfrArray(NA, precBits = 106, dim = c(R, length(psi1)))
  
  omega_hat <- list()
  
  pb = txtProgressBar(min = 0, max = R, initial = 0, style = 3) 
  
  for (i in 1:R) {
    
    setTxtProgressBar(pb, i)
    
    omega_hat[[i]] <- auglag(x0 = u[i,],
                             fn = function(omega) -sum(u[i,]*log(omega)),
                             heq = function(omega) c(sum(omega) - 1, g(omega) - psi_hat),
                             lower = rep(0, m))$par
    
    theta_hat <- omega_hat[[i]]
    
    for (j in N_lower:1) {
      
      theta_hat <- auglag(x0 = theta_hat,
                          fn = function(theta) -sum(omega_hat[[i]]*log(theta)),
                          heq = function(theta) c(sum(theta) - 1, g(theta) - psi1[j]),
                          lower = rep(0, m))$par
      
      L[i, j] <- likelihood(theta_hat)
    }
    
    theta_hat <- omega_hat[[i]]
    
    for (j in (N_lower + 1):N) {
      
      theta_hat <- auglag(x0 = theta_hat,
                          fn = function(theta) -sum(omega_hat[[i]]*log(theta)),
                          heq = function(theta) c(sum(theta) - 1, g(theta) - psi1[j]),
                          lower = rep(0, m))$par
      
      L[i, j] <- likelihood(theta_hat)
    }
  }
  
  return(L)
}

data <- c(1, 1, 2, 4, 7, 10)

sims <- rmultinom(2, length(data), data)

test <- sims %>% 
  t() %>% 
  list() %>% 
  lapply(get_multinom_entropy_IL, 0.05, 250)



