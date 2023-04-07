true_beta1 = 5
true_beta2 = 8
true_sigma = 10
n = 101

# create independent x-values
set.seed(1492)
x1 =(-(n-1)/2):((n-1)/2)
x2 = sample(x1)*0.05 + 0.5

# create dependent y-values according to beta1*x1 + beta2*x2 + N(0, sigma)
y = true_beta1*x1 + true_beta2*x2 + rnorm(n = n, mean = 0, sd = true_sigma)

plot(x1, y, main = "X1 vs Y")
plot(x2, y, main = "X2 vs Y")

likelihood = function(param) {
  
  beta1 = param[1]
  beta2 = param[2]
  sigma = param[3]
  
  pred = beta1*x1 + beta2*x2
  single_likelihoods = dnorm(y, mean = pred, sd = sigma, log = T)
  sum_ll = sum(single_likelihoods)
  return(sum_ll)
}

# Plot the likelihood profile of beta1
beta1_values = function(x) {return(likelihood(c(x, true_beta2, true_sigma)))}
beta1_likelihoods = lapply(seq(3, 7, by = .05), beta1_values)
plot(seq(3, 7, by = .05), beta1_likelihoods , type = "l", 
     xlab = "Values of parameter beta1", ylab = "Log likelihood")

# Plot the likelihood profile of beta2
beta2_values = function(x) {return(likelihood(c(true_beta1, x, true_sigma)))}
beta2_likelihoods = lapply(seq(6, 10, by = .05), beta2_values)
plot(seq(6, 10, by = .05), beta2_likelihoods , type = "l", 
     xlab = "Values of parameter beta2", ylab = "Log likelihood")

# Prior distribution
prior = function(param) {
  
  beta1 = param[1]
  beta2 = param[2]
  sigma = param[3]
  
  beta1_prior = dunif(beta1, min = 0, max = 10, log = TRUE)
  beta2_prior = dunif(beta1, min = 0, max = 10, log = TRUE)
  sigma_prior = dunif(sigma, min = 0, max = 30, log = TRUE)
  return(beta1_prior + beta2_prior + sigma_prior)
}

# Posterior distribution
posterior = function(param) {return(likelihood(param) + prior(param))
}

######## Metropolis algorithm ################
propose = function(param) {return(rnorm(3, mean = param, sd = c(0.1, 1, 1)))
}

run_metropolis_MCMC = function(start_value, iterations) {
  
  chain = array(dim = c(iterations + 1, 3))
  chain[1,] = start_value
  
  for(i in 1:iterations){
    proposal = propose(chain[i,])
    prob = exp(posterior(proposal) - posterior(chain[i,]))
    
    if(runif(1) < prob) {
      
      chain[i+1,] = proposal
    }
    else {
      
      chain[i+1,] = chain[i,]
    }
  }
  return(chain)
}

start_value = c(7, 10, 30)
chain = run_metropolis_MCMC(start_value, 10^5)

burnIn = 10000
acceptance = 1 - mean(duplicated(chain[-(1:burnIn),]))

### Summary: #######################
# Beta1 posterior
hist(chain[-(1:burnIn), 1], nclass = 30, , main = "Posterior of beta1", xlab = "True value = red line" )
abline(v = mean(chain[-(1:burnIn), 1]))
abline(v = true_beta1, col = "red" )

# Beta2 posterior 
hist(chain[-(1:burnIn), 2], nclass = 30, , main = "Posterior of beta2", xlab = "True value = red line" )
abline(v = mean(chain[-(1:burnIn), 2]))
abline(v = true_beta2, col = "red" )

# Sigma posterior
hist(chain[-(1:burnIn), 3], nclass = 30, main = "Posterior of sigma", xlab = "True value = red line")
abline(v = mean(chain[-(1:burnIn), 3]))
abline(v = true_sigma, col = "red" )

# Beta1 chain values
plot(chain[-(1:burnIn), 1], type = "l", xlab = "True value = red line" , main = "Chain values of beta1", )
abline(h = true_beta1, col = "red" )

# Beta2 chain values
plot(chain[-(1:burnIn), 2], type = "l", xlab = "True value = red line" , main = "Chain values of beta2", )
abline(h = true_beta2, col = "red" )

# Sigma chain values
plot(chain[-(1:burnIn), 3], type = "l", xlab = "True value = red line" , main = "Chain values of sigma", )
abline(h = true_sigma, col = "red" )

# for comparison:
summary(lm(y ~ -1 + x1 + x2))



