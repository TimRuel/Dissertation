true_beta = 5
true_sigma = 10
n = 101

# create independent x-values
x =(-(n-1)/2):((n-1)/2)

# create dependent y-values according to beta*x + N(0, sigma)
y = true_beta * x + rnorm(n = n, mean = 0, sd = true_sigma)

plot(x, y, main = "Test Data")


likelihood = function(param) {
  
  beta = param[1]
  sigma = param[2]
  
  pred = beta*x
  single_likelihoods = dnorm(y, mean = pred, sd = sigma, log = T)
  sum_ll = sum(single_likelihoods)
  return(sum_ll)
}

# Example: plot the likelihood profile of the slope beta
slope_values = function(x) {return(likelihood(c(x, true_sigma)))}
slope_likelihoods = lapply(seq(3, 7, by = .05), slope_values)
plot(seq(3, 7, by = .05), slope_likelihoods , type = "l", xlab = "values of slope parameter beta", ylab = "Log likelihood")


# Prior distribution
prior = function(param) {
  
  beta = param[1]
  sigma = param[2]
  
  beta_prior = dunif(beta, min = 0, max = 10, log = T)
  sigma_prior = dunif(sigma, min = 0, max = 30, log = T)
  return(beta_prior + sigma_prior)
}

# Posterior distribution
posterior = function(param) {return(likelihood(param) + prior(param))
}

######## Metropolis algorithm ################
propose = function(param) {return(rnorm(2, mean = param, sd = c(0.1, 0.3)))
}

run_metropolis_MCMC = function(start_value, iterations) {
  
  chain = array(dim = c(iterations + 1, 2))
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

start_value = c(7, 30)
chain = run_metropolis_MCMC(start_value, 100000)

burnIn = 10000
acceptance = 1 - mean(duplicated(chain[-(1:burnIn),]))

### Summary: #######################
par(mfrow = c(2, 2))
hist(chain[-(1:burnIn), 1], nclass = 30, , main = "Posterior of beta", xlab = "True value = red line" )
abline(v = mean(chain[-(1:burnIn), 1]))
abline(v = true_beta, col = "red" )
hist(chain[-(1:burnIn), 2], nclass = 30, main = "Posterior of sigma", xlab = "True value = red line")
abline(v = mean(chain[-(1:burnIn), 2]))
abline(v = true_sigma, col = "red" )
plot(chain[-(1:burnIn), 1], type = "l", xlab = "True value = red line" , main = "Chain values of beta", )
abline(h = true_beta, col = "red" )
plot(chain[-(1:burnIn), 2], type = "l", xlab = "True value = red line" , main = "Chain values of sigma", )
abline(h = true_sigma, col = "red" )

# for comparison:
summary(lm(y ~ -1 + x))

