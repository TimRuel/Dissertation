library(plyr)
library(tidyverse)
library(LaplacesDemon)
library(doParallel)
cl <- makeCluster(6)
registerDoParallel(cl)

source("utils.R")

data <- c(1, 1, 2, 4, 7, 10)

step_size <- 0.01

psi_grid <- data |> 
  length() |> 
  log() |> 
  round_any(step_size, ceiling) |> 
  seq(0, to = _, step_size)

R <- 250

u <- rdirichlet(R, rep(1, length(data))) 









x <- iris[which(iris[,5] != "setosa"), c(1,5)]
trials <- 10000
ptime <- system.time({
  r <- foreach(icount(trials), .combine=cbind) %dopar% {
    ind <- sample(100, 100, replace=TRUE)
    result1 <- glm(x[ind,2]~x[ind,1], family=binomial(logit))
    coefficients(result1)
    }
  })[3]
ptime

stime <- system.time({
  r <- foreach(icount(trials), .combine=cbind) %do% {
    ind <- sample(100, 100, replace=TRUE)
    result1 <- glm(x[ind,2]~x[ind,1], family=binomial(logit))
    coefficients(result1)
    }
  })[3]
stime

stopCluster(cl)





cl <- makeCluster(6, type="SOCK")
registerDoSNOW(cl)
m <- matrix(rnorm(90000), 30000, 3)
foreach(i=1:nrow(m), .combine=rbind) %dopar% (m[i,] / mean(m[i,]))
stopCluster(cl)




library(doSNOW)
library(itertools)
cl <- makeCluster(6, type="SOCK")
registerDoSNOW(cl)
m <- matrix(rnorm(90000), 30000, 3)

foreach(s=isplitRows(m, chunks=getDoParWorkers()), .combine='rbind') %dopar% {
  t(apply(s, 1, function(x) x / mean(x)))
}