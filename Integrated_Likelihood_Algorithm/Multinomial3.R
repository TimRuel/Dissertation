library(LaplacesDemon)

n <- c(1, 1, 2, 4, 7, 10)

m <- length(n)

u <- rdirichlet(1, rep(1, m)) 

ddirichlet(u, rep(1, m))
