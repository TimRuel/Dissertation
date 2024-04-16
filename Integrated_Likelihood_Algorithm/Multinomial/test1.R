library(doSNOW)
library(progress)

test <- function(data, psi_grid) {
  
  theta_MLE <- data / sum(data)
  
  psi_MLE <- PoI_fn(theta_MLE)
  
  u <- LaplacesDemon::rdirichlet(1, rep(1, length(data)))
  
  omega_hat <- get_omega_hat(u, psi_MLE)
  
  L <- psi_grid |> 
    purrr::map(get_theta_hat, omega_hat) |> 
    sapply(likelihood, data)
  
  return(L)
}

iterations <- 250
pb <- progress_bar$new(
  format = "letter = :letter [:bar] :elapsed | eta: :eta",
  total = iterations,    # 100 
  width = 60)

# progress_letter <- rep(LETTERS[1:10], 10)  # token reported in progress bar
# 
# # allowing progress bar to be used in foreach -----------------------------
progress <- function(n){
  pb$tick()
}



cl <- makeCluster(12)
registerDoSNOW(cl)

opts <- list(progress = progress)
a <- Sys.time()
result <- foreach(i = 1:iterations, .combine = rbind, 
                  .options.snow = opts) %dopar%
  {
    s <- test(data, psi_grid)
    return(s)
  }
close(pb)
stopCluster(cl)

b <- Sys.time()
b-a