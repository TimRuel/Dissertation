################################################################################
#################################### GENERAL ###################################
################################################################################

log_likelihood <- function(theta, data) sum(data * log(theta), na.rm = TRUE)

neg_log_likelihood <- function(theta, data) -log_likelihood(theta, data)

likelihood <- function(theta, data) exp(log_likelihood(theta, data))

entropy <- function(theta) -sum(theta * log(theta), na.rm = TRUE)

euclidean_distance <- function(omega, u) sqrt(sum((omega - u)^2))

get_psi_grid <- function(data, step_size, num_std_errors, split = FALSE) {
  
  n <- sum(data)
  
  m <- length(data)
  
  theta_MLE <- data / n
  
  psi_MLE <- entropy(theta_MLE)
  
  sigma <- theta_MLE*diag(m) - matrix(theta_MLE) %*% theta_MLE
  
  psi_MLE_SE <- sqrt(sum(matrix(1 + log(theta_MLE)) %*% (1 + log(theta_MLE)) * sigma, na.rm = TRUE) / n)
  
  MoE <- num_std_errors * psi_MLE_SE
  
  psi_grid <- (psi_MLE + MoE*c(-1, 1)) |> 
    (\(x) c(max(0, x[1]), min(log(m), x[2])))() |> 
    plyr::round_any(step_size, floor) |> 
    (\(x) seq(x[1], x[2], step_size))()
  
  if (split) {
    
    psi_grid_list <- psi_grid |> 
      split(factor(psi_grid > psi_MLE)) |> 
      purrr::modify_in(1, rev) |> 
      unname()
    
    return(psi_grid_list)
  }
  
  return(psi_grid)
}

get_omega_hat <- function(objective_fn, psi_MLE, u_params, tol, return_u = FALSE) {
  
  u <- LaplacesDemon::rdirichlet(1, u_params) |> 
    as.double()
  
  omega_hat <- nloptr::auglag(x0 = rep(1 / length(u), length(u)),
                              fn = function(omega) objective_fn(omega, u),
                              heq = function(omega) c(sum(omega) - 1, entropy(omega) - psi_MLE),
                              lower = rep(0, length(u)),
                              localsolver = "LBFGS")$par
  
  if (abs(entropy(omega_hat) - psi_MLE) < tol) {
    
    if (return_u) {
      
      return(rbind(u, omega_hat))
    }
    
    return(omega_hat)
  }
  
  else {
    
    get_omega_hat(objective_fn, psi_MLE, u_params, tol, return_u)
  }
}

get_theta_hat <- function(init_guess, psi, omega_hat) {
  
  fn <- function(theta) neg_log_likelihood(theta, omega_hat)
  gr <- function(theta) nloptr::nl.grad(theta, fn)
  heq <- function(theta) c(sum(theta) - 1, entropy(theta) - psi)
  heqjac <- function(theta) nloptr::nl.jacobian(theta, heq)
  
  tictoc::tic()
  
  out <- nloptr::auglag(x0 = init_guess,
                        fn = fn,
                        gr = gr,
                        heq = heq,
                        heqjac = heqjac,
                        lower = rep(0, length(omega_hat)),
                        localsolver = "LBFGS")
  
  if (out$convergence < 0) return(init_guess)
  
  theta_hat <- out$par
  
  return(theta_hat)
}

################################################################################
############################ INTEGRATED LIKELIHOOD ############################# 
################################################################################

custom_combine_IL <- function(...) {
  
  arglist <- list(...)
  
  data <- arglist[[1]]$data
  
  preallocations <- matrix(data = NA, nrow = length(arglist), ncol = length(data))
  
  for (result in arglist) {preallocations[result$i,] <- result$preallocations}
  
  return(list(data = data, preallocations = preallocations))
}

get_IL_preallocations <- function(data_sims, u_params, R, tol, chunk_size) {
  
  p <- progressr::progressor(steps = R * length(data_sims))
  
  IL_preallocations <-
    
    foreach(
      data = data_sims,
      .combine = "list",
      .multicombine = TRUE,
      .maxcombine = length(data_sims),
      .options.future = list(seed = TRUE,
                             chunk.size = chunk_size)
      
    ) %:%
    
    foreach(
      i = 1:R,
      .combine = custom_combine_IL,
      .multicombine = TRUE,
      .maxcombine = R,
      .options.future = list(seed = TRUE)
      
    ) %dofuture% {
      
      p()
      
      psi_MLE <- entropy(data / sum(data))
      
      omega_hat <- get_omega_hat(neg_log_likelihood, psi_MLE, u_params, tol, return_u = FALSE)
      
      return(list(data = data, preallocations = omega_hat, i = i))
    }
  
  return(IL_preallocations)
}

get_integrated_log_likelihood_vals.aux <- function(omega_hat, data, psi_grid_list) {
  
  l1 <- psi_grid_list |> 
    purrr::map(
      \(psi_grid) psi_grid |> 
        purrr::accumulate(
          \(acc, nxt) get_theta_hat(acc, nxt, omega_hat), 
          .init = omega_hat
          ) |> 
        magrittr::extract(-1) |> 
        purrr::map_dbl(log_likelihood, data)
      ) |> 
    purrr::modify_in(1, rev) |> 
    unlist()
  
  return(l1)
}

get_integrated_log_likelihood_vals <- function(data, step_size, num_std_errors, u_params, R = 250, tol = 0.0001, chunk_size) {
  
  p <- progressr::progressor(steps = R)
  
  psi_MLE <- entropy(data / sum(data))
  
  psi_grid_list <- data |> 
    get_psi_grid(step_size, num_std_errors, split = TRUE)
  
  foreach(
    
    i = 1:R,
    .combine = "rbind",
    .multicombine = TRUE,
    .maxcombine = R,
    .options.future = list(seed = TRUE,
                           chunk.size = chunk_size)
    
  ) %dofuture% {
    
    p()
    
    get_omega_hat(neg_log_likelihood, psi_MLE, u_params, tol, return_u = FALSE) |>
      get_integrated_log_likelihood_vals.aux(data, psi_grid_list)
  } |> 
    matrixStats::colLogSumExps() |>
    (`-`)(log(R))
}

get_integrated_log_likelihood_sims <- function(preallocations, step_size = 0.01, num_std_errors = 4, chunk_size) {
  
  R <- preallocations[[1]]$preallocations |> 
    nrow()
  
  p <- progressr::progressor(steps = R * length(preallocations))
  
  foreach(
    preallocation = preallocations,
    .combine = "list",
    .multicombine = TRUE,
    .maxcombine = length(preallocations),
    .options.future = list(seed = TRUE,
                           chunk.size = chunk_size)
    
  ) %:%
    
    foreach(
      i = 1:R,
      .combine = "rbind",
      .multicombine = TRUE,
      .maxcombine = R,
      .options.future = list(seed = TRUE)
      
    ) %dofuture% {
      
      p()
      
      data <- preallocation$data
      
      psi_grid_list <- get_psi_grid(data, step_size, num_std_errors, split = TRUE)
      
      omega_hat <- preallocation$preallocations[i,]
      
      get_integrated_log_likelihood_vals.aux(omega_hat, data, psi_grid_list)
    } |>
    map(\(x) x |> 
          matrixStats::colLogSumExps() |>
          (`-`)(log(R))
    )
}

################################################################################
######################### MODIFIED INTEGRATED LIKELIHOOD ####################### 
################################################################################

custom_combine_mod_IL <- function(...) {
  
  arglist <- list(...)
  
  preallocations <- list()
  
  for (result in arglist) {preallocations[[result$i]] <- result$preallocations}

  data <- arglist[[1]]$data

  return(list(data = data, preallocations = preallocations))
}

get_mod_IL_preallocations <- function(data_sims, Q, alpha, R, tol, chunk_size) {
  
  p <- progressr::progressor(steps = R * length(data_sims))
  
  mod_IL_preallocations <- 
    
    foreach(
      data = data_sims,
      .combine = "list",
      .multicombine = TRUE,
      .maxcombine = num_sims,
      .options.future = list(seed = TRUE,
                             chunk.size = chunk_size)
      
    ) %:%
    
    foreach(
      i = 1:R,
      .combine = custom_combine_mod_IL,
      .multicombine = TRUE,
      .maxcombine = R,
      .options.future = list(seed = TRUE)
      
    ) %dofuture% {
      
      p()
      
      psi_MLE <- entropy(data / sum(data))
      
      u_params <- data + alpha
      
      preallocations <- get_omega_hat(Q, psi_MLE, u_params, tol, return_u = TRUE)
      
      return(list(data = data, preallocations = preallocations, i = i))
    }
  
  return(mod_IL_preallocations)
}

get_mod_integrated_log_likelihood_vals <- function(data, Q, step_size, num_std_errors, u_params, R = 250, tol = 0.0001, chunk_size) {
  
  p <- progressr::progressor(steps = R)
  
  psi_MLE <- entropy(data / sum(data))
  
  psi_grid_list <- data |> 
    get_psi_grid(step_size, num_std_errors, split = TRUE)
  
  foreach(
    
    i = 1:R,
    .combine = "rbind",
    .multicombine = TRUE,
    .maxcombine = R,
    .options.future = list(seed = TRUE,
                           chunk.size = chunk_size)
    
  ) %dofuture% {
    
    p()
    
    mat <- get_omega_hat(Q, psi_MLE, u_params, tol, return_u = TRUE)
    
    u <- mat["u",]
    
    omega_hat <- mat["omega_hat",]
    
    log_like_u <- log_likelihood(u, data)
    
    omega_hat |> 
      get_integrated_log_likelihood_vals.aux(data, psi_grid_list) |> 
      (`-`)(log_like_u)
  } |> 
    matrixStats::colLogSumExps() |>
    (`-`)(log(R))
}

get_mod_integrated_log_likelihood_sims <- function(preallocations, step_size = 0.01, num_std_errors = 4, chunk_size) {
  
  R <- preallocations[[1]]$preallocations |> 
    length()
  
  p <- progressr::progressor(steps = R * length(preallocations))
  
  mod_integrated_log_likelihood_sims <-
    
    foreach(
      preallocation = preallocations,
      .combine = "list",
      .multicombine = TRUE,
      .maxcombine = length(preallocations),
      .options.future = list(seed = TRUE,
                             chunk.size = chunk_size)
      
    ) %:%
    
    foreach(
      i = 1:R,
      .combine = "rbind",
      .multicombine = TRUE,
      .maxcombine = R,
      .options.future = list(seed = TRUE)
      
    ) %dofuture% {
      
      p()
      
      data <- preallocation$data
      
      psi_grid_list <- get_psi_grid(data, step_size, num_std_errors, split = TRUE)
      
      psi_MLE <- entropy(data / sum(data))
      
      u <- preallocation$preallocations[[i]]["u",]
      
      omega_hat <- preallocation$preallocations[[i]]["omega_hat",]
      
      log_like_u <- log_likelihood(u, data)
      
      omega_hat |>
        get_integrated_log_likelihood_vals.aux(data, psi_grid_list) |>
        (`-`)(log_like_u)
    } |>
    map(\(x) x |>
          matrixStats::colLogSumExps() |>
          (`-`)(log(R)))
}

################################################################################
############################## PROFILE LIKELIHOOD ############################## 
################################################################################

get_profile_log_likelihood <- function(data, psi_grid_list) {
  
  theta_MLE <- data / sum(data)
  
  l_p <- psi_grid_list |> 
    purrr::map(
      \(x) purrr::accumulate(
        x,
        \(acc, nxt) get_theta_hat(acc, nxt, theta_MLE), 
        .init = theta_MLE
      ) |> 
        magrittr::extract(-1) |> 
        purrr::map_dbl(likelihood, data)
    ) |> 
    purrr::modify_in(1, rev) |> 
    unlist() |> 
    log()
  
  return(l_p)
}


