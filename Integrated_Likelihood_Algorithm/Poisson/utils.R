################################################################################
#################################### GENERAL ###################################
################################################################################

log_likelihood <- function(theta, data) {
  
  data |> 
    map_dbl(sum) |> 
    (`*`)(log(theta)) |> 
    (`-`)(map_dbl(data, length)*theta) |> 
    sum(na.rm = TRUE)
}

neg_log_likelihood <- function(theta, data) -log_likelihood(theta, data)

likelihood <- function(theta, data) exp(log_likelihood(theta, data))

weighted_sum <- function(theta, weights) sum(theta * weights, na.rm = TRUE)

distance <- function(a, b) sum((a - b)^2)

get_psi_grid <- function(data, weights, step_size, num_std_errors, split = FALSE) {
  
  psi_MLE <- data |> 
    map_dbl(mean) |> 
    weighted_sum(weights)
  
  psi_MLE_SE <- data |> 
    map_dbl(\(x) sum(x) / length(x)^2) |> 
    (`*`)(weights^2) |> 
    sum() |> 
    sqrt()
  
  MoE <- num_std_errors * psi_MLE_SE
  
  psi_grid <- (psi_MLE + MoE*c(-1, 1)) |> 
    (\(x) c(max(0, x[1]), x[2]))() |> 
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

# get_omega_hat <- function(objective_fn, data, weights, u_params, tol, return_u = FALSE) {
#   
#   psi_MLE <- data |> 
#     map_dbl(mean) |> 
#     weighted_sum(weights)
#   
#   alpha <- u_params$alpha
#   
#   beta <- u_params$beta
#   
#   u <- rgamma(length(data), alpha, beta)
#   
#   f <- function(omega) objective_fn(omega, u)
#   f.gr <- function(omega) nloptr::nl.grad(omega, f)
#   fcon <- function(omega) weighted_sum(omega, weights) - psi_MLE
#   fcon.jac <- function(omega) nloptr::nl.jacobian(omega, fcon)
#   
#   init_guess <- data |> 
#     map_dbl(mean) |>
#     (`+`)(rnorm(length(data), sd = psi_MLE/weights)) |> 
#     abs()
#   
#   omega_hat <- nloptr::auglag(x0 = init_guess,
#                               fn = f,
#                               gr = f.gr,
#                               heq = fcon,
#                               heqjac = fcon.jac,
#                               lower = rep(1e-6, length(data)),
#                               localsolver = "LBFGS")$par
#   
#   if (abs(weighted_sum(omega_hat, weights) - psi_MLE) < tol) {
#     
#     if (return_u) {
#       
#       return(rbind(u, omega_hat))
#     }
#     
#     return(omega_hat)
#   }
#   
#   else {
#     get_omega_hat(objective_fn, data, weights, u_params, tol, return_u)
#   }
# }

# get_omega_hat <- function(data, weights) {
#   
#     psi_MLE <- data |>
#       map_dbl(mean) |>
#       weighted_sum(weights)
#     
#     remainder_index <- sample(length(data), 1)
#     
#     omega_hat_minus_one <- runif(length(data) - 1, 0, psi_MLE / weights[-remainder_index])
#     
#     remainder <- (psi_MLE - weighted_sum(omega_hat_minus_one, weights[-remainder_index])) / weights[remainder_index]
#     
#     omega_hat <- append(omega_hat_minus_one, remainder, after = remainder_index - 1)
#     
#     return(omega_hat)
# }

get_omega_hat <- function(data, weights, dist, dist_params, return_u = FALSE) {
  
  psi_MLE <- data |>
    map_dbl(mean) |>
    weighted_sum(weights)
  
  u <- data |> 
    length() |> 
    c(dist_params) |> 
    do.call(dist, args = _)
  
  omega_hat <- u |> 
    (`/`)(sum(u * weights)) |> 
    (`*`)(psi_MLE)
  
  if (return_u) return(rbind(u, omega_hat))
  
  return(omega_hat)
}

get_theta_hat <- function(lambda, m, omega_hat, weights) (m * omega_hat) / (m + lambda * weights)

get_lambda <- function(init_guess, psi, m, omega_hat, weights) {
  
  f <- function(lambda) {
    
    lambda |> 
      get_theta_hat(m, omega_hat, weights) |> 
      weighted_sum(weights) |> 
      distance(psi)
  }
  
  f.gr <- function(lambda) nloptr::nl.grad(lambda, f)
  
  lambda <- nloptr::auglag(x0 = init_guess,
                           fn = f,
                           gr = f.gr,
                           lower = -min(m / weights),
                           # upper = min(m / weights),
                           localsolver = "LBFGS")$par
  
  return(lambda)
}

# get_lambda <- function(init_guess, psi, m, omega_hat, weights) {
#   
#   objective <- function(lambda) {
#     
#     lambda |> 
#       get_theta_hat(m, omega_hat, weights) |> 
#       weighted_sum(weights) |> 
#       distance(psi)
#   }
#   
#   out <- nlm(f = objective, 
#              p = init_guess,
#              fscale = 0,
#              iterlim = 10000)
#   
#   lambda <- out$estimate
#   
#   if (lambda <= -min(m / weights)) get_lambda(lambda + runif(1, 0, 2), psi, m, omega_hat, weights)
#   
#   return(lambda)
# }

# get_theta_hat <- function(init_guess, psi, omega_hat, weights) {
#   
#   f <- function(theta) neg_log_likelihood(theta, omega_hat)
#   f.gr <- function(theta) nloptr::nl.grad(theta, f)
#   fcon <- function(theta) weighted_sum(theta, weights) - psi
#   fcon.jac <- function(theta) nloptr::nl.jacobian(theta, fcon)
#   
#   theta_hat <- nloptr::auglag(x0 = init_guess,
#                               fn = f,
#                               gr = f.gr,
#                               heq = fcon,
#                               heqjac = fcon.jac,
#                               lower = rep(1e-10, length(omega_hat)),
#                               localsolver = "LBFGS")$par
#   
#   return(theta_hat)
# }

# get_theta_hat_list <- function(init_guess, data, weights, psi_grid_list) {
# 
#   psi_grid_list |>
#     purrr::map(
#       \(x) {
#         p <- progressr::progressor(along = x)
#         purrr::accumulate(
#           x,
#           \(acc, nxt) {
#             p()
#             get_theta_hat(acc, nxt, data, weights)
#             },
#           .init = init_guess
#           )
#         } |>
#         magrittr::extract(-1)
#       ) |>
#     purrr::modify_in(1, rev) |>
#     unlist(recursive = FALSE)
#   }

get_theta_hat_list <- function(init_guess, data, weights, psi_grid_list) {
  
  p <- progressr::progressor(steps = 2*length(psi_grid_list[[1]]))
  
  foreach(
    
    psi_grid = psi_grid_list,
    .options.future = list(seed = TRUE)
    
  ) %dofuture% {
    
    psi_grid |> purrr::accumulate(
      \(acc, nxt) {
        p()
        get_theta_hat(acc, nxt, data, weights)
      },
      .init = init_guess
    ) |> 
      magrittr::extract(-1) 
  } |> 
    purrr::modify_in(1, rev) |>
    unlist(recursive = FALSE)
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

get_IL_preallocations <- function(data_sims, weights, u_params, R, tol, chunk_size) {
  
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
      
      omega_hat <- get_omega_hat(neg_log_likelihood, data, weights, u_params, tol, return_u = FALSE)
      
      return(list(data = data, preallocations = omega_hat, i = i))
    }
  
  return(IL_preallocations)
}

# get_integrated_log_likelihood_vals.aux <- function(omega_hat, data, weights, psi_grid_list) {
# 
#   psi_grid_list |>
#     purrr::map(
#       \(psi_grid) psi_grid |>
#         purrr::accumulate(
#           \(acc, nxt) get_theta_hat(acc, nxt, omega_hat, weights),
#           .init = omega_hat
#         ) |>
#         magrittr::extract(-1) |>
#         purrr::map_dbl(log_likelihood, data)
#     ) |>
#     purrr::modify_in(1, rev) |>
#     unlist()
#   }
# 
# get_integrated_log_likelihood_vals.aux <- function(omega_hat, data, weights, psi_grid_list) {
#   
#   init_guess <- data |> 
#     map_dbl(mean) |> 
#     (`+`)(0.5)
#   
#   out <- get_theta_hat_list(init_guess, omega_hat, weights, psi_grid_list) |> 
#     sapply(log_likelihood, data)
#   
#   return(out)
# }
# 
# get_integrated_log_likelihood_vals.aux <- function(omega_hat, data, weights, profile_theta_hat_matrix, num_workers) {
#   
#   psi_grid <- profile_theta_hat_matrix |> 
#     rownames() |> 
#     as.double()
#   
#   chunk_size <- ceiling(length(psi_grid) / num_workers)
#   
#   p <- progressr::progressor(along = psi_grid)
#   
#   foreach(
#     
#     psi = psi_grid,
#     .combine = "c",
#     .multicombine = TRUE,
#     .maxcombine = length(psi_grid),
#     .options.future = list(seed = TRUE,
#                            chunk.size = chunk_size)
#     
#   ) %dofuture% {
#     
#     p()
#     
#     profile_theta_hat_matrix[as.character(psi),] |> 
#       get_theta_hat(psi, omega_hat, weights) |> 
#       log_likelihood(data)
#   }
# }

custom_combine <- function(...) {
  
  vals <- list(...) |>    
    purrr::modify_in(1, rev) |>
    unlist()
  
  return(vals)
}

# Left off here 9/23. Need to rework get_integrated_log_likelihood_vals to accomodate updated form of get_theta_hat. 
# Replace accumulate function with purrr::map as initial guesses are no longer required. 
# Compare new way to old.
# Other things to explore: using nlm to search for CI endpoints more efficiently, replacing psi_grid approach


get_integrated_log_likelihood_vals <- function(data, weights, step_size, num_std_errors, dist, dist_params, R = 250, chunk_size) {
  
  psi_grid <- get_psi_grid(data, weights, step_size, num_std_errors, split = FALSE)
  
  p <- progressr::progressor(steps = R)
  
  m <- data |> 
    map_dbl(length)
  
  l_diff <- foreach(
    
    i = 1:R,
    .combine = "rbind",
    .multicombine = TRUE,
    .maxcombine = R,
    .options.future = list(seed = TRUE,
                           chunk.size = chunk_size)
    
  ) %dofuture% {
    
    p()
    
    omega_hat <- get_omega_hat(data, weights, dist, dist_params)
    
    init_guess <- 0
    
    psi_grid |>
      purrr::accumulate(
        \(acc, nxt) get_lambda(acc, nxt, m, omega_hat, weights),
        .init = init_guess) |>
      magrittr::extract(-1) |>
      purrr::map_dbl(\(lambda) {
        lambda |>
          get_theta_hat(m, omega_hat, weights) |> 
          log_likelihood(data)
      })
  }
  
  log_I_hat <- l_diff |>
    matrixStats::colLogSumExps() |>
    (`-`)(log(R))
  
  I_hat <- exp(log_I_hat)
  
  L_ratio <- exp(l_diff)
  
  s_squared <- L_ratio |> 
    apply(2, \(x) var(x) * (R-1) / R)
  
  s <- sqrt(s_squared)
  
  cv <- s / I_hat / sqrt(R)
  
  w_star_squared_sum <- L_ratio |> 
    sweep(2, colSums(L_ratio), FUN = '/') |>
    apply(2, \(x) x^2) |> 
    colSums()
  
  return(list(vals = log_I_hat, 
              L_ratio = L_ratio,
              s_squared = s_squared,
              cv = cv,
              w_star_squared_sum = w_star_squared_sum))
  }

# get_integrated_log_likelihood_vals <- function(data, weights, step_size, num_std_errors, dist, dist_params, R = 250, chunk_size) {
# 
#   psi_grid_list <- get_psi_grid(data, weights, step_size, num_std_errors, split = TRUE)
#   
#   p <- progressr::progressor(steps = 2 * length(psi_grid_list[[1]]) * R)
#   
#   init_guess <- data |> 
#     map_dbl(mean) |> 
#     (`+`)(0.5)
#   
#   omega_hat <- R |> 
#     replicate(get_omega_hat(data, weights, dist, dist_params), simplify = FALSE) 
#   
#   foreach(
#     
#     i = 1:R,
#     .combine = "rbind",
#     .multicombine = TRUE,
#     .maxcombine = R,
#     .options.future = list(seed = TRUE,
#                            chunk.size = chunk_size)
#     
#   ) %:%
#     
#     foreach(
#       
#       psi_grid = psi_grid_list,
#       .combine = custom_combine,
#       .options.future = list(seed = TRUE)
#       
#     )  %dofuture% {
#       
#       psi_grid |> 
#         purrr::accumulate(
#           \(acc, nxt) {
#             p()
#             get_theta_hat(acc, nxt, omega_hat[[i]], weights)
#             },
#           .init = init_guess
#           ) |> 
#         magrittr::extract(-1) |> 
#         purrr::map_dbl(log_likelihood, data)
#       } |>
#     matrixStats::colLogSumExps() |>
#     (`-`)(log(R))
#   }

get_integrated_log_likelihood_sims <- function(preallocations, weights, step_size, num_std_errors, chunk_size) {
  
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
      
      psi_grid_list <- get_psi_grid(data, weights, step_size, num_std_errors, split = TRUE)
      
      omega_hat <- preallocation$preallocations[i,]
      
      get_integrated_log_likelihood_vals.aux(omega_hat, data, weights, psi_grid_list)
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

get_mod_IL_preallocations <- function(data_sims, weights, Q, prior, R, tol, chunk_size) {
  
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
      
      u_params <- list(alpha = data + prior$alpha,
                       beta = 1 + prior$beta)
      
      preallocations <- get_omega_hat(Q, data, u_params, tol, return_u = TRUE)
      
      return(list(data = data, preallocations = preallocations, i = i))
    }
  
  return(mod_IL_preallocations)
}

# get_mod_integrated_log_likelihood_vals <- function(data, weights, Q, step_size, num_std_errors, u_params, R = 250, tol = 0.0001, chunk_size) {
#   
#   p <- progressr::progressor(steps = R)
#   
#   psi_grid_list <- get_psi_grid(data, weights, step_size, num_std_errors, split = TRUE)
#   
#   foreach(
#     
#     i = 1:R,
#     .combine = "rbind",
#     .multicombine = TRUE,
#     .maxcombine = R,
#     .options.future = list(seed = TRUE,
#                            chunk.size = chunk_size)
#     
#   ) %dofuture% {
#     
#     p()
#     
#     alpha <- u_params$alpha
#     
#     beta <- u_params$beta
#     
#     mat <- get_omega_hat(Q, data, weights, u_params, tol, return_u = TRUE)
#     
#     u <- mat["u",]
#     
#     omega_hat <- mat["omega_hat",]
#     
#     log_like_u <- log_likelihood(u, data)
#     
#     omega_hat |> 
#       get_integrated_log_likelihood_vals.aux(data, weights, psi_grid_list) |> 
#       (`-`)(log_like_u)
#   } |> 
#     matrixStats::colLogSumExps() |>
#     (`-`)(log(R))
# }

get_mod_integrated_log_likelihood_vals <- function(data, weights, step_size, num_std_errors, dist, dist_params, R = 250, chunk_size) {
  
  psi_grid <- get_psi_grid(data, weights, step_size, num_std_errors, split = FALSE)
  
  p <- progressr::progressor(steps = R)
  
  m <- data |> 
    map_dbl(length)
  
  l_diff <- foreach(
    
    i = 1:R,
    .combine = "rbind",
    .multicombine = TRUE,
    .maxcombine = R,
    .options.future = list(seed = TRUE,
                           chunk.size = chunk_size)
    
  ) %dofuture% {
    
    p()
    
    omega_hat_mat <- get_omega_hat(data, weights, dist, dist_params, return_u = TRUE)
    u <- omega_hat_mat["u",]
    omega_hat <- omega_hat_mat["omega_hat",]
    log_like_u <- log_likelihood(u, data)
    
    init_guess <- 0
    
    psi_grid |>
      purrr::accumulate(
        \(acc, nxt) get_lambda(acc, nxt, m, omega_hat, weights),
        .init = init_guess) |>
      magrittr::extract(-1) |>
      purrr::map_dbl(\(lambda) {
        lambda |>
          get_theta_hat(m, omega_hat, weights) |> 
          log_likelihood(data) |> 
          (`-`)(log_like_u)
      })
  } 
  
  log_I_hat <- l_diff |>
    matrixStats::colLogSumExps() |>
    (`-`)(log(R))
  
  I_hat <- exp(log_I_hat)
  
  L_ratio <- exp(l_diff)
  
  s_squared <- L_ratio |> 
    apply(2, \(x) var(x) * (R-1) / R)
  
  s <- sqrt(s_squared)
  
  cv <- s / I_hat / sqrt(R)
  
  w_star_squared_sum <- L_ratio |> 
    sweep(2, colSums(L_ratio), FUN = '/') |>
    apply(2, \(x) x^2) |> 
    colSums()
  
  return(list(vals = log_I_hat, 
              L_ratio = L_ratio,
              s_squared = s_squared,
              cv = cv,
              w_star_squared_sum = w_star_squared_sum))
}

# get_mod_integrated_log_likelihood_vals <- function(data, weights, step_size, num_std_errors, dist, dist_params, R = 250, chunk_size) {
#   
#   p <- progressr::progressor(steps = R)
#   
#   psi_grid_list <- get_psi_grid(data, weights, step_size, num_std_errors, split = TRUE)
#   
#   foreach(
#     
#     i = 1:R,
#     .combine = "rbind",
#     .multicombine = TRUE,
#     .maxcombine = R,
#     .options.future = list(seed = TRUE,
#                            chunk.size = chunk_size)
#     
#   ) %dofuture% {
#     
#     p()
#     
#     mat <- get_omega_hat(data, weights, dist, dist_params, return_u = TRUE)
#     
#     u <- mat["u",]
#     
#     omega_hat <- mat["omega_hat",]
#     
#     log_like_u <- log_likelihood(u, data)
#     
#     omega_hat |> 
#       get_integrated_log_likelihood_vals.aux(data, weights, psi_grid_list) |> 
#       (`-`)(log_like_u)
#   } |> 
#     matrixStats::colLogSumExps() |>
#     (`-`)(log(R))
# }

get_mod_integrated_log_likelihood_sims <- function(preallocations, weights, step_size, num_std_errors, chunk_size) {
  
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
      
      psi_grid_list <- get_psi_grid(data, weights, step_size, num_std_errors, split = TRUE)
      
      u <- preallocation$preallocations[[i]]["u",]
      
      omega_hat <- preallocation$preallocations[[i]]["omega_hat",]
      
      log_like_u <- log_likelihood(u, data)
      
      omega_hat |>
        get_integrated_log_likelihood_vals.aux(data, weights, psi_grid_list) |>
        (`-`)(log_like_u)
    } |>
    map(\(x) x |>
          matrixStats::colLogSumExps() |>
          (`-`)(log(R)))
}

################################################################################
############################## PROFILE LIKELIHOOD ############################## 
################################################################################

# get_profile_log_likelihood <- function(data, weights, psi_grid_list) {
#   
#   l_p <- psi_grid_list |> 
#     purrr::map(
#       \(x) purrr::accumulate(
#         x,
#         \(acc, nxt) get_theta_hat(acc, nxt, data, weights), 
#         .init = map_dbl(data, mean)
#       ) |> 
#         magrittr::extract(-1) |> 
#         sapply(log_likelihood, data)
#     ) |>
#     purrr::modify_in(1, rev) |>
#     unlist()
#   
#   return(l_p)
# }

get_profile_log_likelihood <- function(data, weights, step_size, num_std_errors) {
  
  psi_grid_list <- get_psi_grid(data, weights, step_size, num_std_errors, split = TRUE)
  
  theta_mle <- data |> 
    map_dbl(mean)
  
  m <- data |> 
    map_dbl(length)
  
  p <- progressr::progressor(along = psi_grid_list)
  
  foreach(
    
    psi_grid = psi_grid_list,
    .options.future = list(seed = TRUE)
    
  ) %dofuture% {
    
    p()
    
    psi_grid |> 
      purrr::accumulate(
        \(acc, nxt) get_lambda(acc, nxt, m, theta_mle, weights),
        .init = 0) |> 
      magrittr::extract(-1) |> 
      purrr::map_dbl(\(lambda) {
        lambda |>
          get_theta_hat(m, theta_mle, weights) |> 
          log_likelihood(data)
      })
  } |> 
    purrr::modify_in(1, rev) |>
    unlist(recursive = FALSE)
}
