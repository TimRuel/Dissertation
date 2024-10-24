################################################################################
#################################### GENERAL ###################################
################################################################################
library(tidyverse)

log_likelihood <- function(theta, data) {
  
  data |> 
    dpois(theta, log = TRUE) |> 
    sum(na.rm = TRUE)
}

likelihood <- function(theta, data) {
  
  data |> 
    dpois(theta) |> 
    prod()
}

neg_log_likelihood <- function(theta, data) -log_likelihood(theta, data)

dot_product <- function(x, y) sum(x * y, na.rm = TRUE)

distance <- function(a, b) sum((a - b)^2)

get_psi_hat_se <- function(data, weights) sqrt(sum(data*(weights^2)))

get_psi_grid <- function(data, weights, step_size, num_std_errors, split = FALSE) {
  
  psi_hat <- dot_product(data, weights)
  
  psi_hat_se <- get_psi_hat_se(data, weights)
  
  MoE <- num_std_errors * psi_hat_se
  
  psi_grid <- (psi_hat + MoE*c(-1, 1)) |> 
    (\(x) c(max(0, x[1]), x[2]))() |> 
    plyr::round_any(step_size, floor) |> 
    (\(x) seq(x[1], x[2], step_size))()
  
  if (split) {
    
    psi_grid_list <- psi_grid |> 
      split(factor(psi_grid > psi_hat)) |> 
      purrr::modify_in(1, rev) |> 
      unname()
    
    return(psi_grid_list)
  }
  
  return(psi_grid)
}

get_dist_list <- function(MC_params, dist_type = "importance") {
  
  switch(MC_params$method,
         vanilla_MC = MC_params$nominal,
         self_norm_IS = ,
         regression_IS = ,
         basic_IS = switch(dist_type, 
                           importance = MC_params$importance,
                           nominal = MC_params$nominal,
                           stop("Invalid distribution type")),
         stop("Invalid Monte Carlo method")
  )
}

dist_sampler <- function(MC_params, R, simplify = FALSE) {
  
  dist_list <- get_dist_list(MC_params)
  
  rng <- dist_list$rng
  
  rng_params <- dist_list$rng_params
  
  rng |> 
    do.call(args = rng_params) |> 
    replicate(R, expr = _, simplify = simplify)
}

get_importance_weights <- function(u_list, MC_params) {
  
  MC_params$method |> 
    
    switch(
      
      vanilla_MC = 1,
      
      self_norm_IS = ,
      
      regression_IS = ,
      
      basic_IS = {
        
        w <- u_list |>
          map_dbl(\(u) {

            p <- MC_params$nominal$density |>
              do.call(args = c(x = list(u), MC_params$nominal$density_params)) |>
              prod()

            q <- MC_params$importance$density |>
              do.call(args = c(x = list(u), MC_params$importance$density_params)) |>
              prod()

            p / q
          }
          )

        if (MC_params$method == "self_norm_IS") w <- w / mean(w)

        w
      }
    )
}

Q <- function(u, psi_hat, weights) u / dot_product(u, weights) * psi_hat

get_omega_hat_list <- function(u_list, psi_hat, weights) purrr::map(u_list, \(u) Q(u, psi_hat, weights))

get_theta_hat <- function(lambda, omega_hat, weights) omega_hat / (1 + lambda * weights)

get_lambda <- function(init_lambda, psi, omega_hat, weights) {
  
  psi_dist <- function(lambda) {
    
    lambda |> 
      get_theta_hat(omega_hat, weights) |> 
      dot_product(weights) |> 
      distance(psi)
  }
  
  psi_dist.gr <- function(lambda) nloptr::nl.grad(lambda, psi_dist)
  
  nloptr::auglag(x0 = init_lambda,
                 fn = psi_dist,
                 gr = psi_dist.gr,
                 lower = -min(1 / weights),
                 localsolver = "LBFGS")$par
}

get_lambdas_accumulate <- function(psi_grid, omega_hat, weights, init_lambda) {
  
  psi_grid |>
    purrr::accumulate(
      \(acc, nxt) get_lambda(acc, nxt, omega_hat, weights),
      .init = init_lambda) |>
    magrittr::extract(-1)
}

get_lambdas_map <- function(psi_grid, omega_hat, weights, init_lambda) {
  
  psi_grid |>
    purrr::map_dbl(\(psi) get_lambda(init_lambda, psi, omega_hat, weights))
}

get_L_tilde <- function(psi_grid, omega_hat, weights, init_lambda, lambda_method = "accumulate") {
  
  lambda_method |> 
    
    switch(
      
      accumulate = get_lambdas_accumulate(psi_grid, omega_hat, weights, init_lambda),
      
      map = get_lambdas_map(psi_grid, omega_hat, weights, init_lambda)
      
    ) |> 
    
    purrr::map_dbl(\(lambda) {
      lambda |>
        get_theta_hat(omega_hat, weights) |> 
        likelihood(data)
    }
    )
}

################################################################################
############################ INTEGRATED LIKELIHOOD ############################# 
################################################################################

get_L_tilde_mat <- function(psi_grid, omega_hat_list, chunk_size, lambda_method, init_lambda = 0) {
  
  p <- progressr::progressor(along = omega_hat_list)
  
  foreach(
    
    omega_hat = omega_hat_list,
    .combine = "rbind",
    .multicombine = TRUE,
    .maxcombine = length(omega_hat_list),
    .options.future = list(seed = TRUE,
                           chunk.size = chunk_size)
    
  ) %dofuture% {
    
    p()
    
    get_L_tilde(psi_grid, omega_hat, weights, init_lambda, lambda_method)
  }
}

get_L_bar <- function(L_tilde_mat, w, MC_params) {
  
  MC_params$method |> 
    
    switch(
      
      vanilla_MC = ,
      
      self_norm_IS = ,
      
      basic_IS = {
        
        weighted_vals <- L_tilde_mat |> 
          sweep(1, w, '*')
      },
      
      regression_IS = {
        
        w_bar <- mean(w)
        
        w_centered <- w - w_bar
        
        denom <- sum(w_centered^2)
        
        num <- L_tilde_mat |> 
          sweep(1, w*w_centered, '*') |> 
          colSums()
        
        beta_hat <- num / denom
        
        weighted_vals <- sweep(L_tilde_mat, 1, w, '*') - outer(w - 1, beta_hat, '*')
      }
    )
  
  estimate <- colMeans(weighted_vals)
  
  var_estimate <- matrixStats::colVars(weighted_vals)
  
  return(list(weighted_vals = weighted_vals,
              estimate = estimate, 
              var_estimate = var_estimate))
}

get_integrated_likelihood <- function(data, 
                                      weights, 
                                      psi_grid, 
                                      R,
                                      MC_params,
                                      lambda_method,
                                      init_lambda,
                                      chunk_size) {
  
  psi_hat <- dot_product(data, weights)
  
  u_list <- dist_sampler(MC_params, R)
  
  omega_hat_list <- get_omega_hat_list(u_list, psi_hat, weights)
  
  L_tilde_mat <- get_L_tilde_mat(psi_grid, omega_hat_list, chunk_size, lambda_method, init_lambda)
  
  importance_weights <- get_importance_weights(u_list, MC_params)
  
  L_bar <- get_L_bar(L_tilde_mat, importance_weights, MC_params)
  
  l_bar <- log(L_bar$estimate)

  return(list(L_bar = L_bar,
              l_bar = l_bar,
              L_tilde_mat = L_tilde_mat,
              importance_weights = importance_weights,
              omega_hat_list = omega_hat_list,
              u_list = u_list, 
              psi_hat = psi_hat, 
              MC_params = MC_params,
              data = data,
              psi_grid = psi_grid))
}

get_importance_sampling_diagnostics <- function(integrated_log_likelihood) {
  
  MC_params <- integrated_log_likelihood$MC_params
  
  data <- integrated_log_likelihood$data
  
  weights <- integrated_log_likelihood$importance_weights
  
  L_tilde_mat <- integrated_log_likelihood$L_tilde_mat
  
  w_bar <- mean(weights)
  
  CV_global <- sd(weights) / w_bar
  
  ESS_global <- sum(weights)^2 / sum(weights^2)
  
  num <- L_tilde_mat |> 
    sweep(1, weights, '*') 
  
  w_tilde <- num |> 
    sweep(2, colSums(num), FUN = '/')
  
  CV_by_psi <- sqrt(length(weights) * colSums(w_tilde^2) - 1)
  
  ESS_by_psi <- 1 / colSums(w_tilde^2)
  
  density_args <- data.frame(count = factor(data),
                             nominal_shape = MC_params$nominal$density_params$shape,
                             nominal_rate = MC_params$nominal$density_params$rate,
                             importance_shape = MC_params$importance$density_params$shape,
                             importance_rate = MC_params$importance$density_params$rate) |> 
    distinct() |> 
    arrange(count)
  
  grid <- seq(0, 10, 0.1)
  
  gamma_dens <- ddply(density_args, "count", function(df) {
    
    data.frame(
      
      x = grid,
      nominal = dgamma(grid, mean(df$nominal_shape), mean(df$nominal_rate)),
      importance = dgamma(grid, mean(df$importance_shape), mean(df$importance_rate))
    )
  }) |> 
    tidyr::pivot_longer(cols = c("nominal", "importance"),
                        values_to = "density")
  
  density_plots <- gamma_dens |> 
    ggplot() +
    geom_line(aes(x = x,
                  y = density,
                  color = name)) +
    facet_wrap(~ count) + 
    scale_color_discrete(name = "Density Type") +
    theme_minimal() + 
    theme(legend.position = c(0.65, 0.1),
          legend.direction = "horizontal") + 
    guides(col = guide_legend(title.position = "top",
                              title.hjust = 0.5))
  
  importance_weights_hist <- weights |> 
    data.frame() |> 
    ggplot() + 
    geom_histogram(aes(x = weights))
  
  variance_plot <- data.frame(psi = integrated_log_likelihood$psi_grid,
                              variance = integrated_log_likelihood$L_bar$var_estimate) |>  
    ggplot() +
    geom_point(aes(x = psi, y = variance),
               size = 0.1) +
    ylab("Variance") +
    xlab(expression(psi)) +
    theme_minimal() +
    theme(axis.line = element_line())
  
  return(list(w_bar = w_bar,
              CV_global = CV_global,
              ESS_global = ESS_global,
              CV_by_psi = CV_by_psi,
              ESS_by_psi = ESS_by_psi,
              density_plots = density_plots,
              importance_weights_hist = importance_weights_hist,
              variance_plot = variance_plot))
}

get_profile_log_likelihood <- function(data, weights, step_size, num_std_errors) {
  
  psi_grid_list <- get_psi_grid(data, weights, step_size, num_std_errors, split = TRUE)
  
  theta_mle <- data
  
  p <- progressr::progressor(along = psi_grid_list)
  
  foreach(
    
    psi_grid = psi_grid_list,
    .options.future = list(seed = TRUE)
    
  ) %dofuture% {
    
    p()
    
    psi_grid |> 
      purrr::accumulate(
        \(acc, nxt) get_lambda(acc, nxt, theta_mle, weights),
        .init = 0) |> 
      magrittr::extract(-1) |> 
      purrr::map_dbl(\(lambda) {
        lambda |>
          get_theta_hat(theta_mle, weights) |> 
          log_likelihood(data)
      })
  } |> 
    purrr::modify_in(1, rev) |>
    unlist(recursive = FALSE)
}
