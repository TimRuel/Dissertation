p <- Vectorize(dgamma, vectorize.args = "x")

q <- Vectorize(dchisq, vectorize.args = "x")

# q <- Vectorize(dlnorm, vectorize.args = "x")

# Prior
alpha_0 <- 0.5

beta_0 <- 0.5

# Likelihood
x <- 7

# Nominal
alpha_1 <- alpha_0 + x

beta_1 <- beta_0 + 1

# Importance
beta_2 <- 1.5
  
alpha_2 <- beta_2 / beta_1 * (alpha_1 - 1) + 1

xlim <- c(0, 10)

ggplot() +
  stat_function(fun = p, 
                args = list(shape = alpha_1, rate = beta_1),
                xlim = xlim,
                aes(color = "nominal")) +
  stat_function(fun = p, 
                args = list(shape = alpha_2, rate = beta_2),
                xlim = xlim,
                aes(color = "importance")) +
  scale_color_discrete(name = "Density") +
  theme_minimal() + 
  theme(legend.position = c(0.8, 0.8),
        legend.background = element_rect())


alpha_prior <- 2

beta_prior <- 2

alpha_posterior <- alpha_prior + data

beta_posterior <- beta_prior + rep(1, length(data))

rng <- rgamma

density <- dgamma

nominal_rng_params <- list(n = n, shape = alpha_posterior, rate = beta_posterior)

nominal_density_params <- list(shape = alpha_posterior, rate = beta_posterior)

importance_rng_params <- list(n = n, shape = 3, rate = 4)

importance_density_params <- list(shape = 3, rate = 4)

MC_params <- list(method = "importance", 
                  importance = list(rng = rng, 
                                    rng_params = importance_rng_params,
                                    density = density, 
                                    density_params = importance_density_params), 
                  nominal = list(rng = rng, 
                                 rng_params = nominal_rng_params,
                                 density = density, 
                                 density_params = nominal_density_params))

rng |> 
  do.call(args = nominal_rng_params)

u_list |> 
  map_dbl(\(u) sum(log(do.call(density, args = c(x = list(u), nominal_density_params)))))

get_importance_weights(u_list, MC_params)
  

MC_params <- list(method = "vanilla", 
                  dist_list = list(rng = rng, dist_params = nominal_dist_params))

dist_sampler(MC_params, R)

step_size <- 0.1

num_std_errors <- 1

psi_grid <- get_psi_grid(data, weights, step_size, num_std_errors, split = FALSE)

psi_hat <- dot_product(data, weights)

R <- 10

u_list <- dist_sampler(MC_params, R)

omega_hat_list <- get_omega_hat_list(u_list, psi_hat, weights)

lambda_method <- "accumulate"

chunk_size <- 1

plan(multisession, workers = I(10))

l_tilde_mat <- get_l_tilde_mat(psi_grid, omega_hat_list, chunk_size, lambda_method)

integrated_log_likelihood$importance_weights |> 
  data.frame() |> 
  ggplot() + 
  geom_histogram(aes(x = integrated_log_likelihood.importance_weights))

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
  pivot_longer(cols = c("nominal", "importance"),
               values_to = "density")

gamma_dens |> 
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


get_importance_weights <- function(u_list, MC_params) {
  
  MC_params$method |> 
    
    switch(
      
      vanilla = 0,
      
      importance = {
        
        u_list |> 
          map(\(u) {
            
            p <- MC_params$nominal$density |> 
              do.call(args = c(x = list(u), MC_params$nominal$density_params)) 
            
            q <- MC_params$importance$density |> 
              do.call(args = c(x = list(u), MC_params$importance$density_params)) 
            
            p / q
          }
          )
      }
    )
}


R <- 250

u_list <- dist_sampler(MC_params, R)

w <- get_importance_weights(u_list, MC_params)

for (i in 1:18) {
  
  w |> 
    map_dbl(\(x) x[i]) |> 
    hist(main = i)
}

for (i in 1:18) {
  
  w |> 
    map_dbl(\(x) x[i]) |> 
    hist(main = i)
}

mat <- matrix(unlist(w), ncol = 18, byrow = TRUE)

mat |> 
  sweep(2, colSums(mat), '/') 

plot(psi_grid, integrated_log_likelihood$l_bar)



basic_IS_diagnostics$variance_plot

self_norm_IS_diagnostics$variance_plot

regression_IS_diagnostics$variance_plot

methods = c("vanilla_MC", "basic_IS", "self_norm_IS", "regression_IS")

data.frame(psi = psi_grid,
           vanilla_MC = log(integrated_likelihood_vanilla_MC$L_bar$var_estimate),
           basic_IS = log(integrated_likelihood_basic_IS$L_bar$estimate),
           self_norm_IS = log(integrated_likelihood_self_norm_IS$L_bar$estimate),
           regression_IS = log(integrated_likelihood_regression_IS$L_bar$estimate)) |> 
  tidyr::pivot_longer(cols = all_of(methods),
                      names_to = "Method",
                      values_to = "Log_Variance") |>
  mutate(Method = Method |>
           as_factor() |>
           fct_inorder()) |> 
  ggplot() +
  geom_point(aes(x = psi, y = Log_Variance, color = Method),
             size = 0.1) +
  ylab("Log Variance") +
  scale_color_brewer(palette = "Set1") +
  xlab(expression(psi)) +
  theme_minimal() +
  theme(axis.line = element_line())

basic_IS_diagnostics

self_norm_IS_diagnostics

regression_IS_diagnostics




get_L_bar < function(psi, data, weights, MC_params, R) {
  
  psi_hat <- dot_product(data, weights)
  
  u_list <- dist_sampler(MC_params, R)
  
  omega_hat_list <- get_omega_hat_list(u_list, psi_hat, weights)
  
  omega_hat_list |> 
    map_dblget_L_tilde(psi, omega_hat, weights, 0, "map")
}

crit <- qchisq(0.95, 1) / 2






