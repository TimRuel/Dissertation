# Config -----------------------------------------------------------------

inject_shared_args <- function(cfg, shared_name = "shared_args") {
  
  shared <- cfg[[shared_name]]
  
  inject <- function(x) {
    if (is.list(x)) {
      if (!is.null(x$args)) {
        
        x$args <- modifyList(shared, x$args)
      }
      x <- lapply(x, inject)
    }
    x
  }
  
  cfg[names(cfg) != shared_name] <- lapply(cfg[names(cfg) != shared_name], inject)
  
  cfg$shared_args <- NULL
  
  return(cfg)
}

load_config <- function(population) {
  
  base_cfg <- config::get(file = "config/base.yml")
  
  pop_cfg <- population |> 
    sprintf("config/population_%s.yml", ... = _) |> 
    config::get(file = _) |> 
    inject_shared_args()
  
  cfg <- modifyList(base_cfg, pop_cfg)
  
  pop_result_base <- file.path(cfg$output_dir, sprintf("population_%s", population))
  dir.create(pop_result_base, recursive = TRUE, showWarnings = FALSE)
  
  existing_runs <- list.dirs(pop_result_base, full.names = FALSE, recursive = FALSE)
  existing_runs <- existing_runs[grepl("^run_\\d{3}$", existing_runs)]
  
  if (length(existing_runs) == 0) {
    
    next_run_num <- 1
  } else {
    
    run_nums <- as.integer(sub("run_", "", existing_runs))
    next_run_num <- max(run_nums) + 1
  }
  
  run_id <- sprintf("run_%03d", next_run_num)
  result_dir <- file.path(pop_result_base, run_id)
  
  cfg$run_id <- run_id
  cfg$result_dir <- result_dir
  
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  
  return(cfg)
}


# Data --------------------------------------------------------------------

get_theta_0 <- function(cfg, envir = environment()) {
  
  list2env(cfg$theta_0$dist, envir)
  
  do.call(name, args) |>
    (\(x) if (avg) x / mean(x) * avg else x)() |> 
    round(1)
}

get_Y <- function(theta_0) rpois(length(theta_0), theta_0)

get_weights <- function(cfg) {
  
  list2env(cfg$weights$dist, envir = environment())
  
  do.call(name, args) |> 
    (\(x) if (avg) x / mean(x) * avg else x)() 
}

# Pseudolikelihood ------------------------------------------------

log_likelihood <- function(theta, Y) sum(Y * log(theta) - theta, na.rm = TRUE)

likelihood <- function(theta, Y) exp(log_likelihood(theta, Y))

dot_product <- function(x, y) sum(x * y, na.rm = TRUE)

get_theta_hat <- function(lambda, omega_hat, weights) omega_hat / (1 - lambda * weights)

get_lambda_grid <- function(cfg, envir = environment()) {
  
  list2env(cfg$lambda_grid, envir)
  
  a <- interval$start
  
  b <- interval$stop
  
  lambda_neg <- 1 / seq(1/a, -1/epsilon, -step_size)
  
  lambda_pos <- 1 / seq(1/epsilon, 1/b, -step_size)
  
  lambda_grid <- c(lambda_neg, 0, lambda_pos)
  
  return(lambda_grid)
}

get_branch_df <- function(omega_hat, weights, lambda_grid) {
  
  theta_hat_list <- lambda_grid |> 
    purrr::map(\(lambda) get_theta_hat(lambda, omega_hat, weights))
  
  psi_grid <- theta_hat_list |> 
    purrr::map_dbl(\(theta_hat) sum(theta_hat * weights))
  
  log_L_tilde <- theta_hat_list |> 
    purrr::map_dbl(\(theta_hat) log_likelihood(theta_hat, Y))
  
  branch_df <- data.frame(psi = psi_grid, 
                          Integrated = log_L_tilde)
  
  return(branch_df)
}

make_plot <- function(df) {
  
  p <- withr:: with_package("ggplot2", {
    df |> 
      ggplot(aes(x = psi, y = .data[[names(df)[2]]])) + 
      geom_point(color = "cyan", size = 3, alpha = 0.7) + 
      theme_minimal(base_size = 15) +  # Minimal theme with a larger base font size
      theme(
        plot.background = element_rect(fill = "#2E2E2E", color = NA),  # Dark background for the whole plot
        panel.background = element_rect(fill = "#3A3A3F", color = "#1A1A1A", linewidth = 2),  # Lighter panel with a border
        axis.text = element_text(color = "white"),  # White axis labels
        axis.title = element_text(color = "white"),  # White axis titles
        plot.title = element_text(color = "white", size = 18, face = "bold"),  # White title
        plot.caption = element_text(color = "gray", size = 10),  # Gray caption
        panel.grid = element_line(color = "gray30", linetype = "dashed")  # Subtle grid lines
      ) + 
      labs(
        x = "\u03C8",
        y = paste(names(df)[[2]], "Log-Likelihood")
      )
  })
  
  return(p)
}

# Integrated Likelihood ---------------------------------------------------

get_psi_grid <- function(cfg, envir = environment()) {
  
  list2env(cfg$data, envir)
  list2env(cfg$weights, envir)
  list2env(cfg$psi_grid, envir)
  
  psi_hat <- dot_product(Y, weights)
  
  psi_hat_se <- sqrt(sum(Y*(weights^2)))
  
  MoE <- num_std_errors * psi_hat_se
  
  psi_grid <- (psi_hat + MoE) |>
    (\(x) c(max(0, x[1]), x[2]))() |>
    plyr::round_any(step_size, floor) |>
    (\(x) seq(x[1], x[2], step_size))()
  
  return(psi_grid)
}

sample_from_omega_psi_hat <- function(psi_hat, weights) {
  
  n <- length(weights)

  constraints <- hitandrun::simplexConstraints(n)
  constraints$constr[1,] <- weights
  constraints$rhs[[1]] <- psi_hat
  
  omega_hat <- hitandrun::hitandrun(constraints, 1) |> c()
  
  return(omega_hat)
}

generate_branches <- function(cfg, envir = environment()) {
  
  list2env(cfg$data, envir)
  list2env(cfg$weights, envir)
  list2env(cfg$parallel_specs, envir)
  
  num_branches <- num_workers * chunk_size
  
  lambda_grid <- get_lambda_grid(cfg)
  
  psi_hat <- dot_product(Y, weights)
  
  future::plan(future::multisession, workers = I(num_workers))
  
  `%dofuture%` <- doFuture::`%dofuture%`
  
  result <- foreach::foreach(
    
    i = 1:num_branches,
    .combine = "list",
    .multicombine = TRUE,
    .maxcombine = num_branches,
    .errorhandling = "remove",
    .options.future = list(seed = TRUE,
                           chunk.size = chunk_size)
    
  ) %dofuture% {
    
    omega_hat <- sample_from_omega_psi_hat(psi_hat, weights)
    
    branch_df <- get_branch_df(omega_hat, weights, lambda_grid)
    
    list(omega_hat = omega_hat,
         branch_df = branch_df)
  }
  
  future::plan(future::sequential)
  
  omega_hat <- lapply(result, `[[`, 1)
  branch_dfs <- lapply(result, `[[`, 2)
  
  list(omega_hat = omega_hat,
       branch_dfs = branch_dfs)
}

get_log_L_bar <- function(branches, psi_grid) {
  
  df_list <- branches$branch_dfs
  
  branch_interp_list <- df_list |> 
    lapply(\(df) {
      approx(df$psi, df$Integrated, xout = psi_grid, rule = 2)$y
      })
  
  branches_matrix <- do.call(rbind, branch_interp_list)
  
  log_R <- branches_matrix |> 
    nrow() |> 
    log()
  
  log_L_bar <- matrixStats::colLogSumExps(branches_matrix, na.rm = TRUE) - log_R
  
  log_L_bar_df <- data.frame(psi = psi_grid,
                             Integrated = log_L_bar)
  
  return(list(df = log_L_bar_df,
              branches_matrix = branches_matrix))
}

get_integrated_LL <- function(cfg) {
  
  branches <- generate_branches(cfg)
  
  psi_grid <- get_psi_grid(cfg)
  
  log_L_bar <- get_log_L_bar(branches, psi_grid)
  
  log_L_bar_plot <- make_plot(log_L_bar$df)
  
  print(log_L_bar_plot)
  
  output <- list(theta_0 = cfg$theta_0$values,
                 Y = cfg$data$Y,
                 weights = cfg$weights$values,
                 branches = branches,
                 psi_grid = psi_grid,
                 log_L_bar = log_L_bar,
                 plot = log_L_bar_plot)
  
  saveRDS(output, file.path(cfg$result_dir, "IL_object.rds"))
  
  return(output)
}

# Profile Likelihood ------------------------------------------------------

get_profile_LL <- function(cfg, envir = environment()) {
  
  list2env(cfg$data, envir)
  
  lambda_grid <- get_lambda_grid(cfg)
  
  psi_grid <- get_psi_grid(cfg)
  
  weights <- cfg$weights$values
  
  log_L_p_df <- get_branch_df(Y, weights, lambda_grid) |> 
    dplyr::rename(Profile = Integrated) |> 
    dplyr::filter(dplyr::between(psi, head(psi_grid, 1), tail(psi_grid, 1)))
  
  log_L_p_plot <- make_plot(log_L_p_df)
  
  print(log_L_p_plot)
  
  output <- list(theta_0 = cfg$theta_0$values,
                 Y = Y,
                 weights = weights,
                 log_L_p_df = log_L_p_df,
                 plot = log_L_p_plot)
  
  saveRDS(output, file.path(cfg$result_dir, "PL_object.rds"))
  
  return(output)
}

# Comparison --------------------------------------------------------------

get_PLL_names <- function(PLL_df_merged) {
  
  PLL_df |> 
    names() |> 
    setdiff("psi")
}

lengthen <- function(PLL_df_merged) {
  
  PLL_names <- get_PLL_names(PLL_df_merged)
  
  PLL_df_merged |>
    tidyr::pivot_longer(cols = dplyr::all_of(PLL_names),
                        names_to = "PLL",
                        values_to = "value") |>
    dplyr::mutate(PLL = forcats::fct_inorder(PLL)) |> 
    tidyr::drop_na()
}

fit_spline_models <- function(PLL_df_merged_long) {
  
  PLL_names <- PLL_df_merged_long |>
    dplyr::pull(PLL) |> 
    levels()
  
  PLL_df_merged_long |>
    dplyr::group_split(PLL) |> 
    purrr::map(\(df) with(df, smooth.spline(psi, value))) |>
    purrr::set_names(PLL_names)
}

get_PLL <- function(model, psi) predict(model, psi)$y

find_MLE <- function(PLL_df_long, model) {
  
  objective <- function(psi) get_PLL(model, psi)
  
  lower <- df_long |> 
    dplyr::select(psi) |> 
    min()
  
  upper <- df_long |> 
    dplyr::select(psi) |> 
    max()
  
  optimize(objective, 
           lower = lower, 
           upper = upper, 
           maximum = TRUE)
}

get_MLE_data <- function(PLL_df_merged_long, spline_fitted_models) {
  
  spline_fitted_models |>
    purrr::map(\(model) find_MLE(PLL_df_merged_long, model)) |>
    do.call(rbind.data.frame, args = _) |> 
    dplyr::rename(MLE = maximum,
                  Maximum = objective) |>
    dplyr::mutate(MLE_label = c("hat(psi)[I]", "hat(psi)[p]")) |>
    tibble::rownames_to_column("PLL")
}

subtract_max <- function(model, maximum) {
  
  function(psi) get_PLL(model, psi) - maximum
}

make_base_plot <- function() {
  
  base_plot <- withr:: with_package("ggplot2", {
    ggplot() +   
      theme_minimal() +
      theme(
        panel.background = element_rect(fill = "#444444", color = NA),  
        plot.background = element_rect(fill = "#2A2A2A", color = NA),  
        panel.grid.major = element_line(color = "#4C4C4C"),  
        panel.grid.minor = element_line(color = "#333333"),  
        axis.ticks = element_line(color = "white"),
        axis.text = element_text(color = "white"),
        axis.title = element_text(color = "white"),
        strip.text = element_text(color = "white"),
        plot.title = element_text(color = "white", face = "bold"),
        legend.background = element_rect(fill = "#444444"),
        legend.text = element_text(color = "white"),
        legend.title = element_text(color = "white")
        )
  })
  
  return(base_plot)
}

plot_PLL_points <- function(PLL_df_merged_long) {
  
  base_plot <- make_base_plot()
  
  PLL_points <- withr:: with_package("ggplot2", {
    base_plot +
      geom_point(aes(x = psi, y = value, color = PLL),
                 data = PLL_df_long,
                 size = 3) +
      ylab("Pseudo-Log-Likelihood") +
      scale_color_brewer(palette = "Set1") +
      xlab(expression(psi))
  })
  
  return(PLL_points)
}

find_root <- function(objective, interval) {
  
  tryCatch(
    
    uniroot(objective,
            interval)$root,
    
    error = function(e) return(NA)
  )
}

find_PLL_CI_endpoints <- function(PLLs, PLL_df_long, alpha) {
  
  PLL_types <- names(PLLs)
  
  search_intervals <- purrr::map2(PLLs, 
                                  PLL_types, 
                                  \(PLL, type) {
                                    PLL_df_long |> 
                                      dplyr::filter(PLL == type) |> 
                                      dplyr::select(psi) |> 
                                      range()
                                    })
  
  endpoints <- purrr::map2(PLLs, 
                           search_intervals,
                           \(PLL, interval) {
                             
                             
                           })
  
  critical_value <- qchisq(1 - alpha, 1) / 2
  
  objective <- function(psi) PLL(psi) + critical_value
  
  lower_endpoint <- find_root(objective, lower_interval)
  
  upper_endpoint <- find_root(objective, upper_interval)
  
  endpoints <- c(lower_endpoint, upper_endpoint) |> 
    round(3)
  
  return(endpoints)
}

compare_PLLs <- function(PLL_df_merged, alpha = 0.05) {
  
  PLL_df_merged_long <- lengthen(PLL_df_merged)
  
  PLL_points <- plot_PLL_points(PLL_df_merged_long)
  
  print(PLL_points)
  
  spline_fitted_models <- fit_spline_models(PLL_df_merged_long)
  
  MLE_data <- get_MLE_data(PLL_df_merged_long, spline_fitted_models)
  
  PLLs <- purrr::map2(spline_fitted_models, MLE_data$Maximum, subtract_max)
  
  critical_val <- qchisq(1 - alpha, 1) / 2
  
  psi_min <- min(PLL_df_merged_long$psi)
  
  psi_max <- max(PLL_df_merged_long$psi)
  
  PLL_CIs <- purrr::map2(PLLs, )
    
    
    
    PLLs |> 
    map2(MLE_data$MLE,
         function(curve, MLE) {
           
           lower_bound <- tryCatch(
             
             uniroot(function(psi) curve(psi) + crit,
                     interval = c(min_psi, MLE))$root,
             
             error = function(e) return(NA)
           ) |> 
             round(3)
           
           upper_bound <- tryCatch(
             
             uniroot(function(psi) curve(psi) + crit,
                     interval = c(MLE, max_psi))$root,
             
             error = function(e) return(NA)
           ) |>
             round(3)
           
           return(c(lower_bound, upper_bound))
         }
    )
}


