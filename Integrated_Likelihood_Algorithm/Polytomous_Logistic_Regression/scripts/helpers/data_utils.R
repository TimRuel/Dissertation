# Data Generation ---------------------------------------------------------

get_X1 <- function(X1_level_specs) {
  
  X1_levels <- names(X1_level_specs)
  
  X1_ref_level <- X1_level_specs |> 
    purrr::map_lgl(\(x) x$ref_level) |> 
    which() |> 
    names()
  
  m <- purrr::map_dbl(X1_level_specs, \(x) x$num_obs)
  
  X1_levels |>
    rep(times = m) |>
    factor() |>
    relevel(X1_ref_level)
}

# get_X1 <- function(cfg, envir = environment()) {
#   
#   list2env(cfg$X1, envir)
#   
#   X1_levels <- names(levels)
#   
#   X1_ref_level <- levels |> 
#     purrr::map_lgl(\(x) x$ref_level) |> 
#     which() |> 
#     names()
#   
#   m <- purrr::map_dbl(levels, \(x) x$num_obs)
#   
#   X1_levels |>
#     rep(times = m) |>
#     factor() |>
#     relevel(X1_ref_level)
# }

get_X2_within_X1_level <- function(X1_level_spec) {
  
  list2env(X2_specs, environment())
  
  dist |> 
    do.call(dist_params) |> 
    (\(X2) interval[1] + diff(interval) * X2)()
}

get_X2 <- function(X1_level_specs) {
  
  X1_level_specs |> 
    purrr::imap(\(level_spec, level) {
      
      list2env(level_spec$X2, envir)
      
      X2_vals <- if (is.null(values)) {
        list2env(specs, envir)
        args <- as.list(c(num_obs, shape_params))
        samples <- do.call(rbeta, args)
        interval[1] + diff(interval) * samples
        
      } else {
        values
      }
      
      set_names(X2_vals, rep(level, num_obs))
    }) |> 
    unname() |>
    unlist()
}

get_X2 <- function(X1_specs, envir = environment()) {
  
  list2env(X1_specs, envir)
  
  X1_levels |> 
    purrr::imap(\(level, X1_level) {
      
      list2env(level$X2, envir)
      
      X2_vals <- if (is.null(values)) {
        list2env(specs, envir)
        args <- as.list(c(num_obs, shape_params))
        samples <- do.call(rbeta, args)
        interval[1] + diff(interval) * samples
        
      } else {
        values
      }
      
      set_names(X2_vals, rep(X1_level, num_obs))
    }) |> 
    unname() |>
    unlist()
}

get_Y_probs <- function(cfg, envir = environment()) {
  
  list2env(cfg$X1, envir)
  
  X1_ref_level <- levels |> 
    purrr::map_lgl(\(x) x$ref_level) |> 
    which() |> 
    names()
  
  m <- levels |> 
    purrr::map_dbl(\(x) x$num_obs)
  
  X1 <- levels |>
    names() |> 
    rep(times = m) |>
    factor() |>
    relevel(X1_ref_level)
  
  X2 <- get_X2(cfg, envir)
  
  X_design <- model.matrix(~ X1*X2 - 1)
  
  Beta_0 <- cfg$Beta_0$values
  
  X_design %*% cbind(0, Beta_0) |>
    apply(1, softmax) |>
    t() |>
    data.frame() |>
    rename_with( ~ paste0("Y", 1:(ncol(Beta_0) + 1))) |>
    mutate(X1 = X1,
           X2 = X2) |>
    select(X1, X2, everything())
}

get_Y <- function(cfg, envir = environment()) {
  
  
}