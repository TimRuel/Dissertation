get_X1 <- function(X1_levels) {
  
  X1_level_names <- names(X1_levels)
  
  X1_ref_level <- X1_levels |> 
    map_lgl(\(x) x$ref_level) |> 
    which() |> 
    names()
  
  m <- map_dbl(X1_levels, \(x) x$m)
  
  X1_level_names |>
    rep(times = m) |>
    factor() |>
    relevel(X1_ref_level)
}

get_X2 <- function(X1_levels) {
  
  X1_levels |> 
    imap(\(level, h) {
      m <- level$m
      list2env(level$X2, environment())
      args <- as.list(c(m, dist$params))
      samples <- do.call(dist$name, args)
      X2_obs <- support[1] + diff(support) * samples
      set_names(X2_obs, rep(h, m))
    }) |> 
    unname() |>
    unlist()
}

get_design_matrix <- function(formula, ...) {
  
  formula <- as.formula(formula)

  df <- data.frame(...)
  mf <- model.frame(formula, data = df)
  mm <- model.matrix(attr(mf, "terms"), data = mf)
  
  attr(mm, "original_model_frame") <- mf
  attr(mm, "terms") <- terms(formula, data = df)
  attr(mm, "formula") <- formula
  attr(mm, "contrasts") <- attr(mm, "contrasts")
  
  return(mm)
}

recover_original_data <- function(X_design) {
  
  mf <- attr(X_design, "original_model_frame")
  
  if (is.null(mf)) {
    stop("No model frame stored in the design matrix.")
  }
  
  as.data.frame(mf)
}

get_Y_probs <- function(X_design, Beta_0) {
  
  df <- recover_original_data(X_design)
  
  Y_probs <- X_design %*% cbind(0, Beta_0) |>
    apply(1, softmax) |>
    t() |>
    data.frame() |>
    rename_with( ~ paste0("Y", 1:(ncol(Beta_0) + 1)))
  
  cbind(df, Y_probs)
}

get_Y <- function(Y_probs) {
  
  df <- Y_probs |>
    select(starts_with("Y"))
  
  J <- ncol(df)
  
  df |>
    apply(1, \(prob) sample(1:J, size = 1, prob = prob)) |>
    unlist() |>
    unname() |>
    factor(levels = 1:J)
}

get_data <- function(X1_levels, formula, Beta_0) {
  
  X1 <- get_X1(X1_levels)
  
  X2 <- get_X2(X1_levels)
  
  X_design <- get_design_matrix(formula, X1, X2)
  
  Y_probs <- get_Y_probs(X_design, Beta_0)
  
  Y <- get_Y(Y_probs)
  
  model_df <- data.frame(X1 = X1,
                   X2 = X2,
                   Y)
  
  data <- list(model_df = model_df,
               X_design = X_design,
               Y_probs = Y_probs)
  
  return(data)
}
