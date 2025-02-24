
# The below three functions all assume X_design is an object created by model.matrix( ~ X1*X2 - 1), 
# where X1 is a factor and X2 is continuous

get_X1_levels <- function(X_design) {
  
  is_X1_main_effect <- attr(X_design, "assign") == 1
  X1_main_effects <- colnames(X_design)[is_X1_main_effect]
  X1_levels <- substr(X1_main_effects, nchar(X1_main_effects), nchar(X1_main_effects))
  return(X1_levels)
}

get_X1_ref_level <- function(X_design) {
  
  X1_levels <- get_X1_levels(X_design)
  
  pattern <- paste0("[", paste(X1_levels, collapse = ""), "]")
  
  is_interaction <- attr(X_design, "assign") == 3
  
  interaction_terms <- colnames(X_design)[is_interaction]
  
  X1_interactions <- pattern |> 
    gregexpr(interaction_terms) |> 
    regmatches(interaction_terms, m = _) |> 
    unlist()
  
  X1_ref_level <- setdiff(X1_levels, X1_interactions)
  
  return(X1_ref_level)
}  

get_X2_samples <- function(X1, sigma) {
  
  m <- table(X1)
  
  X2_samples <- c()
  
  for (h in names(m)) {
    
    mu <- log(as.numeric(X1[X1==h]) + 1) |> 
      unique()
    
    X2_sample <- rnorm(m[h], mean = mu, sd = sigma)
    
    X2_samples <- c(X2_samples, X2_sample)
  }
  
  names(X2_samples) <- rep(names(m), times = m)
  
  return(X2_samples)
}

extract_X_h <- function(X_design, h, drop_zero_cols = FALSE) {
  
  X1_levels <- get_X1_levels(X_design)
  
  X1_ref_level <- get_X1_ref_level(X_design)
  
  X1_main_effect_names <- paste0("X1", X1_levels)
  
  X1X2_interaction_names <- paste0("X1", setdiff(X1_levels, X1_ref_level)) |> 
    paste0(":X2")
  
  colnames(X_design) <- c(X1_main_effect_names, "X2", X1X2_interaction_names)
  
  rows_to_keep <- X_design[, paste0("X1", h)] == 1
  
  X_h <- X_design[rows_to_keep,]
  
  if (drop_zero_cols) {
    
    cols_to_keep <- grepl(paste0("^X1", h, "$|^X2$|^X1", h, ":X2$"), colnames(X_design))
    
    X_h <- X_h[, cols_to_keep] |> 
      as.matrix(ncol = length(cols_to_keep))
    
    colnames(X_h) <- colnames(X_design)[cols_to_keep]
  }
  
  return(X_h)
}

get_entropy_ranges <- function(num_classes, levels) {
  
  num_ranges <- length(levels)
  
  entropy_ranges <- seq(0, log(num_classes), length.out = num_ranges + 1) |> 
    (\(x) mapply(c, x[-length(x)], x[-1], SIMPLIFY = FALSE))() |> 
    purrr::map(\(x) {
      
      midpoint <- mean(x)
      desired_length <- (x[2] - x[1]) * 0.5
      return(midpoint + c(-1, 1) * desired_length / 2)}
      )
  
  names(entropy_ranges) <- levels
  
  return(entropy_ranges)
}

get_Beta_0_h_entropy <- function(Beta_0_h, h, X1, X2) {
  
  X_design <- model.matrix(~ X1 * X2 - 1)
  
  X_design_h <- extract_X_h(X_design, h, drop_zero_cols = TRUE)
  
  X_design_h %*% cbind(0, Beta_0_h) |> 
    apply(1, softmax) |> 
    t() |> 
    colMeans() |> 
    entropy()
}

make_Beta_0_h_con_fn <- function(X_design, h, Beta2, X2) {
  
  X1_levels <- get_X1_levels(X_design)
  
  X1_ref_level <- get_X1_ref_level(X_design)
  
  X1 <- h |> 
    rep(length(X2)) |> 
    factor(levels = X1_levels) |> 
    relevel(X1_ref_level)
  
  J <- length(Beta2) + 1

  range <- get_entropy_ranges(J, X1_levels)[[h]]

  if (h == X1_ref_level) {

    Beta_0_h_con_fn <- function(vec) {
      
      Beta_0_h <- rbind(vec, Beta2) |>
        unname()
      
      Beta_0_h_entropy <- get_Beta_0_h_entropy(Beta_0_h, h, X1, X2)

      return(c(range[1] - Beta_0_h_entropy, Beta_0_h_entropy - range[2]))
    }
  }

  else {

    Beta_0_h_con_fn <- function(vec) {
      
      Beta_0_h <- rbind(head(vec, J - 1),
                        Beta2,
                        tail(vec, J - 1)) |> 
        unname()

      Beta_0_h_entropy <- get_Beta_0_h_entropy(Beta_0_h, h, X1, X2)

      return(c(range[1] - Beta_0_h_entropy, Beta_0_h_entropy - range[2]))
    }
  }

  return(Beta_0_h_con_fn)
}

get_Beta_0_h <- function(X_design, h, Beta2, X2) {
  
  X_h <- extract_X_h(X_design, h, drop_zero_cols = TRUE)
  
  q <- ncol(X_h) - 1
  
  J <- length(Beta2) + 1
  
  Beta_0_h_con_fn <- make_Beta_0_h_con_fn(X_design, h, Beta2, X2)
  
  init_guess <- rnorm(q * (J - 1))
  
  result <- nloptr::auglag(x0 = init_guess,
                           fn = function(vec) 0,
                           hin = Beta_0_h_con_fn,
                           localsolver = "LBFGS",
                           deprecatedBehavior = FALSE)$par
  
  Beta_0_h <- result |> 
    matrix(nrow = q,
           ncol = J - 1,
           byrow = TRUE)
}

get_Beta_0 <- function(X_design, Beta2, n_samples, sigma) {
  
  p <- ncol(X_design)
  
  J <- length(Beta2) + 1
  
  X1_levels <- get_X1_levels(X_design)
  
  X1_ref_level <- get_X1_ref_level(X_design)
  
  Beta_0 <- matrix(NA,
                   nrow = p,
                   ncol = J - 1)
  
  rownames(Beta_0) <- colnames(X_design)
  
  Beta_0["X2",] <- Beta2
  
  X1 <- X1_levels |> 
    rep(each = n_samples) |> 
    factor(levels = X1_levels) |> 
    relevel(X1_ref_level)
  
  X2_samples <- get_X2_samples(X1, sigma)
  
  for (i in seq_along(X1_levels)) {
    
    h <- X1_levels[i]
    
    X2 <- X2_samples[names(X2_samples) == h]
    
    Beta_0_h <- get_Beta_0_h(X_design, h, Beta2, X2)
    
    Beta_0[grepl(h, rownames(Beta_0)), ] <- Beta_0_h
  }
  
  return(Beta_0)
}

Beta2 <- rnorm(J - 1, sd = 0.3)

Beta_0 <- get_Beta_0(X_design, Beta2, n_samples, sigma)



