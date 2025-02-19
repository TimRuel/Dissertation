get_entropy_ranges <- function(num_classes, levels) {
  
  num_ranges <- length(levels) + 2
  
  entropy_ranges <- seq(0, log(num_classes), length.out = num_ranges + 1) |> 
    (\(x) mapply(c, x[-length(x)], x[-1], SIMPLIFY = FALSE))() |> 
    tail(-1) |> 
    head(-1)
  
  names(entropy_ranges) <- levels
  
  return(entropy_ranges)
}

get_Beta_0_h_entropy <- function(Beta_0_h, h) {
  
  X_h <- extract_X_h(X, h, drop_zero_cols = TRUE)
  
  X_h %*% cbind(0, Beta_0_h) |> 
    apply(1, softmax) |> 
    t() |> 
    colMeans() |> 
    entropy()
}

make_Beta_0_h_con_fn <- function(h, X2_coef) {
  
  X_h <- extract_X_h(X, h, drop_zero_cols = TRUE)
  
  range <- get_entropy_ranges(J, X1_levels)[[h]]
  
  if (h == X1_ref_level) {
    
    Beta_0_h_con_fn <- function(vec) {
      
      entropy <- rbind(vec, X2_coef) |> 
        unname() |> 
        get_Beta_0_h_entropy(h)
      
      return(c(range[1] - entropy, entropy - range[2]))
    }
  }
  
  else {
    
    Beta_0_h_con_fn <- function(vec) {
      
      entropy <- rbind(head(vec, J - 1), 
                       X2_coef,
                       tail(vec, J - 1)) |> 
        unname() |> 
        get_Beta_0_h_entropy(h)
      
      return(c(range[1] - entropy, entropy - range[2]))
    }
  }
  
  return(Beta_0_h_con_fn)
}

get_Beta_0 <- function(X2_coef) {
  
  Beta_0 <- matrix(NA,
                   nrow = p,
                   ncol = J - 1)
  
  rownames(Beta_0) <- colnames(X)
  
  Beta_0["X2",] <- X2_coef
  
  for (h in X1_levels) {
    
    X_h <- extract_X_h(X, h, drop_zero_cols = TRUE)
    
    q <- ncol(X_h) - 1
    
    Beta_0_h_con_fn <- make_Beta_0_h_con_fn(h, X2_coef)
    
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
    
    Beta_0[grepl(h, rownames(Beta_0)), ] <- Beta_0_h
  }
  
  true_probs <- X %*% cbind(0, Beta_0) |> 
    apply(1, softmax) |> 
    t()
  
  ranges <- get_entropy_ranges(J, X1_levels)
  
  done <- TRUE
  
  for (i in 1:length(m)) {
    
    done1 <- true_probs[(1 + (i-1) * m[i]):(m[i]*i),] |> 
      colMeans() |> 
      entropy() |> 
      between(ranges[[i]][1], ranges[[i]][2]) |> 
      all()
    
    done <- c(done, done1)
    
    # Y <- true_probs |>
    #   apply(1, \(prob) sample(1:J, size = 1, prob = prob)) |> 
    #   unlist() |>
    #   unname() |> 
    #   factor(levels = 1:J)
    # 
    # done2 <- min(table(Y)) >= 3
    # 
    # done <- all(done1, done2)
  }
  
  if (all(done)) return(Beta_0)
  
  # if (done) return(Beta_0)
  
  else get_Beta_0(X2_coef)
}
