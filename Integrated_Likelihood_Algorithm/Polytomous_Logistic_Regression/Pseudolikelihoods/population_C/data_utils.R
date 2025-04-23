get_entropy_ranges <- function(num_classes, levels) {
  
  num_ranges <- length(levels) + 2
  
  entropy_ranges <- seq(0, log(num_classes), length.out = num_ranges + 1) |> 
    (\(x) mapply(c, x[-length(x)], x[-1], SIMPLIFY = FALSE))() |> 
    tail(-1) |> 
    head(-1)
  
  names(entropy_ranges) <- levels
  
  return(entropy_ranges)
}

get_Beta_0_h_entropy <- function(Beta_0_h) {
  
  cbind(0, Beta_0_h) |> 
    softmax() |> 
    entropy()
}

make_Beta_0_h_con_fn <- function(h) {
  
  range <- get_entropy_ranges(J, X1_levels)[[h]]
  
  Beta_0_h_con_fn <- function(Beta_0_h) {
    
    Beta_0_h <- Beta_0_h |> 
      matrix(nrow = 1,
             ncol = J - 1,
             byrow = TRUE)
    
    entropy <- get_Beta_0_h_entropy(Beta_0_h)
    
    return(c(range[1] - entropy, entropy - range[2]))
  }
  
  return(Beta_0_h_con_fn)
}

get_Beta_0 <- function() {
  
  Beta_0 <- matrix(NA,
                   nrow = p,
                   ncol = J - 1)
  
  rownames(Beta_0) <- colnames(X)
  
  for (h in X1_levels) {
    
    Beta_0_h_con_fn <- make_Beta_0_h_con_fn(h)
    
    init_guess <- rnorm(J - 1)
    
    result <- nloptr::auglag(x0 = init_guess,
                             fn = function(vec) 0,
                             hin = Beta_0_h_con_fn,
                             localsolver = "LBFGS",
                             deprecatedBehavior = FALSE)$par
    
    Beta_0_h <- result |> 
      matrix(nrow = 1,
             ncol = J - 1,
             byrow = TRUE)
    
    Beta_0[grepl(h, rownames(Beta_0)), ] <- Beta_0_h
  }
  
  true_probs <- X %*% cbind(0, Beta_0) |> 
    apply(1, softmax) |> 
    t() |> 
    unique()
  
  ranges <- get_entropy_ranges(J, X1_levels)
    
  done <- true_probs |> 
    apply(1, entropy) |> 
    purrr::map2_lgl(ranges, \(num, range) between(num, range[1], range[2])) |> 
    all()
  
  if (done) return(Beta_0)
  
  else get_Beta_0()
}

