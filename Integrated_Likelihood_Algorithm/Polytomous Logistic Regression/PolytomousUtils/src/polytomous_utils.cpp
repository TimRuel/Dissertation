// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>

using namespace Rcpp;

#include <Rcpp.h>
using namespace Rcpp;

// ----------------------------
// rowLogSumExps Function 
// ----------------------------
// [[Rcpp::export]]
NumericVector rowLogSumExps_rcpp(NumericMatrix Y_hat) {
  int n = Y_hat.nrow();
  int Jm1 = Y_hat.ncol();
  int J = Jm1 + 1; // Account for extra column
  
  // Create a new matrix with an additional first column of zeros
  NumericMatrix Y_hat_extended(n, J);
  
  for (int i = 0; i < n; i++) {
    Y_hat_extended(i, 0) = 0.0; // First column of zeros
    for (int j = 0; j < Jm1; j++) {
      Y_hat_extended(i, j + 1) = Y_hat(i, j);
    }
  }
  
  // Compute row-wise logSumExp
  NumericVector logsumexp(n);
  
  for (int i = 0; i < n; i++) {
    double max_val = Y_hat_extended(i, 0);
    
    for (int j = 1; j < J; j++) {
      if (Y_hat_extended(i, j) > max_val) {
        max_val = Y_hat_extended(i, j);
      }
    }
    
    double exp_sum = 0.0;
    for (int j = 0; j < J; j++) {
      exp_sum += std::exp(Y_hat_extended(i, j) - max_val);
    }
    
    logsumexp[i] = max_val + std::log(exp_sum);
  }
  
  return logsumexp;
}

// ----------------------------
// Log-Likelihood Function 
// ----------------------------
// [[Rcpp::export]]
double log_likelihood_rcpp(NumericVector Beta, NumericMatrix X_one_hot, NumericMatrix Y_one_hot,
                           int Jm1, int p, int n) {
  
  // Convert Beta vector into a matrix
  NumericMatrix Beta_mat(p, Jm1);
  std::memcpy(Beta_mat.begin(), Beta.begin(), Beta.size() * sizeof(double));
  
  // Compute Y_hat = X_one_hot %*% Beta_mat
  NumericMatrix Y_hat(n, Jm1);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < Jm1; j++) {
      double sum = 0.0;
      for (int k = 0; k < p; k++) {
        sum += X_one_hot(i, k) * Beta_mat(k, j);
      }
      Y_hat(i, j) = sum;
    }
  }
  
  // Compute logsumexp using our function
  NumericVector logsumexp = rowLogSumExps_rcpp(Y_hat);
  
  // Compute log-likelihood
  double log_likelihood = 0.0;
  for (int i = 0; i < n; i++) {
    double row_sum_Y_hat = 0.0;
    for (int j = 0; j < Jm1; j++) {
      row_sum_Y_hat += Y_one_hot(i, j) * Y_hat(i, j);
    }
    log_likelihood += row_sum_Y_hat - logsumexp[i];
  }
  
  return log_likelihood;
}

// -------------------------------------
// Adjusted Softmax Function for Vectors
// -------------------------------------
// [[Rcpp::export]]
NumericVector softmax_adj_vector_rcpp(NumericVector x) {
  int n = x.size();
  NumericVector exp_x(n);
  double sum_exp = 0;
  
  // Compute exponentials and sum
  for (int i = 0; i < n; i++) {
    exp_x[i] = std::exp(x[i]);
    sum_exp += exp_x[i];
  }
  
  // Normalize by sum_exp + 1
  for (int i = 0; i < n; i++) {
    exp_x[i] /= (1 + sum_exp);
  }
  
  return exp_x;
}

// --------------------------------------
// Adjusted Softmax Function for Matrices
// --------------------------------------
// [[Rcpp::export]]
NumericMatrix softmax_adj_matrix_rcpp(NumericMatrix x) {
  int rows = x.nrow();
  int cols = x.ncol();
  NumericMatrix exp_x(rows, cols);
  NumericVector row_sums(rows);
  
  // Compute exponentials and row sums
  for (int i = 0; i < rows; i++) {
    double sum_exp = 0;
    for (int j = 0; j < cols; j++) {
      exp_x(i, j) = std::exp(x(i, j));
      sum_exp += exp_x(i, j);
    }
    row_sums[i] = sum_exp;
  }
  
  // Normalize each row
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      exp_x(i, j) /= (1 + row_sums[i]);
    }
  }
  
  return exp_x;
}

// ----------------------------
// Softmax Function for Vectors
// ----------------------------
// [[Rcpp::export]]
NumericVector softmax_vector_rcpp(NumericVector x) {
  int n = x.size();
  NumericVector exp_x(n);
  double sum_exp = 0;
  
  // Compute exponentials and sum
  for (int i = 0; i < n; i++) {
    exp_x[i] = std::exp(x[i]);
    sum_exp += exp_x[i];
  }
  
  // Normalize
  for (int i = 0; i < n; i++) {
    exp_x[i] /= sum_exp;
  }
  
  return exp_x;
}

// -----------------------------
// Softmax Function for Matrices
// -----------------------------
// [[Rcpp::export]]
NumericMatrix softmax_matrix_rcpp(NumericMatrix x) {
  int rows = x.nrow();
  int cols = x.ncol();
  NumericMatrix exp_x(rows, cols);
  NumericVector row_sums(rows);
  
  // Compute exponentials and row sums
  for (int i = 0; i < rows; i++) {
    double sum_exp = 0;
    for (int j = 0; j < cols; j++) {
      exp_x(i, j) = std::exp(x(i, j));
      sum_exp += exp_x(i, j);
    }
    row_sums[i] = sum_exp;
  }
  
  // Normalize each row
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      exp_x(i, j) /= row_sums[i];
    }
  }
  
  return exp_x;
}

// ----------------------------
// Entropy Calculation
// ----------------------------
// // [[Rcpp::export]]
// double PoI_fn_rcpp(NumericMatrix Beta, NumericVector X_h_one_hot) {
//   int p = Beta.nrow();
//   int Jm1 = Beta.ncol();
//   
//   // Compute X_h_one_hot %*% Beta_extended without explicitly creating Beta_extended
//   NumericVector X_h_Beta(Jm1 + 1, 0.0); // Initialize with zeros (first column is implicitly zero)
//   
//   for (int j = 0; j < Jm1; j++) {
//     double dot_product = 0.0;
//     for (int i = 0; i < p; i++) {
//       dot_product += X_h_one_hot[i] * Beta(i, j);
//     }
//     X_h_Beta[j + 1] = dot_product; // Shift results to simulate Beta_extended
//   }
//   
//   // Compute softmax in place
//   double max_val = max(X_h_Beta);
//   NumericVector exp_vals(Jm1 + 1);
//   double sum_exp = 0.0;
//   
//   for (int j = 0; j < Jm1 + 1; j++) {
//     exp_vals[j] = std::exp(X_h_Beta[j] - max_val); // Avoid overflow issues
//     sum_exp += exp_vals[j];
//   }
//   
//   NumericVector probs(Jm1 + 1);
//   for (int j = 0; j < Jm1 + 1; j++) {
//     probs[j] = exp_vals[j] / sum_exp;
//   }
//   
//   // Compute entropy: H(p) = -sum(p * log(p))
//   double entropy = 0.0;
//   for (int j = 0; j < Jm1 + 1; j++) {
//     if (probs[j] > 0) { // Avoid log(0)
//       entropy -= probs[j] * std::log(probs[j]);
//     }
//   }
//   
//   return entropy;
// }

// [[Rcpp::export]]
double PoI_fn_rcpp(NumericMatrix Beta, NumericMatrix X_h_one_hot) {
  int p = Beta.nrow();
  int Jm1 = Beta.ncol();
  int n = X_h_one_hot.nrow();  // Number of rows in X_h_one_hot
  
  NumericVector avg_probs(Jm1 + 1, 0.0); // Initialize mean probability vector
  
  // Iterate over each row of X_h_one_hot
  for (int r = 0; r < n; r++) {
    NumericVector X_h_Beta(Jm1 + 1, 0.0); // First column implicitly zero
    
    // Compute X_h_one_hot[r, ] %*% Beta (shifting Beta values)
    for (int j = 0; j < Jm1; j++) {
      double dot_product = 0.0;
      for (int i = 0; i < p; i++) {
        dot_product += X_h_one_hot(r, i) * Beta(i, j);
      }
      X_h_Beta[j + 1] = dot_product;
    }
    
    // Compute softmax
    double max_val = max(X_h_Beta); // Stabilization
    double sum_exp = 0.0;
    NumericVector exp_vals(Jm1 + 1);
    
    for (int j = 0; j < Jm1 + 1; j++) {
      exp_vals[j] = std::exp(X_h_Beta[j] - max_val);
      sum_exp += exp_vals[j];
    }
    
    for (int j = 0; j < Jm1 + 1; j++) {
      avg_probs[j] += exp_vals[j] / sum_exp; // Directly accumulate into avg_probs
    }
  }
  
  // Finalize mean probabilities
  for (int j = 0; j < Jm1 + 1; j++) {
    avg_probs[j] /= n;
  }
  
  // Compute entropy of the mean probability vector
  double entropy = 0.0;
  for (int j = 0; j < Jm1 + 1; j++) {
    if (avg_probs[j] > 0) { // Avoid log(0)
      entropy -= avg_probs[j] * std::log(avg_probs[j]);
    }
  }
  
  return entropy;
}

// ----------------------------
// Beta_hat Objective Function
// ----------------------------
// [[Rcpp::export]]
double Beta_hat_obj_fn_rcpp(NumericVector Beta, NumericMatrix X_one_hot, 
                            NumericMatrix omega_hat, int Jm1, int p, int n) {
  
  // Reshape Beta into a matrix
  NumericMatrix Beta_mat(p, Jm1);
  for (int j = 0; j < Jm1; j++) {
    for (int i = 0; i < p; i++) {
      Beta_mat(i, j) = Beta[j * p + i];
    }
  }
  
  // Compute Y_hat = X_one_hot %*% Beta_mat
  NumericMatrix Y_hat(n, Jm1);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < Jm1; j++) {
      double sum = 0.0;
      for (int k = 0; k < p; k++) {
        sum += X_one_hot(i, k) * Beta_mat(k, j);
      }
      Y_hat(i, j) = sum;
    }
  }
  
  // Compute X_omega = X_one_hot %*% omega_hat
  NumericMatrix X_omega(n, Jm1);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < Jm1; j++) {
      double sum = 0.0;
      for (int k = 0; k < p; k++) {
        sum += X_one_hot(i, k) * omega_hat(k, j);
      }
      X_omega(i, j) = sum;
    }
  }
  
  // Compute softmax probabilities
  NumericMatrix probs = softmax_adj_matrix_rcpp(X_omega);
  
  // Compute logsumexp using our function
  NumericVector logsumexp = rowLogSumExps_rcpp(Y_hat);
  
  // Compute log likelihood term
  double log_L = 0.0;
  for (int i = 0; i < n; i++) {
    double row_sum_Y_hat = 0.0;
    for (int j = 0; j < Jm1; j++) {
      row_sum_Y_hat += probs(i, j) * Y_hat(i, j);
    }
    log_L += row_sum_Y_hat - logsumexp[i];
  }
  
  return -log_L;
}

// ----------------------------
// Beta_hat Constraint Function
// ----------------------------
// [[Rcpp::export]]
double Beta_hat_con_fn_rcpp(NumericVector Beta, NumericMatrix X_h_one_hot, double psi, int Jm1, int p) {
  // Reshape Beta into a matrix
  NumericMatrix Beta_mat(p, Jm1);
  for (int j = 0; j < Jm1; j++) {
    for (int i = 0; i < p; i++) {
      Beta_mat(i, j) = Beta[j * p + i];
    }
  }

  // Compute PoI function
  double entropy = PoI_fn_rcpp(Beta_mat, X_h_one_hot);

  return entropy - psi;
}

// ----------------------------
// Omega_hat Objective Function
// ----------------------------
// [[Rcpp::export]]
double omega_hat_obj_fn_rcpp(NumericVector Beta, NumericMatrix X_h_one_hot, int Jm1, int p, double psi_hat) {
  // Reshape Beta into a matrix
  NumericMatrix Beta_mat(p, Jm1);
  for (int j = 0; j < Jm1; j++) {
    for (int i = 0; i < p; i++) {
      Beta_mat(i, j) = Beta[j * p + i];
    }
  }

  // Compute PoI function
  double entropy = PoI_fn_rcpp(Beta_mat, X_h_one_hot);

  return std::abs(entropy - psi_hat);
}

// --------------------------------------
// Omega_hat Equality Constraint Function
// --------------------------------------
// [[Rcpp::export]]
double omega_hat_eq_con_fn_rcpp(NumericVector Beta, NumericMatrix X_h_one_hot, int Jm1, int p, double psi_hat) {
  // Reshape Beta into a matrix
  NumericMatrix Beta_mat(p, Jm1);
  for (int j = 0; j < Jm1; j++) {
    for (int i = 0; i < p; i++) {
      Beta_mat(i, j) = Beta[j * p + i];
    }
  }
  
  // Compute PoI function
  double entropy = PoI_fn_rcpp(Beta_mat, X_h_one_hot);
  
  return entropy - psi_hat;
}

// ----------------------------------------
// Omega_hat Inequality Constraint Function
// ----------------------------------------
// [[Rcpp::export]]
double omega_hat_ineq_con_fn_rcpp(NumericVector Beta, NumericMatrix X_one_hot, NumericMatrix Y_one_hot, int Jm1, int p, int n, double threshold) {
  // Compute log-likelihood
  double log_likelihood = log_likelihood_rcpp(Beta, X_one_hot, Y_one_hot, Jm1, p, n);

  return -log_likelihood - threshold;
}

