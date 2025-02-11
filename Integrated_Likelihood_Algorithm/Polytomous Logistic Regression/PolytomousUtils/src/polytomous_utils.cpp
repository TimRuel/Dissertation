// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>

using namespace Rcpp;

// ----------------------------
// Log-Likelihood Function 
// ----------------------------
// [[Rcpp::export]]
double log_likelihood_rcpp(NumericVector Beta, NumericMatrix X_one_hot, NumericMatrix Y_one_hot,
                           int Jm1, int p, int n) {
  
  // Convert Beta vector into a matrix
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
  
  // Compute row-wise exponential sum for normalization
  NumericVector row_sums(n, 1.0); // Start with 1 for log(1 + sum(exp(Y_hat)))
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < Jm1; j++) {
      row_sums[i] += std::exp(Y_hat(i, j));
    }
  }
  
  // Compute log likelihood = sum(rowSums(Y_one_hot * Y_hat)) - sum(log(1 + row_sums))
  double log_likelihood = 0.0;
  for (int i = 0; i < n; i++) {
    double row_sum_Y_hat = 0.0;
    for (int j = 0; j < Jm1; j++) {
      row_sum_Y_hat += Y_one_hot(i, j) * Y_hat(i, j);
    }
    log_likelihood += row_sum_Y_hat - std::log(row_sums[i]);
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
// [[Rcpp::export]]
double entropy_rcpp(NumericVector p) {
  double entropy = 0.0;
  for (int i = 0; i < p.size(); i++) {
    if (p[i] > 0) {
      entropy -= p[i] * std::log(p[i]);
    }
  }
  return entropy;
}

// ----------------------------
// PoI Function
// ----------------------------
// [[Rcpp::export]]
double PoI_fn_rcpp(NumericMatrix Beta, NumericVector X_h_one_hot) {
  int p = Beta.nrow();
  int Jm1 = Beta.ncol();
  
  // Extend Beta by adding a column of zeros at the beginning
  NumericMatrix Beta_extended(p, Jm1 + 1);
  for (int i = 0; i < p; i++) {
    for (int j = 0; j < Jm1; j++) {
      Beta_extended(i, j + 1) = Beta(i, j); // Shift existing columns right
    }
  }
  
  // Compute X_h_one_hot %*% Beta_extended
  NumericVector X_h_Beta(Jm1 + 1);
  for (int j = 0; j < Jm1 + 1; j++) {
    double dot_product = 0.0;
    for (int i = 0; i < p; i++) {
      dot_product += X_h_one_hot[i] * Beta_extended(i, j);
    }
    X_h_Beta[j] = dot_product;
  }
  
  // Compute softmax probabilities
  NumericVector probs = softmax_vector_rcpp(X_h_Beta);
  
  // Compute entropy
  double entropy = entropy_rcpp(probs);
  
  return entropy;
}

// ----------------------------
// Beta_hat Objective Function
// ----------------------------
// [[Rcpp::export]]
double Beta_hat_obj_fn_rcpp(NumericVector Beta, NumericMatrix X_one_hot, 
                            NumericMatrix omega_hat, double lambda,
                            int Jm1, int p, int n) {
  
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
  
  // Compute row-wise exponential sum for normalization
  NumericVector row_sums(n, 1.0); // Start with 1 for log(1 + sum(exp(Y_hat)))
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < Jm1; j++) {
      row_sums[i] += std::exp(Y_hat(i, j));
    }
  }
  
  // Compute log likelihood term
  double log_L = 0.0;
  for (int i = 0; i < n; i++) {
    double row_sum = 0.0;
    for (int j = 0; j < Jm1; j++) {
      row_sum += probs(i, j) * Y_hat(i, j);
    }
    log_L += row_sum - std::log(row_sums[i]);
  }
  
  // Compute regularization term
  double reg_term = 0.0;
  for (int i = 0; i < Beta.size(); i++) {
    reg_term += Beta[i] * Beta[i];
  }
  
  return -log_L + lambda * reg_term;
}

// ----------------------------
// Beta_hat Constraint Function
// ----------------------------
// [[Rcpp::export]]
double Beta_hat_con_fn_rcpp(NumericVector Beta, NumericVector X_h_one_hot, double psi, int Jm1, int p) {
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
double omega_hat_obj_fn_rcpp(NumericVector Beta, NumericVector X_h_one_hot, int Jm1, int p, double psi_hat) {
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

// -----------------------------
// Omega_hat Constraint Function
// -----------------------------
// [[Rcpp::export]]
double omega_hat_con_fn_rcpp(NumericVector Beta, NumericMatrix X_one_hot, NumericMatrix Y_one_hot, int Jm1, int p, int n, double threshold) {
  // Compute log-likelihood
  double log_likelihood = log_likelihood_rcpp(Beta, X_one_hot, Y_one_hot, Jm1, p, n);

  return -log_likelihood - threshold;
}

