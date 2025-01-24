#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector get_Beta_hat(
    NumericVector init_guess, 
    double psi, 
    NumericMatrix omega_hat, 
    NumericMatrix X_one_hot, 
    NumericMatrix X_h_one_hot
) {
  // Softmax function to be applied on probabilities
  auto adj_softmax = [](NumericVector x) {
    NumericVector result(x.size());
    double exp_sum = sum(exp(x));
    for (int i = 0; i < x.size(); i++) {
      result[i] = exp(x[i]) / (1 + exp_sum);
    }
    return result;
  };
  
  // Compute probabilities using the softmax
  NumericMatrix probs = X_one_hot * omega_hat;  // Matrix multiplication for X_one_hot and omega_hat
  for (int i = 0; i < probs.nrow(); i++) {
    probs(i, _) = adj_softmax(probs(i, _));  // Apply softmax to each row
  }
  
  // Objective function to minimize (f)
  Function f = [=](NumericVector Beta) {
    // Convert Beta into a matrix form for proper multiplication
    NumericMatrix Beta_mat(Beta.begin(), omega_hat.nrow(), omega_hat.ncol(), false);  // Convert Beta to matrix
    NumericMatrix Y_hat = X_one_hot * Beta_mat;  // Matrix multiplication for X_one_hot and Beta_mat
    NumericVector row_sums = rowSums(probs * Y_hat); // Sum across rows
    return -sum(row_sums - log(1 + rowSums(exp(Y_hat))));  // Calculate the negative log-likelihood
  };
  
  // Gradient of the objective function
  Function f_gr = [=](NumericVector Beta) {
    return nloptr::nl.grad(Beta, f);  // Gradient using nloptr
  };
  
  // Constraint function (entropy constraint)
  Function fcon = [=](NumericVector Beta) {
    // Convert Beta into matrix form
    NumericMatrix Beta_mat(Beta.begin(), omega_hat.nrow(), omega_hat.ncol(), false);  // Convert Beta to matrix
    NumericVector entropy = X_h_one_hot * cbind(0, Beta_mat);  // Compute entropy from X_h_one_hot and Beta_mat
    return get_entropy(entropy) - psi;  // Ensure entropy matches psi
  };
  
  // Jacobian of the constraint function
  Function fcon_jac = [=](NumericVector Beta) {
    return nloptr::nl.jacobian(Beta, fcon);  // Jacobian using nloptr
  };
  
  // Call the auglag function for constrained optimization
  List auglag_result = auglag(
    _["x0"] = init_guess,
    _["fn"] = f,
    _["heq"] = fcon,
    _["gr"] = f_gr,
    _["heqjac"] = fcon_jac,
    _["localsolver"] = "LBFGS"
  );
  
  // Extract the optimized Beta from the result
  NumericVector result = auglag_result["par"];
  
  return result;
}
