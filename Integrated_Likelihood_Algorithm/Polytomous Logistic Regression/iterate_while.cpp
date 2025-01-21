#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix iterate_while(
    double psi_max,               // Starting point (argmax of function)
    double branch_max,            // Maximum branch value
    NumericVector initial_guess,  // Initial guess for the iterative process
    Function get_Beta_hat,        // Function to compute Beta_hat
    Function get_log_likelihood,  // Function to compute log-likelihood
    NumericMatrix omega_hat,      // Input matrix for entropy calculation
    NumericMatrix X_one_hot,      // Input matrix for likelihood calculation
    NumericVector X_h_one_hot,    // Input vector for entropy calculation
    NumericMatrix Y_one_hot,      // Input vector for likelihood calculation
    double crit,                  // Criterion for stopping condition
    double step_size              // Step size for incrementing psi
) {
  // Initialize variables
  double stopping_branch_val = branch_max - crit;
  double current_branch_val_forward = branch_max;
  double current_branch_val_backward = branch_max;
  
  double current_psi_val_forward = std::ceil(psi_max / step_size) * step_size;
  double current_psi_val_backward = std::floor(psi_max / step_size) * step_size;
  
  std::vector<NumericVector> results;  // Store results temporarily
  std::vector<double> psi_values;     // Store psi values
  
  // Forward loop
  NumericVector current_guess = clone(initial_guess);
  while (current_branch_val_forward > stopping_branch_val) {
    NumericVector Beta_hat = as<NumericVector>(
      get_Beta_hat(current_guess, current_psi_val_forward, omega_hat, X_one_hot, X_h_one_hot)
    );
    results.push_back(Beta_hat);
    psi_values.push_back(current_psi_val_forward);
    
    // Update log-likelihood and psi value
    current_branch_val_forward = as<double>(
      get_log_likelihood(Beta_hat, X_one_hot, Y_one_hot)
    );
    current_psi_val_forward += step_size;
    current_guess = Beta_hat;
  }
  
  // Backward loop
  current_guess = clone(initial_guess);
  while (current_branch_val_backward > stopping_branch_val) {
    NumericVector Beta_hat = as<NumericVector>(
      get_Beta_hat(current_guess, current_psi_val_backward, omega_hat, X_one_hot, X_h_one_hot)
    );
    results.insert(results.begin(), Beta_hat);  // Insert at the beginning
    psi_values.insert(psi_values.begin(), current_psi_val_backward);
    
    // Update log-likelihood and psi value
    current_branch_val_backward = as<double>(
      get_log_likelihood(Beta_hat, X_one_hot, Y_one_hot)
    );
    current_psi_val_backward -= step_size;
    current_guess = Beta_hat;
  }
  
  // Create the output matrix
  int n_rows = results.size();
  int vec_size = initial_guess.size();
  NumericMatrix output(n_rows, vec_size);
  
  for (int i = 0; i < n_rows; i++) {
    output(i, _) = results[i];
  }
  
  // Add row names as psi values
  CharacterVector row_names(n_rows);
  for (int i = 0; i < n_rows; i++) {
    row_names[i] = std::to_string(psi_values[i]);
  }
  rownames(output) = row_names;
  
  return output;
}
