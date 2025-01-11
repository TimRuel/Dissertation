#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix iterate_over(
    NumericVector iteration_values, 
    NumericVector initial_guess, 
    Function custom_function, 
    List args, 
    bool fill_bottom_to_top = false
) {
  int n = iteration_values.size(); 
  int vec_size = initial_guess.size();
  NumericMatrix results(n, vec_size);  // Initialize matrix with n rows.
  
  NumericVector current_value = initial_guess;  // Start with initial guess.
  
  // Fill results based on the order
  for (int i = 0; i < n; i++) {
    // Perform the calculation
    current_value = as<NumericVector>(
      custom_function(current_value, iteration_values[i], args)
    );
    
    // Determine where to store the result
    int row_index = fill_bottom_to_top ? (n - 1 - i) : i;
    results(row_index, _) = current_value;
  }
  
  return results;
}
