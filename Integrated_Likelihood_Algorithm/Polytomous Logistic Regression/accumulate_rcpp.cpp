#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List accumulate_rcpp(NumericVector x, Function f, NumericVector init) {
  int n = x.size();
  
  // Initialize the result list
  List result(n);
  
  // Start with the initial state
  NumericVector current = init;
  
  for (int i = 0; i < n; i++) {
    // Apply the function directly to 'current' and 'x[i]'
    current = as<NumericVector>(f(current, x[i]));
    
    // Store the result in the list
    result[i] = NumericVector(current.begin(), current.end());
  }
  
  return result;
}