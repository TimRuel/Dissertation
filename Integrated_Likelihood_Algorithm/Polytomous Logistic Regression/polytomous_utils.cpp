// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double log_likelihood_rcpp(NumericVector Beta, NumericMatrix X_one_hot, NumericMatrix Y_one_hot, int Jm1, int p) {
  mat Beta_mat(Beta.begin(), p, Jm1, false);
  mat Y_hat = as<mat>(X_one_hot) * Beta_mat;
  double log_likelihood = accu(sum(as<mat>(Y_one_hot) % Y_hat, 1) - log(1 + sum(exp(Y_hat), 1)));
  return log_likelihood;
}

// [[Rcpp::export]]
mat softmax_adj_mat_rcpp(const mat& x) {
  mat exp_x = exp(x); 
  mat row_sums = sum(exp_x, 1); 
  mat row_sums_expanded = repmat(row_sums, 1, exp_x.n_cols); 
  return exp_x / (1 + row_sums_expanded); 
}

// [[Rcpp::export]]
rowvec softmax_adj_rowvec_rcpp(const rowvec& x) {
  rowvec exp_x = exp(x); // Exponentiate each element
  double row_sum = sum(exp_x); // Sum of exponentials in the row
  return exp_x / (1 + row_sum); // Normalize by row_sum + 1
}

// [[Rcpp::export]]
mat softmax_mat_rcpp(const mat& x) {
  mat exp_x = exp(x);
  mat row_sums = sum(exp_x, 1);
  mat row_sums_expanded = repmat(row_sums, 1, exp_x.n_cols);
  return exp_x / row_sums_expanded;
}

// [[Rcpp::export]]
rowvec softmax_rowvec_rcpp(const rowvec& x) {
  rowvec exp_x = exp(x); // Exponentiate each element
  double row_sum = sum(exp_x); // Sum of exponentials in the row
  return exp_x / row_sum; // Normalize by row_sum + 1
}

// [[Rcpp::export]]
double entropy_rcpp(const rowvec& p) {
  rowvec log_p = log(p);
  log_p.replace(-datum::inf, 0);
  return -accu(p % log_p);
}

// [[Rcpp::export]]
double PoI_fn_rcpp(const mat& Beta, const mat& X_h_one_hot) {
  mat Beta_extended = join_horiz(zeros<mat>(Beta.n_rows, 1), Beta);
  mat exp_X_h_Beta = X_h_one_hot * Beta_extended;
  rowvec exp_X_h_Beta_row = exp_X_h_Beta.row(0);
  rowvec softmax_probs = softmax_rowvec_rcpp(exp_X_h_Beta);
  double entropy = entropy_rcpp(softmax_probs);
  return entropy;
}

// [[Rcpp::export]]
double Beta_hat_obj_fn_rcpp(NumericVector Beta, NumericMatrix X_one_hot, NumericMatrix omega_hat, double lambda) {
  int p = omega_hat.nrow(), Jm1 = omega_hat.ncol();
  mat Beta_mat(Beta.begin(), p, Jm1, false);
  mat Y_hat = as<mat>(X_one_hot) * Beta_mat;
  mat probs = softmax_adj_mat_rcpp(as<mat>(X_one_hot) * as<mat>(omega_hat));
  double obj = -accu(sum(probs % Y_hat, 1) - log(1 + sum(exp(Y_hat), 1))) + lambda * accu(Beta_mat % Beta_mat);
  return obj;
}

// [[Rcpp::export]]
double Beta_hat_con_fn_rcpp(NumericVector Beta, NumericMatrix X_h_one_hot, double psi, int Jm1, int p) {
  mat Beta_mat(Beta.begin(), p, Jm1, false);
  return PoI_fn_rcpp(Beta_mat, as<mat>(X_h_one_hot)) - psi;
}

// [[Rcpp::export]]
double omega_hat_obj_fn_rcpp(NumericVector Beta, NumericMatrix X_h_one_hot, int Jm1, int p, double psi_hat) {
  mat Beta_mat(Beta.begin(), p, Jm1, false);
  double entropy = PoI_fn_rcpp(Beta_mat, as<mat>(X_h_one_hot));
  return std::abs(entropy - psi_hat);
}

// [[Rcpp::export]]
double omega_hat_con_fn_rcpp(NumericVector Beta, NumericMatrix X_one_hot, NumericMatrix Y_one_hot, int Jm1, int p, double threshold) {
  return -log_likelihood_rcpp(Beta, X_one_hot, Y_one_hot, Jm1, p) - threshold;
}











