// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// iterate_while
NumericMatrix iterate_while(double psi_max, double branch_max, NumericVector initial_guess, Function get_Beta_hat, Function get_log_likelihood, NumericMatrix omega_hat, NumericMatrix X_one_hot, NumericVector X_h_one_hot, NumericMatrix Y_one_hot, double crit, double step_size, double noise_sd, bool verbose);
RcppExport SEXP _iterate_iterate_while(SEXP psi_maxSEXP, SEXP branch_maxSEXP, SEXP initial_guessSEXP, SEXP get_Beta_hatSEXP, SEXP get_log_likelihoodSEXP, SEXP omega_hatSEXP, SEXP X_one_hotSEXP, SEXP X_h_one_hotSEXP, SEXP Y_one_hotSEXP, SEXP critSEXP, SEXP step_sizeSEXP, SEXP noise_sdSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type psi_max(psi_maxSEXP);
    Rcpp::traits::input_parameter< double >::type branch_max(branch_maxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type initial_guess(initial_guessSEXP);
    Rcpp::traits::input_parameter< Function >::type get_Beta_hat(get_Beta_hatSEXP);
    Rcpp::traits::input_parameter< Function >::type get_log_likelihood(get_log_likelihoodSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type omega_hat(omega_hatSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X_one_hot(X_one_hotSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type X_h_one_hot(X_h_one_hotSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Y_one_hot(Y_one_hotSEXP);
    Rcpp::traits::input_parameter< double >::type crit(critSEXP);
    Rcpp::traits::input_parameter< double >::type step_size(step_sizeSEXP);
    Rcpp::traits::input_parameter< double >::type noise_sd(noise_sdSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(iterate_while(psi_max, branch_max, initial_guess, get_Beta_hat, get_log_likelihood, omega_hat, X_one_hot, X_h_one_hot, Y_one_hot, crit, step_size, noise_sd, verbose));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_iterate_iterate_while", (DL_FUNC) &_iterate_iterate_while, 13},
    {NULL, NULL, 0}
};

RcppExport void R_init_iterate(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}