
setwd("C:/Northwestern/Dissertation/Integrated_Likelihood_Algorithm/Polytomous Logistic Regression")

# setwd("/home/tbr0780/Dissertation/Integrated_Likelihood_Algorithm/Polytomous Logistic Regression")

Rcpp::Rcpp.package.skeleton(name = "PolytomousUtils", cpp_files = "polytomous_utils.cpp", example_code = FALSE) 

install.packages("PolytomousUtils", repos = NULL, type = "source")

