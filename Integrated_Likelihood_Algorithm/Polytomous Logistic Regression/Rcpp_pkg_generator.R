
setwd("C:/Northwestern/Dissertation/Integrated_Likelihood_Algorithm/Polytomous Logistic Regression")

Rcpp::Rcpp.package.skeleton(name = "iterate", cpp_files = c("iterate_while_1.cpp", "iterate_while_2.cpp"), example_code = FALSE) 

install.packages("iterate", repos = NULL, type = "source")



setwd("C:/Northwestern/Dissertation/Integrated_Likelihood_Algorithm/Polytomous Logistic Regression")

RcppArmadillo::RcppArmadillo.package.skeleton(name = "PolytomousUtils", code_files = "polytomous_utils.cpp", example_code = FALSE) 

install.packages("PolytomousUtils", repos = NULL, type = "source")


