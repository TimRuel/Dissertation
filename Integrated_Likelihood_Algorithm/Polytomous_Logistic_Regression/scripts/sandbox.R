library(tidyverse)
library(doFuture)
library(here)
library(yaml)
library(fs)
library(Rcpp)
library(nloptr)
library(PolytomousUtils)

setwd("C:/Northwestern/Dissertation/Integrated_Likelihood_Algorithm/Polytomous_Logistic_Regression/experiments/exp_v1.0.2/simulations/sim_04/iter_0005")
config <- yaml::read_yaml("config_snapshot.yml")
model_df <- readRDS("data/model_df.rds")
X_design <- readRDS("data/X_design.rds")

miceadds::source.all("../../../../../scripts/helpers/", print.source = FALSE)

list2env(config$model_specs, environment())
list2env(config$optimization_specs, environment())

X1_levels <- config$X1_levels
formula <- as.formula(formula)
ml_model <- fit_multinomial_logistic_model(model_df, formula)
Y_design <- get_Y_design(model_df)
Beta_MLE <- get_Beta_MLE(ml_model)
threshold <- get_threshold(Beta_MLE, X_design, Y_design, IL$threshold_offset)
h <- get_X1_level_of_interest(X1_levels)
X_h_design <- get_X_h_design(X_design, X1_levels)
psi_hat <- get_psi_hat(ml_model, X1_levels)
n_h <- nrow(X_h_design)
psi_endpoints_PL <- get_psi_endpoints(psi_hat, Beta_MLE, X_h_design, PL$num_std_errors, J, n_h)
psi_grid <- get_psi_grid(psi_endpoints_PL, PL$step_size, J)
Jm1 <- J - 1
init_guess_sd <- 1
maxtime <- -1
lambda <- 0.1
max_retries <- 10

omega_hat_eq_con_fn <- function(Beta) omega_hat_eq_con_fn_rcpp(Beta, X_h_design, Jm1, p, psi_hat)

omega_hat_ineq_con_fn <- function(Beta) omega_hat_ineq_con_fn_rcpp(Beta, X_design, Y_design, Jm1, p, n, threshold)

omega_hat <- get_omega_hat(omega_hat_eq_con_fn, omega_hat_ineq_con_fn, Jm1, p, 10)

make_omega_hat_branch_fn <- function(omega_hat,
                                     X_design,
                                     Y_design,
                                     X_h_design,
                                     Jm1,
                                     p,
                                     n) {
  
  Beta_hat_obj_fn <- function(Beta) Beta_hat_obj_fn_rcpp(Beta, X_design, omega_hat, Jm1, p, n)
  
  get_Beta_hat <- function(con_fn, init_guess) get_Beta_hat_template(Beta_hat_obj_fn, con_fn, init_guess)
  
  function(psi) {
    
    con_fn <- function(Beta) Beta_hat_con_fn_rcpp(Beta, X_h_design, psi, Jm1, p)
    
    Beta_hat <- get_Beta_hat(con_fn, c(omega_hat))
    
    return(log_likelihood_rcpp(Beta_hat, X_design, Y_design, Jm1, p, n))
  }
}

omega_hat_branch_fn <- make_omega_hat_branch_fn(omega_hat,
                                                X_design,
                                                Y_design,
                                                X_h_design,
                                                Jm1,
                                                p,
                                                n)

optimize_omega_hat_branch <- function(omega_hat_branch_fn, J) {
  
  optimize(omega_hat_branch_fn,
           interval = c(0, log(J)),
           maximum = TRUE,
           tol = 0.1)$maximum
}

omega_hat_branch_mode <- optimize_omega_hat_branch(omega_hat_branch_fn, J)

num_std_errors <- 3

psi_endpoints <- get_psi_endpoints(psi_hat, Beta_MLE, X_h_design, num_std_errors, J, n_h)

get_psi_grid(psi_endpoints, 0.05, J)

# Example matrix: rows are curves, columns are x-values
mat <- matrix(rnorm(50), nrow = 5)  # 5 rows (curves), 10 columns (x-points)
x_vals <- seq(0, 1, length.out = ncol(mat))

# Plot each row as a separate line
matplot(x_vals, t(mat), type = "l", lty = 1, col = 1:nrow(mat),
        xlab = "X", ylab = "Y", main = "Multiple Lines from Matrix Rows")
legend("topright", legend = paste("Curve", 1:nrow(mat)), col = 1:nrow(mat), lty = 1)

