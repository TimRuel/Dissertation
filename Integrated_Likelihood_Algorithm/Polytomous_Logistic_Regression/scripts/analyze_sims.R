# scripts/analyze_sims.R

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
})

proj_subdir <- here("Integrated_Likelihood_Algorithm", "Polytomous_Logistic_Regression")
proj_path <- function(...) here(proj_subdir, ...)

# Load helpers
miceadds::source.all(proj_path("scripts", "helpers"), print.source = FALSE)

