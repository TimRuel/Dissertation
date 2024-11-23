# library(future)
# library(doFuture)
# library(zeallot)
# library(purrr)
# library(dplyr)
# library(plyr)
# library(stringr)
# library(progressr)
# library(tictoc)
library(tidyverse)
library(purrr)

# handlers(global = TRUE)
# handlers("cli")

# num_cores <- Sys.getenv("SLURM_NPROCS") |>
#   as.numeric()

# num_cores <- availableCores() |>
#   as.numeric()

# num_cores <- parallel::detectCores() |>
#   as.numeric()

print("Choose your R script that generates your sample data.")
data_file_path <- file.choose() 

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source(data_file_path)
source("utils.R")

set.seed(seed)

