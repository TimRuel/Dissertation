
print("Choose your population parameters R script.")
population_params_file_path <- file.choose() 

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source(population_params_file_path)
source("../../utils2.R")

set.seed(seed)
         
theta_0 <- n |> 
 c(theta_0_range) |> 
 as.list() |> 
 do.call(runif, args = _) |> 
 (\(x) if (theta_0_sum) x / sum(x) * theta_0_sum else x)() 

data <- n |> 
  rpois(theta_0) |> 
  as.numeric()

weights <- n |>
 c(weights_range) |>
 as.list() |>
 do.call(runif, args = _) |>
 (\(x) x / mean(x) * weights_mean)()



