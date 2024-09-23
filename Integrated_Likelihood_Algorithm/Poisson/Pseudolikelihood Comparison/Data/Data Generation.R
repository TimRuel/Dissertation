
print("Choose your population parameters R script.")
population_params_file_path <- file.choose() 

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source(population_params_file_path)
source("../../utils.R")

set.seed(seed)
         
theta_0 <- d |> 
 c(theta_0_range) |> 
 as.list() |> 
 do.call(runif, args = _) |> 
 (\(x) if (theta_0_sum) x / sum(x) * theta_0_sum else x)() 

data <- n |> 
 map2(
   theta_0,
   \(x, y) rpois(x, y) |> 
     as.numeric())

weights <- d |>
 c(weights_range) |>
 as.list() |>
 do.call(runif, args = _) |>
 (\(x) x / mean(x) * weights_mean)()



