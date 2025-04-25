library(foreach)
library(doFuture)
library(future)

TASK_ID <- Sys.getenv("TASK_ID")
num_cores <- Sys.getenv("MC_CORES") |> 
  as.integer()

plan(multisession, workers = I(num_cores)) 

results <- foreach(i = 1:num_cores, 
                   .combine = c) %dofuture% 
  {
  paste("Running on worker", i, "- PID:", Sys.getpid())
    }

# Print results
message <- paste(results, collapse = "\n")

saveRDS(list(message = message,
             num_cores = num_cores), paste0("out_", TASK_ID, ".rda"))

