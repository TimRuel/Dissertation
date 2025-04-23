

TASK_ID <- Sys.getenv("TASK_ID")
MC_CORES <- Sys.getenv("MC_CORES")

num_cores <- parallel::detectCores()

saveRDS(list(MC_CORES = MC_CORES, num_cores = num_cores), paste0("test", TASK_ID, ".rda"))

test0 <- readRDS("test0.rda")

