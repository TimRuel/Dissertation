inject_shared_args <- function(cfg, shared_name = "shared_args") {
  
  shared <- cfg[[shared_name]]
  
  inject <- function(x) {
    if (is.list(x)) {
      if (!is.null(x$args)) {
        
        x$args <- modifyList(shared, x$args)
      }
      x <- lapply(x, inject)
    }
    x
  }
  
  cfg[names(cfg) != shared_name] <- lapply(cfg[names(cfg) != shared_name], inject)
  
  cfg$shared_args <- NULL
  
  return(cfg)
}

load_config <- function(population) {
  
  base_cfg <- config::get(file = "config/base.yml")
  
  pop_cfg <- population |> 
    sprintf("config/population_%s.yml", ... = _) |> 
    config::get(file = _) |> 
    inject_shared_args()
  
  cfg <- modifyList(base_cfg, pop_cfg)
  
  cfg$X1_specs$X1_levels$base_X2 <- NULL
  
  pop_result_base <- file.path(cfg$output_dir, sprintf("population_%s", population))
  dir.create(pop_result_base, recursive = TRUE, showWarnings = FALSE)
  
  existing_runs <- list.dirs(pop_result_base, full.names = FALSE, recursive = FALSE)
  existing_runs <- existing_runs[grepl("^run_\\d{3}$", existing_runs)]
  
  if (length(existing_runs) == 0) {
    
    next_run_num <- 1
  } else {
    
    run_nums <- as.integer(sub("run_", "", existing_runs))
    next_run_num <- max(run_nums) + 1
  }
  
  run_id <- sprintf("run_%03d", next_run_num)
  result_dir <- file.path(pop_result_base, run_id)
  
  cfg$run_id <- run_id
  cfg$result_dir <- result_dir
  
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  
  return(cfg)
}

