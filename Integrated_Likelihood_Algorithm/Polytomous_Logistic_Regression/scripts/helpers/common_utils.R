softmax <- function(x) exp(x) / sum(exp(x))

softmax_adj <- function(x) exp(x) / (1 + sum(exp(x)))

entropy <- function(p) -sum(p * log(p), na.rm = TRUE)

`%||%` <- function(a, b) if (!is.null(a)) a else b

save_list_objects <- function(object_list, dir_path) {
  
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  for (object in names(object_list)) saveRDS(object_list[[object]], file = file.path(dir_path, paste0(object, ".rds")))
}

save_list_plots <- function(plots_list, dir_path) {
  
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  for (plot in names(plots_list)) {
    
    ggsave(filename = here(dir_path, paste0(plot, ".png")), 
           plot = plots_list[[plot]],
           width = 20,
           height = 9,
           dpi = "retina",           
           units = "in")
  }
}


