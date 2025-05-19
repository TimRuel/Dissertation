with_file_lock <- function(lock_path, timeout = 60, stale_after = 5 * 60, verbose = TRUE, code) {
  start_time <- Sys.time()
  
  # Create lock directory if needed
  dir_create(path_dir(lock_path))
  
  repeat {
    if (!file_exists(lock_path)) {
      # Create lock file
      write_lines(as.character(Sys.getpid()), lock_path)
      if (verbose) message("[Lock] Acquired lock: ", lock_path)
      break
    } else {
      # Check if lock is stale
      lock_age <- difftime(Sys.time(), file_info(lock_path)$modification_time, units = "secs")
      if (lock_age > stale_after) {
        if (verbose) message("[Lock] Stale lock detected. Removing: ", lock_path)
        file_delete(lock_path)
        next
      }
      
      if (verbose) message("[Lock] Waiting for lock to release: ", lock_path)
      Sys.sleep(0.5)
      
      if (difftime(Sys.time(), start_time, units = "secs") > timeout) {
        stop("[Lock] Timeout reached while waiting for lock at: ", lock_path)
      }
    }
  }
  
  on.exit({
    if (file_exists(lock_path)) {
      file_delete(lock_path)
      if (verbose) message("[Lock] Released lock: ", lock_path)
    }
  })
  
  force(code)
}
