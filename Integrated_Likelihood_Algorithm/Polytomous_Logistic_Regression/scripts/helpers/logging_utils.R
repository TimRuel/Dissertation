# scripts/helpers/logging_utils.R
save_run_metadata <- function(metadata, output_dir) {
  write_yaml(metadata, file.path(output_dir, "metadata.yml"))
}

append_run_log <- function(metadata, log_path) {
  log_entry <- tibble(
    run_id = metadata$run_id,
    mode = metadata$mode,
    iteration_id = metadata$iteration_id,
    slurm_array_id = metadata$slurm_array_id,
    start_time = metadata$timestamp_start,
    end_time = metadata$timestamp_end,
    elapsed_seconds = metadata$elapsed_seconds,
    git_commit = metadata$git_commit
  )
  
  if (fs::file_exists(log_path)) {
    readr::read_csv(log_path, show_col_types = FALSE) %>%
      bind_rows(log_entry) %>%
      readr::write_csv(log_path)
  } else {
    readr::write_csv(log_entry, log_path)
  }
}
