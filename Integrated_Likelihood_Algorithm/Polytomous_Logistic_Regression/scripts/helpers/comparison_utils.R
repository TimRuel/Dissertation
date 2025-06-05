
get_LL_df_long <- function(LL_df) {
  
  LL_df |>
    pivot_longer(cols = -psi,
                 names_to = "pseudolikelihood",
                 values_to = "value") |>
    mutate(pseudolikelihood = pseudolikelihood |>
             as_factor())
}

get_spline_models <- function(LL_df_long) {
  
  LL_df_long |>
    drop_na(value) |>
    group_by(pseudolikelihood) |>
    group_map(~ smooth.spline(.x$psi, .x$value, spar = 0.3)) |>
    set_names(levels(LL_df_long$pseudolikelihood))
}

get_MLE_data <- function(spline_models) {
  
  spline_models |>
    imap(
      \(mod, pseudolikelihood) {
        optimize(\(psi) predict(mod, psi)$y,
                 lower = LL_df_long |>
                   filter(pseudolikelihood == pseudolikelihood) |> 
                   select(psi) |>
                   min(),
                 upper = LL_df_long |>
                   filter(pseudolikelihood == pseudolikelihood) |> 
                   select(psi) |>
                   max(),
                 maximum = TRUE) |>
          data.frame() |>
          mutate(MLE = as.numeric(maximum),
                 Maximum = as.numeric(objective)) |>
          select(MLE, Maximum)
      }) |>
    do.call(rbind, args = _) |>
    mutate(MLE_label = c("hat(psi)[IL]", "hat(psi)[PL]")) |>
    rownames_to_column("Source")
}

get_pseudolikelihoods <- function(spline_models, MLE_data) {
  
  spline_models |> 
    map2(MLE_data$Maximum,
         \(mod, maximum) \(psi) predict(mod, psi)$y - maximum)
}

get_confidence_intervals <- function(pseudolikelihoods, 
                                     LL_df_long,
                                     MLE_data, 
                                     alpha, 
                                     J) {
  
  crit <- qchisq(1 - alpha, 1) / 2
  
  list(pseudolikelihoods,
       names(pseudolikelihoods),
       MLE_data$MLE) |> 
    pmap(
      
      \(pseudolikelihood, name, MLE) {
        
        psi_endpoints <- LL_df_long |> 
          filter(pseudolikelihood == name) |> 
          drop_na() |> 
          pull(psi) |> 
          (\(x) c(head(x, 1), tail(x, 1)))()
        
        lower_bound <- tryCatch(
          uniroot(\(psi) pseudolikelihood(psi) + crit,
                  interval = c(psi_endpoints[1], MLE))$root,
          error = function(e) return(0)
          ) |> 
          round(3)
        
        upper_bound <- tryCatch( 
          uniroot(\(psi) pseudolikelihood(psi) + crit,
                  interval = c(MLE, psi_endpoints[2]))$root,
          error = function(e) return(log(J))
          ) |>
          round(3)
        
        return(c(lower_bound, upper_bound))
        }
    )
}

get_stat_fns <- function(pseudolikelihoods, LL_df) {
  
  psi_endpoints <- LL_df |> 
    pull(psi) |> 
    (\(x) c(head(x, 1), tail(x, 1)))()
  
  pseudolikelihoods |> 
    imap(
      \(pseudolikelihood, name) {
        
        stat_function(fun = pseudolikelihood,
                      geom = "line",
                      aes(color = name),
                      linewidth = 1.5,
                      show.legend = FALSE,
                      xlim = psi_endpoints)
        }
    )
}
