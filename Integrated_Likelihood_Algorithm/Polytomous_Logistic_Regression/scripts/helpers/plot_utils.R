plot_base <- function() {
  
  ggplot() +   
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "#444444", color = NA),  # Dark panel
      plot.background = element_rect(fill = "#2A2A2A", color = NA),  # Dark background
      panel.grid.major = element_line(color = "#4C4C4C"),  # Darker grid
      panel.grid.minor = element_line(color = "#333333"),  
      axis.ticks = element_line(color = "white"),
      axis.text = element_text(color = "white"),
      axis.title = element_text(color = "white"),
      strip.text = element_text(color = "white"),
      plot.title = element_text(color = "white", face = "bold"),
      legend.background = element_rect(fill = "#444444"),
      legend.text = element_text(color = "white"),
      legend.title = element_text(color = "white"))
}

get_entropy_df <- function(Y_probs) {
  
  Y_probs |> 
    mutate(X2 = plyr::round_any(X2, 0.01)) |> 
    group_by(X1) |> 
    distinct(X2, .keep_all = TRUE) |> 
    arrange(X2) |> 
    rowwise() |>
    mutate(entropy = entropy(c_across(-X2))) |> 
    select(X1, X2, entropy) |> 
    data.frame()
}

plot_entropy <- function(X1_levels, Y_probs, title) {
  
  entropy_ranges <- map(X1_levels, \(level) level$X2$entropy_range)
  
  df <- get_entropy_df(Y_probs)
  
  plot_base() + 
    geom_line(aes(x = X2, y = entropy, color = X1),
              linewidth = 1,
              data = df) +
    geom_hline(yintercept = unlist(entropy_ranges), color = "black", linewidth = 0.6, linetype = "dashed") +
    labs(title = title, x = "X2", y = "Entropy", color = "X1 Level")
}

plot_Y_probs <- function(Y_probs) {
  
  Y_probs_long <- Y_probs |> 
    pivot_longer(cols = starts_with("Y"),
                 names_to = "Y_class",
                 names_prefix = "Y",
                 values_to = "Probability")
  
  plot_base() + 
    geom_point(aes(x = X2, y = Probability, color = Y_class),
               size = 1,
               data = Y_probs_long) +
    facet_wrap(~X1, 
               ncol = 1,
               labeller = labeller(X1 = function(x) paste("X1 =", x)),
               scales = "free_x") +
    scale_x_continuous(breaks = seq(0, 4, 0.5), expand = c(0.05, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0.1, 0)) +
    labs(title = "True Probabilities of Y Classes by X1 and X2", 
         x = "X2", 
         y = "Probability", 
         color = "Class of Y")
}

plot_data <- function(df) {
  
  plot_base() + 
    geom_histogram(aes(x = X2, fill = Y),
                   color = "black",
                   binwidth = 0.1, 
                   position = "stack",
                   data = df) +
    facet_wrap(~X1, 
               ncol = 1,
               labeller = labeller(X1 = function(x) paste("X1 =", x)),
               scales = "free_x") +
    scale_x_continuous(breaks = seq(-2, 4, 0.5)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(title = "Observed Y Values by X1 and X2", x = "X2", y = "Count", fill = "Observed Y")
}

get_theoretical_entropy_plot <- function(X1_levels, pY_0) {
  
  title <- "Theoretical Entropy of Probabilities of Y Classes by X1 and X2"

  entropy_plot <- plot_entropy(X1_levels, 
                               pY_0, 
                               title)
  
  return(entropy_plot)
}

get_observed_plots <- function(X1_levels, Y_probs, df) {
  
  title <- "Observed Entropy of Probabilities of Y Classes by X1 and X2"
  
  observed_entropy_plot <- plot_entropy(X1_levels, 
                                        Y_probs, 
                                        title)
  
  Y_probs_plot <- plot_Y_probs(Y_probs)
  
  data_plot <- plot_data(df)
  
  plots <- list(observed_entropy_plot = observed_entropy_plot,
                Y_probs_plot = Y_probs_plot,
                data_plot = data_plot)
  
  return(plots)
}

get_LL_plot <- function(df) {
  
  df |>
    ggplot(aes(x = psi, y = .data[[names(df)[2]]])) +
    geom_point(color = "cyan", size = 3, alpha = 0.7) +
    theme_minimal(base_size = 15) +  # Minimal theme with a larger base font size
    theme(
      plot.background = element_rect(fill = "#2E2E2E", color = NA),  # Dark background for the whole plot
      panel.background = element_rect(fill = "#3A3A3F", color = "#1A1A1A", linewidth = 2),  # Lighter panel with a border
      axis.text = element_text(color = "white"),  # White axis labels
      axis.title = element_text(color = "white"),  # White axis titles
      plot.title = element_text(color = "white", size = 18, face = "bold"),  # White title
      plot.caption = element_text(color = "gray", size = 10),  # Gray caption
      panel.grid = element_line(color = "gray30", linetype = "dashed")  # Subtle grid lines
    ) +
    labs(
      x = "\u03C8",
      y = paste(names(df)[[2]], "Log-Likelihood")
    )
}

get_branches_plot <- function(mat) {

  df <- as.data.frame(mat)
  df$CurveID <- paste0("Curve_", seq_len(nrow(df)))
  
  df_long <- df |> 
    pivot_longer(-CurveID, names_to = "X", values_to = "Y") |> 
    mutate(X = as.numeric(X))
  
  plot_base() + 
    theme(legend.position = "none") +
    geom_line(data = df_long, aes(x = X, y = Y, group = CurveID, color = CurveID), linewidth = 1) +
    labs(title = "Integrated Log-Likelihood Branches", x = "\u03C8", y = "Log-Likelihood")
}

get_LL_comparison_table <- function(MLE_data, 
                                    conf_ints, 
                                    psi_0,
                                    alpha) {
  
  MLE_data |> 
    select(Source, MLE) |> 
    mutate(
      MLE = sprintf("%.3f", MLE),
      conf_int = conf_ints |> 
        map_chr(\(x) sprintf("(%.3f, %.3f)", x[1], x[2])),
      length = conf_ints |> 
        map_dbl(diff) |> 
        (\(x) sprintf("%.3f", x))()
    ) |> 
    arrange(as.numeric(length)) |> 
    add_row(Source = "Truth",
            MLE = sprintf("%.3f", psi_0),
            conf_int = "--",
            length = "--") |> 
    kbl(col.names = c("Source", 
                      "MLE",
                      paste(sprintf("%1.0f%%", 100*(1-alpha)), "Confidence Interval"),
                      "CI Length"),   
        align = "c") |> 
    kable_material_dark(full_width = FALSE) |> 
    row_spec(3, color = "green", bold = TRUE) |> 
    column_spec(c(3, 4), color = "white")
}

get_LL_comparison_plot <- function(stat_fns, 
                                   LL_df_long,
                                   MLE_data,
                                   alpha) {
  
  c(stat_fn_IL, stat_fn_PL) %<-% stat_fns
  
  crit <- qchisq(1 - alpha, 1) / 2
  
  psi_endpoints <- LL_df_long |> 
    pull(psi) |> 
    (\(x) c(head(x, 1), tail(x, 1)))()  
  
  y_min <- -crit - 0.5
  
  label_data <- MLE_data |>
    add_row(Source = "Truth",
            MLE = psi_0,
            MLE_label = "psi[0]")
  
  plot_base() +
    stat_fn_IL +
    stat_fn_PL +
    geom_hline(yintercept = 0,
               linetype = 5) +
    geom_hline(aes(yintercept = -crit),
               linewidth = 1.5) +
    geom_text(aes(1.5, -crit, label = paste(sprintf("%1.0f%%", 100*(1-alpha)), "CI Line"), vjust = 1.3)) +
    geom_vline(aes(xintercept = MLE,
                   color = Source),
               data = label_data,
               show.legend = FALSE) +
    geom_label_repel(aes(x = MLE,
                         y = y_min / 2,
                         label = MLE_label,
                         color = Source),
                     data = label_data,
                     direction = "y",
                     parse = TRUE,
                     show.legend = FALSE,
                     seed = 7835) +
    geom_point(aes(x = psi, y = value - max(value, na.rm = TRUE)),
               data = LL_df_long |> 
                 filter(pseudolikelihood == "Integrated"),
               size = 1) +
    geom_point(aes(x = psi, y = value - max(value, na.rm = TRUE)),
               data = LL_df_long |> 
                 filter(pseudolikelihood == "Profile"),
               size = 1) +
    ylab("Log-Likelihood") +
    scale_x_continuous(expand = c(0, 0),
                       limits = psi_endpoints) +
    scale_y_continuous(expand = c(0, 0),
                       limits = c(y_min, 0.1)) +
    scale_color_brewer(palette = "Set1") +
    xlab(expression(psi))
}

