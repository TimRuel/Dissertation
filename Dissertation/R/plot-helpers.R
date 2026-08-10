# =========================================================================
# plot-helpers.R
#
# Shared plotting vocabulary for the dissertation's figures.
#
# This is the layer that IS common across experiments. Nothing here knows
# which experiment it is drawing — the functions take a tidy data frame in
# the column vocabulary produced by R/load-results.R (estimand, model,
# site, N, method, conf_level, estimate, lower, upper) and return a plot.
#
# Per-experiment figure scripts live in Dissertation/figures/ and compose
# these however that experiment needs. Add a function here only when a
# second figure script would use it; one-off plotting code belongs in the
# figure script that needs it.
# =========================================================================

library(ggplot2)
library(dplyr)
library(forcats)
library(scales)

# Output directory for figures used in the dissertation
FIG_DIR <- "Dissertation/files/Images"

# -------------------------------------------------------------------------
# Write a figure, refusing to write an empty one
#
# PL/IL results come from real experiments, so a figure whose experiment
# has not been run or downloaded yet must be SKIPPED — leaving any
# previously committed PNG intact — rather than overwritten with an
# axis-only panel. Every skip is reported so a missing figure is never
# silent.
# -------------------------------------------------------------------------
save_figure <- function(
  filename,
  df,
  build,
  width = 6.5,
  height = 8,
  dpi = 300,
  dir = FIG_DIR
) {
  if (nrow(df) == 0L) {
    message("⏭  Skipping ", filename, " — no data for this combination")
    return(invisible(NULL))
  }

  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }

  p <- build(df)

  ggsave(file.path(dir, filename), p, width = width, height = height,
    dpi = dpi)
  message("✔ Wrote ", filename)

  invisible(p)
}

# -------------------------------------------------------------------------
# A comparison figure needs both models present; one model alone would
# render a single line under a title promising a comparison.
# -------------------------------------------------------------------------
require_both_models <- function(df) {
  if (n_distinct(df$model) >= 2L) df else df[0, ]
}

# -------------------------------------------------------------------------
# Full comparison: all methods, faceted by confidence level. Sites ordered
# by N (ascending) so sparsity-driven differences (e.g. vs. Chao) are
# visually apparent left-to-right rather than requiring cross-reference
# to a table.
# -------------------------------------------------------------------------
plot_interval_comparison <- function(df, y_label, method_colors = NULL) {
  df <- df |>
    mutate(site_ord = fct_reorder(factor(site), N))

  p <- df |>
    ggplot(aes(x = site_ord, y = estimate, color = method)) +
    geom_pointrange(
      aes(ymin = lower, ymax = upper),
      position = position_dodge(width = 0.7),
      size = 0.35,
      linewidth = 0.6
    ) +
    facet_wrap(~ scales::percent(conf_level), ncol = 1) +
    labs(x = "Site (ordered by increasing N)", y = y_label, color = "Method") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())

  if (!is.null(method_colors)) {
    p <- p + scale_color_manual(values = method_colors)
  }
  p
}

# -------------------------------------------------------------------------
# Compact supplement: interval width only, easier to compare precision
# -------------------------------------------------------------------------
plot_interval_width <- function(df, y_label, group_var = "method") {
  df <- df |>
    mutate(
      width = upper - lower,
      site_ord = fct_reorder(factor(site), N),
      group_val = .data[[group_var]]
    )

  df |>
    ggplot(aes(x = site_ord, y = width, color = group_val, group = group_val)) +
    geom_line(aes(group = group_val), alpha = 0.4) +
    geom_point(size = 1.6) +
    facet_wrap(~ scales::percent(conf_level), ncol = 1, scales = "free_y") +
    labs(
      x = "Site (ordered by increasing N)",
      y = paste(y_label, "interval width"),
      color = if (group_var == "model") "Model" else "Method"
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())
}

# -------------------------------------------------------------------------
# Pseudolikelihood curve overlay for one site
#
# Uses the psi/rel_loglik grids exported by likelyr-simulations'
# analyze_app.R (via as_curve_df()). rel_loglik puts the profile and
# integrated curves on a common vertical scale with each maximum at 0,
# which is what makes their shapes comparable; psi_mle is drawn as a
# reference line.
# -------------------------------------------------------------------------
plot_curve_overlay <- function(df, x_label, ncol = 4, trim = TRUE) {
  if (isTRUE(trim) && "above_crit" %in% names(df)) {
    df <- df |> filter(above_crit)
  }

  df |>
    mutate(site_label = paste0("Site ", site, "  (N = ", N, ")")) |>
    ggplot(aes(x = psi, y = rel_loglik, color = method)) +
    geom_line(linewidth = 0.6) +
    geom_vline(
      aes(xintercept = psi_mle),
      linetype = "dashed",
      color = "grey50",
      linewidth = 0.4
    ) +
    facet_wrap(~site_label, ncol = ncol, scales = "free_x") +
    labs(
      x = x_label,
      y = "Relative log-pseudolikelihood",
      color = "Method"
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank())
}
