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
# Method palette — EMPHASIS, not categorical
#
# The argument these figures serve is about the integrated likelihood, with
# the profile likelihood as its direct rival and the Wald-type benchmarks as
# context. That is the "emphasis" form: one accent hue, one rival hue, and
# the reference methods in de-emphasis gray. Four equal-weight categorical
# hues would bury the series the chapter is actually about.
#
# Validated with the data-viz six checks (light surface, --pairs all): the
# accent pair passes every check, worst-case CVD ΔE 31.1 (protan) and 35.7
# (normal vision). A four-hue all-categorical alternative was tried and
# rejected — the best four-hue set still WARNs at ΔE 6.2 (protan) between
# its 3rd and 4th slots, which is legal only with secondary encoding.
#
# METHOD_SHAPES supplies that secondary encoding regardless, so the figures
# survive being printed in black and white — identity is never carried by
# color alone.
#
# Chao and the Tiffeau-Mayer unbiased estimator share a point estimator on
# this data (max abs difference 4.4e-4, correlation 0.99997) and differ only
# in variance estimate, so they are deliberately two steps of the same gray
# rather than two hues.
# -------------------------------------------------------------------------
METHOD_COLORS <- c(
  "IL" = "#1F6FEB",
  "PL" = "#D4761A",
  "Chao" = "#8B949E",
  "Unbiased (TM)" = "#57606A"
)

METHOD_SHAPES <- c(
  "IL" = 16, # filled circle  — the accent series
  "PL" = 17, # filled triangle
  "Chao" = 1, # open circle
  "Unbiased (TM)" = 4 # cross
)

# Legend order: the series the argument is about comes first, benchmarks last.
METHOD_ORDER <- c("IL", "PL", "Unbiased (TM)", "Chao")

order_methods <- function(x) {
  factor(x, levels = intersect(METHOD_ORDER, unique(x)))
}

# ggplot does not wrap captions, so a long one silently runs off the right
# edge of the saved PNG — invisible unless you open the file. Wrap at a
# character count matched to the figure's width instead of trusting eyeballs:
# ~14 chars per inch at the caption's rendered size.
wrap_caption <- function(x, fig_width = 6.5) {
  paste(strwrap(x, width = floor(fig_width * 14)), collapse = "\n")
}

# Subtitles render larger than captions, so they overflow sooner — a
# subtitle stating the figure's finding is exactly the string you least want
# truncated. Same idea, tighter budget.
wrap_subtitle <- function(x, fig_width = 6.5) {
  paste(strwrap(x, width = floor(fig_width * 11)), collapse = "\n")
}

# Recessive theme shared by the diagnostic figures: hairline solid grid, no
# minor grid, generous panel spacing.
theme_diagnostic <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, color = "grey88"),
      panel.spacing = unit(1, "lines"),
      plot.title = element_text(face = "plain", size = rel(1.05)),
      plot.subtitle = element_text(color = "grey30", size = rel(0.9)),
      plot.caption = element_text(color = "grey45", size = rel(0.78), hjust = 0),
      legend.position = "top",
      legend.title = element_blank()
    )
}

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
plot_interval_width <- function(df, y_label, group_var = "method",
                                method_colors = NULL) {
  df <- df |>
    mutate(
      width = upper - lower,
      site_ord = fct_reorder(factor(site), N),
      group_val = .data[[group_var]]
    )

  p <- df |>
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

  # Only meaningful when the grouping IS method — a model-grouped panel has
  # "No effects"/"Fixed effects" levels that METHOD_COLORS knows nothing
  # about, and scale_color_manual would drop every series to NA.
  if (!is.null(method_colors) && group_var == "method") {
    p <- p + scale_color_manual(values = method_colors)
  }
  p
}

# -------------------------------------------------------------------------
# Pseudolikelihood curve overlay for one site
#
# Uses the psi/rel_loglik grids exported by likelyr-simulations'
# analyze_app.R (via as_curve_df()). rel_loglik puts the profile and
# integrated curves on a common vertical scale with each maximum at 0,
# which is what makes their shapes comparable; psi_mle is drawn as a
# reference line.
#
# TRIMMING: do NOT use the bundle's `above_crit` column. analyze_app.R
# computes it for the integrated pseudolikelihood only and leaves it NA for
# every Profile row (verified in both ne_entropy/exp_v13 and
# logit_simpson/exp_v6). Since dplyr::filter() drops NA, `filter(above_crit)`
# silently deleted every PL curve while still drawing a "Method" legend —
# the figure looked like a method comparison but showed one method.
#
# Instead derive the cutoff locally from rel_loglik, which exists for both
# methods. That also keeps the trim IDENTICAL across methods; trimming one
# and not the other would distort the very shape comparison the figure is
# for. The default keeps everything inside the widest (99%) likelihood-ratio
# cutoff plus a margin, so no curve is cropped before its interval ends.
# -------------------------------------------------------------------------
plot_curve_overlay <- function(
  df,
  x_label,
  ncol = 4,
  trim = TRUE,
  crit_level = 0.99,
  margin = 1
) {
  if (isTRUE(trim)) {
    if (!"rel_loglik" %in% names(df)) {
      stop("plot_curve_overlay(trim = TRUE) needs a rel_loglik column.",
        call. = FALSE)
    }
    cutoff <- -(stats::qchisq(crit_level, df = 1) / 2 + margin)
    df <- df |> filter(rel_loglik >= cutoff)
  }

  # A method silently vanishing is the failure mode this function already had
  # once. Refuse to do it quietly a second time.
  present <- intersect(METHOD_ORDER, unique(as.character(df$method)))

  if (length(present) < 2L) {
    warning(
      "plot_curve_overlay(): only ", paste(present, collapse = ", "),
      " has curve data — this is a single-method figure, not a comparison.",
      call. = FALSE
    )
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
    scale_color_manual(values = METHOD_COLORS, drop = FALSE) +
    labs(
      x = x_label,
      y = "Relative log-pseudolikelihood",
      color = "Method"
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank())
}

# =========================================================================
# IL/PL DIAGNOSTIC PLOTS
#
# Promoted here from figures/ once a second estimand needed them. Each takes
# an already-assembled frame plus text, because WHICH bias-corrected
# estimator plays the reference role is estimand-specific (Tiffeau-Mayer for
# Simpson, Miller-Madow for Shannon) and belongs in the figure script.
# =========================================================================

# -------------------------------------------------------------------------
# Where the pseudolikelihoods put the point estimate
#
# pts: N, delta, series — series a factor whose third level is the predicted
# plug-in bias. Difference scale, because on the raw scale the sparse sites
# stretch the axis and compress everything else into a band.
# -------------------------------------------------------------------------
plot_bias_anatomy <- function(pts, y_expr, pred_label, title, subtitle,
                              caption) {
  lvls <- levels(pts$series)

  pal <- c(METHOD_COLORS[["IL"]], METHOD_COLORS[["PL"]], "#8B949E")
  shp <- c(16, 17, 5)
  szs <- c(2.1, 2.1, 3.4)
  names(pal) <- names(shp) <- names(szs) <- lvls

  # Direct-label the two data series only. At the largest N the prediction
  # coincides with PL, so a third label collides rather than clarifies.
  edge <- pts |>
    filter(series != pred_label) |>
    group_by(series) |>
    slice_max(N, n = 1, with_ties = FALSE) |>
    ungroup()

  ggplot(pts, aes(x = N, y = delta)) +
    geom_hline(yintercept = 0, linewidth = 0.5, color = "grey35") +
    geom_point(aes(color = series, shape = series, size = series)) +
    geom_text(
      data = edge,
      aes(label = series, color = series),
      hjust = 0,
      nudge_x = 0.9,
      size = 3.1,
      show.legend = FALSE
    ) +
    scale_color_manual(values = pal, breaks = lvls) +
    scale_shape_manual(values = shp, breaks = lvls) +
    scale_size_manual(values = szs, breaks = lvls) +
    scale_x_continuous(expand = expansion(mult = c(0.04, 0.12))) +
    labs(
      x = "Sample size N",
      y = y_expr,
      title = title,
      subtitle = subtitle,
      caption = caption
    ) +
    theme_diagnostic()
}

# -------------------------------------------------------------------------
# Interval shape: lower arm against upper arm, with the identity line
#
# df: lo_arm, up_arm, method, level_lab, breach.
#
# Log-log with coord_fixed() so equal decades make the identity line read as
# a true 45 degrees. Wald intervals land exactly on it by construction;
# likelihood intervals sit above it.
# -------------------------------------------------------------------------
plot_arm_shape <- function(df, title, subtitle, caption) {
  brk <- c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0)

  ggplot(df, aes(x = lo_arm, y = up_arm)) +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.4, color = "grey55") +
    # A SQUARE, not a ring: most breaching points are Chao, whose own mark is
    # an open circle, and a ring around a circle reads as one doubled circle.
    geom_point(
      data = df |> filter(breach),
      shape = 0,
      size = 4.2,
      color = "grey25",
      stroke = 0.5
    ) +
    geom_point(aes(color = method, shape = method), size = 1.9) +
    scale_color_manual(values = METHOD_COLORS, drop = FALSE) +
    scale_shape_manual(values = METHOD_SHAPES, drop = FALSE) +
    scale_x_log10(breaks = brk, labels = label_number()) +
    scale_y_log10(breaks = brk, labels = label_number()) +
    coord_fixed() +
    facet_wrap(~level_lab, nrow = 1) +
    labs(
      x = expression("Lower arm  " * hat(psi) - "lower"),
      y = expression("Upper arm  upper" - hat(psi)),
      title = title,
      subtitle = subtitle,
      caption = caption
    ) +
    theme_diagnostic()
}

# -------------------------------------------------------------------------
# Per-site forest plot against a reference estimate
#
# df: site_lab (factor, ordered by N), estimate, lower, upper, method.
# ref: site_lab, psi_ref — drawn as a vertical tick.
#
# Horizontal with sites on y: rows of two intervals read cleanly in
# portrait, where dodged groups of four do not. Log x because a single wide
# site otherwise squeezes the high-N rows, where the containment question is
# tightest, into a sliver.
# -------------------------------------------------------------------------
plot_forest <- function(df, ref, x_label, title, subtitle, caption,
                        floor_line = NULL, breaks = waiver(),
                        log_x = TRUE) {
  p <- ggplot(df, aes(y = site_lab, x = estimate, color = method,
    shape = method))

  if (!is.null(floor_line)) {
    p <- p + geom_vline(
      xintercept = floor_line,
      linetype = "22",
      linewidth = 0.4,
      color = "grey60"
    )
  }

  p +
    geom_pointrange(
      aes(xmin = lower, xmax = upper),
      orientation = "y",
      position = position_dodge(width = 0.62),
      size = 0.32,
      linewidth = 0.6
    ) +
    geom_point(
      data = ref,
      aes(y = site_lab, x = psi_ref),
      inherit.aes = FALSE,
      shape = 124,
      size = 3.1,
      color = "grey20"
    ) +
    scale_color_manual(values = METHOD_COLORS, drop = FALSE) +
    scale_shape_manual(values = METHOD_SHAPES, drop = FALSE) +
    (if (isTRUE(log_x)) {
      scale_x_log10(breaks = breaks, labels = label_number(accuracy = 0.01))
    } else {
      scale_x_continuous(breaks = breaks)
    }) +
    labs(
      x = x_label,
      y = NULL,
      title = title,
      subtitle = subtitle,
      caption = caption
    ) +
    theme_diagnostic() +
    theme(panel.grid.major.y = element_line(linewidth = 0.2,
      color = "grey92"))
}
