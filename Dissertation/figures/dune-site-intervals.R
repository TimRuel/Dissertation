# =========================================================================
# figures/dune-site-intervals.R
#
# Per-site interval comparisons on the dune meadow dataset: profile and
# integrated likelihood against the Chao and Tiffeau-Mayer benchmarks,
# for Shannon entropy and Simpson's index under the no-effects and
# fixed-effects models.
#
# This is ONE figure script. Each experiment (or family of related
# figures) gets its own script in this directory, declaring at the top
# which experiments it reads. Nothing here is generic — the shared parts
# live in ../R/, and this file is free to do whatever this particular set
# of figures needs.
#
# Expects to be sourced by ../data-viz.R, which loads:
#   ../R/load-results.R      load_applications(), load_bundle()
#   ../R/plot-helpers.R      save_figure(), plot_interval_comparison(), ...
#   ../R/dune-comparators.R  chao_*, unbiased_*, conf_levels
# =========================================================================

# -------------------------------------------------------------------------
# Which experiments supply the PL/IL results for this figure set
#
# This figure genuinely spans several experiments — one per (estimand,
# model) cell — so it uses a registry. That is a property of THIS figure,
# not of the architecture: a figure script that reads a single experiment
# just calls load_bundle() directly.
#
# Add a row when an experiment becomes available. Cells with no row, or
# whose bundle hasn't been downloaded yet, are skipped and their figures
# are not written — no placeholder data is substituted, so a figure on
# disk always reflects real output.
#
# `app` is the directory under experiments/<family>/, which is NOT always
# the spec directory name: the no-effects Simpson experiment uses the
# ne_simpson specs but lives under config/ and experiments/ as
# `logit_simpson`. bundle_path() only knows the experiments/ layout, so
# `app` must be the latter.
#
# Pending: fixed-effects entropy (applications/multinom/fe_entropy) has
# spec files but no experiment config under config/multinom/ yet.
# -------------------------------------------------------------------------

dune_applications <- tribble(
  ~estimand, ~model,       ~family,    ~app,            ~version,
  "Shannon", "No effects", "multinom", "ne_entropy",    "exp_v13",
  "Simpson", "No effects", "multinom", "logit_simpson", "exp_v6"
)

pl_il_intervals <- load_applications(dune_applications, what = "intervals")

if (is.null(pl_il_intervals)) {
  message(
    "\n⚠ No application results available — only the comparator-only ",
    "figures will be written.\n"
  )
  pl_il_intervals <- tibble(
    estimand = character(),
    model = character(),
    site = integer(),
    N = numeric(),
    method = character(),
    conf_level = numeric(),
    estimate = numeric(),
    lower = numeric(),
    upper = numeric()
  )
}

# =========================================================================
# COMBINED FRAMES — comparators plus whatever PL/IL is available
#
# The Simpson comparators (Tiffeau-Mayer, Chao) and the Shannon comparator
# (Chao) are computed locally from `dune`; PL and IL come from the
# experiments above. The Simpson frame therefore always has content even
# with no PL/IL results, while the Shannon frame has Chao plus whichever
# models have been run.
# =========================================================================

simpson_interval_df <- bind_rows(
  unbiased_interval_df,
  chao_interval_df,
  pl_il_intervals |> filter(estimand == "Simpson")
) |>
  arrange(site, model, method, conf_level)

entropy_interval_df <- bind_rows(
  chao_shannon_df,
  pl_il_intervals |> filter(estimand == "Shannon")
) |>
  arrange(site, model, method, conf_level)

# =========================================================================
# FIGURES
# =========================================================================

# ---- No-effects: Shannon entropy (Chao + PL/IL) ----
p_ne_shannon <- save_figure(
  "ne-shannon-intervals.png",
  entropy_interval_df |> filter(model == "No effects"),
  \(df) plot_interval_comparison(df, "Shannon entropy", METHOD_COLORS)
)

# ---- No-effects: Simpson's index, all available methods ----
p_ne_simpson <- save_figure(
  "ne-simpson-intervals.png",
  simpson_interval_df |> filter(model == "No effects"),
  \(df) plot_interval_comparison(df, "Simpson's index", METHOD_COLORS)
)

# ---- No-effects: Simpson's index, interval width supplement ----
p_ne_simpson_width <- save_figure(
  "ne-simpson-width.png",
  simpson_interval_df |> filter(model == "No effects"),
  \(df) plot_interval_width(df, "Simpson's index", method_colors = METHOD_COLORS)
)

# ---- Fixed-effects: Shannon entropy (PL/IL only) ----
p_fe_shannon <- save_figure(
  "fe-shannon-intervals.png",
  entropy_interval_df |> filter(model == "Fixed effects"),
  \(df) plot_interval_comparison(df, "Shannon entropy", METHOD_COLORS)
)

# ---- Fixed-effects: Simpson's index (PL/IL only) ----
p_fe_simpson <- save_figure(
  "fe-simpson-intervals.png",
  simpson_interval_df |>
    filter(model == "Fixed effects", method %in% c("PL", "IL")),
  \(df) plot_interval_comparison(df, "Simpson's index", METHOD_COLORS)
)

# ---- Model comparison: no-effects vs. fixed-effects width, Shannon ----
p_model_width_shannon <- save_figure(
  "model-comparison-width-shannon.png",
  entropy_interval_df |>
    filter(method %in% c("PL", "IL")) |>
    require_both_models(),
  \(df) plot_interval_width(df, "Shannon entropy", group_var = "model")
)

# ---- Model comparison: no-effects vs. fixed-effects width, Simpson ----
p_model_width_simpson <- save_figure(
  "model-comparison-width-simpson.png",
  simpson_interval_df |>
    filter(method %in% c("PL", "IL")) |>
    require_both_models(),
  \(df) plot_interval_width(df, "Simpson's index", group_var = "model")
)

# =========================================================================
# FIG. 3 MIMICRY — Chao vs. Tiffeau-Mayer only
#
# Reproductions of the Tiffeau-Mayer paper's Figure 3 layout (thick Chao
# band behind thin TM pointranges), kept for reference. Defined but not
# written to disk, matching their previous status in data-viz.R; call
# either one at the console to inspect it.
# =========================================================================

plot_quantile_intervals <- function(df, level = 0.90) {
  df <- df |>
    filter(conf_level == level, method %in% c("Unbiased (TM)", "Chao")) |>
    mutate(method = factor(method, levels = c("Unbiased (TM)", "Chao")))

  ggplot() +
    geom_linerange(
      data = df |> filter(method == "Chao"),
      aes(x = site, ymin = lower, ymax = upper),
      color = "#F4A96C",
      linewidth = 4.5,
      alpha = 0.9
    ) +
    geom_pointrange(
      data = df |> filter(method == "Unbiased (TM)"),
      aes(x = site, y = estimate, ymin = lower, ymax = upper, color = method),
      size = 0.5,
      linewidth = 0.7
    ) +
    scale_color_manual(values = c("Unbiased (TM)" = "#2C4A6E"), name = NULL) +
    scale_x_continuous(breaks = seq(1, 19, by = 2)) +
    labs(
      x = "Site",
      y = expression("Simpson index" ~ hat(p)[C]),
      title = paste0(scales::percent(level), " confidence intervals")
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = c(0.85, 0.85)
    ) +
    guides(color = guide_legend(override.aes = list(linewidth = 0.7)))
}

plot_se_mimicry <- function(df) {
  df <- df |>
    mutate(method = factor(method, levels = c("Unbiased (TM)", "Chao")))

  ggplot() +
    geom_linerange(
      data = df |> filter(method == "Chao"),
      aes(x = site, ymin = lower, ymax = upper),
      color = "#F4A96C",
      linewidth = 4.5,
      alpha = 0.9
    ) +
    geom_pointrange(
      data = df |> filter(method == "Unbiased (TM)"),
      aes(x = site, y = estimate, ymin = lower, ymax = upper, color = method),
      size = 0.5,
      linewidth = 0.7
    ) +
    scale_color_manual(values = c("Unbiased (TM)" = "#2C4A6E"), name = NULL) +
    scale_x_continuous(breaks = seq(1, 19, by = 2)) +
    labs(x = "Site", y = expression("Simpson index" ~ hat(p)[C])) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = c(0.85, 0.85)
    ) +
    guides(color = guide_legend(override.aes = list(linewidth = 0.7)))
}
