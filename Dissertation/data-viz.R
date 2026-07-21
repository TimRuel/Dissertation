library(tidyverse)
library(vegan)
library(iNEXT)
library(scales)

data(dune)

# Output directory for figures used in the dissertation
fig_dir <- "Dissertation/files/Images"
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# =========================================================================
# SIMPSON'S INDEX — shared estimating functions
# =========================================================================

simpson_ub <- function(x) {
  N <- sum(x)
  sum(x * (x - 1)) / (N * (N - 1))
}

simpson_ub_var <- function(x) {
  N <- sum(x)
  a <- (4 * (N - 2)) / (N * (N - 1))
  b <- (2 * (2 * N - 3)) / (N * (N - 1))
  c <- 2 / (N * (N - 1))

  p_hat_C <- sum(x * (x - 1)) / (N * (N - 1))
  p_hat_T <- sum(x * (x - 1) * (x - 2)) / (N * (N - 1) * (N - 2))

  a / (1 - b) * p_hat_T - b / (1 - b) * p_hat_C^2 + c / (1 - b) * p_hat_C
}

simpson_poisson_var <- function(x) {
  N <- sum(x)
  p_hat_C <- sum(x * (x - 1)) / (N * (N - 1))
  2 / (N * (N - 1)) * p_hat_C
}

z_from_conf <- function(cl) qnorm(1 - (1 - cl) / 2)

conf_levels <- tibble(
  conf_level = c(0.90, 0.95, 0.99),
  z = z_from_conf(conf_level)
)

site_estimates_simpson <- dune |>
  as_tibble() |>
  mutate(site = row_number(), .before = 1) |>
  rowwise() |>
  mutate(
    N = sum(c_across(-site)),
    estimate = simpson_ub(c_across(-c(site, N))),
    var_ub = simpson_ub_var(c_across(-c(site, N, estimate))),
    var_poisson = simpson_poisson_var(c_across(-c(site, N, estimate, var_ub))),
    var_hat = max(var_ub, var_poisson),
    se = sqrt(var_hat)
  ) |>
  ungroup() |>
  select(site, N, estimate, se)

# ---- Unbiased (Tiffeau-Mayer) ----
unbiased_interval_df <- site_estimates_simpson |>
  cross_join(conf_levels) |>
  mutate(
    estimand = "Simpson",
    model = "No effects",
    method = "Unbiased (TM)",
    lower = estimate - z * se,
    upper = estimate + z * se
  ) |>
  select(estimand, model, site, N, method, conf_level, estimate, lower, upper)

# ---- Chao (Simpson) ----
chao_interval_df <- conf_levels$conf_level |>
  map(function(cl) {
    ChaoSimpson(t(dune), conf = cl) |>
      as_tibble() |>
      mutate(
        site = row_number(),
        estimand = "Simpson",
        model = "No effects",
        method = "Chao",
        conf_level = cl,
        estimate = 1 - Estimator,
        se = Est_s.e.,
        lower = estimate - z_from_conf(cl) * se,
        upper = estimate + z_from_conf(cl) * se
      ) |>
      select(estimand, model, site, method, conf_level, estimate, lower, upper)
  }) |>
  list_rbind() |>
  left_join(site_estimates_simpson |> select(site, N), by = "site")

# =========================================================================
# DUMMY PL / IL DATA (placeholder — replace with real likelihood output)
# =========================================================================

generate_dummy_pl_il <- function(
  estimand_name,
  model_name,
  site_estimates,
  methods = c("PL", "IL"),
  width_multiplier = 1.0,
  seed = 42
) {
  set.seed(seed)
  site_estimates |>
    select(site, N, estimate) |>
    cross_join(conf_levels) |>
    cross_join(tibble(method = methods)) |>
    mutate(
      estimand = estimand_name,
      model = model_name,
      jitter = rnorm(n(), 0, 0.01),
      width_scale = if_else(method == "IL", 0.85, 1.0) * width_multiplier,
      se_placeholder = 0.03 + 0.002 * N,
      half_width = z * se_placeholder * width_scale,
      estimate = estimate + jitter,
      lower = estimate - half_width,
      upper = estimate + half_width
    ) |>
    select(estimand, model, site, N, method, conf_level, estimate, lower, upper)
}

pl_il_simpson_ne <- generate_dummy_pl_il(
  "Simpson", "No effects", site_estimates_simpson
)

# Dummy fixed-effects placeholder — narrower than no-effects, illustrating
# the borrowing-strength effect this section's prose discusses. Replace
# with real fixed-effects likelihood output once available.
pl_il_simpson_fe <- generate_dummy_pl_il(
  "Simpson", "Fixed effects", site_estimates_simpson,
  width_multiplier = 0.75, seed = 43
)

# =========================================================================
# COMBINED SIMPSON DATAFRAME — Unbiased, Chao, PL, IL, both models
# =========================================================================

simpson_interval_df <- bind_rows(
  unbiased_interval_df,
  chao_interval_df,
  pl_il_simpson_ne,
  pl_il_simpson_fe
) |>
  arrange(site, model, method, conf_level)

# =========================================================================
# SHANNON ENTROPY — Chao only (TM paper doesn't cover entropy) + PL/IL
# =========================================================================

site_estimates_shannon <- dune |>
  as_tibble() |>
  mutate(site = row_number(), .before = 1) |>
  rowwise() |>
  mutate(N = sum(c_across(-site))) |>
  ungroup() |>
  select(site, N)

chao_shannon_df <- conf_levels$conf_level |>
  map(function(cl) {
    ChaoShannon(t(dune), conf = cl) |>
      as_tibble() |>
      mutate(
        site = row_number(),
        estimand = "Shannon",
        model = "No effects",
        method = "Chao",
        conf_level = cl,
        estimate = Estimator,
        se = Est_s.e,
        lower = estimate - z_from_conf(cl) * se,
        upper = estimate + z_from_conf(cl) * se
      ) |>
      select(estimand, model, site, method, conf_level, estimate, lower, upper)
  }) |>
  list_rbind() |>
  left_join(site_estimates_shannon, by = "site")

pl_il_shannon_ne <- generate_dummy_pl_il(
  "Shannon", "No effects",
  chao_shannon_df |> filter(conf_level == 0.90) |> select(site, N, estimate)
)

pl_il_shannon_fe <- generate_dummy_pl_il(
  "Shannon", "Fixed effects",
  chao_shannon_df |> filter(conf_level == 0.90) |> select(site, N, estimate),
  width_multiplier = 0.75, seed = 43
)

entropy_interval_df <- bind_rows(
  chao_shannon_df,
  pl_il_shannon_ne,
  pl_il_shannon_fe
) |>
  arrange(site, model, method, conf_level)

# =========================================================================
# GENERIC PLOTTING FUNCTIONS
# =========================================================================

# Full comparison: all methods, faceted by confidence level. Sites ordered
# by N (ascending) so sparsity-driven differences (e.g. vs. Chao) are
# visually apparent left-to-right rather than requiring cross-reference
# to a table.
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

# Compact supplement: interval width only, easier to compare precision
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

# =========================================================================
# ORIGINAL TWO-METHOD PLOTS (quantile-based and +/-1 SE mimicry of Fig. 3)
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

unbiased_interval_df_se <- site_estimates_simpson |>
  mutate(
    method = "Unbiased (TM)",
    lower = estimate - se,
    upper = estimate + se
  ) |>
  select(site, N, method, estimate, lower, upper)

chao_interval_df_se <- ChaoSimpson(t(dune), conf = 0.6827) |>
  as_tibble() |>
  mutate(
    site = row_number(),
    method = "Chao",
    estimate = 1 - Estimator,
    se = Est_s.e.,
    lower = estimate - se,
    upper = estimate + se
  ) |>
  select(site, method, estimate, lower, upper) |>
  left_join(site_estimates_simpson |> select(site, N), by = "site")

simpson_interval_df_se <- bind_rows(unbiased_interval_df_se, chao_interval_df_se) |>
  arrange(site, method)

plot_se_mimicry <- function(df) {
  df <- df |> mutate(method = factor(method, levels = c("Unbiased (TM)", "Chao")))

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

# =========================================================================
# GENERATE, SAVE, AND ASSIGN FIGURES FOR THE DISSERTATION
# =========================================================================

# ---- No-effects: Shannon entropy (PL/IL only) ----
p_ne_shannon <- plot_interval_comparison(
  entropy_interval_df |> filter(model == "No effects"),
  "Shannon entropy"
)
ggsave(file.path(fig_dir, "ne-shannon-intervals.png"), p_ne_shannon,
  width = 6.5, height = 8, dpi = 300)

# ---- No-effects: Simpson's index, all four methods ----
p_ne_simpson <- plot_interval_comparison(
  simpson_interval_df |> filter(model == "No effects"),
  "Simpson's index"
)
ggsave(file.path(fig_dir, "ne-simpson-intervals.png"), p_ne_simpson,
  width = 6.5, height = 8, dpi = 300)

# ---- No-effects: Simpson's index, interval width supplement ----
p_ne_simpson_width <- plot_interval_width(
  simpson_interval_df |> filter(model == "No effects"),
  "Simpson's index"
)
ggsave(file.path(fig_dir, "ne-simpson-width.png"), p_ne_simpson_width,
  width = 6.5, height = 8, dpi = 300)

# ---- Fixed-effects: Shannon entropy (PL/IL only) ----
p_fe_shannon <- plot_interval_comparison(
  entropy_interval_df |> filter(model == "Fixed effects"),
  "Shannon entropy"
)
ggsave(file.path(fig_dir, "fe-shannon-intervals.png"), p_fe_shannon,
  width = 6.5, height = 8, dpi = 300)

# ---- Fixed-effects: Simpson's index (PL/IL only) ----
p_fe_simpson <- plot_interval_comparison(
  simpson_interval_df |> filter(model == "Fixed effects", method %in% c("PL", "IL")),
  "Simpson's index"
)
ggsave(file.path(fig_dir, "fe-simpson-intervals.png"), p_fe_simpson,
  width = 6.5, height = 8, dpi = 300)

# ---- Model comparison: no-effects vs. fixed-effects width, Shannon ----
p_model_width_shannon <- plot_interval_width(
  entropy_interval_df |> filter(method %in% c("PL", "IL")),
  "Shannon entropy",
  group_var = "model"
)
ggsave(file.path(fig_dir, "model-comparison-width-shannon.png"), p_model_width_shannon,
  width = 6.5, height = 8, dpi = 300)

# ---- Model comparison: no-effects vs. fixed-effects width, Simpson ----
p_model_width_simpson <- plot_interval_width(
  simpson_interval_df |> filter(method %in% c("PL", "IL")),
  "Simpson's index",
  group_var = "model"
)
ggsave(file.path(fig_dir, "model-comparison-width-simpson.png"), p_model_width_simpson,
  width = 6.5, height = 8, dpi = 300)