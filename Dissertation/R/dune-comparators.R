# =========================================================================
# dune-comparators.R
#
# Competing (non-likelihood) interval estimators for the dune meadow
# dataset, computed locally from vegan::dune.
#
# These are the benchmarks the profile and integrated likelihoods are
# compared against. They are computed here rather than on Quest because
# they are closed-form, cost nothing, and depend on packages (iNEXT,
# vegan) that have no business being installed in the methods repo.
#
# Provides:
#   conf_levels              the confidence levels used throughout
#   site_estimates_simpson   per-site Simpson estimate + SE
#   unbiased_interval_df     Simpson, Tiffeau-Mayer unbiased estimator
#   chao_interval_df         Simpson, Chao estimator (iNEXT)
#   chao_shannon_df          Shannon entropy, Chao estimator (iNEXT)
#   simpson_interval_df_se   +/-1 SE variants, for the Fig. 3 mimicry plots
#
# All interval frames use the same column vocabulary as
# R/load-results.R's as_interval_df(), so they bind_rows() directly with
# PL/IL results.
# =========================================================================

library(tidyverse)
library(vegan)
library(iNEXT)

data(dune)

z_from_conf <- function(cl) qnorm(1 - (1 - cl) / 2)

conf_levels <- tibble(
  conf_level = c(0.90, 0.95, 0.99),
  z = z_from_conf(conf_level)
)

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
# SHANNON ENTROPY — Chao only (the TM paper doesn't cover entropy)
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

# =========================================================================
# +/-1 SE VARIANTS — inputs to the Fig. 3 mimicry plots
# =========================================================================

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

simpson_interval_df_se <- bind_rows(
  unbiased_interval_df_se,
  chao_interval_df_se
) |>
  arrange(site, method)
