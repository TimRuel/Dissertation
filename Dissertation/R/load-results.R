# =========================================================================
# load-results.R
#
# Reads experiment results produced by likelyr-simulations and reshapes
# them into the column names this repo's plotting functions expect.
#
# Where results live
# ------------------
# The .rds results are NOT committed here. They live in the
# likelyr-simulations checkout, downloaded from Quest by:
#
#   Quest:  make results  EXP_CONFIG=config/<family>/<app>/<exp_vX>/<exp_vX>.yml
#   Local:  make download EXP=<family>/<app>/<exp_vX>
#
# which lands one file at:
#
#   <LIKELYR_SIMS_DIR>/experiments/<family>/<app>/<exp_vX>/analysis/bundle.rds
#
# One copy of every result, no binary churn in this repo's history. What
# IS committed here is the rendered PNGs under files/Images, so the
# dissertation renders from a clean clone without needing Quest access,
# a VPN, or the likelyr-simulations checkout. Only *regenerating* the
# figures needs those.
#
# Override the source location with the LIKELYR_SIMS_DIR environment
# variable (e.g. when regenerating figures on Quest, where the Windows
# default path does not exist).
# =========================================================================

library(dplyr)
library(fs)
library(tibble)

LIKELYR_SIMS_DIR <- Sys.getenv(
  "LIKELYR_SIMS_DIR",
  "C:/Northwestern/likelyr-simulations"
)

# -------------------------------------------------------------------------
# Path to one experiment's bundle
# -------------------------------------------------------------------------
bundle_path <- function(family, app, version) {
  path(
    LIKELYR_SIMS_DIR,
    "experiments",
    family,
    app,
    version,
    "analysis",
    "bundle.rds"
  )
}

# -------------------------------------------------------------------------
# Is an experiment's bundle available locally?
#
# Used by the figure scripts to skip figures whose experiment hasn't been
# run or downloaded yet, rather than substituting placeholder data.
# -------------------------------------------------------------------------
has_bundle <- function(family, app, version) {
  file_exists(bundle_path(family, app, version))
}

# -------------------------------------------------------------------------
# Load one experiment bundle
#
# Returns the named list written by likelyr-simulations/R/bundle_exp.R:
#   $meta plus, depending on experiment kind,
#   application: $estimates $curves $context
#   simulation:  $point_metrics $interval_metrics $invalid_ci_index
# -------------------------------------------------------------------------
load_bundle <- function(family, app, version) {
  p <- bundle_path(family, app, version)

  if (!file_exists(p)) {
    stop(
      "No results bundle at:\n  ",
      p,
      "\n\nOn Quest:  make results  EXP_CONFIG=config/",
      family,
      "/",
      app,
      "/",
      version,
      "/",
      version,
      ".yml",
      "\nThen local: make download EXP=",
      family,
      "/",
      app,
      "/",
      version,
      "\n\n(Or set LIKELYR_SIMS_DIR if the checkout is elsewhere; it is ",
      "currently '",
      LIKELYR_SIMS_DIR,
      "'.)",
      call. = FALSE
    )
  }

  readRDS(p)
}

# -------------------------------------------------------------------------
# Site index from a sim id
#
# Application experiments encode dataset identity in the sim index —
# sim_04 is dune site 4 (see applications/multinom/ne_entropy/data.R,
# which uses the sim index as the row index into vegan::dune).
# -------------------------------------------------------------------------
site_from_sim_id <- function(simulation) {
  as.integer(sub("^sim_0*", "", simulation))
}

# -------------------------------------------------------------------------
# Application interval estimates, in this repo's plotting vocabulary
#
# Column mapping from the analyzer's output to what
# plot_interval_comparison() / plot_interval_width() expect:
#
#   simulation "sim_04"              -> site       4
#   pseudolikelihood "Profile"       -> method     "PL"
#   pseudolikelihood "Integrated"    -> method     "IL"
#   alpha 0.10                       -> conf_level 0.90
#   psi_hat                          -> estimate
#   lower / upper                    -> unchanged
#   context$n_obs                    -> N
#
# conf_level is derived from alpha rather than parsing the "90%" level
# string, so it stays exactly numeric and joins cleanly against the
# locally-computed comparators (Chao, Unbiased TM).
#
# Rows with a non-finite endpoint are dropped: an invalid CI cannot be
# drawn as a pointrange, and silently plotting a half-open interval would
# misrepresent it. The count of dropped rows is reported so a systematic
# failure doesn't pass unnoticed.
# -------------------------------------------------------------------------
as_interval_df <- function(bundle, estimand, model) {
  if (is.null(bundle$estimates)) {
    stop(
      "Bundle for ",
      bundle$meta$experiment,
      " has no $estimates table — kind is '",
      paste(bundle$meta$kind, collapse = ", "),
      "', expected 'application'.",
      call. = FALSE
    )
  }

  context <- bundle$context |>
    select(simulation, n_obs) |>
    distinct()

  df <- bundle$estimates |>
    left_join(context, by = "simulation") |>
    mutate(
      # .env$ so these resolve to the function arguments even if a future
      # analyzer version starts writing estimand/model columns of its own.
      estimand = .env$estimand,
      model = .env$model,
      site = site_from_sim_id(simulation),
      method = case_when(
        pseudolikelihood == "Profile" ~ "PL",
        pseudolikelihood == "Integrated" ~ "IL",
        TRUE ~ pseudolikelihood
      ),
      conf_level = 1 - alpha,
      estimate = psi_hat,
      N = n_obs
    )

  n_invalid <- sum(!is.finite(df$lower) | !is.finite(df$upper))

  if (n_invalid > 0L) {
    message(
      "ℹ ",
      bundle$meta$experiment,
      " (",
      estimand,
      ", ",
      model,
      "): dropping ",
      n_invalid,
      " of ",
      nrow(df),
      " rows with a non-finite CI endpoint."
    )
  }

  df |>
    filter(is.finite(lower), is.finite(upper)) |>
    select(
      estimand,
      model,
      site,
      N,
      method,
      conf_level,
      estimate,
      lower,
      upper
    ) |>
    arrange(site, model, method, conf_level)
}

# -------------------------------------------------------------------------
# Application likelihood curves, in plotting vocabulary
#
# rel_loglik is the curve to plot: it is the log-likelihood shifted so its
# maximum is 0, which is what puts the profile and integrated curves on a
# common vertical scale. above_crit marks the grid region inside the
# widest CI cutoff, useful for trimming long uninformative tails.
# -------------------------------------------------------------------------
as_curve_df <- function(bundle, estimand, model) {
  if (is.null(bundle$curves) || nrow(bundle$curves) == 0L) {
    stop(
      "Bundle for ",
      bundle$meta$experiment,
      " has no $curves table.",
      call. = FALSE
    )
  }

  context <- bundle$context |>
    select(simulation, n_obs, psi_mle) |>
    distinct()

  bundle$curves |>
    left_join(context, by = "simulation") |>
    mutate(
      estimand = .env$estimand,
      model = .env$model,
      site = site_from_sim_id(simulation),
      method = case_when(
        pseudolikelihood == "Profile" ~ "PL",
        pseudolikelihood == "Integrated" ~ "IL",
        TRUE ~ pseudolikelihood
      ),
      N = n_obs
    ) |>
    select(
      estimand,
      model,
      site,
      N,
      method,
      psi,
      loglik,
      any_of(c("rel_loglik", "above_crit")),
      psi_mle
    ) |>
    arrange(site, method, psi)
}

# -------------------------------------------------------------------------
# Load a registry of applications into one long data frame
#
# registry is a data frame with columns estimand, model, family, app,
# version. Entries whose bundle is not present locally are skipped with a
# message rather than erroring, so a figure script can be written once for
# the full set of planned experiments and produce whatever is ready.
#
# Returns NULL if nothing is available, so callers can guard on that
# instead of on an empty-tibble edge case.
# -------------------------------------------------------------------------
load_applications <- function(registry, what = c("intervals", "curves")) {
  what <- match.arg(what)
  reshape <- if (what == "intervals") as_interval_df else as_curve_df

  parts <- list()

  for (i in seq_len(nrow(registry))) {
    r <- registry[i, ]

    if (!has_bundle(r$family, r$app, r$version)) {
      message(
        "⏭  Skipping ",
        r$estimand,
        " / ",
        r$model,
        " — no bundle for ",
        r$family,
        "/",
        r$app,
        "/",
        r$version
      )
      next
    }

    bundle <- load_bundle(r$family, r$app, r$version)

    parts[[length(parts) + 1]] <- reshape(bundle, r$estimand, r$model)
  }

  if (length(parts) == 0L) {
    return(NULL)
  }

  bind_rows(parts)
}
