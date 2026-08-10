# =========================================================================
# explore-results.R — interactive scratchpad for experiment results
#
# Run this at the console (not as part of a render) to get every
# downloaded result into tidy data frames you can plot however you like:
#
#   source("Dissertation/explore-results.R")
#
# Leaves these in your session:
#
#   available   one row per downloaded bundle — what you have locally
#   intervals   point estimates + CI endpoints, one row per
#               site x method x confidence level
#   curves      the psi / log-pseudolikelihood grids behind those CIs
#   context     per-site metadata (N, J, psi_mle, support bounds)
#
# plus every plot helper (plot_interval_comparison, plot_interval_width,
# plot_curve_overlay) and the dune benchmarks (chao_interval_df,
# chao_shannon_df, unbiased_interval_df).
#
# This is deliberately separate from data-viz.R. That script is the
# committed, reproducible figure build and writes PNGs; this one writes
# nothing and exists purely for poking around.
#
# Nothing here is a dependency of the dissertation render.
# =========================================================================

library(tidyverse)

source("Dissertation/R/load-results.R")
source("Dissertation/R/plot-helpers.R")
source("Dissertation/R/dune-comparators.R")

# -------------------------------------------------------------------------
# What's downloaded?
# -------------------------------------------------------------------------
available <- list_bundles()

if (nrow(available) == 0L) {
  message(
    "\nNothing to explore yet. From the likelyr-simulations repo, with the\n",
    "Northwestern VPN active:\n\n",
    "    make download EXP=multinom/ne_entropy/exp_v13\n\n",
    "That copies one file (~100 KB) — the analysis bundle only, never the\n",
    "experiment output.\n"
  )
} else {
  print(available)
}

# -------------------------------------------------------------------------
# Load everything available
#
# load_all_applications() discovers what's present and derives the
# estimand / model labels from the bundle itself, so there's no registry to
# keep in sync while exploring. Simulation bundles are skipped — they hold
# per-iteration metrics rather than per-dataset estimates, and want the
# summarisers in likelyr-simulations/R/data_viz_utils.R instead.
# -------------------------------------------------------------------------
intervals <- load_all_applications("intervals")
curves <- load_all_applications("curves", quiet = TRUE)

context <- if (nrow(available) > 0L) {
  available |>
    filter(str_detect(tables, "context")) |>
    pmap(function(family, app, version, ...) {
      load_bundle(family, app, version)$context |>
        mutate(app = app, version = version, .before = 1)
    }) |>
    list_rbind()
} else {
  NULL
}

# -------------------------------------------------------------------------
# Summary
# -------------------------------------------------------------------------
if (!is.null(intervals)) {
  message(
    "\n✔ intervals: ", nrow(intervals), " rows across ",
    n_distinct(intervals$site), " site(s), ",
    n_distinct(paste(intervals$estimand, intervals$model)), " estimand/model combo(s)"
  )
  print(intervals |> count(app, version, estimand, model, method))
}

if (!is.null(curves)) {
  message("✔ curves:    ", nrow(curves), " grid points")
}

if (!is.null(context)) {
  message("✔ context:   ", nrow(context), " site-level rows")
}

message(
  "\nTry:\n",
  '  intervals |> filter(conf_level == 0.90) |> plot_interval_comparison("Shannon entropy")\n',
  '  curves |> filter(site <= 4) |> plot_curve_overlay("Shannon entropy")\n',
  "  intervals |> mutate(width = upper - lower) |>\n",
  "    select(site, N, method, conf_level, width) |> arrange(site)\n",
  "\nTo compare against the dune benchmarks, bind_rows() with\n",
  "chao_shannon_df / chao_interval_df / unbiased_interval_df — they already\n",
  "share the same column vocabulary.\n"
)
