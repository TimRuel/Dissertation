# =========================================================================
# data-viz.R — regenerate the dissertation's data figures
#
# Run this from the project root:
#   Rscript Dissertation/data-viz.R
#
# Architecture
# ------------
# Figures differ from experiment to experiment, so there is no single
# script that draws them all. The work is split into three layers:
#
#   1. R/load-results.R      Experiment-agnostic. Finds and reads the
#                            bundle.rds an experiment produced, and
#                            reshapes it into this repo's column
#                            vocabulary. Never changes when a new
#                            experiment is added.
#
#   2. R/plot-helpers.R      The shared plotting vocabulary — the plot
#      R/dune-comparators.R  types and benchmark estimators that more
#                            than one figure script uses.
#
#   3. figures/*.R           One script per experiment or family of
#                            related figures. These are SUPPOSED to
#                            differ. Each declares which experiments it
#                            reads, composes layer 2 however it needs,
#                            and saves under its own filenames.
#
# To add figures for a new experiment: write a new figures/*.R. It is
# picked up automatically below. Nothing in layers 1 or 2 needs to change
# unless you find yourself writing the same plot twice.
#
# Where the data comes from
# -------------------------
# Results are NOT stored in this repo. They live in the
# likelyr-simulations checkout, downloaded from Quest. See
# R/load-results.R for the path convention and the refresh commands.
# Figures are committed as PNGs, so the dissertation renders from a clean
# clone without needing Quest access — only regenerating them needs the
# results.
#
# Each figure script is sourced in its own error handler, so one script
# failing (a missing bundle, a package that isn't installed) does not stop
# the others from producing their figures.
# =========================================================================

library(tidyverse)

source("Dissertation/R/load-results.R")
source("Dissertation/R/plot-helpers.R")
source("Dissertation/R/dune-comparators.R")

figure_scripts <- sort(
  list.files("Dissertation/figures", pattern = "\\.R$", full.names = TRUE)
)

if (length(figure_scripts) == 0L) {
  stop("No figure scripts found in Dissertation/figures/", call. = FALSE)
}

failed <- character(0)

for (script in figure_scripts) {
  message("\n═══ ", basename(script), " ", strrep("═", max(0, 50 - nchar(basename(script)))))

  ok <- tryCatch(
    {
      # Each script runs in its own environment so two figure scripts can
      # both use a name like `pl_il_intervals` without colliding. The
      # parent is set explicitly to globalenv() — new.env()'s default of
      # parent.frame() would resolve to a tryCatch internal frame here,
      # and the script would not see the helpers sourced above.
      source(script, local = new.env(parent = globalenv()))
      TRUE
    },
    error = function(e) {
      message("✗ ", basename(script), " failed: ", conditionMessage(e))
      FALSE
    }
  )

  if (!ok) {
    failed <- c(failed, basename(script))
  }
}

message("\n", strrep("═", 56))

if (length(failed) > 0L) {
  message(
    "⚠ ",
    length(failed),
    " of ",
    length(figure_scripts),
    " figure script(s) failed: ",
    paste(failed, collapse = ", ")
  )
} else {
  message("✔ All ", length(figure_scripts), " figure script(s) completed")
}

# Note: fig-simplex-level-sets.R is a standalone illustrative figure with
# no experiment dependency and is not run from here. Move it into
# figures/ if you want it regenerated as part of this sweep.