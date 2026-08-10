# =========================================================================
# figures/ne-entropy-curves.R
#
# Profile vs. integrated log-pseudolikelihood curves for Shannon entropy
# under the no-effects multinomial model, one panel per dune site.
#
# Contrast with dune-site-intervals.R: that figure spans several
# experiments and so keeps a registry. This one reads a SINGLE experiment,
# so it just calls load_bundle() directly. Both are normal — the figure
# decides, not the loader.
#
# Uses the psi/loglik grids exported by likelyr-simulations'
# R/analyze_app.R. These live inside model.rds and never used to leave
# Quest, which is why curve figures weren't previously possible here.
# =========================================================================

FAMILY <- "multinom"
APP <- "ne_entropy"
VERSION <- "exp_v13"

if (!has_bundle(FAMILY, APP, VERSION)) {
  message(
    "⏭  Skipping entropy curve figures — no bundle for ",
    FAMILY, "/", APP, "/", VERSION
  )
} else {
  bundle <- load_bundle(FAMILY, APP, VERSION)

  curves <- as_curve_df(bundle, "Shannon", "No effects")

  # ---- All sites, trimmed to the region inside the widest CI cutoff ----
  save_figure(
    "ne-shannon-curves.png",
    curves,
    \(df) plot_curve_overlay(df, "Shannon entropy"),
    width = 9,
    height = 7
  )

  # ---- Untrimmed, to show the full traversed grid including the tails ----
  save_figure(
    "ne-shannon-curves-full.png",
    curves,
    \(df) plot_curve_overlay(df, "Shannon entropy", trim = FALSE),
    width = 9,
    height = 7
  )
}
