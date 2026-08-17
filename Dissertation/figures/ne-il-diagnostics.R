# =========================================================================
# figures/ne-il-diagnostics.R
#
# Diagnostic figures for the INTEGRATED likelihood under the no-effects
# multinomial model on the dune meadow data, for BOTH estimands.
#
# These are not more per-site interval plots. The per-site comparison
# already exists (figures/dune-site-intervals.R); what it cannot show is
# *why* PL and IL differ, because near-duplicate confidence-level facets and
# four dodged pointranges hide the two things that separate them:
#
#   1. WHERE the point estimate sits. PL is anchored at the plug-in MLE, so
#      it inherits that estimator's exactly-known bias. IL moves off it, in
#      the direction that corrects the plug-in's bias.
#
#   2. The SHAPE of the interval. PL and IL are asymmetric because they
#      invert a likelihood ratio; the Wald benchmarks are exactly symmetric
#      by construction, and so can leave the parameter space.
#
# The two estimands are a genuine pair, not a copy-paste: the plug-in
# estimator is biased UPWARD for Simpson (sum of squares) and DOWNWARD for
# Shannon (entropy), so the reference correction and the sign of PL's offset
# both flip. Each estimand supplies its own bias-corrected reference:
#
#   Simpson  Tiffeau-Mayer unbiased estimator,  predicted bias  +(1-psi)/N
#   Shannon  Miller-Madow corrected entropy,    predicted bias  -(S-1)/(2N)
#
# EVERY NUMBER QUOTED IN A CAPTION IS COMPUTED FROM THE DATA BELOW. Nothing
# is hard-coded, so re-running an experiment cannot leave a figure asserting
# a statistic it no longer shows.
#
# Coverage is NOT addressed here: this is real data with no psi_0, so
# nothing in this file speaks to whether IL's width is justified. That
# argument belongs with the simulation experiments; these figures establish
# the geometry the coverage discussion refers back to.
#
# Expects to be sourced by ../data-viz.R, which loads:
#   ../R/load-results.R      load_applications()
#   ../R/plot-helpers.R      save_figure(), plot_bias_anatomy(), ...
#   ../R/dune-comparators.R  unbiased_interval_df, chao_*, dune
# =========================================================================

MAIN_LEVEL <- 0.90 # psi_hat is level-invariant, so this matters only to the
# forest plot; verified one distinct estimate per site x method.

# Per-site quantities the references are built from. `dune` is attached by
# dune-comparators.R; S is the number of species actually observed at a
# site, which is what the Miller-Madow correction depends on (not J).
dune_site_stats <- tibble(
  site = seq_len(nrow(dune)),
  N = rowSums(dune),
  S = rowSums(dune > 0),
  H_plugin = apply(dune, 1, function(x) {
    p <- x / sum(x)
    p <- p[p > 0]
    -sum(p * log(p))
  })
)

# -------------------------------------------------------------------------
# One spec per estimand
#
# ref_fn returns site, N, psi_ref, pred_bias:
#   psi_ref    the bias-corrected estimate the figures measure against
#   pred_bias  the plug-in MLE's predicted bias RELATIVE to psi_ref, signed
# -------------------------------------------------------------------------
ne_specs <- list(
  list(
    slug = "simpson",
    estimand = "Simpson",
    app = "logit_simpson",
    version = "exp_v6",
    x_label = "Simpson's index",
    ref_name = "Tiffeau-Mayer unbiased estimate",
    pred_label = "Predicted MLE bias  (1 - ψ)/N",
    floor = 1 / 30,
    floor_note = paste(
      "parameter-space floor ψ ≥ 1/30, since J is fixed at 30 for",
      "every site"
    ),
    forest_breaks = c(0.03, 0.05, 0.1, 0.2, 0.4),
    forest_log = TRUE,
    ref_fn = function() {
      unbiased_interval_df |>
        filter(conf_level == MAIN_LEVEL) |>
        select(site, N, psi_ref = estimate) |>
        mutate(pred_bias = (1 - psi_ref) / N)
    },
    comparators = function() bind_rows(unbiased_interval_df, chao_interval_df),
    breach_fn = function(lower, upper) lower < 1 / 30
  ),
  list(
    slug = "shannon",
    estimand = "Shannon",
    app = "ne_entropy",
    version = "exp_v13",
    x_label = "Shannon entropy",
    ref_name = "Miller-Madow corrected estimate",
    pred_label = "Predicted MLE bias  -(S - 1)/2N",
    floor = NULL,
    floor_note = NULL,
    forest_breaks = c(1.0, 1.5, 2.0, 2.5, 3.0),
    forest_log = FALSE,
    ref_fn = function() {
      dune_site_stats |>
        transmute(
          site,
          N,
          psi_ref = H_plugin + (S - 1) / (2 * N),
          pred_bias = -(S - 1) / (2 * N)
        )
    },
    comparators = function() chao_shannon_df,
    # Entropy is bounded in [0, log J]; no site comes near either end here,
    # but the check is kept so the figure reports honestly if one ever does.
    breach_fn = function(lower, upper) lower < 0 | upper > log(30)
  )
)

# =========================================================================
# BUILD
# =========================================================================

for (spec in ne_specs) {
  pl_il <- load_applications(
    tribble(
      ~estimand, ~model, ~family, ~app, ~version,
      spec$estimand, "No effects", "multinom", spec$app, spec$version
    ),
    what = "intervals"
  )

  if (is.null(pl_il)) {
    message(
      "⏭  No ", spec$estimand, " results — IL diagnostics skipped."
    )
    next
  }

  ref <- spec$ref_fn()
  sites_present <- sort(unique(pl_il$site))
  n_present <- length(sites_present)
  n_total <- nrow(dune_site_stats)

  # An experiment that analysed fewer sims than it planned must say so on
  # the figure itself; ne_entropy/exp_v13 analysed 17 of 20.
  coverage_note <- if (n_present < n_total) {
    paste0(
      " Based on the ", n_present, " of ", n_total,
      " sites the experiment analysed (missing: ",
      paste(setdiff(seq_len(n_total), sites_present), collapse = ", "), ")."
    )
  } else {
    ""
  }

  # ---- point estimates, one row per site x method (level-invariant) ----
  pe <- pl_il |>
    filter(conf_level == MAIN_LEVEL) |>
    select(site, N, method, estimate) |>
    left_join(ref |> select(site, psi_ref, pred_bias), by = "site") |>
    mutate(delta = estimate - psi_ref)

  wide <- pe |>
    select(site, N, method, delta, pred_bias) |>
    pivot_wider(names_from = method, values_from = delta)

  # ---------------------------------------------------------------------
  # FIGURE 1 — bias anatomy
  # ---------------------------------------------------------------------
  ratio <- wide$PL / wide$pred_bias
  # Overshoot vs undershoot is read off the data, not assumed: if IL ends up
  # on the opposite side of the reference from PL it crossed past it.
  crossed <- sign(median(wide$IL)) != sign(median(wide$PL))

  pts <- bind_rows(
    ref |>
      filter(site %in% sites_present) |>
      transmute(N, delta = pred_bias, series = spec$pred_label),
    pe |> transmute(N, delta, series = as.character(method))
  ) |>
    mutate(series = factor(series, levels = c("IL", "PL", spec$pred_label)))

  p_bias <- save_figure(
    paste0("ne-", spec$slug, "-bias-anatomy.png"),
    pe,
    function(df) {
      plot_bias_anatomy(
        pts,
        y_expr = bquote(hat(psi) - hat(psi)[.(if (spec$slug == "simpson")
          "unbiased" else "MM")]),
        pred_label = spec$pred_label,
        title = "Where each pseudolikelihood puts the point estimate",
        subtitle = wrap_subtitle(paste0(
          "PL sits at the plug-in MLE's known bias; IL moves off it ",
          if (crossed) {
            "and overshoots the corrected value."
          } else {
            "toward the corrected value but stops short."
          }
        )),
        caption = wrap_caption(
          paste0(
            "Difference from the ", spec$ref_name, ", which is the zero ",
            "line. Observed PL offset ÷ predicted bias = ",
            sprintf("%.3f–%.3f", min(ratio), max(ratio)),
            " (median ", sprintf("%.3f", median(ratio)), ") over the ",
            n_present, " sites. Median IL offset ",
            sprintf("%+.4f", median(wide$IL)), ".", coverage_note
          ),
          fig_width = 6.5
        )
      )
    },
    width = 6.5,
    height = 4.6
  )

  # ---------------------------------------------------------------------
  # FIGURE 2 — interval shape
  # ---------------------------------------------------------------------
  arm <- bind_rows(spec$comparators(), pl_il) |>
    filter(is.finite(lower), is.finite(upper)) |>
    mutate(
      lo_arm = estimate - lower,
      up_arm = upper - estimate,
      breach = spec$breach_fn(lower, upper),
      method = order_methods(method),
      level_lab = scales::percent(conf_level, accuracy = 1)
    ) |>
    filter(lo_arm > 0, up_arm > 0)

  arm_stats <- arm |>
    filter(conf_level == MAIN_LEVEL) |>
    group_by(method) |>
    summarise(r = median(up_arm / lo_arm), .groups = "drop")

  n_breach <- arm |>
    filter(breach) |>
    count(method, name = "n")

  breach_txt <- if (nrow(n_breach) == 0L) {
    paste0(
      "No interval from any method leaves the parameter space here",
      if (!is.null(spec$floor_note)) paste0(" (", spec$floor_note, ")") else "",
      "."
    )
  } else {
    paste0(
      "Boxed points leave the parameter space (",
      paste(n_breach$method, n_breach$n, sep = ": ", collapse = ", "),
      ")."
    )
  }

  p_arms <- save_figure(
    paste0("ne-", spec$slug, "-interval-shape.png"),
    arm,
    function(df) {
      plot_arm_shape(
        df,
        title = paste0(
          "Interval shape: likelihood intervals are asymmetric, ",
          "Wald intervals are not"
        ),
        subtitle = paste(
          "The Wald benchmarks fall exactly on the identity line by",
          "construction; PL and IL sit above it."
        ),
        caption = wrap_caption(
          paste0(
            "Log-log, equal decades, so the line is a true 45°. ",
            breach_txt, " Median arm ratio at ",
            scales::percent(MAIN_LEVEL, accuracy = 1), ": ",
            paste(
              arm_stats$method, sprintf("%.2f", arm_stats$r),
              sep = " ", collapse = ", "
            ), ".", coverage_note
          ),
          fig_width = 8.0
        )
      )
    },
    width = 8.0,
    height = 4.0
  )

  # ---------------------------------------------------------------------
  # FIGURE 3 — forest plot against the corrected reference
  # ---------------------------------------------------------------------
  site_label <- function(s, n) paste0("Site ", s, "  (N = ", n, ")")

  forest <- pl_il |>
    filter(conf_level == MAIN_LEVEL) |>
    mutate(method = order_methods(method), site_lab = site_label(site, N)) |>
    mutate(site_lab = fct_reorder(site_lab, N))

  ref_forest <- ref |>
    filter(site %in% sites_present) |>
    mutate(site_lab = factor(site_label(site, N),
      levels = levels(forest$site_lab)))

  reach <- forest |>
    left_join(ref |> select(site, psi_ref), by = "site") |>
    group_by(method) |>
    summarise(k = sum(psi_ref >= lower & psi_ref <= upper), .groups = "drop")

  p_forest <- save_figure(
    paste0("ne-", spec$slug, "-forest.png"),
    forest,
    function(df) {
      plot_forest(
        df,
        ref = ref_forest,
        x_label = paste0(
          spec$x_label, if (spec$forest_log) "  (log scale)" else ""
        ),
        title = paste0(
          scales::percent(MAIN_LEVEL, accuracy = 1),
          " intervals against the ", spec$ref_name
        ),
        subtitle = paste0(
          "Reaches the reference at ",
          paste(reach$method, paste0(reach$k, "/", n_present),
            sep = " ", collapse = ", "),
          " sites. Sites ordered by N."
        ),
        caption = wrap_caption(
          paste0(
            "Gray tick: ", spec$ref_name, ". ",
            if (!is.null(spec$floor_note)) {
              paste0("Dashed line: ", spec$floor_note, ". ")
            } else {
              ""
            },
            "Benchmark intervals are omitted — see the interval-shape ",
            "figure.", coverage_note
          ),
          fig_width = 6.5
        ),
        floor_line = spec$floor,
        breaks = spec$forest_breaks,
        log_x = spec$forest_log
      )
    },
    width = 6.5,
    height = if (n_present >= 20) 7.0 else 6.2
  )
}
