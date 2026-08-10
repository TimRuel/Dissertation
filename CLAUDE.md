# CLAUDE.md

Guidance for Claude Code / Positron Assistant when working in this repo.

## Environment notes

### Quarto rendering (Positron, Windows) — resolved 2026-08-06

**Symptom:** Ctrl+Shift+K in Positron on `Dissertation/Dissertation.qmd` failed with
`bash: quarto: command not found`, even though Quarto is installed at
`C:\Users\Tim\AppData\Local\Programs\Quarto\bin\quarto.exe` and is on the system PATH.

**Root cause:** Positron's *own process* environment was captured at launch time (inherited
from Explorer), before the Quarto install directory was on PATH. Every integrated terminal
Positron spawns — regardless of which shell profile — inherits that stale copy. Fully closing
and reopening the Positron *app* does not fix this, because Explorer itself caches an old
environment block and keeps handing it to whatever gets relaunched. Only a fresh process
(a newly opened terminal window, or a full sign-out/reboot) re-reads the registry PATH.

A contributing factor: the default integrated terminal profile was set to **"Rtools Bash"**
(`C:/rtools45/usr/bin/bash.exe`), a shell meant for compiling R packages from source, not
general use. This was not the actual root cause (a manually-spawned Rtools Bash login shell
did have Quarto on PATH), but it was a confusing red herring and an odd choice of default
for everyday work.

**Fix applied:**
1. Changed `terminal.integrated.defaultProfile.windows` from `"Rtools Bash"` to `"Git Bash"`
   in `%APPDATA%\Positron\User\settings.json` (the "Rtools Bash" profile definition itself
   was left in place, just no longer the default).
2. Workaround for the stale-PATH problem when it recurs: fully quit Positron (check Task
   Manager for lingering `Positron.exe` processes), open a **brand-new** terminal window,
   verify `quarto --version` works there, then launch Positron from that fresh window
   (`& "C:\Users\Tim\AppData\Local\Programs\Positron\Positron.exe" c:\Northwestern\Dissertation`).
3. Permanent fix if the stale-PATH issue reappears after some future PATH change: sign out
   and back in (or reboot) to refresh Explorer's cached environment block.

**Relevant settings** (`%APPDATA%\Positron\User\settings.json`):
- `quarto.path`: `C:/Users/Tim/AppData/Local/Programs/Quarto/bin/quarto.exe`
- `terminal.integrated.defaultProfile.windows`: `Git Bash`
- `terminal.integrated.profiles.windows.Rtools Bash`: still defined, available via the
  terminal profile picker for R package compilation if needed.

If `quarto: command not found` shows up again in an integrated terminal, check whether
Positron itself was launched from a stale-PATH parent process before assuming a Quarto
reinstall or config change is needed.

## Experiment results pipeline — built 2026-08-07

How PL/IL results get from Quest into this repo's figures. Spans two repos:

- **`c:\Northwestern\Dissertation`** (this one) — prose, figure code, committed PNGs
- **`c:\Northwestern\likelyr-simulations`** — experiment configs, analyzers, and the
  local `experiments/` results tree
- **`c:\Northwestern\likelyr`** — the package supplying `infer()` / `compare()`

### Design principle: analyze remote, ship only the analysis

A single `model.rds` is 0.5–1 MB, so a simulation experiment (48 sims × 1100 iterations)
is tens of GB. Analysis therefore runs on Quest and only derived tables come down.
The experiment folder is never downloaded.

Results live **only** in `likelyr-simulations/experiments/` (gitignored there). They are
deliberately *not* copied into this repo — one copy of the data, no binary churn in the
dissertation's history. What *is* committed here is the rendered PNGs under
`Dissertation/files/Images`, so the dissertation renders from a clean clone with no VPN,
no Quest account, and no likelyr-simulations checkout. Only *regenerating* figures needs
those.

### Two experiment kinds, same naming convention

Both use `exp_vX` / `sim_XX` and the same on-disk layout, so nothing about the path
distinguishes them. `experiment.kind` in the exp yaml selects the analyzer:

| kind | analyzer | what it computes |
|---|---|---|
| `simulation` (default) | `R/analyze_sim.R` | bias, squared error, coverage — frequency properties aggregated **across iterations** |
| `application` | `R/analyze_app.R` | point estimates, CI endpoints, likelihood curves per real dataset — the interesting axis is **across sims** |

`bin/analyze_sim.sh` reads `kind` from the sim yaml, falls back to the sibling `exp_*.yml`,
then defaults to `simulation`. That fallback means sim yamls generated before `kind:`
existed still dispatch correctly — **no `make gen` re-run needed**.

**Why applications needed a separate analyzer:** `analyze_sim.R` drops `lower`/`upper`
before saving (correct for simulations — only coverage and length matter). For an
application those endpoints *are* the result. Also, `psi_0` doesn't exist for real data;
`applications/multinom/ne_entropy/estimand.R` supplies a placeholder `log(J)/2` only
because `likelyr::estimand_spec()` requires one, and `analyze_app.R` deliberately drops
every column derived from it (`error`, `contains_truth`, and hence bias / sq_error /
covered).

### Commands: where to run what

`$HOME` holds the code; `/projects/p32397` holds the data. This needs no special handling
because **every script reads `exp_dir` from the yaml**, never from its own location. Always
run Make from the checkout and let the config locate the data.

```bash
# ── On Quest, from the code checkout ──
cd /home/tbr0780/likelyr-simulations

# Application (~20 sims × 1 iteration — light, login node is fine)
make results EXP_CONFIG=config/multinom/ne_entropy/exp_v13/exp_v13.yml

# Simulation (48 sims × 1100 iterations ≈ 53k models — MUST use Slurm,
# not a login node: one array task per sim + a dependent bundle job)
make submit-analysis EXP_CONFIG=config/multinom/logit_simpson/exp_v3/exp_v3.yml
```

Either path writes `<exp_dir>/analysis/bundle.rds`. Then locally, VPN active:

```bash
cd /c/Northwestern/likelyr-simulations
make download EXP=multinom/ne_entropy/exp_v13     # MODE=tree for per-sim files
```

One file: ~100 KB for an application, ~5 MB for a simulation. Simulation bundles keep
per-iteration rows rather than pre-aggregating, because `get_coverage_summary()` in
`R/data_viz_utils.R` joins point to interval per iteration for the blended /
regime-stratified analysis.

Finally, in this repo: `Rscript Dissertation/data-viz.R`.

### Figure architecture (three layers)

Figures differ per experiment, so there is no one-size-fits-all script. `data-viz.R` is
just a runner: it sources layers 1–2, then every `figures/*.R` in its own environment
inside its own error handler (one failing script doesn't stop the rest).

| Layer | Files | Changes when? |
|---|---|---|
| 1. Load | `Dissertation/R/load-results.R` | Never — experiment-agnostic |
| 2. Vocabulary | `Dissertation/R/plot-helpers.R`, `R/dune-comparators.R` | Only when you'd write the same plot twice |
| 3. Figures | `Dissertation/figures/*.R` — one per experiment | **Supposed to differ** |

**To add figures for a new experiment: write one new `figures/*.R`.** It's picked up
automatically; nothing in layers 1–2 changes.

Two existing scripts show both patterns:
- `figures/dune-site-intervals.R` — spans several experiments, so it keeps a local
  registry tribble. The registry is that figure's concern, not the architecture's.
- `figures/ne-entropy-curves.R` — reads a single experiment, so it calls `load_bundle()`
  directly. Plots the `psi_loglik_df` grids, which used to never leave Quest.

`save_figure()` refuses to write an empty plot: a figure whose experiment isn't available
yet is skipped (leaving any committed PNG intact) rather than overwritten with an
axis-only panel. **No placeholder data is ever substituted** — the old
`generate_dummy_pl_il()` is gone, so a PNG on disk always reflects real output.

### Interactive exploration

`Dissertation/explore-results.R` is the console scratchpad — `source()` it to get
`available` / `intervals` / `curves` / `context` into your session, plus every plot helper
and the dune benchmarks. It writes nothing and is not a render dependency.

It needs no registry: `list_bundles()` discovers whatever has been downloaded, and
`derive_labels()` recovers the estimand from `context$estimand_name` (recorded by
`analyze_app.R`) and the model from the application directory prefix (`ne_` / `fe_` / `re_`).
Simulation bundles are listed but skipped for reshaping — they hold per-iteration metrics
and want the summarisers in `likelyr-simulations/R/data_viz_utils.R` instead.

The figure scripts still pass estimand/model labels explicitly, because a published
figure's legend shouldn't depend on a string parsed out of a spec file.

Results are read via `LIKELYR_SIMS_DIR` (default `C:/Northwestern/likelyr-simulations`).
Override it when regenerating figures somewhere the Windows path doesn't exist.

### Column vocabulary

`load-results.R` maps the analyzer's output to what the plot functions expect:

| bundle column | figure column |
|---|---|
| `simulation` `"sim_04"` | `site` `4` (the sim index *is* the dune row index) |
| `pseudolikelihood` `"Profile"` / `"Integrated"` | `method` `"PL"` / `"IL"` |
| `alpha` `0.10` | `conf_level` `0.90` |
| `psi_hat` | `estimate` |
| `context$n_obs` | `N` (figures order sites by it) |

Shipping `n_obs` in `app_context.rds` is why the figure scripts never reload `vegan::dune`
just to recover a sample size.

### Open items

- **`ne_simpson` and `fe_entropy` have spec files** under
  `likelyr-simulations/applications/multinom/` but **no experiment configs** under
  `config/multinom/`. Until those run, the fixed-effects and Simpson figure cells are
  skipped.
- **The placeholder `psi_0` is still in the estimand spec.** Nothing downstream is wrong
  (`analyze_app.R` filters it at the boundary), but making it genuinely `NULL` is a
  likelyr-level change: `get_point_estimate_df()` in `R/infer-point_estimate.R` computes
  `error = psi_hat - psi_0` unconditionally.
- **`figures/ne-entropy-curves.R` writes two figures no chapter references yet**
  (`ne-shannon-curves.png`, `-full.png`). Delete the file if they aren't wanted.
- **`fig-simplex-level-sets.R` is standalone** and not run by the `data-viz.R` sweep.
  Move it into `figures/` to include it.
- `likelyr-simulations` has **no CLAUDE.md**; this pipeline knowledge lives only here.

### Verified 2026-08-07 (not just written)

Run against real models in `likelyr-simulations/dev/exp_v13`:
- `analyze_app.R` on dune site 1: `n_obs=18`, `J=30`, `psi_upper=log(30)=3.401`,
  Profile `psi_hat` ≈ `psi_mle` (the profile mode is anchored at the MLE by construction),
  nested CIs, no truth-derived columns.
- `bundle_exp.R` both kinds — application 4 sims → 23.5 KB; simulation `exp_v3` 48 sims →
  4.6 MB, correctly stamping `simulation` onto `invalid_ci_index`, the one table
  `analyze_sim.R` writes without that column.
- Kind dispatch across all three yaml cases (in sim yaml / in exp yaml only / absent).
- Figure runner: helpers visible inside figure scripts, per-script variables do *not* leak
  to global, and a deliberately-failing script was contained while the others completed.
