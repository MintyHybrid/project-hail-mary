# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo is

A **reproducible research compendium** (not a conventional software project) built as a
Quarto **book**. It documents a preprint on DNAzyme-like motifs, horizontal gene transfer,
and short-linear-motif conservation in *Potyviridae*. Three things coexist here:

1. **The book** — Quarto `.qmd` chapters in `notebooks/` plus `index.qmd`, rendered to `_book/`.
2. **`potytools`** — an installable, tested R package in `packages/potytools/` holding the
   reusable analysis functions, being prepared for rOpenSci peer review.
3. **The `targets` pipeline** (`_targets.R`) — orchestrates the expensive steps (NCBI fetch,
   MAFFT alignment, IQ-TREE inference) that feed the chapters.

The book is **preprint / work-in-progress and not peer reviewed**; several core claims are
framed as explicit testable hypotheses, not established results. Preserve that hedging when
editing prose.

## Environment

- Dependencies are pinned with **renv** (`renv.lock`). Run `renv::restore()` to reproduce the
  library. `.Rprofile` auto-activates renv.
- The user develops on **Windows + WSL**; the book preview server (`.claude/launch.json`) runs
  under WSL. Prefer running R and Quarto in the environment the user is currently using.

## Common commands

**Render / preview the book** (from repo root):
```r
quarto::quarto_render()          # or `quarto render` / `quarto preview` from a terminal
```
Chapters use `execute: freeze: auto` — rendered results are cached in `_freeze/` (which **is**
git-tracked) and only re-run when a chapter's code changes.

**Run the pipeline:**
```r
targets::tar_make()              # skips up-to-date targets
targets::tar_visnetwork()        # inspect the DAG
targets::tar_read(target_name)   # load a computed result
```
External CLI tools must be on `PATH` for a full pipeline run: `mafft`, `iqtree2`, `hmmscan`.

**Develop `potytools`** (run from `packages/potytools/`, e.g. `devtools::test()` with that as
the working directory — CI checks the package with `working-directory: packages/potytools`):
```r
devtools::load_all()             # load package for interactive use
devtools::test()                 # run all tests
testthat::test_file("tests/testthat/test-codon_analysis.R")   # single test file
devtools::document()             # regenerate NAMESPACE + man/ from roxygen
devtools::check()                # full R CMD check (must stay clean)
lintr::lint_package()            # config in .lintr
```

**Deploy the book:** `quarto publish gh-pages . --no-browser` from a local machine (renders with
your own Quarto/R/renv versions). The `publish.yml` GitHub Action is a manual (`workflow_dispatch`)
fallback only — deploys are intentionally **not** on push so CI's Quarto version can't diverge
from what was tested locally.

## Architecture notes

- **`potytools` is the function library; chapters and scripts are consumers.** Reusable, tested
  logic (GenBank parsing, codon-usage bias RSCU/ENC/CAI/GC, host mono/dicot classification,
  motif-flank extraction, phylogenetics, sequence logos, plotting themes) belongs in
  `packages/potytools/R/` with a corresponding test in `tests/testthat/`. Chapters should call
  `potytools::` functions, not redefine analysis logic inline. When adding an exported function,
  run `devtools::document()` so `NAMESPACE`/`man/` stay in sync.
- **`scripts/`** holds one-off orchestration/data-prep (fetch external data, build an offline
  fallback tree, run ESMFold/HyPhy/ELM, etc.) — run occasionally to regenerate a `data/`
  artifact. Not a function library and not imported by the package.
- **Not all analysis flows through `_targets.R`.** The ELM short-linear-motif search runs via a
  separate gget-based path (`scripts/prep_elm_regions.R` → `scripts/gget_elm_run.py`, see
  notebook 04), and several scripts are invoked manually. Don't assume a result comes from the
  targets DAG.
- **`data/` and `embeds/` are git-ignored** (large or regenerable, possibly unverified). Do not
  commit them or assume their contents are present; regenerate via the pipeline/scripts.
- **Chapters are numbered to match rendered order** and grouped into Parts in `_quarto.yml`. Keep
  the `_quarto.yml` chapter list, file numbering, and cross-references consistent when adding or
  reordering chapters.
- **HTML-only output by design.** PDF is deliberately disabled (commented out in `_quarto.yml`)
  because the interactive widgets (visNetwork, plotly, r3dmol, DT, Observable) have no PDF form.
  Don't re-enable PDF without adding per-widget static fallbacks.
- **Design system lives in `assets/`, not in the SCSS.** `theme-dark.scss` (primary) and
  `theme-light.scss` define *only* palette tokens as `--hm-*` custom properties plus Bootstrap
  variable overrides; both files use identical token names. Every component rule lives once in
  `assets/hail-mary.css` and reads those tokens. Never hard-code a colour in that CSS, and never
  add a rule to one theme file only — that's how the two silently diverge.
- **The site's navigation layer is plain JS driven by one data file.** `assets/hail-mary-data.js`
  holds the chapter list, the reasoning-graph edges, the per-chapter stepper text and the package
  list; `network-nav.js` (landing graph + package graph) and `workbench.js` (per-chapter stepper,
  downloads) both read it. Add a chapter in `_quarto.yml` *and* in that data file.
- **The chapter workbench is injected client-side on purpose.** Writing it into the `.qmd` files
  would invalidate `_freeze/` and force a re-run of the whole expensive pipeline. Keep it that way
  unless you're prepared to re-render everything.
- **Fonts are self-hosted, not CDN-loaded.** `scripts/fetch_fonts.sh` regenerates
  `assets/fonts.css` + `assets/fonts/` (woff2 + OFL licence texts) from Google's CSS2 response, so
  the `unicode-range` values are authoritative — the chapters use Φ (U+03A6) in YXXΦ, which needs
  the greek subset correctly range-gated. The rendered book makes no third-party requests; don't
  reintroduce a `fonts.googleapis.com` link.
- **Figure backgrounds are the known rough edge.** The transparent-background work in
  `potytools::theme_manuscript()` / `hailmary_palette()` is currently reverted in the working tree,
  so the cached `_freeze/` PNGs are opaque light. `assets/hail-mary.css` compensates with a
  light "plate" behind figures in dark mode — a clearly-marked interim block to delete once the
  chapters are re-rendered with transparent backgrounds. Keep new figures/widgets
  transparent-background so they work in both themes.

## Conventions

- The `potytools` package targets **rOpenSci standards** — keep it documented (roxygen2 markdown),
  tested (testthat edition 3), lint-clean, and `R CMD check`-clean. There are `ropensci-skills`
  available (package-standards, package-review, package-release, etc.) for this work.
- Bibliography lives in `references.bib` (manual) and `packages.bib` (package citations); the
  author/ORCID metadata is Christina Muedsam, ORCID 0000-0001-8488-3442.
- Before deleting any `_freeze/` or knitr `_cache/` directory, commit first — re-rendering is
  expensive and the cache is the record of prior runs.
