# Project Hail Mary

> DNAzyme-like motifs, horizontal gene transfer, and short linear motif conservation in Potyviridae

This repository is a reproducible research compendium documenting the discovery of a
DNAzyme-like catalytic domain sequence conserved across several potyvirus species, with
flanking regions matching eukaryotic genes across plant, fungal, and animal genomes —
suggestive of cross-kingdom horizontal gene transfer. It also characterises a conserved
motif in the potyviral CI (cylindrical inclusion) helicase implicated in host-sorting
and immune-trafficking interactions (potentially via the AP-4 adaptor complex).

**This is a preprint / work in progress and has not been peer reviewed.** Several of the
central ideas are presented as explicit, testable hypotheses rather than established
results. See the book's [preprint notice](index.qmd) for the full caveat.

## Structure

| Path | Contents |
|---|---|
| `notebooks/` | The book's chapters (Quarto `.qmd`), numbered to match their rendered order |
| `packages/potytools/` | An R package of the reusable, tested analysis functions used across chapters — GenBank parsing, codon usage bias (RSCU/ENC/CAI), phylogenetics, sequence logos, host classification, motif-flank extraction. Being prepared for rOpenSci peer review |
| `scripts/` | One-off orchestration/data-preparation scripts (fetch external data, build an offline fallback tree, etc.) — run occasionally to regenerate a `data/` artifact, not a function library |
| `shiny/bionexus/` | A companion Shiny app for live GBIF + GLOBI host–vector network exploration |
| `data/` | Sequence data, alignments, and cached lookups (not tracked in git; regenerate via the pipeline below) |
| `_targets.R` | The `targets` pipeline orchestrating the computationally intensive steps (alignment, phylogenetic inference, NCBI/ELM fetching) |

## Reading the book

The rendered book isn't currently published (see [Status](#status) below). To read it
locally:

```r
install.packages(c("quarto"))
quarto::quarto_render()  # or: quarto render, from a terminal
```

## Reproducing the analysis

```r
# Install the pipeline dependencies once
install.packages(c("targets", "tarchetypes"))
remotes::install_local("packages/potytools")

# Run the full pipeline (skips up-to-date steps automatically)
targets::tar_make()
```

External tools required on `PATH` for the full pipeline: `mafft`, `iqtree2`, `hmmscan`.
Package dependencies are pinned via [`renv`](https://rstudio.github.io/renv/); run
`renv::restore()` to install the exact versions used.

## potytools

Analysis functions are collected in the `potytools` R package
([`packages/potytools/`](packages/potytools/)), which is rOpenSci-standards-compliant
(documented, tested, CI-checked) and being prepared for rOpenSci peer review.

```r
remotes::install_local("packages/potytools")
library(potytools)
```

## Status

- Book: actively evolving; content and interpretations may change between renders.
- `potytools`: tests passing, R CMD check clean, not yet submitted for review.
- Not yet published to GitHub Pages.

## License

Code in `packages/potytools/` is MIT-licensed (see
[`packages/potytools/LICENSE`](packages/potytools/LICENSE)). Sequence data and other
third-party inputs under `data/` retain their original source licensing/terms — see
[the data-sources appendix](notebooks/A1_data_sources.qmd) for provenance.
