# potytools

<!-- badges: start -->
[![R-CMD-check](https://github.com/badges/workflow/R-CMD-check/badge.svg)](https://github.com/actions)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![rOpenSci peer-review](https://badges.ropensci.org/XXX_status.svg)](https://github.com/ropensci/software-review/issues/XXX)
<!-- badges: end -->

**potytools** provides a reproducible toolkit for potyvirus genomic analyses,
built around the research in
[Project Hail Mary](https://github.com/).

## What it does

| Module | Functions | Description |
|---|---|---|
| GenBank parsing | `parse_genbank_file()`, `load_genbank_folder()`, `features_as_df()` | Parse `.gb` / `.gbk` files with no Bioconductor dependency |
| Codon usage | `calculate_rscu()`, `calculate_enc()`, `calculate_gc_content()`, `compare_codon_usage()` | RSCU, ENC (Wright 1990), GC by codon position, Fisher-exact comparison |
| Host detection | `detect_host_from_name()`, `classify_host_type()`, `create_host_classification_table()` | Classify potyvirus isolates as monocot- or dicot-infecting |
| ELM motifs | `elm_batch_search()`, `elm_search_api()` | Batch query ELM REST API for short linear motif hits |
| Motif flanks | `extract_motif_flanks()`, `build_concat()` | Extract ± N aa windows around motif hits; build concatenated fragments |
| NCBI fetch | `fetch_custom_sequences()` | Retrieve sequence fragments by accession + coordinates via `rentrez` |

## Installation

```r
# Development version from GitHub
remotes::install_github("yourname/potytools")

# Or from the local packages/ directory in Project Hail Mary
remotes::install_local("packages/potytools")
```

## Quick start

```r
library(potytools)

# Parse a folder of GenBank files
records <- load_genbank_folder("data/genbank/")

# Codon usage of CI protein CDS
cds <- Biostrings::readDNAStringSet("data/ci_cds.fasta")
rscu <- calculate_rscu(cds)
enc  <- sapply(as.character(cds), calculate_enc)

# Host classification
hosts <- create_host_classification_table(names(cds))

# ELM motif search
elm_results <- elm_batch_search("data/ci_flanks_concat.fasta")
```

## rOpenSci

This package is being prepared for [rOpenSci peer review](https://ropensci.org/software-review/).
Contributions and issues welcome.

## Citation

```r
citation("potytools")
```

## License

MIT © chrissi m.
