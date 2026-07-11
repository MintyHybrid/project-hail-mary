# potytools

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

**potytools** provides a reproducible toolkit for potyvirus genomic analyses,
built around the research in
[Project Hail Mary](https://github.com/MintyHybrid/project-hail-mary).

## What it does

| Module | Functions | Description |
|---|---|---|
| GenBank parsing | `parse_genbank_file()`, `load_genbank_folder()`, `features_as_df()` | Parse `.gb` / `.gbk` files with no Bioconductor dependency |
| Codon usage | `calculate_rscu()`, `calculate_enc()`, `calculate_gc_content()`, `compare_codon_usage()` | RSCU, ENC (Wright 1990), GC by codon position, Fisher-exact comparison |
| Host detection | `detect_host_from_name()`, `classify_host_type()`, `create_host_classification_table()` | Classify potyvirus isolates as monocot- or dicot-infecting |
| Motif flanks | `extract_motif_with_flanks()`, `find_motif_homologs()`, `build_concat()` | Extract ± N aa windows around motif hits; concatenate fragments (with a position map) for short-linear-motif analysis |
| NCBI fetch | `fetch_custom_sequences()` | Retrieve sequence fragments by accession + coordinates via `rentrez` |

## Installation

```r
# Development version from GitHub (package lives in a repo subdirectory)
remotes::install_github("MintyHybrid/project-hail-mary", subdir = "packages/potytools")

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

# Concatenate motif-flank fragments (with a position map) for SLM analysis
concat <- build_concat(c(PVY = "MDFDYSKQ", TuMV = "LKPTGGVE"))
```

## rOpenSci

This package is being prepared for [rOpenSci peer review](https://ropensci.org/software-review/).
Contributions and issues welcome.

## Citation

```r
citation("potytools")
```

## Contributing

Contributions are welcome. Please see
[CONTRIBUTING.md](.github/CONTRIBUTING.md) and note that this project is
released with a [Contributor Code of Conduct](CODE_OF_CONDUCT.md); by
contributing you agree to abide by its terms.

## License

MIT © Christina Muedsam
