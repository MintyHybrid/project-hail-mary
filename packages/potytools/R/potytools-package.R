#' potytools: Analysis Tools for Potyvirus Genomics and Short Linear Motifs
#'
#' Functions for potyvirus genomic analysis: GenBank parsing, codon-usage
#' bias, host classification, motif-flank extraction and concatenation, and
#' phylogenetic utilities.
#'
#' @keywords internal
#' @import ggplot2
#' @importFrom dplyr mutate filter select arrange group_by summarise left_join n n_distinct bind_rows desc
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom stringr str_split str_detect str_extract
#' @importFrom purrr map map_dfr keep map_chr map_int map2_lgl
#' @importFrom rentrez entrez_fetch entrez_summary
#' @importFrom Biostrings readAAStringSet readDNAStringSet DNAStringSet AAStringSet writeXStringSet subseq translate width
#' @importFrom seqinr uco s2c GC
#' @importFrom ape read.tree write.tree
#' @importFrom cowplot theme_cowplot ggdraw draw_label plot_grid
#' @importFrom ggseqlogo ggseqlogo theme_logo
#' @importFrom stats fisher.test t.test p.adjust chisq.test sd var setNames
#' @importFrom utils write.csv read.csv read.table head tail
"_PACKAGE"

# Non-standard-evaluation column names used inside dplyr/ggplot2 pipelines.
# Declared to satisfy R CMD check ("no visible binding for global variable").
utils::globalVariables(c(
  "position", "information_content", "Letter", "Score"
))
