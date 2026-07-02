#!/usr/bin/env Rscript
# Build a CI-region phylogeny locally (offline), as a stand-in for the IQ-TREE
# `potyvirus_CI.treefile` that the codon pipeline in notebook 02 would produce.
#
# The codon tree needs rentrez CDS retrieval + IQ-TREE (network + heavy). Here we
# instead build a neighbour-joining tree from the CI region of the LOCAL protein
# alignment, whose tip labels (abbr_accession) match data/CI_motif_states.csv so
# the trait-mapping plot in notebook 02 can annotate it.
#
# Output: data/potyvirus_CI.treefile (Newick)

suppressPackageStartupMessages({
  library(Biostrings)
  library(phangorn)
  library(ape)
})

root      <- "C:/Users/chris/R_projects/project-hail-mary"
prot_file <- file.path(root, "data",
                       "OPSR.Poty.Fig3.v16_align_withblasthits_protmafftadd.txt")
out_file  <- file.path(root, "data", "potyvirus_CI.treefile")

aln <- readAAStringSet(prot_file)

# CI region: same fallback window notebook 02 uses when GenBank CI annotation is
# unavailable (aa 2/7 .. 3/7 of the polyprotein alignment).
aa_len   <- unique(width(aln))[1]
ci_start <- floor(aa_len * 2 / 7)
ci_end   <- ceiling(aa_len * 3 / 7)
ci       <- subseq(aln, ci_start, ci_end)

# Drop sequences that are essentially all-gap in this window
frac_gap <- vapply(as.character(ci),
                   function(s) mean(strsplit(s, "")[[1]] %in% c("-", "X")),
                   numeric(1))
ci <- ci[frac_gap < 0.9]

# Character matrix -> phyDat (amino acids) -> ML distance -> NJ
m <- do.call(rbind, strsplit(as.character(ci), ""))
rownames(m) <- names(ci)
pd <- phangorn::phyDat(m, type = "AA")

message(sprintf("Building NJ tree from %d taxa x %d CI columns (aa %d-%d)...",
                nrow(m), ncol(m), ci_start, ci_end))
dm  <- phangorn::dist.ml(pd)
tr  <- ape::ladderize(ape::nj(dm))
tr  <- ape::multi2di(tr)                 # ensure binary
tr$edge.length[tr$edge.length < 0] <- 0  # NJ can yield small negatives

ape::write.tree(tr, out_file)
message(sprintf("Wrote %s (%d tips)", out_file, ape::Ntip(tr)))
